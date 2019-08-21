// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>
#include <memory>
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/MaterialMapping/MaterialMapping.hpp"
#include "ACTFW/MaterialMapping/MaterialMappingOptions.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/Json/JsonGeometryConverter.hpp"
#include "ACTFW/Plugins/Json/JsonMaterialWriter.hpp"
#include "ACTFW/Plugins/Root/RootMaterialTrackReader.hpp"
#include "ACTFW/Plugins/Root/RootMaterialWriter.hpp"
#include "ACTFW/Propagation/PropagationOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

namespace po = boost::program_options;

/// @brief The material validation example, it runs a propagation
/// and then writes out the material information
///
/// @tparam option_setup_t Type of the option setter
/// @tparam geometry_setup_t Type of the geometry setter
///
/// @param argc the number of argumetns of the call
/// @param atgv the argument list
/// @param optionsSetup is the access struct to the additional options
/// @param geometrySetup is the access struct for the trackingGeometry
///
template <typename options_setup_t, typename geometry_setup_t>
int
materialMappingExample(int              argc,
                       char*            argv[],
                       options_setup_t  optionsSetup,
                       geometry_setup_t geometrySetup)
{

  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addMaterialMappingOptions(desc);
  FW::Options::addPropagationOptions(desc);
  FW::Options::addInputOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  optionsSetup(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // Get the log level
  auto logLevel = FW::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry  = FW::Geometry::build(vm, geometrySetup);
  auto tGeometry = geometry.first;

  /// Default contexts
  Acts::GeometryContext      geoContext;
  Acts::MagneticFieldContext mfContext;

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  // Straight line stepper
  using SlStepper  = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;
  // Make stepper and propagator
  SlStepper  stepper;
  Propagator propagator(std::move(stepper), std::move(navigator));

  auto matCollection = vm["mat-mapping-collection"].template as<std::string>();

  // ---------------------------------------------------------------------------------
  // Input directory & input file handling
  std::string intputDir   = vm["input-dir"].template as<std::string>();
  auto        intputFiles = vm["input-files"].template as<read_strings>();

  if (vm["input-root"].template as<bool>()) {
    // Read the material step information from a ROOT TTree
    FW::Root::RootMaterialTrackReader::Config matTrackReaderRootConfig;
    if (not matCollection.empty()) {
      matTrackReaderRootConfig.collection = matCollection;
    }
    matTrackReaderRootConfig.fileList = intputFiles;
    auto matTrackReaderRoot
        = std::make_shared<FW::Root::RootMaterialTrackReader>(
            matTrackReaderRootConfig);
    sequencer.addReader(matTrackReaderRoot);
  }

  /// The material mapper
  Acts::SurfaceMaterialMapper::Config smmConfig;
  auto smm = std::make_shared<Acts::SurfaceMaterialMapper>(
      smmConfig,
      std::move(propagator),
      Acts::getDefaultLogger("SurfaceMaterialMapper", logLevel));

  /// The material mapping algorithm
  FW::MaterialMapping::Config mmAlgConfig(geoContext, mfContext);
  mmAlgConfig.materialMapper   = smm;
  mmAlgConfig.trackingGeometry = tGeometry;

  // Get the file name from the options
  std::string materialFileName = vm["mat-output-file"].as<std::string>();

  if (!materialFileName.empty() and vm["output-root"].template as<bool>()) {

    // The writer of the indexed material
    FW::Root::RootMaterialWriter::Config rmwConfig("MaterialWriter");
    rmwConfig.fileName = materialFileName + ".root";
    FW::Root::RootMaterialWriter rmwImpl(rmwConfig);
    // Fullfill the IMaterialWriter interface
    using RootWriter = FW::MaterialWriterT<FW::Root::RootMaterialWriter>;
    mmAlgConfig.materialWriters.push_back(
        std::make_shared<RootWriter>(std::move(rmwImpl)));
  }

  if (!materialFileName.empty() and vm["output-json"].template as<bool>()) {
    /// The name of the output file
    std::string fileName = vm["mat-output-file"].template as<std::string>();
    // the material writer
    FW::Json::JsonGeometryConverter::Config jmConverterCfg(
        "JsonGeometryConverter", Acts::Logging::INFO);
    jmConverterCfg.processSensitives
        = vm["mat-output-sensitives"].template as<bool>();
    jmConverterCfg.processApproaches
        = vm["mat-output-approaches"].template as<bool>();
    jmConverterCfg.processRepresenting
        = vm["mat-output-representing"].template as<bool>();
    jmConverterCfg.processBoundaries
        = vm["mat-output-boundaries"].template as<bool>();
    jmConverterCfg.processVolumes
        = vm["mat-output-volumes"].template as<bool>();
    jmConverterCfg.writeData = vm["mat-output-data"].template as<bool>();
    // The writer
    FW::Json::JsonMaterialWriter jmwImpl(jmConverterCfg,
                                         materialFileName + ".json");
    // Fullfill the IMaterialWriter interface
    using JsonWriter = FW::MaterialWriterT<FW::Json::JsonMaterialWriter>;
    mmAlgConfig.materialWriters.push_back(
        std::make_shared<JsonWriter>(std::move(jmwImpl)));
  }

  // Create the material mapping
  auto mmAlg = std::make_shared<FW::MaterialMapping>(mmAlgConfig);

  // Append the Algorithm
  sequencer.addAlgorithm(mmAlg);

  // Initiate the run
  sequencer.run();
  // Return success code
  return 0;
}
