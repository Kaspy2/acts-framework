// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <boost/program_options.hpp>

#include "ACTFW/Digitization/DigitizationOptions.hpp"
#include "ACTFW/EventData/Barcode.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"

#include "FatrasDigitizationBase.hpp"
#include "FatrasEvgenBase.hpp"
#include "FatrasSimulationBase.hpp"

/// @brief The Fatras example
///
/// This instantiates the geometry and runs fast track simultion
///
/// @tparam options_setup_t are the callable example options
/// @tparam geometry_setup_t Type of the geometry getter struct
///
/// @param argc the number of argumetns of the call
/// @param argv the argument list
/// @param optionsSetup is a callable options struct
/// @param geometrySetup is a callable geometry getter
template <typename options_setup_t, typename geometry_setup_t>
int
fatrasExample(int               argc,
              char*             argv[],
              options_setup_t&  optionsSetup,
              geometry_setup_t& geometrySetup)
{
  using boost::program_options::value;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addParticleGunOptions(desc);
  FW::Options::addPythia8Options(desc);
  FW::Options::addRandomNumbersOptions(desc);
  FW::Options::addBFieldOptions(desc);
  FW::Options::addFatrasOptions(desc);
  FW::Options::addDigitizationOptions(desc);
  FW::Options::addOutputOptions(desc);
  desc.add_options()("evg-input-type",
                     value<std::string>()->default_value("pythia8"),
                     "Type of evgen input 'gun', 'pythia8'");
  // Add specific options for this geometry
  optionsSetup(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  auto logLevel = FW::Options::readLogLevel(vm);

  // Create the random number engine
  auto randomNumberSvcCfg = FW::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc
      = std::make_shared<FW::RandomNumbers>(randomNumberSvcCfg);
  // Create the barcode service
  auto barcodeSvc = std::make_shared<FW::BarcodeSvc>(FW::BarcodeSvc::Config());

  // The geometry, material and decoration
  auto geometry          = FW::Geometry::build(vm, geometrySetup);
  auto tGeometry         = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) { sequencer.addContextDecorator(cdr); }

  // (A) EVGEN
  // Setup the evgen input to the simulation
  setupEvgenInput(vm, sequencer, barcodeSvc, randomNumberSvc);

  // (B) SIMULATION
  // Setup the simulation
  setupSimulation(vm, sequencer, tGeometry, barcodeSvc, randomNumberSvc);

  // (C) DIGITIZATION
  // Setup the digitization
  setupDigitization(vm, sequencer, barcodeSvc, randomNumberSvc);

  // (D) TRUTH TRACKING

  // (E) PATTERN RECOGNITION

  return sequencer.run();
}
