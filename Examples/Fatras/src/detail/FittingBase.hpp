// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Fitting/FittingAlgorithm.hpp"
#include "ACTFW/Fitting/FittingOptions.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/Root/RootPerformanceValidation.hpp"
#include "ACTFW/Plugins/Root/RootTrackWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GeometryID.hpp"

/// @brief Fitting setup
///
/// @tparam bfield_t Type of the bfield for the fitting to be set up
///
/// @param fieldMap The field map for the fitting setup
/// @param sequencer The framework sequencer
/// @param vm The boost variable map to resolve
/// @param tGeometry The TrackingGeometry for the tracking setup
/// @param barcodesSvc The barcode service to be used for the fitting
/// @param randomNumberSvc The random number service to be used for the
/// fitting
template <typename bfield_t>
void
setupFittingAlgorithm(bfield_t                                      fieldMap,
                      FW::Sequencer&                                sequencer,
                      po::variables_map&                            vm,
                      std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
                      std::shared_ptr<FW::BarcodeSvc>               barcodeSvc,
                      std::shared_ptr<FW::RandomNumbersSvc> randomNumberSvc)
{

  // Read the log level
  Acts::Logging::Level logLevel = FW::Options::readLogLevel(vm);

  // Create a navigator for this tracking geometry
  Acts::Navigator cNavigator(tGeometry);
  cNavigator.resolvePassive   = false;
  cNavigator.resolveMaterial  = true;
  cNavigator.resolveSensitive = true;

  using ChargedStepper    = Acts::EigenStepper<bfield_t>;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Acts::Navigator>;

  ChargedStepper    cStepper(std::move(fieldMap));
  ChargedPropagator cPropagator(std::move(cStepper), std::move(cNavigator));

  using Updator  = Acts::GainMatrixUpdator<Acts::BoundParameters>;
  using Smoother = Acts::GainMatrixSmoother<Acts::BoundParameters>;

  using KalmanFitter = Acts::KalmanFitter<ChargedPropagator, Updator, Smoother>;

  KalmanFitter kFitter(cPropagator);

  using FittingAlgorithm = FW::FittingAlgorithm<KalmanFitter>;

  typename FittingAlgorithm::Config fittingConfig
      = FW::Options::readFittingConfig<po::variables_map, KalmanFitter>(
          vm, kFitter);
  fittingConfig.randomNumberSvc = randomNumberSvc;

  // Finally the fitting algorithm
  sequencer.addAlgorithm(
      std::make_shared<FittingAlgorithm>(fittingConfig, logLevel));

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // Write fitted tracks as ROOT files
  if (vm["output-root"].template as<bool>()) {
    FW::Root::RootTrackWriter::Config tWriterRootConfig;
    tWriterRootConfig.trackCollection = fittingConfig.trackCollection;
    tWriterRootConfig.simulatedEventCollection
        = fittingConfig.simulatedEventCollection;
    tWriterRootConfig.simulatedHitCollection
        = fittingConfig.simulatedHitCollection;
    tWriterRootConfig.filePath
        = FW::joinPaths(outputDir, fittingConfig.trackCollection + ".root");
    tWriterRootConfig.treeName = fittingConfig.trackCollection;
    sequencer.addWriter(
        std::make_shared<FW::Root::RootTrackWriter>(tWriterRootConfig));
  }

  // Write performance plots as ROOT files
  if (vm["output-root"].template as<bool>()) {
    FW::ResPlotTool::Config                     resPlotToolConfig;
    FW::Root::RootPerformanceValidation::Config perfValidationConfig;
    perfValidationConfig.resPlotToolConfig = resPlotToolConfig;
    perfValidationConfig.trackCollection   = fittingConfig.trackCollection;
    perfValidationConfig.simulatedEventCollection
        = fittingConfig.simulatedEventCollection;
    perfValidationConfig.simulatedHitCollection
        = fittingConfig.simulatedHitCollection;
    perfValidationConfig.filePath = FW::joinPaths(
        outputDir, fittingConfig.trackCollection + "_performance.root");
    sequencer.addWriter(std::make_shared<FW::Root::RootPerformanceValidation>(
        perfValidationConfig));
  }
}

/// @brief Fitting setup
///
/// @tparam bfield_t Type of the bfield for the fitting to be set up
///
/// @param fieldMap The field map for the fitting setup
/// @param sequencer The framework sequencer
/// @param vm The boost variable map to resolve
/// @param tGeometry The TrackingGeometry for the tracking setup
/// @param barcodesSvc The barcode service to be used for the fitting
/// @param randomNumberSvc The random number service to be used for the
/// fitting
template <typename vmap_t>
void
setupFitting(vmap_t&                                       vm,
             FW::Sequencer&                                sequencer,
             std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
             std::shared_ptr<FW::BarcodeSvc>               barcodeSvc,
             std::shared_ptr<FW::RandomNumbersSvc>         randomNumberSvc)
{
  // create BField service
  auto bField  = FW::Options::readBField<po::variables_map>(vm);
  auto field2D = std::get<std::shared_ptr<InterpolatedBFieldMap2D>>(bField);
  auto field3D = std::get<std::shared_ptr<InterpolatedBFieldMap3D>>(bField);

  if (field2D) {
    // Define the interpolated b-field
    using BField = Acts::SharedBField<InterpolatedBFieldMap2D>;
    BField fieldMap(field2D);
    // now setup of the fitting and append it to the sequencer
    setupFittingAlgorithm(std::move(fieldMap),
                          sequencer,
                          vm,
                          tGeometry,
                          barcodeSvc,
                          randomNumberSvc);
  } else if (field3D) {
    // Define the interpolated b-field
    using BField = Acts::SharedBField<InterpolatedBFieldMap3D>;
    BField fieldMap(field3D);
    // now setup of the fitting and append it to the sequencer
    setupFittingAlgorithm(std::move(fieldMap),
                          sequencer,
                          vm,
                          tGeometry,
                          barcodeSvc,
                          randomNumberSvc);
  } else if (vm["bf-context-scalable"].template as<bool>()) {
    using SField = FW::BField::ScalableBField;
    SField fieldMap(*std::get<std::shared_ptr<SField>>(bField));
    // now setup of the fitting and append it to the sequencer
    setupFittingAlgorithm(std::move(fieldMap),
                          sequencer,
                          vm,
                          tGeometry,
                          barcodeSvc,
                          randomNumberSvc);
  } else {
    // Create the constant  field
    using CField = Acts::ConstantBField;
    CField fieldMap(*std::get<std::shared_ptr<CField>>(bField));
    // now setup of the fitting and append it to the sequencer
    setupFittingAlgorithm(std::move(fieldMap),
                          sequencer,
                          vm,
                          tGeometry,
                          barcodeSvc,
                          randomNumberSvc);
  }
}
