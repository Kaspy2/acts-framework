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
template <typename bfield_t>
void
setupFittingAlgorithm(bfield_t                                      fieldMap,
                      FW::Sequencer&                                sequencer,
                      po::variables_map&                            vm,
                      std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
                      std::shared_ptr<FW::BarcodeSvc>               barcodeSvc)
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
  fittingConfig.simulatedHitCollection
      = vm["fatras-sim-hits"].template as<std::string>();
  fittingConfig.simulatedEventCollection
      = vm["fatras-sim-particles"].template as<std::string>();
  fittingConfig.trackCollection
      = vm["fitted-tracks"].template as<std::string>();

  // Finally the fitting algorithm
  auto fittingAlgorithm
      = std::make_shared<FittingAlgorithm>(fittingConfig, logLevel);

  // Finalize the squencer setting and run
  sequencer.appendEventAlgorithms({fittingAlgorithm});

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // Write fitted tracks as ROOT files
  std::shared_ptr<FW::Root::RootTrackWriter> tWriterRoot = nullptr;
  if (vm["output-root"].template as<bool>()) {
    FW::Root::RootTrackWriter::Config tWriterRootConfig;
    tWriterRootConfig.trackCollection = fittingConfig.trackCollection;
    tWriterRootConfig.simulatedEventCollection
        = fittingConfig.simulatedEventCollection;
    tWriterRootConfig.filePath
        = FW::joinPaths(outputDir, fittingConfig.trackCollection + ".root");
    tWriterRootConfig.treeName = fittingConfig.trackCollection;
    auto tWriterRoot
        = std::make_shared<FW::Root::RootTrackWriter>(tWriterRootConfig);

    sequencer.addWriters({tWriterRoot});
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
template <typename vmap_t>
void
setupFitting(vmap_t&                                       vm,
             FW::Sequencer&                                sequencer,
             std::shared_ptr<const Acts::TrackingGeometry> tGeometry,
             std::shared_ptr<FW::BarcodeSvc>               barcodeSvc)
{
  // create BField service
  auto bField = FW::Options::readBField<vmap_t>(vm);
  // a charged propagator
  if (bField.first) {
    // create the shared field
    using BField = Acts::SharedBField<Acts::InterpolatedBFieldMap>;
    BField fieldMap(bField.first);
    // now setup of the fitting algorithm and append it to the sequencer
    setupFittingAlgorithm(
        std::move(fieldMap), sequencer, vm, tGeometry, barcodeSvc);
  } else {
    // create the shared field
    using CField = Acts::ConstantBField;
    CField fieldMap(*bField.second);
    // now setup of the fitting algorithm and append it to the sequencer
    setupFittingAlgorithm(
        std::move(fieldMap), sequencer, vm, tGeometry, barcodeSvc);
  }
}
