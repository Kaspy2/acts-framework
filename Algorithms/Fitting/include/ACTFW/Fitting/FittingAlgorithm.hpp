// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>

#include "ACTFW/Barcode/Barcode.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace FW {

template <typename kalman_Fitter_t>
class FittingAlgorithm : public BareAlgorithm
{
public:
  // A few initialisations and definitionas
  using Identifier  = Acts::GeometryID;
  using TrackState  = Acts::TrackState<Identifier, Acts::BoundParameters>;
  using TrackMap    = std::map<barcode_type, std::vector<TrackState>>;
  using ParticleMap = std::map<barcode_type, Data::SimParticle>;

  /// Nested configuration struct
  struct Config
  {

    Config(kalman_Fitter_t fitter) : kFitter(std::move(fitter)) {}
    /// input hit collection
    std::string simulatedHitCollection = "";
    /// input event collection
    std::string simulatedEventCollection = "";
    /// output track collection
    std::string trackCollection = "";
    /// kalmanFitter instance
    kalman_Fitter_t kFitter;
    /// FW random number service
    std::shared_ptr<RandomNumbersSvc> randomNumberSvc = nullptr;
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  FittingAlgorithm(const Config&        cfg,
                   Acts::Logging::Level level = Acts::Logging::INFO);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  FW::ProcessCode
  execute(const FW::AlgorithmContext& ctx) const final override;

private:
  Config m_cfg;  /// config struct
};

}  // namespace FW

#include "FittingAlgorithm.ipp"
