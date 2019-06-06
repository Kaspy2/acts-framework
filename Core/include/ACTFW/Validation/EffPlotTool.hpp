// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <string>
#include "ACTFW/Barcode/Barcode.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Helpers.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace FW {

// Tools to make efficiency plots to show tracking efficiency.
// For the moment, the efficiency is taken as the fraction of successfully
// smoothed track over all tracks
class EffPlotTool
{
  using Identifier  = Acts::GeometryID;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;
  using Track
      = std::vector<Acts::TrackState<Identifier, Acts::BoundParameters>>;

public:
  /// @brief The nested configuration struct
  struct Config
  {
    std::map<std::string, PlotHelpers::Binning> varBinning
        = {{"Eta", PlotHelpers::Binning("#eta", 50, -4, 4)},
           {"Phi", PlotHelpers::Binning("#phi", 100, -3.15, 3.15)},
           {"Pt", PlotHelpers::Binning("pT [GeV/c]", 100, 0, 1000)}};
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  EffPlotTool(const Config&        cfg,
              Acts::Logging::Level level = Acts::Logging::INFO);

  /// Destructor
  ~EffPlotTool();

  /// @brief fill efficiency plots
  /// @param track a vector of trackStates
  /// @param truthParticle a vector of truth hits
  void
  fill(const Track& track, const Data::SimParticle& truthParticle);

  /// @brief write the efficiency plots to file
  void
  write();

private:
  TEfficiency* m_trackeff_vs_eta;  ///< Tracking efficiency vs eta
  TEfficiency* m_trackeff_vs_phi;  ///< Tracking efficiency vs phi
  TEfficiency* m_trackeff_vs_pT;   ///< Tracking efficiency vs pT

  Config                              m_cfg;     ///< The Config class
  std::unique_ptr<const Acts::Logger> m_logger;  ///< The logging instance

  /// The logger
  const Acts::Logger&
  logger() const
  {
    return *m_logger;
  }
};

}  // namespace FW
