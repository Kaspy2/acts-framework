// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

class TFile;
class TTree;

namespace FW {
namespace Root {

  using Identifier  = Acts::GeometryID;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;
  using TrackState        = Acts::TrackState<Identifier, Acts::BoundParameters>;
  using TrackMap          = std::map<barcode_type, std::vector<TrackState>>;
  using SimParticleVector = std::vector<Data::SimHit<Data::SimParticle>>;

  /// @class RootPerformanceValidation
  ///
  /// Write out the residual and pull of track parameters and efficiency, i.e.
  /// fraction of smoothed track
  //  into output file
  class RootPerformanceValidation final : public WriterT<TrackMap>
  {
  public:
    using Base = WriterT<TrackMap>;
    /// @brief The nested configuration struct
    struct Config
    {
      std::string trackCollection;            ///< track collection to write
      std::string simulatedEventCollection;   ///< truth particle collection
      std::string simulatedHitCollection;     ///< truth hit collection
      std::string filePath;                   ///< path of the output file
      std::string fileMode = "RECREATE";      ///< file access mode
      TFile*      rootFile = nullptr;         ///< common root file
      ResPlotTool::Config resPlotToolConfig;  ///< ResPlotTool config
      EffPlotTool::Config effPlotToolConfig;  ///< EffPlotTool config
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    RootPerformanceValidation(const Config&        cfg,
                              Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootPerformanceValidation() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// @brief Write method called by the base class
    /// @param [in] ctx is the algorithm context for event information
    /// @param [in] tracks are what to be written out
    ProcessCode
    writeT(const AlgorithmContext& ctx, const TrackMap& tracks) final override;

  private:
    Config       m_cfg;                   ///< The config class
    TFile*       m_outputFile{nullptr};   ///< The output file
    int          m_eventNr{0};            ///< The event number
    ResPlotTool* m_resPlotTool{nullptr};  ///< The plot tool for residual&pull
    EffPlotTool* m_effPlotTool{nullptr};  ///< The plot tool for efficiency
  };

}  // namespace Root
}  // namespace FW
