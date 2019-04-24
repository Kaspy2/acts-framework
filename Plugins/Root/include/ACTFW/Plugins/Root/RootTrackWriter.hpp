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
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "TEfficiency.h"

class TFile;
class TTree;

namespace FW {
namespace Root {

  using Identifier  = Acts::GeometryID;
  using Measurement = Acts::
      Measurement<Identifier, Acts::ParDef::eLOC_0, Acts::ParDef::eLOC_1>;
  using TrackState = Acts::TrackState<Identifier, Acts::BoundParameters>;
  using TrackMap   = std::map<barcode_type, std::vector<TrackState>>;
  /// Write out a track (i.e. a vector of trackState at the moment) into a TTree
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  ///
  /// Each entry in the TTree corresponds to one track for optimum writing
  /// speed. The event number is part of the written data.
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  class RootTrackWriter final : public WriterT<TrackMap>
  {
  public:
    using Base = WriterT<TrackMap>;
    /// @brief The nested configuration struct
    struct Config
    {
      std::string trackCollection;           ///< track collection to write
      std::string simulatedEventCollection;  ///< truth particle collection
      std::string filePath;                  ///< path of the output file
      std::string fileMode = "RECREATE";     ///< file access mode
      std::string treeName = "tracks";       ///< name of the output tree
      TFile*      rootFile = nullptr;        ///< common root file
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    RootTrackWriter(const Config&        cfg,
                    Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootTrackWriter() override;

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
    Config     m_cfg;         ///< The config class
    std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
    TFile*     m_outputFile{nullptr};  ///< The output file
    TTree*     m_outputTree{nullptr};  ///< The output tree
    int        m_eventNr{0};           ///< the event number

    unsigned long m_t_barcode{0};  ///< Truth particle barcode
    float         m_t_charge;      ///< Truth particle charge
    float         m_t_vx;          ///< Truth particle vertex x
    float         m_t_vy;          ///< Truth particle vertex y
    float         m_t_vz;          ///< Truth particle vertex z
    float         m_t_ipx;         ///< Truth particle initial momentum px
    float         m_t_ipy;         ///< Truth particle initial momentum py
    float         m_t_ipz;         ///< Truth particle initial momentum pz
    float         m_t_itheta;      ///< Truth particle initial momentum theta
    float         m_t_iphi;        ///< Truth particle initial momentum phi
    float         m_t_ieta;        ///< Truth particle initial momentum eta
    float         m_t_ipT;         ///< Truth particle initial momentum pT

    int                m_nStates{0};  ///< number of states
    std::vector<int>   m_volumeID;    ///< volume identifier
    std::vector<int>   m_layerID;     ///< layer identifier
    std::vector<int>   m_moduleID;    ///< surface identifier
    std::vector<float> m_lx_uncalib;  ///< uncalibrated measurement local x
    std::vector<float> m_ly_uncalib;  ///< uncalibrated measurement local y
    std::vector<float> m_x_uncalib;   ///< uncalibrated measurement global x
    std::vector<float> m_y_uncalib;   ///< uncalibrated measurement global y
    std::vector<float> m_z_uncalib;   ///< uncalibrated measurement global y

    int m_nPredicted{0};          ///< number of states with predicted parameter
    std::vector<bool>  m_prt;     ///< predicted status
    std::vector<float> m_lx_prt;  ///< predicted local x
    std::vector<float> m_ly_prt;  ///< predicted local y
    std::vector<float> m_resid_x_prt;  ///< residual x from predicted
    std::vector<float> m_resid_y_prt;  ///< residual y from predicted
    std::vector<float> m_pull_x_prt;   ///< pull x from predicted
    std::vector<float> m_pull_y_prt;   ///< pull y from predicted
    std::vector<float> m_x_prt;        ///< predicted global x
    std::vector<float> m_y_prt;        ///< predicted global y
    std::vector<float> m_z_prt;        ///< predicted global z
    std::vector<float> m_px_prt;       ///< predicted momentum px
    std::vector<float> m_py_prt;       ///< predicted momentum py
    std::vector<float> m_pz_prt;       ///< predicted momentum pz
    std::vector<float> m_theta_prt;    ///< predicted momentum theta
    std::vector<float> m_eta_prt;      ///< predicted momentum eta
    std::vector<float> m_phi_prt;      ///< predicted momentum phi
    std::vector<float> m_pT_prt;       ///< predicted momentum pT

    int m_nFiltered{0};           ///< number of states with filtered parameter
    std::vector<bool>  m_flt;     ///< filtered status
    std::vector<float> m_lx_flt;  ///< filtered local x
    std::vector<float> m_ly_flt;  ///< filtered local y
    std::vector<float> m_resid_x_flt;  ///< residual x from filtered
    std::vector<float> m_resid_y_flt;  ///< residual y from filtered
    std::vector<float> m_pull_x_flt;   ///< pull x from filtered
    std::vector<float> m_pull_y_flt;   ///< pull y from filtered
    std::vector<float> m_x_flt;        ///< filtered global x
    std::vector<float> m_y_flt;        ///< filtered global y
    std::vector<float> m_z_flt;        ///< filtered global z
    std::vector<float> m_px_flt;       ///< filtered momentum px
    std::vector<float> m_py_flt;       ///< filtered momentum py
    std::vector<float> m_pz_flt;       ///< filtered momentum pz
    std::vector<float> m_theta_flt;    ///< filtered momentum theta
    std::vector<float> m_eta_flt;      ///< filtered momentum eta
    std::vector<float> m_phi_flt;      ///< filtered momentum phi
    std::vector<float> m_pT_flt;       ///< filtered momentum pT

    int m_nSmoothed{0};           ///< number of states with smoothed parameter
    std::vector<bool>  m_smt;     ///< smoothed status
    std::vector<float> m_lx_smt;  ///< smoothed local x
    std::vector<float> m_ly_smt;  ///< smoothed local y
    std::vector<float> m_resid_x_smt;  ///< residual x from smoothed
    std::vector<float> m_resid_y_smt;  ///< residual y from smoothed
    std::vector<float> m_pull_x_smt;   ///< pull x from filtered
    std::vector<float> m_pull_y_smt;   ///< pull y from filtered
    std::vector<float> m_x_smt;        ///< smoothed global x
    std::vector<float> m_y_smt;        ///< smoothed global y
    std::vector<float> m_z_smt;        ///< smoothed global z
    std::vector<float> m_px_smt;       ///< smoothed momentum px
    std::vector<float> m_py_smt;       ///< smoothed momentum py
    std::vector<float> m_pz_smt;       ///< smoothed momentum pz
    std::vector<float> m_theta_smt;    ///< smoothed momentum theta
    std::vector<float> m_eta_smt;      ///< smoothed momentum eta
    std::vector<float> m_phi_smt;      ///< smoothed momentum phi
    std::vector<float> m_pT_smt;       ///< smoothed momentum pT
  };

}  // namespace Root
}  // namespace FW
