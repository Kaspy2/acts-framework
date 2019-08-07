// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"

class TFile;
class TTree;

namespace FW {
namespace Root {

  using VertexAndTracksWriter = WriterT<std::vector<VertexAndTracks>>;

  /// Write out vertices together with associated tracks into a TTree
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  ///
  /// A common file can be provided for to the writer to attach his TTree,
  /// this is done by setting the Config::rootFile pointer to an existing file
  ///
  /// Safe to use from multiple writer threads - uses a std::mutex lock.
  class RootVertexAndTracksWriter final : public VertexAndTracksWriter
  {
  public:
    /// @brief The nested configuration struct
    struct Config
    {
      std::string collection;             ///< particle collection to write
      std::string filePath;               ///< path of the output file
      std::string fileMode = "RECREATE";  ///< file access mode
      std::string treeName = "event";     ///< name of the output tree
      std::shared_ptr<FW::BarcodeSvc>
             barcodeSvc;          ///< the barcode service to decode (optional)
      TFile* rootFile = nullptr;  ///< common root file
    };

    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    RootVertexAndTracksWriter(const Config&        cfg,
                              Acts::Logging::Level level = Acts::Logging::INFO);

    /// Virtual destructor
    ~RootVertexAndTracksWriter() override;

    /// End-of-run hook
    ProcessCode
    endRun() final override;

  protected:
    /// @brief Write method called by the base class
    /// @param [in] context is the algorithm context for event information
    /// @param [in] vertexAndTracksCollection is the VertexAndTracks collection
    ProcessCode
    writeT(const AlgorithmContext&             context,
           const std::vector<VertexAndTracks>& vertexAndTracksCollection)
        final override;

  private:
    Config     m_cfg;         ///< The config class
    std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
    TFile*     m_outputFile{nullptr};  ///< The output file
    TTree*     m_outputTree{nullptr};  ///< The output tree
    int        m_eventNr{0};           ///< the event number of

    /// The vertex positions
    std::vector<float> m_vx;
    std::vector<float> m_vy;
    std::vector<float> m_vz;

    /// The track information
    std::vector<float> m_d0;
    std::vector<float> m_z0;
    std::vector<float> m_phi;
    std::vector<float> m_theta;
    std::vector<float> m_qp;
    std::vector<float> m_time;
    std::vector<int> m_vtxID;

    /// Pointers to the vectors
    std::vector<float>* m_ptrVx    = &m_vx;
    std::vector<float>* m_ptrVy    = &m_vy;
    std::vector<float>* m_ptrVz    = &m_vz;
    std::vector<float>* m_ptrD0    = &m_d0;
    std::vector<float>* m_ptrZ0    = &m_z0;
    std::vector<float>* m_ptrPhi   = &m_phi;
    std::vector<float>* m_ptrTheta = &m_theta;
    std::vector<float>* m_ptrQP    = &m_qp;
    std::vector<float>* m_ptrTime  = &m_time;
    std::vector<int>* m_ptrVtxID = &m_vtxID;

    /// @brief Clears all vectors
    void
    ClearAll();
  };

}  // namespace Root
}  // namespace FW
