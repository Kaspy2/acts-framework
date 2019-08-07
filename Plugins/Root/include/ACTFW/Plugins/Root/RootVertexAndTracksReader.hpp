// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include <vector>
#include "ACTFW/Framework/IReader.hpp"
#include "ACTFW/Framework/IService.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

class TChain;

namespace FW {

namespace Root {

  /// @class RootVertexAndTracksReader
  ///
  /// @brief Reads in vertex and tracks information from a root file
  /// and fills it into a format to be understood by the vertexing algorithms
  class RootVertexAndTracksReader : public IReader
  {
  public:
    /// @brief The nested configuration struct
    struct Config
    {
      std::string outputCollection = "vertexAndTracksCollection";
      std::string treeName = "event";  ///< name of the output tree
      std::vector<std::string> fileList;         ///< The name of the input file

      unsigned int batchSize = 1;  ///!< Batch

      /// The default logger
      std::shared_ptr<const Acts::Logger> logger;

      /// The name of the service
      std::string name;

      /// Constructor
      /// @param lname The name of the Material reader
      /// @parqam lvl The log level for the logger
      Config(const std::string&   lname = "VertexAndTracksReader",
             Acts::Logging::Level lvl   = Acts::Logging::INFO)
        : logger(Acts::getDefaultLogger(lname, lvl)), name(lname)
      {
      }
    };

    /// Constructor
    /// @param cfg The Configuration struct
    RootVertexAndTracksReader(const Config& cfg);

    /// Destructor
    ~RootVertexAndTracksReader();

    /// Framework name() method
    std::string
    name() const final override;

    /// Skip a few events in the IO stream
    /// @param [in] nEvents is the number of skipped events
    ProcessCode
    skip(size_t nEvents) final override;

    /// Read out data from the input stream
    ///
    /// @param context The algorithm context
    ProcessCode
    read(const FW::AlgorithmContext& context) final override;

    /// Return the number of events
    virtual size_t
    numEvents() const final override;

  private:
    /// Private access to the logging instance
    const Acts::Logger&
    logger() const
    {
      return *m_cfg.logger;
    }

    /// The config class
    Config m_cfg;

    /// mutex used to protect multi-threaded reads
    std::mutex m_read_mutex;

    /// The number of events
    size_t m_events = 0;

    /// The input tree name
    TChain* m_inputChain = nullptr;

    int m_eventNr = 0;

    std::vector<float>* m_ptrVx    = new std::vector<float>;
    std::vector<float>* m_ptrVy    = new std::vector<float>;
    std::vector<float>* m_ptrVz    = new std::vector<float>;
    std::vector<float>* m_ptrD0    = new std::vector<float>;
    std::vector<float>* m_ptrZ0    = new std::vector<float>;
    std::vector<float>* m_ptrPhi   = new std::vector<float>;
    std::vector<float>* m_ptrTheta = new std::vector<float>;
    std::vector<float>* m_ptrQP    = new std::vector<float>;
    std::vector<float>* m_ptrTime  = new std::vector<float>;
    std::vector<int>* m_ptrVtxID = new std::vector<int>;

  };

  inline std::string
  RootVertexAndTracksReader::name() const
  {
    return m_cfg.name;
  }

  inline size_t
  RootVertexAndTracksReader::numEvents() const
  {
    return m_events;
  }

  inline ProcessCode RootVertexAndTracksReader::skip(size_t /*nEvents*/)
  {
    return ProcessCode::SUCCESS;
  }

}  // namespace Root
}  // namespace FW
