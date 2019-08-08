// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootVertexAndTracksWriter.hpp"
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

FW::Root::RootVertexAndTracksWriter::RootVertexAndTracksWriter(
    const FW::Root::RootVertexAndTracksWriter::Config& cfg,
    Acts::Logging::Level                               level)
  : VertexAndTracksWriter(cfg.collection, "RootVertexAndTracksWriter", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("vx", &m_ptrVx);
    m_outputTree->Branch("vy", &m_ptrVy);
    m_outputTree->Branch("vz", &m_ptrVz);

    m_outputTree->Branch("d0", &m_ptrD0);
    m_outputTree->Branch("z0", &m_ptrZ0);
    m_outputTree->Branch("phi", &m_ptrPhi);
    m_outputTree->Branch("theta", &m_ptrTheta);
    m_outputTree->Branch("qp", &m_ptrQP);
    m_outputTree->Branch("time", &m_ptrTime);
    m_outputTree->Branch("vtxID", &m_ptrVtxID);
    m_outputTree->Branch("trkCov", &m_ptrTrkCov);
  }
}

FW::Root::RootVertexAndTracksWriter::~RootVertexAndTracksWriter()
{
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode
FW::Root::RootVertexAndTracksWriter::endRun()
{
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote event to tree '" << m_cfg.treeName << "' in '"
                                      << m_cfg.filePath
                                      << "'");
  }
  return ProcessCode::SUCCESS;
}

void
FW::Root::RootVertexAndTracksWriter::ClearAll()
{
  m_vx.clear();
  m_vy.clear();
  m_vz.clear();
  m_d0.clear();
  m_z0.clear();
  m_phi.clear();
  m_theta.clear();
  m_qp.clear();
  m_time.clear();
  m_vtxID.clear();
  m_trkCov.clear();
}

FW::ProcessCode
FW::Root::RootVertexAndTracksWriter::writeT(
    const AlgorithmContext&             context,
    const std::vector<VertexAndTracks>& vertexAndTracksCollection)
{

  if (m_outputFile == nullptr || vertexAndTracksCollection.empty()) {
    return ProcessCode::SUCCESS;
  }

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  ClearAll();

  // Get the event number
  m_eventNr = context.eventNumber;

  for (auto& vertexAndTracks : vertexAndTracksCollection) {

    // Collect the vertex information
    m_vx.push_back(vertexAndTracks.vertex.position.x());
    m_vy.push_back(vertexAndTracks.vertex.position.y());
    m_vz.push_back(vertexAndTracks.vertex.position.z());

    for (auto& track : vertexAndTracks.tracks) {
      // Collect the track information
      m_d0.push_back(track.parameters()[Acts::ParDef::eLOC_D0]);
      m_z0.push_back(track.parameters()[Acts::ParDef::eLOC_Z0]);
      m_phi.push_back(track.parameters()[Acts::ParDef::ePHI]);
      m_theta.push_back(track.parameters()[Acts::ParDef::eTHETA]);
      m_qp.push_back(track.parameters()[Acts::ParDef::eQOP]);
      m_time.push_back(track.parameters()[Acts::ParDef::eT]);
      // Current vertex index as vertex ID
      m_vtxID.push_back(m_vx.size() - 1);

      // Save track covariance
      Acts::BoundSymMatrix cov = *track.covariance();
      std::vector<double>  vec(36);
      Eigen::Map<Acts::BoundSymMatrix>(vec.data(), cov.rows(), cov.cols())
          = cov;

      for (auto i : vec) {
        std::cout << i << std::endl;
      }
      std::cout << std::endl;

      m_trkCov.push_back(vec);
    }
  }

  m_outputTree->Fill();

  return ProcessCode::SUCCESS;
}
