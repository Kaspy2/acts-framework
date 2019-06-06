// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootPerformanceValidation.hpp"
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;

FW::Root::RootPerformanceValidation::RootPerformanceValidation(
    const FW::Root::RootPerformanceValidation::Config& cfg,
    Acts::Logging::Level                               level)
  : Base(cfg.trackCollection, "RootPerformanceValidation", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // Input track and truth collection name
  if (m_cfg.trackCollection.empty()) {
    throw std::invalid_argument("Missing input track collection");
  } else if (m_cfg.simulatedEventCollection.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  } else if (m_cfg.simulatedHitCollection.empty()) {
    throw std::invalid_argument("Missing input hit collection");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }

  // Initialize the residual and efficiency plots tool
  m_resPlotTool = new FW::ResPlotTool(m_cfg.resPlotToolConfig, level);
  m_effPlotTool = new FW::EffPlotTool(m_cfg.effPlotToolConfig, level);
  if (m_resPlotTool == nullptr) {
    throw std::bad_alloc();
  } else if (m_effPlotTool == nullptr) {
    throw std::bad_alloc();
  }
}

FW::Root::RootPerformanceValidation::~RootPerformanceValidation()
{
  delete m_resPlotTool;
  delete m_effPlotTool;
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode
FW::Root::RootPerformanceValidation::endRun()
{
  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool->write();
    m_effPlotTool->write();
    ACTS_INFO("Write performance plots to '" << m_cfg.filePath << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootPerformanceValidation::writeT(const AlgorithmContext& ctx,
                                            const TrackMap&         tracks)
{
  if (m_outputFile == nullptr) return ProcessCode::SUCCESS;

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // Read truth particles from input collection
  const std::vector<FW::Data::SimVertex<>>* simulatedEvent = nullptr;
  if (ctx.eventStore.get(m_cfg.simulatedEventCollection, simulatedEvent)
      == FW::ProcessCode::ABORT)
    return FW::ProcessCode::ABORT;
  ACTS_DEBUG("Read collection '" << m_cfg.simulatedEventCollection << "' with "
                                 << simulatedEvent->size()
                                 << " vertices");

  // Read truth hits from input collection
  const FW::DetectorData<geo_id_value, Data::SimHit<Data::SimParticle>>* simHits
      = nullptr;
  if (ctx.eventStore.get(m_cfg.simulatedHitCollection, simHits)
      == FW::ProcessCode::ABORT) {
    throw std::ios_base::failure("Retrieve truth hit collection "
                                 + m_cfg.simulatedHitCollection
                                 + " failure!");
  }
  ACTS_DEBUG("Retrieved hit data '" << m_cfg.simulatedHitCollection
                                    << "' from event store.");
  // Get the map of truth particle
  ACTS_DEBUG("Get the truth particles.");
  std::map<barcode_type, Data::SimParticle> particles;
  for (auto& vertex : *simulatedEvent) {
    for (auto& particle : vertex.outgoing()) {
      particles.insert(std::make_pair(particle.barcode(), particle));
    }
  }

  // Get the map of truth hits on a module
  ACTS_DEBUG("Get the truth hits.");
  std::map<Acts::GeometryID, std::vector<Data::SimHit<Data::SimParticle>>>
      hitsOnModule;
  for (auto& vData : (*simHits)) {
    for (auto& lData : vData.second) {
      for (auto& sData : lData.second) {
        for (auto& hit : sData.second) {
          auto geoID = hit.surface->geoID();
          auto hits  = hitsOnModule.find(geoID);
          if (hits == hitsOnModule.end()) {
            hitsOnModule[geoID]
                = std::vector<Data::SimHit<Data::SimParticle>>();
            hits = hitsOnModule.find(geoID);
          }
          (hits->second).push_back(hit);
        }
      }
    }
  }

  // Loop over the tracks
  for (auto& track : tracks) {
    // find the truth Particle for this track
    Data::SimParticle truthParticle;
    if (particles.find(track.first) != particles.end()) {
      ACTS_DEBUG("Find the truth particle with barcode = " << track.first);
      truthParticle = particles.find(track.first)->second;
    } else {
      ACTS_WARNING("Truth particle with barcode = " << track.first
                                                    << "not found.");
    }

    // find the truth hits for this track
    TruthTrack truthTrack;
    for (auto& state : track.second) {
      // get the geometry ID
      auto geoID = state.referenceSurface().geoID();
      // get all truth hits on this module
      auto hitsOnThisModule = hitsOnModule.find(geoID)->second;
      // lambda to find the truth hit belonging to a given truth track
      barcode_type                    t_barcode = track.first;
      Data::SimHit<Data::SimParticle> truthHit;
      auto                            findTruthHit = [&t_barcode, &truthHit](
          std::vector<Data::SimHit<Data::SimParticle>> hits) -> bool {
        for (auto& hit : hits) {
          if (hit.particle.barcode() == t_barcode) {
            truthHit = hit;
            return true;
          }
        }
        return false;
      };

      // get the truth hit corresponding to this trackState
      if (findTruthHit(hitsOnThisModule)) {
        ACTS_DEBUG("Find the truth hit for trackState on"
                   << " : volume = "
                   << geoID.value(Acts::GeometryID::volume_mask)
                   << " : layer = "
                   << geoID.value(Acts::GeometryID::layer_mask)
                   << " : module = "
                   << geoID.value(Acts::GeometryID::sensitive_mask));
        truthTrack.push_back(truthHit);
      } else {
        ACTS_WARNING("Truth hit for trackState on"
                     << " : volume = "
                     << geoID.value(Acts::GeometryID::volume_mask)
                     << " : layer = "
                     << geoID.value(Acts::GeometryID::layer_mask)
                     << " : module = "
                     << geoID.value(Acts::GeometryID::sensitive_mask)
                     << " not found!");
      }
    }  // all states

    // fill the plots
    m_resPlotTool->fill(track.second, truthTrack);
    m_effPlotTool->fill(track.second, truthParticle);

  }  // all tracks

  // fill residual and pull details into additional hists
  m_resPlotTool->refinement();

  return ProcessCode::SUCCESS;
}
