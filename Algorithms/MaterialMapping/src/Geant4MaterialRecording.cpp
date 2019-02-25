// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/MaterialMapping/Geant4MaterialRecording.hpp"
#include <iostream>
#include <stdexcept>
#include "ACTFW/Plugins/Geant4/MMDetectorConstruction.hpp"
#include "ACTFW/Plugins/Geant4/MMEventAction.hpp"
#include "ACTFW/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "ACTFW/Plugins/Geant4/MMRunAction.hpp"
#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include "FTFP_BERT.hh"
#include "G4RunManager.hh"

FW::Geant4MaterialRecording::Geant4MaterialRecording(
    const FW::Geant4MaterialRecording::Config& cnf,
    Acts::Logging::Level                       level)
  : FW::BareAlgorithm("Geant4MaterialRecording", level)
  , m_cfg(cnf)
  , m_runManager(std::make_unique<G4RunManager>())
{

  /// Check if the geometry should be accessed over the geant4 service
  if (m_cfg.geant4Service) {
    m_runManager->SetUserInitialization(m_cfg.geant4Service->geant4Geometry());
  } else if (!m_cfg.gdmlFile.empty()) {
    /// Access the geometry from the gdml file
    ACTS_INFO(
        "received Geant4 geometry from GDML file: " << m_cfg.gdmlFile.c_str());
    FW::Geant4::MMDetectorConstruction* detConstruction
        = new FW::Geant4::MMDetectorConstruction();
    detConstruction->setGdmlInput(m_cfg.gdmlFile.c_str());
    m_runManager->SetUserInitialization(
        detConstruction);  // constructs detector (calls Construct in
                           // Geant4DetectorConstruction)
  } else {
    throw std::invalid_argument("Missing geometry input for Geant4");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new FW::Geant4::MMPrimaryGeneratorAction(
      "geantino", 1000., m_cfg.seed1, m_cfg.seed2));
  FW::Geant4::MMRunAction* runaction = new FW::Geant4::MMRunAction();
  m_runManager->SetUserAction(runaction);
  m_runManager->SetUserAction(new FW::Geant4::MMEventAction());
  m_runManager->SetUserAction(new FW::Geant4::MMSteppingAction());
  m_runManager->Initialize();
}

FW::ProcessCode
FW::Geant4MaterialRecording::execute(FW::AlgorithmContext context) const
{

  /// Begin with the simulation
  m_runManager->BeamOn(m_cfg.tracksPerEvent);
  ///
  std::vector<Acts::RecordedMaterialTrack> recordedMaterialTracks
      = FW::Geant4::MMEventAction::Instance()->MaterialTracks();

  ACTS_INFO("Received " << recordedMaterialTracks.size()
                        << " MaterialTracks. Writing them now onto file...");

  // Write the recorded material to the event store
  if (context.eventStore.add(m_cfg.recordedMaterialCollection,
                             std::move(recordedMaterialTracks))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  return FW::ProcessCode::SUCCESS;
}