// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/OREventAction.hpp"
#include <stdexcept>
#include "ACTFW/Plugins/Geant4/ORPrimaryGeneratorAction.hpp"
#include "ACTFW/Plugins/Geant4/ORSteppingAction.hpp"
#include "G4Event.hh"
#include "G4RunManager.hh"

FW::Geant4::OREventAction* FW::Geant4::OREventAction::fgInstance = nullptr;

FW::Geant4::OREventAction*
FW::Geant4::OREventAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::Geant4::OREventAction::OREventAction() : G4UserEventAction()
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::Geant4::OREventAction::~OREventAction()
{
  fgInstance = nullptr;
}

void
FW::Geant4::OREventAction::BeginOfEventAction(const G4Event*)
{
  // reset the collection of material steps
  m_particles.clear();
  ORSteppingAction::Instance()->Reset();
}

void
FW::Geant4::OREventAction::EndOfEventAction(const G4Event*)
{
	std::vector<ParticleRecord> outgoingParticles;
	for(std::map<int, std::vector<ParticleRecord>>::iterator it = m_particles.begin(); it != m_particles.end(); it++)
		outgoingParticles.push_back(it->second.back());
	m_events.push_back(outgoingParticles);
}

void
FW::Geant4::OREventAction::Reset()
{
	m_particles.clear();
}