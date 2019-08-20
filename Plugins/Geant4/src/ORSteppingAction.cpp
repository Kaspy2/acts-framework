// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/ORSteppingAction.hpp"
#include "ACTFW/Plugins/Geant4/OREventAction.hpp"
#include <stdexcept>
#include "Acts/Utilities/Units.hpp"
#include "G4Material.hh"
#include "G4Step.hh"

FW::Geant4::ORSteppingAction* FW::Geant4::ORSteppingAction::fgInstance
    = nullptr;

FW::Geant4::ORSteppingAction*
FW::Geant4::ORSteppingAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::Geant4::ORSteppingAction::ORSteppingAction()
  : G4UserSteppingAction()
// m_volMgr(MaterialRunAction::Instance()->getGeant4VolumeManager())
{
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::Geant4::ORSteppingAction::ORSteppingAction(OREventAction* eventAction) : m_EventAction(eventAction)
{
}


FW::Geant4::ORSteppingAction::~ORSteppingAction()
{
  fgInstance = nullptr;
}

void
FW::Geant4::ORSteppingAction::UserSteppingAction(const G4Step* step)
{
	ParticleRecord p;
	p.position[0] = step->GetPostStepPoint()->GetPosition().x();
	p.position[1] = step->GetPostStepPoint()->GetPosition().y();
	p.position[2] = step->GetPostStepPoint()->GetPosition().z();
	p.momentum[0] = step->GetPostStepPoint()->GetMomentum().x();
	p.momentum[1] = step->GetPostStepPoint()->GetMomentum().y();
	p.momentum[2] = step->GetPostStepPoint()->GetMomentum().z();
	p.pdg = step->GetTrack()->GetDynamicParticle()->GetPDGcode();
	p.energy = step->GetPostStepPoint()->GetTotalEnergy();
	p.mass = step->GetPostStepPoint()->GetMass();
	p.charge = step->GetPostStepPoint()->GetCharge();
	p.trackid = step->GetTrack()->GetTrackID();
	p.parentid = step->GetTrack()->GetParentID();
	m_EventAction->AddParticle(p);
}

void
FW::Geant4::ORSteppingAction::Reset()
{
}
