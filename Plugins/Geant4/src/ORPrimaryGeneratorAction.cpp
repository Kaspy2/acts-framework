// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/ORPrimaryGeneratorAction.hpp"
#include <stdexcept>
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

FW::Geant4::ORPrimaryGeneratorAction*
    FW::Geant4::ORPrimaryGeneratorAction::fgInstance
    = nullptr;

FW::Geant4::ORPrimaryGeneratorAction::ORPrimaryGeneratorAction(
    const G4String& particleName,
    G4double        energy,
	 G4bool lockAngle,
	 G4double phi,
	 G4double theta,
	 G4bool lockPosition,
	 G4ThreeVector pos,
    G4int           randomSeed1,
    G4int           randomSeed2)
  : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr), m_lockAngle(lockAngle), m_phi(phi), m_theta(theta), m_lockPosition(lockPosition), m_pos(pos)
{
  // configure the run
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
  G4int nofParticles = 1;
  fParticleGun       = std::make_unique<G4ParticleGun>(nofParticles);

  // default particle kinematic
  G4ParticleTable*      particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(energy);
  G4UnitDefinition::PrintUnitsTable();

  // set the random seeds
  CLHEP::HepRandom::getTheEngine()->setSeed(randomSeed1, randomSeed2);
}

FW::Geant4::ORPrimaryGeneratorAction::~ORPrimaryGeneratorAction()
{
  fgInstance = nullptr;
}

FW::Geant4::ORPrimaryGeneratorAction*
FW::Geant4::ORPrimaryGeneratorAction::Instance()
{
  // Static acces function via G4RunManager
  return fgInstance;
}

void
FW::Geant4::ORPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // this function is called at the begining of event
  G4double phi   = m_lockAngle ? m_phi : -M_PI + G4UniformRand() * 2. * M_PI;
  G4double theta = m_lockAngle ? m_theta : G4UniformRand() * M_PI;
  // build a direction
  m_direction
      = G4ThreeVector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  m_position = m_lockPosition ? m_pos : G4ThreeVector(
      0., 0., 0.);  /// @todo make configurable G4RandGauss::shoot(0., 150.));
  // set to the particle gun and
  fParticleGun->SetParticleMomentumDirection(m_direction);
  fParticleGun->SetParticlePosition(m_position);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
