// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// OREventAction.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include <memory>
#include <map>
#include <string>
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "G4UserEventAction.hh"
#include "globals.hh"

/// @namespace FW::Geant4:: Namespace for geant4 material mapping
namespace FW {

namespace Geant4 {

	/// @brief Data container of a particle
	// TODO: maybe replace this by an Acts particle
	struct ParticleRecord
	{
		std::array<double, 3> position, momentum;
		int pdg;
		double energy, mass;
		int charge, trackid, parentid;
		std::string volume, process;
	};
	
  class ORSteppingAction;

  /// @class OREventAction
  ///
  /// @brief Writes out particles record
  ///
  class OREventAction : public G4UserEventAction
  {
  public:
    /// Constructor
    OREventAction();

    /// Virtual destructor
    ~OREventAction() override;

    /// Static access method
    static OREventAction*
    Instance();

    /// Interface method for begin of the event
    /// @param event is the G4Event to be processed
    /// @note resets the material step action
    void
    BeginOfEventAction(const G4Event* event) final override;

    /// Interface method for end of event
    /// @param event is the G4Event to be processed
    /// @note this method is writing out the material track records
    void
    EndOfEventAction(const G4Event* event) final override;

    void AddParticle(ParticleRecord& p)
    {
		m_particles[p.trackid].push_back(p);
	}
	
	std::vector<std::vector<ParticleRecord>>
	outcomingParticles() const
	{
		return m_events;
	}
	
    /// Interface method
    /// @note does nothing
    void
    Reset();


  private:
    /// Instance of the EventAction
    static OREventAction* fgInstance;

	std::map<int, std::vector<ParticleRecord>> m_particles;
	std::vector<std::vector<ParticleRecord>> m_events;
  };  


}  // namespace Geant4
}  // namespace FW
