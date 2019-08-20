// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMMaterialStepAction.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "G4UserSteppingAction.hh"
#include "globals.hh"

namespace FW {
namespace Geant4 {

class OREventAction;
struct ParticleRecord;

  /// @class ORSteppingAction
  ///
  /// @brief Collects the particle
  class ORSteppingAction : public G4UserSteppingAction
  {
  public:
    /// Constructor
    ORSteppingAction();
    
    ORSteppingAction(OREventAction* eventAction);

    /// Destructor
    ~ORSteppingAction() override;

    /// Static access method to the instance
    static ORSteppingAction*
    Instance();

    /// @brief Interface Method doing the step
    /// @note it creates and collects the MaterialInteraction entities
    /// @param step is the Geant4 step of the particle
    void
    UserSteppingAction(const G4Step* step) final override;

    /// Interface reset method
    /// @note it clears the collected step vector
    void
    Reset();

  private:
    /// Instance of the SteppingAction
    static ORSteppingAction* fgInstance;
    
    OREventAction*  m_EventAction;
  };
}  // namespace Geant4
}  // namespace FW
