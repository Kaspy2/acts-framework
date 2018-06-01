// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTFW_ALGORITHMS_DIGITIZATIONALGORITHM_H
#define ACTFW_ALGORITHMS_DIGITIZATIONALGORITHM_H

#include <memory>

#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/ProcessCode.hpp"
#include "ACTS/Utilities/Logger.hpp"

namespace Acts {
class PlanarModuleStepper;
}
namespace FW {

class RandomNumbersSvc;

class DigitizationAlgorithm : public FW::BareAlgorithm
{
public:
  struct Config
  {
    /// input hit collection
    std::string simulatedHitsCollection;
    /// output space point collection
    std::string spacePointsCollection;
    /// output clusters collection
    std::string clustersCollection;
    /// FW random number service
    std::shared_ptr<RandomNumbersSvc> randomNumbers = nullptr;
    /// module stepper
    std::shared_ptr<Acts::PlanarModuleStepper> planarModuleStepper = nullptr;
  };

  DigitizationAlgorithm(const Config&        cnf,
                        Acts::Logging::Level level = Acts::Logging::INFO);

  FW::ProcessCode
  execute(FW::AlgorithmContext ctx) const override;

protected:
  Config m_cfg;
};

}  // namespace FW

#endif
