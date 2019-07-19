// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Simulation/EventToTrackConverter.hpp"

#include <iostream>
#include <stdexcept>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

FW::EventToTrackConverterAlgorithm::EventToTrackConverterAlgorithm(
    const FW::EventToTrackConverterAlgorithm::Config& cfg,
    Acts::Logging::Level                              level)
  : FW::BareAlgorithm("EventToTrackConverterAlgorithm", level), m_cfg(cfg)
{
  if (m_cfg.inputCollection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.outputCollection.empty()) {
    throw std::invalid_argument("Missing output collection");
  } else if (m_cfg.randomNumberSvc == nullptr) {
    throw std::invalid_argument("Missing random number service");
  }
}

FW::ProcessCode
FW::EventToTrackConverterAlgorithm::execute(
    const AlgorithmContext& context) const
{

  const auto& vertexCollection
      = context.eventStore.get<std::vector<FW::Data::SimVertex<>>>(
          m_cfg.inputCollection);

  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface
      = Acts::Surface::makeShared<Acts::PerigeeSurface>(m_cfg.refPosition);

  // Set up stepper
  Acts::EigenStepper<Acts::ConstantBField> stepper(m_cfg.bField);

  // Set up propagator with void navigator
  Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>> propagator(
      stepper);

  // Set up propagator options
  Acts::PropagatorOptions<> pOptions(context.geoContext,
                                     context.magFieldContext);
  pOptions.direction = Acts::backward;

  // Create random number generator and spawn gaussian distribution
  FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  // Vector to store tracks extracted from event
  std::vector<Acts::BoundParameters> trackCollection;

  // Start looping over all vertices in current event
  for (auto& vtx : vertexCollection) {

    // Iterate over all particle emerging from current vertex
    for (auto const& particle : vtx.out) {
      const Acts::Vector3D& ptclMom = particle.momentum();

      // Define start track params
      Acts::CurvilinearParameters start(
          nullptr, particle.position(), ptclMom, particle.q(), 0.);
      // Run propagator
      auto result = propagator.propagate(start, *perigeeSurface, pOptions);
      if (!result.ok()) {
        continue;
      }

      // get perigee parameters
      const auto& perigeeParameters = (*result).endParameters->parameters();

      auto newTrackParams = perigeeParameters;

      if (m_cfg.doSmearing) {

        // Calculate pt-dependent IP resolution
        const double particlePt
            = Acts::VectorHelpers::perp(ptclMom) / Acts::units::_GeV;
        const double ipRes = m_cfg.ipResA * std::exp(-m_cfg.ipResB * particlePt)
            + m_cfg.ipResC;
        // except for IP resolution, following variances are rough guesses
        // Gaussian distribution for IP resolution
        std::normal_distribution<double> gaussDist_IP(0., ipRes);
        // Gaussian distribution for angular resolution
        std::normal_distribution<double> gaussDist_angular(0., 0.1);
        // Gaussian distribution for q/p (momentum) resolution
        std::normal_distribution<double> gaussDist_qp(
            0., 0.1 * perigeeParameters[4]);

        double rn_d0 = gaussDist_IP(rng);
        double rn_z0 = gaussDist_IP(rng);
        double rn_ph = gaussDist_angular(rng);
        double rn_th = gaussDist_angular(rng);
        double rn_qp = gaussDist_qp(rng);

        Acts::TrackParametersBase::ParVector_t smrdParamVec;
        smrdParamVec << rn_d0, rn_z0, rn_ph, rn_th, rn_qp, 0.;

        // Update track parameters
        newTrackParams += smrdParamVec;
        // Correct for phi and theta wrap
        correctPhiThetaPeriodicity(newTrackParams[0], newTrackParams[1]);

        // Update track covariance
        std::unique_ptr<Acts::BoundSymMatrix> covMat
            = std::make_unique<Acts::BoundSymMatrix>();
        covMat->setZero();
        covMat->diagonal() << rn_d0 * rn_d0, rn_z0 * rn_z0, rn_ph * rn_ph,
            rn_th * rn_th, rn_qp * rn_qp, 1.;

        trackCollection.push_back(Acts::BoundParameters(context.geoContext,
                                                        std::move(covMat),
                                                        newTrackParams,
                                                        perigeeSurface));
      } else {
        trackCollection.push_back(Acts::BoundParameters(
            context.geoContext, nullptr, newTrackParams, perigeeSurface));
      }

    }  // end iteration over all particle at vertex
  }    // end iteration over all vertices

  // write the SpacePoints to the EventStore
  context.eventStore.add(m_cfg.outputCollection, std::move(trackCollection));

  return FW::ProcessCode::SUCCESS;
}

void
FW::EventToTrackConverterAlgorithm::correctPhiThetaPeriodicity(
    double& phiIn,
    double& thetaIn) const
{
  double tmpPhi = std::fmod(phiIn, 2 * M_PI);  // temp phi
  if (tmpPhi > M_PI) {
    tmpPhi -= 2 * M_PI;
  }
  if (tmpPhi < -M_PI && tmpPhi > -2 * M_PI) {
    tmpPhi += 2 * M_PI;
  }

  double tmpTht = std::fmod(thetaIn, 2 * M_PI);  // temp theta
  if (tmpTht < -M_PI) {
    tmpTht = std::abs(tmpTht + 2 * M_PI);
  } else if (tmpTht < 0) {
    tmpTht *= -1;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
  }
  if (tmpTht > M_PI) {
    tmpTht = 2 * M_PI - tmpTht;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
  }

  phiIn   = tmpPhi;
  thetaIn = tmpTht;
}
