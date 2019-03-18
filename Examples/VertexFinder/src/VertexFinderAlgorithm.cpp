// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "VertexFinderAlgorithm.hpp"
#include "ACTFW/Random/RandomNumberDistributions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

#include <tbb/tbb.h>

#include "Acts/Vertexing/VertexEventData/LinearizedTrack.hpp"
#include "Acts/Vertexing/VertexEventData/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterUtils/LinearizedTrackFactory.hpp"

FWE::VertexFinderAlgorithm::VertexFinderAlgorithm(const Config&        cfg,
                                                  Acts::Logging::Level level)
  : FW::BareAlgorithm("VertexFinder", level), m_cfg(cfg)
{
}

/// @brief Algorithm that receives an Evgen input event, runs over all
/// vertices and smears corresponding tracks.
/// Track collections belonging to a certain vertex (truth-based vertex finder
/// emuluation)
/// are then passed to vertex fitter to fit vertex position.
FW::ProcessCode
FWE::VertexFinderAlgorithm::execute(FW::AlgorithmContext context) const
{

  std::cout << "START VertexFinderAlgorithm" << std::endl;

  const double eta_cut = 3.0;

  /// Define parameter for pt-dependent IP resolution
  /// of the form sigma_d/z(p_t[GeV]) = A*exp(-B*p_t[GeV]) + C
  const double ipResA = 100.7439 * Acts::units::_um;
  const double ipResB = 0.23055;
  const double ipResC = 20. * Acts::units::_um;

  /// Create and fill input event
  const std::vector<FW::Data::SimVertex<>>* inputEvent = nullptr;
  if (context.eventStore.get(m_cfg.collection, inputEvent)
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  /// Define perigee surface center coordinates
  const Acts::Vector3D surfaceCenter(0., 0., 0.);

  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface
      = Acts::Surface::makeShared<Acts::PerigeeSurface>(surfaceCenter);

  /// Set up stepper
  Acts::EigenStepper<Acts::ConstantBField> stepper(m_cfg.bField);

  /// Set up propagator with void navigator
  Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>> propagator(
      stepper);

  /// Set up propagator options
  Acts::PropagatorOptions<> options;

  /// Create random number generator and spawn gaussian distribution
  FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  /// typedefs for simplicity
  std::vector<Acts::BoundParameters> trackCollection;

  /// Start looping over all vertices in current event
  for (auto& vtx : (*inputEvent)) {

    /// Iterate over all particle emerging from current vertex
    for (auto const& particle : vtx.out) {

      const Acts::Vector3D& ptclMom = particle.momentum();
      /// Calculate pseudo-rapidity
      const double eta = Acts::VectorHelpers::eta(ptclMom);
      /// Only charged particles for |eta| < 2.5
      if (particle.q() != 0 && std::abs(eta) < eta_cut) {
        /// Define start track params
        Acts::CurvilinearParameters start(
            nullptr, particle.position(), ptclMom, particle.q());

        /// Run propagator
        const auto result
            = propagator.propagate(start, *perigeeSurface, options);

        if (result.status == Acts::PropagatorStatus::SUCCESS) {

          // get perigee parameters
          const auto& perigeeParameters = result.endParameters->parameters();

          if (std::abs(perigeeParameters[0]) > 30
              || std::abs(perigeeParameters[1]) > 200) {
            continue;
          }

          /// Calculate pt-dependent IP resolution
          const double pclPt
              = Acts::VectorHelpers::perp(ptclMom) / Acts::units::_GeV;
          const double ipRes = ipResA * std::exp(-ipResB * pclPt) + ipResC;

          /// except for IP resolution, following variances are rough guesses
          /// Gaussian distribution for IP resolution
          FW::GaussDist gaussDist_IP(0., ipRes);
          /// Gaussian distribution for angular resolution
          FW::GaussDist gaussDist_angular(0., 0.1);
          /// Gaussian distribution for q/p (momentum) resolution
          FW::GaussDist gaussDist_qp(0., 0.1 * perigeeParameters[4]);

          double rn_d0 = gaussDist_IP(rng);
          double rn_z0 = gaussDist_IP(rng);
          double rn_ph = gaussDist_angular(rng);
          double rn_th = gaussDist_angular(rng);
          double rn_qp = gaussDist_qp(rng);

          double smrd_d0    = perigeeParameters[0] + rn_d0;
          double smrd_z0    = perigeeParameters[1] + rn_z0;
          double smrd_phi   = perigeeParameters[2] + rn_ph;
          double smrd_theta = perigeeParameters[3] + rn_th;
          double srmd_qp    = perigeeParameters[4] + rn_qp;

          /// smearing can bring theta out of range ->close to beam line ->
          /// discard
          if (smrd_theta < 0 || smrd_theta > M_PI) {
            continue;
          }

          double new_eta = -log(tan(smrd_theta / 2));
          if (std::abs(new_eta) > eta_cut) continue;

          Acts::TrackParametersBase::ParVector_t paramVec;
          paramVec << smrd_d0, smrd_z0, smrd_phi, smrd_theta, srmd_qp;

          /// Fill vector of smeared tracks
          std::unique_ptr<Acts::ActsSymMatrixD<5>> covMat
              = std::make_unique<Acts::ActsSymMatrixD<5>>();
          covMat->setZero();
          (*covMat)(0, 0) = rn_d0 * rn_d0;
          (*covMat)(1, 1) = rn_z0 * rn_z0;
          (*covMat)(2, 2) = rn_ph * rn_ph;
          (*covMat)(3, 3) = rn_th * rn_th;
          (*covMat)(4, 4) = rn_qp * rn_qp;

          Acts::BoundParameters currentBoundParams(
              std::move(covMat), paramVec, perigeeSurface);

          trackCollection.push_back(currentBoundParams);
        }
      }
    }
  }

  VertexFinder::State state;

  // find vertices
  Acts::VertexingStatus status
      = m_cfg.vertexFinder->find(trackCollection, state);

  std::vector<Acts::Vertex<Acts::BoundParameters>> fittedVertices;

  if (status == Acts::VertexingStatus::SUCCESS) {
    fittedVertices = state.vertexCollection;
  }

  int vtxCount = 0;
  for (auto& vtx : fittedVertices) {
    vtxCount++;
    Acts::Vector3D pos = vtx.position();
    std::cout << "Reconstructed vertex " << vtxCount << " at position:"
              << "(" << pos[Acts::eX] << "," << pos[Acts::eY] << ","
              << pos[Acts::eZ] << "). Number of tracks: " << vtx.tracks().size()
              << std::endl;
  }

  if (context.eventStore.add(m_cfg.collectionOut, std::move(fittedVertices))
      != FW::ProcessCode::SUCCESS) {
    return FW::ProcessCode::ABORT;
  }

  return FW::ProcessCode::SUCCESS;
}