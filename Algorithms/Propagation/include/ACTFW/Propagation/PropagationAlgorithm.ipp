// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename propagator_t>
std::unique_ptr<Acts::ActsSymMatrixD<5>>
PropagationAlgorithm<propagator_t>::generateCovariance(
    FW::RandomEngine& rnd,
    FW::GaussDist&    gauss) const
{
  if (m_cfg.covarianceTransport && m_cfg.randomNumbers) {
    // we start from the correlation matrix
    auto newCov = std::make_unique<Acts::ActsSymMatrixD<5>>(m_cfg.correlations);
    // then we draw errors according to the error values
    Acts::ActsVectorD<5> covs_smeared = m_cfg.covariances;
    for (size_t k = 0; k < 5; ++k) covs_smeared[k] *= gauss(rnd);
    // and apply a double loop
    for (size_t i = 0; i < 5; ++i)
      for (size_t j = 0; j < 5; ++j) {
        (*newCov)(i, j) *= covs_smeared[i];
        (*newCov)(i, j) *= covs_smeared[j];
      }
    return std::move(newCov);
  }
  return nullptr;
}

template <typename propagator_t>
PropagationAlgorithm<propagator_t>::PropagationAlgorithm(
    const PropagationAlgorithm<propagator_t>::Config& cfg,
    Acts::Logging::Level                              loglevel)
  : BareAlgorithm("PropagationAlgorithm", loglevel), m_cfg(cfg)
{
}

template <typename propagator_t>
ProcessCode
PropagationAlgorithm<propagator_t>::execute(AlgorithmContext context) const
{
  // Create a random number generator
  FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  // Setup random number distributions for some quantities
  FW::GaussDist   d0Dist(0., m_cfg.d0Sigma);
  FW::GaussDist   z0Dist(0., m_cfg.z0Sigma);
  FW::UniformDist phiDist(m_cfg.phiRange.first, m_cfg.phiRange.second);
  FW::UniformDist etaDist(m_cfg.etaRange.first, m_cfg.etaRange.second);
  FW::UniformDist ptDist(m_cfg.ptRange.first, m_cfg.ptRange.second);
  FW::UniformDist qDist(0., 1.);

  Acts::PerigeeSurface surface({0., 0., 0.});

  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(m_cfg.ntests);

  // loop over number of particles
  for (size_t it = 0; it < m_cfg.ntests; ++it) {
    /// get the d0 and z0
    double d0     = d0Dist(rng);
    double z0     = z0Dist(rng);
    double phi    = phiDist(rng);
    double eta    = etaDist(rng);
    double theta  = 2 * atan(exp(-eta));
    double pt     = ptDist(rng);
    double p      = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop    = charge / p;
    // parameters
    Acts::ActsVectorD<5> pars;
    pars << d0, z0, phi, theta, qop;
    // some screen output
    std::unique_ptr<Acts::ActsSymMatrixD<5>> cov = nullptr;

    // execute the test for charged particles
    std::vector<Acts::detail::Step> testSteps;
    if (charge) {
      // charged extrapolation - with hit recording
      Acts::BoundParameters startParameters(
          std::move(cov), std::move(pars), surface);
      testSteps = executeTest<Acts::TrackParameters>(startParameters);
    } else {
      // execute the test for neeutral particles
      Acts::NeutralBoundParameters neutralParameters(
          std::move(cov), std::move(pars), surface);
      testSteps = executeTest<Acts::NeutralParameters>(neutralParameters);
    }
    propagationSteps.push_back(testSteps);
  }

  // write simulated data to the event store
  // - the simulated particles
  if (context.eventStore.add(m_cfg.propagationStepCollection,
                             std::move(propagationSteps))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  return ProcessCode::SUCCESS;
}
