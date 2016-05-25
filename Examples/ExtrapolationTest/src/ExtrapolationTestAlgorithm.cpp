#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Framework/MsgStreamMacros.hpp"
#include "ACTFW/Random/RandomNumbers.hpp"
#include "ACTFW/Writers/IExtrapolationCellWriter.hpp"
#include "ExtrapolationTestAlgorithm.hpp"
#include <iostream>

FWE::ExtrapolationTestAlgorithm::ExtrapolationTestAlgorithm(const FWE::ExtrapolationTestAlgorithm::Config& cfg) :
    FW::Algorithm(cfg),
    m_cfg(cfg)
{}

FWE::ExtrapolationTestAlgorithm::~ExtrapolationTestAlgorithm()
{}

/** Framework finalize mehtod */
FW::ProcessCode FWE::ExtrapolationTestAlgorithm::initialize(std::shared_ptr<FW::WhiteBoard> eStore,
                                                            std::shared_ptr<FW::WhiteBoard> jStore)
{
    // call the algorithm initialize for setting the stores
    if ( FW::Algorithm::initialize(eStore,jStore) != FW::ProcessCode::SUCCESS){
        MSG_FATAL("Algorithm::initialize() did not succeed!");
        return FW::ProcessCode::SUCCESS;
    }
    MSG_VERBOSE("initialize successful.");
    return FW::ProcessCode::SUCCESS;
}

/** Framework execode method */
FW::ProcessCode FWE::ExtrapolationTestAlgorithm::execute(size_t eventNumber)
{

    // loop
    for (size_t iex = 0; iex < m_cfg.testsPerEvent; ++iex){
        // gaussian d0 and z0
        double d0  = drawGauss(m_cfg.d0Defs);
        double z0  = drawGauss(m_cfg.z0Defs);
        double phi = drawUniform(m_cfg.phiRange);
        double eta = drawUniform(m_cfg.etaRange);
        double theta = 2.*atan(exp(-eta));
        double pt  = drawUniform(m_cfg.ptRange);
        double p   = pt/sin(theta);
        double q   = drawUniform({{0.,1.}}) > 0.5 ? 1. : -1.;
        
        Acts::Vector3D momentum(p*sin(theta)*cos(phi), p*sin(theta)*sin(phi), p*cos(theta));
        std::unique_ptr<Acts::ActsSymMatrixD<5> > cov;
        Acts::ActsVectorD<5> pars; pars << d0, z0, phi, theta, q/p;
        // perigee parameters
        MSG_VERBOSE("Building parameters from Perigee with (" << d0 << ", " << z0 << ", " << phi << ", " << theta << ", " << q/p);
        // charged extrapolation
        Acts::PerigeeSurface pSurface(Acts::Vector3D(0.,0.,0.));
        
        // neutral extrapolation
        if (m_cfg.parameterType){
            Acts::BoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
            if (executeTestT<Acts::TrackParameters>(startParameters) != FW::ProcessCode::SUCCESS)
                MSG_WARNING("Test of parameter extrapolation did not succeed.");
        
        } else {
            // charged extrapolation
            Acts::NeutralBoundParameters startParameters(std::move(cov),std::move(pars),pSurface);
            if (executeTestT<Acts::NeutralParameters>(startParameters) != FW::ProcessCode::SUCCESS)
                MSG_WARNING("Test of parameter extrapolation did not succeed.");
       }
        
    }
    // return SUCCESS to the frameword
    return FW::ProcessCode::SUCCESS;
}

/** Framework finalize mehtod */
FW::ProcessCode FWE::ExtrapolationTestAlgorithm::finalize()
{
    MSG_VERBOSE("initialize successful.");
    return FW::ProcessCode::SUCCESS;
}

double FWE::ExtrapolationTestAlgorithm::drawGauss(const std::array<double,2>& pars) const
{
    double mean  = pars[0];
    double sigma = pars[1];
    return mean+m_cfg.randomNumbers->draw(FW::Distribution::gauss)*sigma;
}

double FWE::ExtrapolationTestAlgorithm::drawUniform(const std::array<double,2>& range) const
{
    double low  = range[0];
    double high = range[1];
    double delta = high-low;
    return low+m_cfg.randomNumbers->draw(FW::Distribution::uniform)*delta;
}


