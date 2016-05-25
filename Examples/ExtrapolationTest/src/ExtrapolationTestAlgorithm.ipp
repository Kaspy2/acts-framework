//
//  ExtrapolationTestAlgorithm.ipp
//  ACTS-Development
//
//  Created by Andreas Salzburger on 19/05/16.
//
//

#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/IExtrapolationEngine.hpp"
#include "ACTS/EventData/ParticleDefinitions.hpp"

template <class T> FW::ProcessCode ExtrapolationTestAlgorithm::executeTestT(const T& startParameters) {


    // setup the extrapolation how you'd like it
    Acts::ExtrapolationCell<T> ecc(startParameters);
    //ecc.setParticleHypothesis(m_cfg.particleType);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::StopAtBoundary);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::FATRAS);
    if (m_cfg.collectSensitive)     ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectSensitive);
    if (m_cfg.collectPassive)       ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectPassive);
    if (m_cfg.collectBoundary)      ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectBoundary);
    if (m_cfg.collectMaterial)      ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectMaterial);
    if (m_cfg.sensitiveCurvilinear) ecc.sensitiveCurvilinear = true;    
    
    // force a stop in the extrapoaltion mode
    if (m_cfg.pathLimit > 0.) {
        ecc.pathLimit = m_cfg.pathLimit;
        ecc.addConfigurationMode(Acts::ExtrapolationMode::StopWithPathLimit);
    }
    // screen output
    MSG_DEBUG("===> forward extrapolation - collecting information <<===");
    // call the extrapolation engine
    Acts::ExtrapolationCode eCode = m_cfg.extrapolationEngine->extrapolate(ecc);
    if (eCode.isFailure()){
        MSG_WARNING("Extrapolation failed.");
        return FW::ProcessCode::ABORT;
    }
    // write out if configured
    if (m_cfg.extrapolationCellWriter)
        m_cfg.extrapolationCellWriter->write(ecc);
    // return success
    return FW::ProcessCode::SUCCESS;

}

