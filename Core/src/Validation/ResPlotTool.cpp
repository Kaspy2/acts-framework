// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Validation/ResPlotTool.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::theta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

FW::ResPlotTool::ResPlotTool(const FW::ResPlotTool::Config& cfg,
                             Acts::Logging::Level           level)
  : m_cfg(cfg), m_logger(Acts::getDefaultLogger("ResPlotTool", level))
{
}

void
FW::ResPlotTool::book(ResPlotTool::ResPlotCache& resPlotCache) const
{
  PlotHelpers::Binning bEta      = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bR        = m_cfg.varBinning.at("R");
  PlotHelpers::Binning bZ        = m_cfg.varBinning.at("Z");
  PlotHelpers::Binning bResidual = m_cfg.varBinning.at("Residual");
  PlotHelpers::Binning bPull     = m_cfg.varBinning.at("Pull");
  ACTS_DEBUG("Initialize the histograms for residual and pull plots");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    // residual distributions
    resPlotCache.res[parName]
        = PlotHelpers::bookHisto(Form("res_%s", parName.c_str()),
                                 Form("Residual of %s", parName.c_str()),
                                 bResidual);
    // residual vs eta scatter plots
    resPlotCache.res_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_eta", parName.c_str()),
                                 Form("Residual of %s vs eta", parName.c_str()),
                                 bEta,
                                 bResidual);
    // residual mean in each eta bin
    resPlotCache.resmean_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_eta", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bEta);
    // residual width in each eta bin
    resPlotCache.reswidth_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_eta", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bEta);
    // residual vs r scatter plots
    resPlotCache.res_vs_r[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_r", parName.c_str()),
                                 Form("Residual of %s vs r", parName.c_str()),
                                 bR,
                                 bResidual);
    // residual mean in each r bin
    resPlotCache.resmean_vs_r[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_r", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bR);
    // residual width in each r bin
    resPlotCache.reswidth_vs_r[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_r", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bR);
    // residual mean vs z scatter plots
    resPlotCache.res_vs_z[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_z", parName.c_str()),
                                 Form("Residual of %s vs z", parName.c_str()),
                                 bZ,
                                 bResidual);
    // residual mean in each z bin
    resPlotCache.resmean_vs_z[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_z", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bZ);
    // residual width in each z bin
    resPlotCache.reswidth_vs_z[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_z", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bZ);

    // pull distritutions
    resPlotCache.pull[parName]
        = PlotHelpers::bookHisto(Form("pull_%s", parName.c_str()),
                                 Form("Pull of %s", parName.c_str()),
                                 bPull);
    // pull vs eta scatter plots
    resPlotCache.pull_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_eta", parName.c_str()),
                                 Form("Pull of %s vs eta", parName.c_str()),
                                 bEta,
                                 bPull);
    // pull mean in each eta bin
    resPlotCache.pullmean_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_eta", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bEta);
    // pull width in each eta bin
    resPlotCache.pullwidth_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_eta", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bEta);
    // pull vs r scatter plots
    resPlotCache.pull_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_r", parName.c_str()),
                                 Form("Pull of %s vs r", parName.c_str()),
                                 bR,
                                 bPull);
    // pull mean in each r bin
    resPlotCache.pullmean_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_r", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bR);
    // pull width in each r bin
    resPlotCache.pullwidth_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_r", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bR);
    // pull mean vs z scatter plots
    resPlotCache.pull_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_z", parName.c_str()),
                                 Form("Pull of %s vs z", parName.c_str()),
                                 bZ,
                                 bPull);
    // pull mean in each z bin
    resPlotCache.pullmean_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_z", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bZ);
    // pull width in each z bin
    resPlotCache.pullwidth_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_z", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bZ);
  }
}

void
FW::ResPlotTool::clear(ResPlotCache& resPlotCache) const
{
  ACTS_DEBUG("Delete the hists.");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    delete resPlotCache.res.at(parName);
    delete resPlotCache.res_vs_eta.at(parName);
    delete resPlotCache.resmean_vs_eta.at(parName);
    delete resPlotCache.reswidth_vs_eta.at(parName);
    delete resPlotCache.res_vs_r.at(parName);
    delete resPlotCache.resmean_vs_r.at(parName);
    delete resPlotCache.reswidth_vs_r.at(parName);
    delete resPlotCache.res_vs_z.at(parName);
    delete resPlotCache.resmean_vs_z.at(parName);
    delete resPlotCache.reswidth_vs_z.at(parName);
    delete resPlotCache.pull.at(parName);
    delete resPlotCache.pull_vs_eta.at(parName);
    delete resPlotCache.pullmean_vs_eta.at(parName);
    delete resPlotCache.pullwidth_vs_eta.at(parName);
    delete resPlotCache.pull_vs_r.at(parName);
    delete resPlotCache.pullmean_vs_r.at(parName);
    delete resPlotCache.pullwidth_vs_r.at(parName);
    delete resPlotCache.pull_vs_z.at(parName);
    delete resPlotCache.pullmean_vs_z.at(parName);
    delete resPlotCache.pullwidth_vs_z.at(parName);
  }
}

void
FW::ResPlotTool::write(const ResPlotTool::ResPlotCache& resPlotCache) const
{
  ACTS_DEBUG("Write the hists to output file.");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    resPlotCache.res.at(parName)->Write();
    resPlotCache.res_vs_eta.at(parName)->Write();
    resPlotCache.resmean_vs_eta.at(parName)->Write();
    resPlotCache.reswidth_vs_eta.at(parName)->Write();
    resPlotCache.res_vs_r.at(parName)->Write();
    resPlotCache.resmean_vs_r.at(parName)->Write();
    resPlotCache.reswidth_vs_r.at(parName)->Write();
    resPlotCache.res_vs_z.at(parName)->Write();
    resPlotCache.resmean_vs_z.at(parName)->Write();
    resPlotCache.reswidth_vs_z.at(parName)->Write();
    resPlotCache.pull.at(parName)->Write();
    resPlotCache.pull_vs_eta.at(parName)->Write();
    resPlotCache.pullmean_vs_eta.at(parName)->Write();
    resPlotCache.pullwidth_vs_eta.at(parName)->Write();
    resPlotCache.pull_vs_r.at(parName)->Write();
    resPlotCache.pullmean_vs_r.at(parName)->Write();
    resPlotCache.pullwidth_vs_r.at(parName)->Write();
    resPlotCache.pull_vs_z.at(parName)->Write();
    resPlotCache.pullmean_vs_z.at(parName)->Write();
    resPlotCache.pullwidth_vs_z.at(parName)->Write();
  }
}

void
FW::ResPlotTool::fill(ResPlotTool::ResPlotCache&   resPlotCache,
                      const Acts::GeometryContext& gctx,
                      const TrackStateVector&      trackStates,
                      const SimParticleVector&     truthParticles) const
{

  // Get the map of truth hits with geoID
  ACTS_DEBUG("Get the truth hits.");
  std::map<Acts::GeometryID, Data::SimHit<Data::SimParticle>> simHits;
  for (auto& hit : truthParticles) {
    auto geoID = hit.surface->geoID();
    simHits.insert(std::make_pair(geoID, hit));
  }

  // get the distribution of residual/pull
  for (auto& state : trackStates) {
    ParVector_t truthParameter;
    float       truthEta, truthR, truthZ;
    auto        geoID = state.referenceSurface().geoID();
    // get truth parameter at a trackState
    if (simHits.find(geoID) != simHits.end()) {
      Data::SimHit<Data::SimParticle> truthHit = simHits.find(geoID)->second;
      Acts::Vector2D                  hitlocal;
      state.referenceSurface().globalToLocal(
          gctx, truthHit.position, truthHit.direction, hitlocal);
      truthParameter[Acts::ParDef::eLOC_0] = hitlocal.x();
      truthParameter[Acts::ParDef::eLOC_1] = hitlocal.y();
      truthParameter[Acts::ParDef::ePHI]   = phi(truthHit.particle.momentum());
      truthParameter[Acts::ParDef::eTHETA]
          = theta(truthHit.particle.momentum());
      truthParameter[Acts::ParDef::eQOP]
          = truthHit.particle.q() / truthHit.particle.momentum().norm();
      truthEta = eta(truthHit.position);
      truthR   = perp(truthHit.position);
      truthZ   = truthHit.position.z();
    } else {
      ACTS_WARNING("Truth hit for state on "
                   << " : volume = "
                   << geoID.value(Acts::GeometryID::volume_mask)
                   << " : layer = "
                   << geoID.value(Acts::GeometryID::layer_mask)
                   << " : module = "
                   << geoID.value(Acts::GeometryID::sensitive_mask)
                   << " not found!");
      truthParameter[Acts::ParDef::eLOC_0] = -99;
      truthParameter[Acts::ParDef::eLOC_1] = -99;
      truthParameter[Acts::ParDef::ePHI]   = -99;
      truthParameter[Acts::ParDef::eTHETA] = -99;
      truthParameter[Acts::ParDef::eQOP]   = -99;
      truthEta                             = -99;
      truthR                               = -99;
      truthZ                               = -99;
    }

    // get the track paramter and error of track parameter at a trackState
    if (state.parameter.smoothed) {
      auto smoothed       = *state.parameter.smoothed;
      auto trackParameter = smoothed.parameters();
      auto covariance     = *smoothed.covariance();
      // fill the histograms for residual and pull
      for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
        std::string parName  = m_cfg.paramNames.at(parID);
        float       residual = trackParameter[parID] - truthParameter[parID];
        PlotHelpers::fillHisto(resPlotCache.res.at(parName), residual);
        PlotHelpers::fillHisto(
            resPlotCache.res_vs_eta.at(parName), truthEta, residual);
        PlotHelpers::fillHisto(
            resPlotCache.res_vs_r.at(parName), truthR, residual);
        PlotHelpers::fillHisto(
            resPlotCache.res_vs_z.at(parName), truthZ, residual);
        if (covariance(parID, parID) > 0) {
          float pull = residual / sqrt(covariance(parID, parID));
          PlotHelpers::fillHisto(resPlotCache.pull[parName], pull);
          PlotHelpers::fillHisto(
              resPlotCache.pull_vs_eta.at(parName), truthEta, pull);
          PlotHelpers::fillHisto(
              resPlotCache.pull_vs_r.at(parName), truthR, pull);
          PlotHelpers::fillHisto(
              resPlotCache.pull_vs_z.at(parName), truthZ, pull);
        } else {
          ACTS_WARNING("Track parameter :" << parName << " has covariance = "
                                           << covariance(parID, parID)
                                           << " which is smaller than 0 !");
        }
      }
    }
  }  // all states
}

// get the mean and width of residual/pull in each eta bin and fill them into
// histograms
void
FW::ResPlotTool::refinement(ResPlotTool::ResPlotCache& resPlotCache) const
{
  PlotHelpers::Binning bEta = m_cfg.varBinning.at("Eta");
  PlotHelpers::Binning bR   = m_cfg.varBinning.at("R");
  PlotHelpers::Binning bZ   = m_cfg.varBinning.at("Z");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    for (int j = 1; j <= bEta.nBins; j++) {
      TH1D* temp_res = resPlotCache.res_vs_eta.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_eta_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_res,
                            j,
                            resPlotCache.resmean_vs_eta.at(parName),
                            resPlotCache.reswidth_vs_eta.at(parName));

      TH1D* temp_pull = resPlotCache.pull_vs_eta.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_pull,
                            j,
                            resPlotCache.pullmean_vs_eta.at(parName),
                            resPlotCache.pullwidth_vs_eta.at(parName));
    }

    for (int j = 1; j <= bR.nBins; j++) {
      TH1D* temp_res = resPlotCache.res_vs_r.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_r_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_res,
                            j,
                            resPlotCache.resmean_vs_r.at(parName),
                            resPlotCache.reswidth_vs_r.at(parName));

      TH1D* temp_pull = resPlotCache.pull_vs_r.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_r_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_pull,
                            j,
                            resPlotCache.pullmean_vs_r.at(parName),
                            resPlotCache.pullwidth_vs_r.at(parName));
    }

    for (int j = 1; j <= bZ.nBins; j++) {
      TH1D* temp_res = resPlotCache.res_vs_z.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_z_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_res,
                            j,
                            resPlotCache.resmean_vs_z.at(parName),
                            resPlotCache.reswidth_vs_z.at(parName));

      TH1D* temp_pull = resPlotCache.pull_vs_z.at(parName)->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_z_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_pull,
                            j,
                            resPlotCache.pullmean_vs_z.at(parName),
                            resPlotCache.pullwidth_vs_z.at(parName));
    }
  }
}
