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
  PlotHelpers::Binning bEta      = m_cfg.varBinning["Eta"];
  PlotHelpers::Binning bR        = m_cfg.varBinning["R"];
  PlotHelpers::Binning bZ        = m_cfg.varBinning["Z"];
  PlotHelpers::Binning bResidual = m_cfg.varBinning["Residual"];
  PlotHelpers::Binning bPull     = m_cfg.varBinning["Pull"];
  ACTS_DEBUG("Initialize the histograms for residual and pull plots");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    // residual distributions
    m_res[parName]
        = PlotHelpers::bookHisto(Form("res_%s", parName.c_str()),
                                 Form("Residual of %s", parName.c_str()),
                                 bResidual);
    // residual vs eta scatter plots
    m_res_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_eta", parName.c_str()),
                                 Form("Residual of %s vs eta", parName.c_str()),
                                 bEta,
                                 bResidual);
    // residual mean in each eta bin
    m_resmean_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_eta", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bEta);
    // residual width in each eta bin
    m_reswidth_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_eta", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bEta);
    // residual vs r scatter plots
    m_res_vs_r[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_r", parName.c_str()),
                                 Form("Residual of %s vs r", parName.c_str()),
                                 bR,
                                 bResidual);
    // residual mean in each r bin
    m_resmean_vs_r[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_r", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bR);
    // residual width in each r bin
    m_reswidth_vs_r[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_r", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bR);
    // residual mean vs z scatter plots
    m_res_vs_z[parName]
        = PlotHelpers::bookHisto(Form("res_%s_vs_z", parName.c_str()),
                                 Form("Residual of %s vs z", parName.c_str()),
                                 bZ,
                                 bResidual);
    // residual mean in each z bin
    m_resmean_vs_z[parName]
        = PlotHelpers::bookHisto(Form("resmean_%s_vs_z", parName.c_str()),
                                 Form("Residual mean of %s", parName.c_str()),
                                 bZ);
    // residual width in each z bin
    m_reswidth_vs_z[parName]
        = PlotHelpers::bookHisto(Form("reswidth_%s_vs_z", parName.c_str()),
                                 Form("Residual width of %s", parName.c_str()),
                                 bZ);

    // pull distritutions
    m_pull[parName]
        = PlotHelpers::bookHisto(Form("pull_%s", parName.c_str()),
                                 Form("Pull of %s", parName.c_str()),
                                 bPull);
    // pull vs eta scatter plots
    m_pull_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_eta", parName.c_str()),
                                 Form("Pull of %s vs eta", parName.c_str()),
                                 bEta,
                                 bPull);
    // pull mean in each eta bin
    m_pullmean_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_eta", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bEta);
    // pull width in each eta bin
    m_pullwidth_vs_eta[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_eta", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bEta);
    // pull vs r scatter plots
    m_pull_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_r", parName.c_str()),
                                 Form("Pull of %s vs r", parName.c_str()),
                                 bR,
                                 bPull);
    // pull mean in each r bin
    m_pullmean_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_r", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bR);
    // pull width in each r bin
    m_pullwidth_vs_r[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_r", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bR);
    // pull mean vs z scatter plots
    m_pull_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pull_%s_vs_z", parName.c_str()),
                                 Form("Pull of %s vs z", parName.c_str()),
                                 bZ,
                                 bPull);
    // pull mean in each z bin
    m_pullmean_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pullmean_%s_vs_z", parName.c_str()),
                                 Form("Pull mean of %s", parName.c_str()),
                                 bZ);
    // pull width in each z bin
    m_pullwidth_vs_z[parName]
        = PlotHelpers::bookHisto(Form("pullwidth_%s_vs_z", parName.c_str()),
                                 Form("Pull width of %s", parName.c_str()),
                                 bZ);
  }
}

FW::ResPlotTool::~ResPlotTool()
{
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    delete m_res[parName];
    delete m_res_vs_eta[parName];
    delete m_resmean_vs_eta[parName];
    delete m_reswidth_vs_eta[parName];
    delete m_res_vs_r[parName];
    delete m_resmean_vs_r[parName];
    delete m_reswidth_vs_r[parName];
    delete m_res_vs_z[parName];
    delete m_resmean_vs_z[parName];
    delete m_reswidth_vs_z[parName];
    delete m_pull[parName];
    delete m_pull_vs_eta[parName];
    delete m_pullmean_vs_eta[parName];
    delete m_pullwidth_vs_eta[parName];
    delete m_pull_vs_r[parName];
    delete m_pullmean_vs_r[parName];
    delete m_pullwidth_vs_r[parName];
    delete m_pull_vs_z[parName];
    delete m_pullmean_vs_z[parName];
    delete m_pullwidth_vs_z[parName];
  }
}

void
FW::ResPlotTool::write()
{
  ACTS_DEBUG("Write the hists to output file.");
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    m_res[parName]->Write();
    m_res_vs_eta[parName]->Write();
    m_resmean_vs_eta[parName]->Write();
    m_reswidth_vs_eta[parName]->Write();
    m_res_vs_r[parName]->Write();
    m_resmean_vs_r[parName]->Write();
    m_reswidth_vs_r[parName]->Write();
    m_res_vs_z[parName]->Write();
    m_resmean_vs_z[parName]->Write();
    m_reswidth_vs_z[parName]->Write();
    m_pull[parName]->Write();
    m_pull_vs_eta[parName]->Write();
    m_pullmean_vs_eta[parName]->Write();
    m_pullwidth_vs_eta[parName]->Write();
    m_pull_vs_r[parName]->Write();
    m_pullmean_vs_r[parName]->Write();
    m_pullwidth_vs_r[parName]->Write();
    m_pull_vs_z[parName]->Write();
    m_pullmean_vs_z[parName]->Write();
    m_pullwidth_vs_z[parName]->Write();
  }
}

void
FW::ResPlotTool::fill(const Track& track, const TruthTrack& truthTrack)
{

  // Get the map of truth hits with geoID
  ACTS_DEBUG("Get the truth hits.");
  std::map<Acts::GeometryID, Data::SimHit<Data::SimParticle>> simHits;
  for (auto& hit : truthTrack) {
    auto geoID = hit.surface->geoID();
    simHits.insert(std::make_pair(geoID, hit));
  }

  // get the distribution of residual/pull
  for (auto& state : track) {
    ParVector_t truthParameter;
    float       truthEta, truthR, truthZ;
    auto        geoID = state.referenceSurface().geoID();
    // get truth parameter at a trackState
    if (simHits.find(geoID) != simHits.end()) {
      Data::SimHit<Data::SimParticle> truthHit = simHits.find(geoID)->second;
      Acts::Vector2D                  hitlocal;
      state.referenceSurface().globalToLocal(
          truthHit.position, truthHit.direction, hitlocal);
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
        PlotHelpers::fillHisto(m_res[parName], residual);
        PlotHelpers::fillHisto(m_res_vs_eta[parName], truthEta, residual);
        PlotHelpers::fillHisto(m_res_vs_r[parName], truthR, residual);
        PlotHelpers::fillHisto(m_res_vs_z[parName], truthZ, residual);
        if (covariance(parID, parID) > 0) {
          float pull = residual / sqrt(covariance(parID, parID));
          PlotHelpers::fillHisto(m_pull[parName], pull);
          PlotHelpers::fillHisto(m_pull_vs_eta[parName], truthEta, pull);
          PlotHelpers::fillHisto(m_pull_vs_r[parName], truthR, pull);
          PlotHelpers::fillHisto(m_pull_vs_z[parName], truthZ, pull);
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
FW::ResPlotTool::refinement()
{
  PlotHelpers::Binning bEta = m_cfg.varBinning["Eta"];
  PlotHelpers::Binning bR   = m_cfg.varBinning["R"];
  PlotHelpers::Binning bZ   = m_cfg.varBinning["Z"];
  for (unsigned int parID = 0; parID < Acts::NGlobalPars; parID++) {
    std::string parName = m_cfg.paramNames.at(parID);
    for (int j = 1; j <= bEta.nBins; j++) {
      TH1D* temp_res = m_res_vs_eta[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_eta_Histo", j), j, j);
      PlotHelpers::anaHisto(
          temp_res, j, m_resmean_vs_eta[parName], m_reswidth_vs_eta[parName]);

      TH1D* temp_pull = m_pull_vs_eta[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      PlotHelpers::anaHisto(temp_pull,
                            j,
                            m_pullmean_vs_eta[parName],
                            m_pullwidth_vs_eta[parName]);
    }

    for (int j = 1; j <= bR.nBins; j++) {
      TH1D* temp_res = m_res_vs_r[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_r_Histo", j), j, j);
      PlotHelpers::anaHisto(
          temp_res, j, m_resmean_vs_r[parName], m_reswidth_vs_r[parName]);

      TH1D* temp_pull = m_pull_vs_r[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_r_Histo", j), j, j);
      PlotHelpers::anaHisto(
          temp_pull, j, m_pullmean_vs_r[parName], m_pullwidth_vs_r[parName]);
    }

    for (int j = 1; j <= bZ.nBins; j++) {
      TH1D* temp_res = m_res_vs_z[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Residual_vs_z_Histo", j), j, j);
      PlotHelpers::anaHisto(
          temp_res, j, m_resmean_vs_z[parName], m_reswidth_vs_z[parName]);

      TH1D* temp_pull = m_pull_vs_z[parName]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_z_Histo", j), j, j);
      PlotHelpers::anaHisto(
          temp_pull, j, m_pullmean_vs_z[parName], m_pullwidth_vs_z[parName]);
    }
  }
}
