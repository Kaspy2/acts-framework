// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Root/RootTrackWriter.hpp"
#include <TFile.h>
#include <TTree.h>
#include <ios>
#include <stdexcept>
#include "Acts/Utilities/Helpers.hpp"
//#include <boost/variant.hpp>
//#include <boost/variant/static_visitor.hpp>
 
using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::theta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

FW::Root::RootTrackWriter::RootTrackWriter(
                                           const FW::Root::RootTrackWriter::Config& cfg,
                                           Acts::Logging::Level                        level)
  : Base(cfg.collection, "RootTrackWriter", level)
  , m_cfg(cfg)
  , m_outputFile(cfg.rootFile)
{
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  if (m_outputFile == nullptr) {
    m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
    if (m_outputFile == nullptr) {
      throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
    }
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr)
    throw std::bad_alloc();
  else {
    // I/O parameters
    m_outputTree->Branch("event_nr", &m_eventNr);
    m_outputTree->Branch("nStates", &m_nStates);
    m_outputTree->Branch("volume_id", &m_volumeID);
    m_outputTree->Branch("layer_id", &m_layerID);
    m_outputTree->Branch("module_id", &m_moduleID);
    m_outputTree->Branch("l_x_uncalib", &m_lx_uncalib);
    m_outputTree->Branch("l_y_uncalib", &m_ly_uncalib);
    m_outputTree->Branch("g_x_uncalib", &m_x_uncalib);
    m_outputTree->Branch("g_y_uncalib", &m_y_uncalib);
    m_outputTree->Branch("g_z_uncalib", &m_z_uncalib);

    m_outputTree->Branch("predicted", &m_prt);
    m_outputTree->Branch("l_x_prt", &m_lx_prt);
    m_outputTree->Branch("l_y_prt", &m_ly_prt);
    m_outputTree->Branch("resid_x_prt", &m_resid_x_prt);
    m_outputTree->Branch("resid_y_prt", &m_resid_y_prt);
    m_outputTree->Branch("g_x_prt", &m_x_prt);
    m_outputTree->Branch("g_y_prt", &m_y_prt);
    m_outputTree->Branch("g_z_prt", &m_z_prt);
    m_outputTree->Branch("px_prt", &m_px_prt);
    m_outputTree->Branch("py_prt", &m_py_prt);
    m_outputTree->Branch("pz_prt", &m_pz_prt);
    m_outputTree->Branch("theta_prt", &m_theta_prt);
    m_outputTree->Branch("eta_prt", &m_eta_prt);
    m_outputTree->Branch("phi_prt", &m_phi_prt);
    m_outputTree->Branch("pT_prt", &m_pT_prt);

    m_outputTree->Branch("smoothed", &m_smt);
    m_outputTree->Branch("l_x_smt", &m_lx_smt);
    m_outputTree->Branch("l_y_smt", &m_ly_smt);
    m_outputTree->Branch("resid_x_smt", &m_resid_x_smt);
    m_outputTree->Branch("resid_y_smt", &m_resid_y_smt);
    m_outputTree->Branch("g_x_smt", &m_x_smt);
    m_outputTree->Branch("g_y_smt", &m_y_smt);
    m_outputTree->Branch("g_z_smt", &m_z_smt);
    m_outputTree->Branch("px_smt", &m_px_smt);
    m_outputTree->Branch("py_smt", &m_py_smt);
    m_outputTree->Branch("pz_smt", &m_pz_smt);
    m_outputTree->Branch("theta_smt", &m_theta_smt);
    m_outputTree->Branch("eta_smt", &m_eta_smt);
    m_outputTree->Branch("phi_smt", &m_phi_smt);
    m_outputTree->Branch("pT_smt", &m_pT_smt);
 
    m_outputTree->Branch("filtered", &m_flt);
    m_outputTree->Branch("l_x_flt", &m_lx_flt);
    m_outputTree->Branch("l_y_flt", &m_ly_flt);
    m_outputTree->Branch("resid_x_flt", &m_resid_x_flt);
    m_outputTree->Branch("resid_y_flt", &m_resid_y_flt);
    m_outputTree->Branch("g_x_flt", &m_x_flt);
    m_outputTree->Branch("g_y_flt", &m_y_flt);
    m_outputTree->Branch("g_z_flt", &m_z_flt);
    m_outputTree->Branch("px_flt", &m_px_flt);
    m_outputTree->Branch("py_flt", &m_py_flt);
    m_outputTree->Branch("pz_flt", &m_pz_flt);
    m_outputTree->Branch("theta_flt", &m_theta_flt);
    m_outputTree->Branch("eta_flt", &m_eta_flt);
    m_outputTree->Branch("phi_flt", &m_phi_flt);
    m_outputTree->Branch("pT_flt", &m_pT_flt);
 
    m_outputTree->Branch("charge", &m_charge);
    m_outputTree->Branch("barcode", &m_barcode, "barcode/l");
  }
}

FW::Root::RootTrackWriter::~RootTrackWriter()
{
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode
FW::Root::RootTrackWriter::endRun()
{
  if (m_outputFile) {
    m_outputFile->cd();
    m_outputTree->Write();
    ACTS_INFO("Wrote particles to tree '" << m_cfg.treeName << "' in '"
              << m_cfg.filePath
              << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode
FW::Root::RootTrackWriter::writeT(
                                  const AlgorithmContext&               ctx,
                                  const TrackMap&                       tracks)
{

  if (m_outputFile == nullptr) return ProcessCode::SUCCESS;

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Get the event number
  m_eventNr = ctx.eventNumber;

  // loop over the process vertices
  for (auto& track : tracks) {
    /// collect the information
    m_barcode = track.first;
    m_nStates = 0;
    for(auto& state: track.second){
      auto geoID = state.referenceSurface().geoID();
      m_volumeID.push_back(geoID.value(Acts::GeometryID::volume_mask)); 
      m_layerID.push_back(geoID.value(Acts::GeometryID::layer_mask)); 
      m_moduleID.push_back(geoID.value(Acts::GeometryID::sensitive_mask)); 
      // get the uncalibrated measurment
      auto meas = boost::get<Measurement>(*state.measurement.uncalibrated);
      // get local position 
      Acts::Vector2D local(meas.parameters()[Acts::ParDef::eLOC_0],
                           meas.parameters()[Acts::ParDef::eLOC_1]);
      // get global position 
      Acts::Vector3D global(0, 0, 0);
      Acts::Vector3D mom(1, 1, 1);
      // transform local into global position information
      meas.referenceSurface().localToGlobal(local, mom, global);
      
      // push the measurment 
      m_lx_uncalib.push_back(local.x());
      m_ly_uncalib.push_back(local.y());
      m_x_uncalib.push_back(global.x());
      m_y_uncalib.push_back(global.y());
      m_z_uncalib.push_back(global.z());

      // get the predicted parameter 
      bool predicted = false;
      if(state.parameter.predicted){
        predicted = true;
        auto parameter = *state.parameter.predicted;
        // push the predicted parameter 
        m_lx_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]),
          m_ly_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_resid_x_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] - local.x());
        m_resid_y_prt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] - local.y());
        m_x_prt.push_back(parameter.position().x());
        m_y_prt.push_back(parameter.position().y());
        m_z_prt.push_back(parameter.position().z());
        m_px_prt.push_back(parameter.momentum().x());
        m_py_prt.push_back(parameter.momentum().y());
        m_pz_prt.push_back(parameter.momentum().z());
        m_theta_prt.push_back(theta(parameter.momentum()));
        m_phi_prt.push_back(phi(parameter.momentum()));
        m_eta_prt.push_back(parameter.eta());
        m_pT_prt.push_back(parameter.pT());
        m_charge.push_back(parameter.charge());
      }

      // get the filtered parameter 
      bool filtered = false;
      if(state.parameter.filtered){
        filtered = true;
        auto parameter = *state.parameter.filtered;
        // push the filtered parameter 
        m_lx_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]),
          m_ly_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_resid_x_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] - local.x());
        m_resid_y_flt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] - local.y());
        m_x_flt.push_back(parameter.position().x());
        m_y_flt.push_back(parameter.position().y());
        m_z_flt.push_back(parameter.position().z());
        m_px_flt.push_back(parameter.momentum().x());
        m_py_flt.push_back(parameter.momentum().y());
        m_pz_flt.push_back(parameter.momentum().z());
        m_theta_flt.push_back(theta(parameter.momentum()));
        m_phi_flt.push_back(phi(parameter.momentum()));
        m_eta_flt.push_back(parameter.eta());
        m_pT_flt.push_back(parameter.pT());
      }

      // get the smoothed parameter 
      bool smoothed = false;
      if(state.parameter.smoothed){
        smoothed = true;
        auto parameter = *state.parameter.smoothed;
        // push the smoothed parameter 
        m_lx_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0]),
          m_ly_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1]);
        m_resid_x_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_0] - local.x());
        m_resid_y_smt.push_back(parameter.parameters()[Acts::ParDef::eLOC_1] - local.y());
        m_x_smt.push_back(parameter.position().x());
        m_y_smt.push_back(parameter.position().y());
        m_z_smt.push_back(parameter.position().z());
        m_px_smt.push_back(parameter.momentum().x());
        m_py_smt.push_back(parameter.momentum().y());
        m_pz_smt.push_back(parameter.momentum().z());
        m_theta_smt.push_back(theta(parameter.momentum()));
        m_phi_smt.push_back(phi(parameter.momentum()));
        m_eta_smt.push_back(parameter.eta());
        m_pT_smt.push_back(parameter.pT());
      }
      m_prt.push_back(predicted);
      m_flt.push_back(filtered);
      m_smt.push_back(smoothed);
      m_nStates ++;
    }
    m_outputTree->Fill();
    // now reset
    m_volumeID.clear();
    m_layerID.clear();
    m_moduleID.clear();
    m_lx_uncalib.clear();
    m_ly_uncalib.clear();
    m_x_uncalib.clear();
    m_y_uncalib.clear();
    m_z_uncalib.clear();

    m_prt.clear();
    m_lx_prt.clear();
    m_ly_prt.clear();
    m_resid_x_prt.clear();
    m_resid_y_prt.clear();
    m_x_prt.clear();
    m_y_prt.clear();
    m_z_prt.clear();
    m_px_prt.clear();
    m_py_prt.clear();
    m_pz_prt.clear();
    m_theta_prt.clear();
    m_phi_prt.clear();
    m_eta_prt.clear();
    m_pT_prt.clear();
    m_charge.clear();

    m_flt.clear();
    m_lx_flt.clear();
    m_ly_flt.clear();
    m_resid_x_flt.clear();
    m_resid_y_flt.clear();
    m_x_flt.clear();
    m_y_flt.clear();
    m_z_flt.clear();
    m_px_flt.clear();
    m_py_flt.clear();
    m_pz_flt.clear();
    m_theta_flt.clear();
    m_phi_flt.clear();
    m_eta_flt.clear();
    m_pT_flt.clear();
    
    m_smt.clear();
    m_lx_smt.clear();
    m_ly_smt.clear();
    m_resid_x_smt.clear();
    m_resid_y_smt.clear();
    m_x_smt.clear();
    m_y_smt.clear();
    m_z_smt.clear();
    m_px_smt.clear();
    m_py_smt.clear();
    m_pz_smt.clear();
    m_theta_smt.clear();
    m_phi_smt.clear();
    m_eta_smt.clear();
    m_pT_smt.clear();
  }

  return ProcessCode::SUCCESS;
}
