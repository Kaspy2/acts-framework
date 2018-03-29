// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Csv/CsvParticleWriter.hpp"
#include <fstream>
#include <ios>
#include <stdexcept>
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "ACTFW/Barcode/BarcodeSvc.hpp"

FW::Csv::CsvParticleWriter::CsvParticleWriter(
    const FW::Csv::CsvParticleWriter::Config& cfg,
    Acts::Logging::Level                      level)
  : Base(cfg.collection, "CsvParticleWriter", level), m_cfg(cfg)
{
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (!m_cfg.barcodeSvc) {
    throw std::invalid_argument("Missing barcode service");
  }
}

FW::ProcessCode
FW::Csv::CsvParticleWriter::writeT(
    const FW::AlgorithmContext&             ctx,
    const std::vector<Acts::ProcessVertex>& vertices)
{
  std::string pathOs
      = perEventFilepath(m_cfg.outputDir, "particles.csv", ctx.eventNumber);
  std::ofstream os(pathOs, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + pathOs + "' to write");
  }

  const std::map<barcode_type,size_t>* hitsPerParticle = nullptr;
  if (m_cfg.hitsPerParticleCollection != "" &&
      ctx.eventStore.get(m_cfg.hitsPerParticleCollection, hitsPerParticle) 
        == ProcessCode::ABORT )
  {
    throw std::ios_base::failure("Could not retrieve hits/particle reference map.");
  }
  
  std::string pathOsn
      = perEventFilepath(m_cfg.outputDir, "sterile-particles.csv", ctx.eventNumber);
  std::ofstream osn(pathOsn, std::ofstream::out | std::ofstream::trunc);
  if (!osn) {
    throw std::ios_base::failure("Could not open '" + pathOsn + "' to write");
  }
  
  std::string pathHpp
      = perEventFilepath(m_cfg.outputDir, "hits-per-particles.csv", ctx.eventNumber);
  std::ofstream hpp(pathHpp, std::ofstream::out | std::ofstream::trunc);
  if (!hpp) {
    throw std::ios_base::failure("Could not open '" + pathHpp + "' to write");
  }
  
  // write csv header
  os << "particle_id,";
  os << "vx,vy,vz,";
  os << "px,py,pz,";
  os << "q\n";

  // write csv header
  osn << "particle_id,";
  osn << "vx,vy,vz,";
  osn << "px,py,pz,";
  osn << "q\n";
  
  // write csv header
  hpp << "particle_id,nhits\n";
  
  // write one line per particle
  os << std::setprecision(m_cfg.outputPrecision);
  for (auto& vertex : vertices) {
    auto& vtx = vertex.position();
    for (auto& particle : vertex.outgoingParticles()) {
      if (hitsPerParticle->find(particle.barcode()) != hitsPerParticle->end()){
        os << particle.barcode() << ",";
        os << vtx.x() << ",";
        os << vtx.y() << ",";
        os << vtx.z() << ",";
        os << particle.momentum().x() << ",";
        os << particle.momentum().y() << ",";
        os << particle.momentum().z() << ",";
        os << particle.charge() << '\n';
      } else {
        osn << particle.barcode() << ",";
        osn << vtx.x() << ",";
        osn << vtx.y() << ",";
        osn << vtx.z() << ",";
        osn << particle.momentum().x() << ",";
        osn << particle.momentum().y() << ",";
        osn << particle.momentum().z() << ",";
        osn << particle.charge() << '\n';
      }
    }
  }
  
  // write the hits per particle 
  for (auto& hppi : (*hitsPerParticle)){
    if (hppi.second)
      hpp << hppi.first << "," << hppi.second << "\n";
  }
  
  return ProcessCode::SUCCESS;
}
