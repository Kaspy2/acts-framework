// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "VertexFindingAlgorithm.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <iostream>

FWE::VertexFindingAlgorithm::VertexFindingAlgorithm(const Config&        cfg,
                                                    Acts::Logging::Level level)
  : FW::BareAlgorithm("VertexFinding", level), m_cfg(cfg)
{
}

/// @brief Algorithm that receives all selected tracks from an event
/// and finds and fits its vertices
FW::ProcessCode
FWE::VertexFindingAlgorithm::execute(const FW::AlgorithmContext& context) const
{

  const auto& inputTrackCollection
      = context.eventStore.get<std::vector<Acts::BoundParameters>>(
          m_cfg.trackCollection);

  // The geometry context
  Acts::GeometryContext geoContext;
  // The magnetic Field context
  Acts::MagneticFieldContext magFieldContext;

  // Vertex finder options
  Acts::VertexFinderOptions<Acts::BoundParameters> vFinderOptions(
      geoContext, magFieldContext);

  // Find vertices
  auto res = m_cfg.vertexFinder->find(inputTrackCollection, vFinderOptions);

  if (res.ok()) {
    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    ACTS_INFO("Found " << vertexCollection.size() << " vertices in event.");

    unsigned int count = 0;
    for (const auto& vtx : vertexCollection) {
      ACTS_INFO(++count << ". vertex at "
                        << "("
                        << vtx.position()[0]
                        << ","
                        << vtx.position()[1]
                        << ","
                        << vtx.position()[2]
                        << ") with "
                        << vtx.tracks().size()
                        << " tracks.");
    }
  } else {
    ACTS_ERROR("Error in vertex finder.");
  }

  return FW::ProcessCode::SUCCESS;
}
