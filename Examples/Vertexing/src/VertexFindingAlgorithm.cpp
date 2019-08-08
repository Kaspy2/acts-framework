// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "VertexFindingAlgorithm.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"

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

  using Propagator = Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>>;
  using VertexFitter = Acts::FullBilloirVertexFitter<Acts::ConstantBField,
                                                     Acts::BoundParameters,
                                                     Propagator>;
  using VertexFinder = Acts::IterativeVertexFinder<Acts::ConstantBField,
                                                   Acts::BoundParameters,
                                                   Propagator,
                                                   VertexFitter>;

  const auto& input = context.eventStore.get<std::vector<FW::VertexAndTracks>>(
      m_cfg.trackCollection);

  // Set up constant B-Field
  Acts::ConstantBField bField(m_cfg.bField);

  // Set up Eigenstepper
  Acts::EigenStepper<Acts::ConstantBField> stepper(bField);
  // Set up propagator with void navigator
  Propagator propagator(stepper);

  // Set up Billoir Vertex Fitter
  VertexFitter::Config vertexFitterCfg(bField, propagator);
  VertexFitter         vertexFitter(vertexFitterCfg);

  // Set up Iterative Vertex Finder
  VertexFinder::Config finderCfg(bField, std::move(vertexFitter), propagator);
  VertexFinder         vertexFinder(finderCfg);

  // The geometry context
  Acts::GeometryContext geoContext;
  // The magnetic Field context
  Acts::MagneticFieldContext magFieldContext;

  // Vertex finder options
  Acts::VertexFinderOptions<Acts::BoundParameters> vFinderOptions(
      geoContext, magFieldContext);

  std::vector<Acts::BoundParameters> inputTrackCollection;

  int counte = 0;
  for (auto& bla : input) {
    counte += bla.tracks.size();
  }

  ACTS_INFO("Truth vertices in event: " << input.size());

  for (auto& vertexAndTracks : input) {
    ACTS_INFO("\t True vertex at (" << vertexAndTracks.vertex.position[0] << ","
                                    << vertexAndTracks.vertex.position[1]
                                    << ","
                                    << vertexAndTracks.vertex.position[2]
                                    << ") with "
                                    << vertexAndTracks.tracks.size()
                                    << " tracks.");
    inputTrackCollection.insert(inputTrackCollection.end(),
                                vertexAndTracks.tracks.begin(),
                                vertexAndTracks.tracks.end());
  }

  // Find vertices
  auto res = vertexFinder.find(inputTrackCollection, vFinderOptions);

  if (res.ok()) {
    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    ACTS_INFO("Found " << vertexCollection.size() << " vertices in event.");

    unsigned int count = 0;
    for (const auto& vtx : vertexCollection) {
      ACTS_INFO("\t" << ++count << ". vertex at "
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
