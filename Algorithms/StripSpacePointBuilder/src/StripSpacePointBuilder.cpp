// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/StripSpacePointBuilder/StripSpacePointBuilder.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include "ACTFW/Framework/WhiteBoard.hpp"

///
/// @note Used abbreviation: "Strip Detector Element" -> SDE
///

FW::StripSpacePointBuilder::StripSpacePointBuilder(const FW::StripSpacePointBuilder::Config& cfg,
                                       Acts::Logging::Level level)
  : FW::BareAlgorithm("StripSpacePointBuilder", level), m_cfg(cfg)
{
  // Check that all mandatory configuration parameters are present
  if (m_cfg.collectionIn.empty()) {
    throw std::invalid_argument("Missing input collection");
  }

  if (m_cfg.collectionOut.empty()) {
    throw std::invalid_argument("Missing output collection");
  }

  // Print chosen configuration
  ACTS_DEBUG("Space Point Finder settings: ");
  ACTS_DEBUG("- intput collection: " << m_cfg.collectionIn);
  ACTS_DEBUG("- output collection: " << m_cfg.collectionOut);
  ACTS_DEBUG("- difference in squared theta: " << m_cfg.diffTheta2);
  ACTS_DEBUG("- difference in squared phi: " << m_cfg.diffPhi2);
  ACTS_DEBUG("- difference in distances: " << m_cfg.diffDist);
  ACTS_DEBUG("- vertex position: " << m_cfg.vertex);
}

FW::ProcessCode
FW::StripSpacePointBuilder::clusterReading(AlgorithmContext& ctx,
                                     const DetData*&   detData) const
{
  // Load hit data from Whiteboard
  if (ctx.eventStore.get(m_cfg.collectionIn, detData)
      != FW::ProcessCode::SUCCESS) {
    ACTS_DEBUG("Unable to receive event data");
    return FW::ProcessCode::ABORT;
  }

  ACTS_DEBUG("Event data successfully received");
  return FW::ProcessCode::SUCCESS;
}

Acts::Vector2D
FW::StripSpacePointBuilder::localCoords(const Acts::PlanarModuleCluster& hit) const
{
  // Local position information
  auto           par = hit.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D
FW::StripSpacePointBuilder::globalCoords(const Acts::PlanarModuleCluster& hit) const
{
  // Receive corresponding surface
  auto& clusterSurface = hit.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(localCoords(hit), mom, pos);

  return pos;
}

double
FW::StripSpacePointBuilder::differenceOfHits(
    const Acts::PlanarModuleCluster& hit1,
    const Acts::PlanarModuleCluster& hit2) const
{
  // Calculate the global position of the hits
  Acts::Vector3D pos1 = globalCoords(hit1);
  Acts::Vector3D pos2 = globalCoords(hit2);

  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > m_cfg.diffDist) return -1.;

  // Calculate the angles of the hits
  double phi1, theta1, phi2, theta2;
  phi1   = (pos1 - m_cfg.vertex).phi();
  theta1 = (pos1 - m_cfg.vertex).theta();
  phi2   = (pos2 - m_cfg.vertex).phi();
  theta2 = (pos2 - m_cfg.vertex).theta();

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > m_cfg.diffTheta2) {
    ACTS_DEBUG(
        "Squared theta angle " << diffTheta2 << " between positions ("
                               << pos1.x()
                               << ", "
                               << pos1.y()
                               << ", "
                               << pos1.z()
                               << ") and ("
                               << pos2.x()
                               << ", "
                               << pos2.y()
                               << ", "
                               << pos2.z()
                               << ") too large - points are not combined");
    return -1.;
  }

  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > m_cfg.diffPhi2) {
    ACTS_DEBUG("Squared phi angle " << diffPhi2 << " between positions ("
                                    << pos1.x()
                                    << ", "
                                    << pos1.y()
                                    << ", "
                                    << pos1.z()
                                    << ") and ("
                                    << pos2.x()
                                    << ", "
                                    << pos2.y()
                                    << ", "
                                    << pos2.z()
                                    << ") too large - points are not combined");
    return -1.;
  }

  // Return the squared distance between both hits
  return diffTheta2 + diffPhi2;
}

void
FW::StripSpacePointBuilder::combineHits(
    const std::vector<Acts::PlanarModuleCluster>&    vec1,
    const std::vector<Acts::PlanarModuleCluster>&    vec2,
    std::vector<FW::StripSpacePointBuilder::CombinedHits>& combHits) const
{
  // TODO: only the closest differences get selected -> some points are not
  // taken into account
  // Declare helper variables
  double                             currentDiff;
  FW::StripSpacePointBuilder::CombinedHits tmpCombHits;
  double                             diffMin;
  unsigned int                       hitMin;

  // Walk through all hits on both surfaces
  for (unsigned int iVec1 = 0; iVec1 < vec1.size(); iVec1++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of hits
    hitMin = vec2.size();
    for (unsigned int iVec2 = 0; iVec2 < vec2.size(); iVec2++) {
      // Calculate the distances between the hits
      currentDiff = differenceOfHits(vec1[iVec1], vec2[iVec2]);
      // Store the closest hits (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        hitMin  = iVec2;
      }
    }
    // Store the best (=closest) result
    if (hitMin < vec2.size()) {
      tmpCombHits.hitModule1 = &(vec1[iVec1]);
      tmpCombHits.hitModule2 = &(vec2[hitMin]);
      tmpCombHits.diff       = diffMin;
      combHits.push_back(tmpCombHits);
    }
  }
}

void
FW::StripSpacePointBuilder::findOverlappingClusters(
    const DetData*&                                               detData,
    std::vector<std::vector<FW::StripSpacePointBuilder::CombinedHits>>& allCombHits)
    const
{
  // Declare temporary storage
  std::vector<std::vector<Acts::PlanarModuleCluster> const*> modules;
  unsigned int                                               index;
  std::vector<FW::StripSpacePointBuilder::CombinedHits>            combHits;

  // Loop over the planar clusters in this event
  for (auto& volumeData : *detData)
    for (auto& layerData : volumeData.second) {
      modules.clear();
      combHits.clear();

      for (auto& moduleData : layerData.second) {
        index = moduleData.first - 1;
        if (moduleData.first + 1 >= modules.size()) modules.resize(index + 1);
        // Store the data
        modules[index] = &(moduleData.second);
      }
      switch (modules.size()) {
      // Easy exit if no surfaces are available
      case 0:
        ACTS_VERBOSE("Layer " << layerData.first << " stored no surfaces");
        break;
      // Easy exit if only a single surface is available
      case 1:
        ACTS_VERBOSE("Layer " << layerData.first << " stored a single surface");
        break;
      // Form combinations between all hits on both surfaces
      case 2:
        combineHits(*(modules[0]), *(modules[1]), combHits);
        break;
      // Try to combine the hits of every surface with each other
      default:
        std::vector<FW::StripSpacePointBuilder::CombinedHits> tmpCombHits;
        for (unsigned int i = 0; i < modules.size(); i++)
          for (unsigned int j = i + 1; j < modules.size(); j++) {
            combineHits(*(modules[i]), *(modules[j]), tmpCombHits);
            combHits.insert(
                combHits.end(), tmpCombHits.begin(), tmpCombHits.end());
          }
      }
      // Store the found combination candidates
      allCombHits.push_back(std::move(combHits));
    }
}

void
FW::StripSpacePointBuilder::filterCombinations(
    std::vector<std::vector<FW::StripSpacePointBuilder::CombinedHits>>& allCombHits)
    const
{
  // Walk through the layers
  for (const auto& layer : allCombHits)
    // Walk through every hit combination of a layer
    for (const auto& combination : layer) {
      for (const auto& combinationCompare : layer)
        // Check, if combinations infere and resolve the problem.
        // Hits on the first module are different from each other by
        // construction but a hit on the second module can appear in multiple
        // combinations.
        if ((combination.hitModule1 != combinationCompare.hitModule1)
            && (combination.hitModule2 == combinationCompare.hitModule2))
          ACTS_INFO(
              "Warning: Multiple candidate selection not implemented yet");
    }
}

std::pair<Acts::Vector3D, Acts::Vector3D>
FW::StripSpacePointBuilder::endsOfStrip(const Acts::PlanarModuleCluster& hit) const
{
  // Calculate the local coordinates of the hit
  const Acts::Vector2D local = localCoords(hit);

  // Receive the binning
  auto& sur        = hit.referenceSurface();
  auto  genDetElem = dynamic_cast<const FWGen::GenericDetectorElement*>(
      sur.associatedDetectorElement());
  auto segment = dynamic_cast<const Acts::CartesianSegmentation*>(
      &(genDetElem->digitizationModule()->segmentation()));
  auto& binData     = segment->binUtility().binningData();
  auto& boundariesX = binData[0].boundaries();
  auto& boundariesY = binData[1].boundaries();

  // Search the x-/y-bin hit
  size_t binX = binData[0].searchLocal(local);
  size_t binY = binData[1].searchLocal(local);

  Acts::Vector2D topLocal, bottomLocal;

  if (boundariesX[binX + 1] - boundariesX[binX]
      < boundariesY[binY + 1] - boundariesY[binY]) {
    // Set the top and bottom end of the strip in local coordinates
    topLocal = {(boundariesX[binX] + boundariesX[binX + 1]) / 2,
                boundariesY[binY + 1]};
    bottomLocal
        = {(boundariesX[binX] + boundariesX[binX + 1]) / 2, boundariesY[binY]};
  } else {
    // Set the top and bottom end of the strip in local coordinates
    topLocal
        = {boundariesX[binX], (boundariesY[binY] + boundariesY[binY + 1]) / 2};
    bottomLocal = {boundariesX[binX + 1],
                   (boundariesY[binY] + boundariesY[binY + 1]) / 2};
  }

  // Calculate the global coordinates of the top and bottom end of the strip
  Acts::Vector3D topGlobal, bottomGlobal, mom;  // mom is a dummy variable
  sur.localToGlobal(topLocal, mom, topGlobal);
  sur.localToGlobal(bottomLocal, mom, bottomGlobal);

  // Return the top and bottom end of the strip in global coordinates
  return std::make_pair(topGlobal, bottomGlobal);
}

bool
FW::StripSpacePointBuilder::recoverSpacePoint(
    FW::StripSpacePointBuilder::SpacePointParameters& spaPoPa) const
{
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (m_cfg.stripLengthGapTolerance <= 0.) return false;
  spaPoPa.qmag = spaPoPa.q.mag();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended
      = spaPoPa.limit + m_cfg.stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) return false;
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.)
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) return false;

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory. The
  /// would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale
      = spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot
        = (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot
        = (mOvershoot > nOvershoot) ? mOvershoot : nOvershoot;
    // Move m and n towards 0
    spaPoPa.m -= biggerOvershoot;
    spaPoPa.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // Check if both overshoots are in the same direction
  if (spaPoPa.m < -1. && spaPoPa.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spaPoPa.m + 1.);
    double nOvershoot
        = -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot
        = (mOvershoot > nOvershoot) ? mOvershoot : nOvershoot;
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

void
FW::StripSpacePointBuilder::storeSpacePoint(const Acts::Vector3D& spacePoint,
                                      const CombinedHits&   hits,
                                      DetData&              stripClusters) const
{
  // Receive the identification of the digitized hits on the first surface
  Identifier       id(hits.hitModule1->identifier());
  Acts::GeometryID geoID(id.value());

  // The covariance is currently set to 0.
  Acts::ActsSymMatrixD<2> cov;
  cov << 0., 0., 0., 0.;

  // Get the local coordinates of the space point
  Acts::Vector2D local;
  hits.hitModule1->referenceSurface().globalToLocal(
      spacePoint, {0., 0., 0.}, local);

  // Build the space point
  Acts::PlanarModuleCluster pCluster(
      hits.hitModule1->referenceSurface(),
      Identifier(geoID.value()),
      std::move(cov),
      local[0],
      local[1],
      std::move(hits.hitModule1->digitizationCells()),
      {hits.hitModule1->truthVertices()});

  // Insert into the cluster map
  FW::Data::insert(stripClusters,
                   geoID.value(Acts::GeometryID::volume_mask),
                   geoID.value(Acts::GeometryID::layer_mask),
                   geoID.value(Acts::GeometryID::sensitive_mask),
                   std::move(pCluster));
}

FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>
FW::StripSpacePointBuilder::calculateSpacePoints(
    std::vector<std::vector<CombinedHits>>& allCombHits,
    const DetData*&                         detData) const
{
  // Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  DetData                                    stripClusters;
  FW::StripSpacePointBuilder::SpacePointParameters spaPoPa;

  // Walk over every found candidate pair
  for (const auto& layers : allCombHits)
    for (const auto& hits : layers) {
      // Calculate the ends of the SDEs
      const auto& ends1 = endsOfStrip(*(hits.hitModule1));
      const auto& ends2 = endsOfStrip(*(hits.hitModule2));

      /// The following algorithm is meant for finding the position on the first
      /// strip if there is a corresponding hit on the second strip. The
      /// resulting point is a point x on the first surfaces. This point is
      /// along a line between the points a (top end of the strip)
      /// and b (bottom end of the strip). The location can be parametrized as
      /// 	2 * x = (1 + m) a + (1 - m) b
      /// as function of the scalar m. m is a parameter in the interval
      /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
      /// from the vertex to the hit on the second strip y is needed to be a
      /// multiple k of the vector from vertex to the hit on the first strip x.
      /// As a consequence of this demand y = k * x needs to be on the
      /// connecting line between the top (c) and bottom (d) end of
      /// the second strip. If both hits correspond to each other, the condition
      /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
      /// needs to be fulfilled. Inserting the first equation into this
      /// equation leads to the condition for m as given in the following
      /// algorithm and therefore to the calculation of x.
      /// The same calculation can be repeated for y. Its corresponding
      /// parameter will be named n.

      spaPoPa.reset();
      spaPoPa.q  = ends1.first - ends1.second;
      spaPoPa.s  = ends1.first + ends1.second - 2 * m_cfg.vertex;
      spaPoPa.r  = ends2.first - ends2.second;
      spaPoPa.t  = ends2.first + ends2.second - 2 * m_cfg.vertex;
      spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
      spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
      spaPoPa.m  = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

      // Set the limit for the parameter
      if (spaPoPa.limit == 1. && m_cfg.stripLengthTolerance != 0.)
        spaPoPa.limit = 1. + m_cfg.stripLengthTolerance;

      // Check if m and n can be resolved in the interval (-1, 1)
      if (fabs(spaPoPa.m) <= spaPoPa.limit
          && fabs(spaPoPa.n
                  = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs))
              <= spaPoPa.limit)
        // Store the space point
        storeSpacePoint(
            0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q),
            hits,
            stripClusters);
      else
          /// If this point is reached then it was not possible to resolve both
          /// points such that they are on their SDEs
          /// The following code treats a possible recovery of points resolved
          /// slightly outside of the SDE.
          /// @note This procedure is an indirect variation of the vertex
          /// position.
          // Check if a recovery the point(s) and store them if successful
          if (recoverSpacePoint(spaPoPa))
        storeSpacePoint(
            0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q),
            hits,
            stripClusters);
    }
  // Return the resolved hits
  return stripClusters;
}

FW::ProcessCode
FW::StripSpacePointBuilder::execute(AlgorithmContext ctx) const
{
  ACTS_DEBUG("::execute() called for event " << ctx.eventNumber);

  // DetData is typename of DetectorData<geo_id_value,
  // Acts::PlanarModuleCluster>
  const DetData*                                               detData;
  std::vector<std::vector<FW::StripSpacePointBuilder::CombinedHits>> allCombHits;
  DetData                                                      stripClusters;

  // Receive all hits from the Whiteboard
  clusterReading(ctx, detData);

  // Extract hits that occur in opposite SCTs
  findOverlappingClusters(detData, allCombHits);

  // Filter entries
  //~ filterCombinations(allCombHits);

  // Calculate the SpacePoints measured by the combination
  stripClusters = calculateSpacePoints(allCombHits, detData);

  // Write to Whiteboard
  if (ctx.eventStore.add(m_cfg.collectionOut, std::move(stripClusters))
      != ProcessCode::SUCCESS)
    return ProcessCode::ABORT;
  return FW::ProcessCode::SUCCESS;
}
