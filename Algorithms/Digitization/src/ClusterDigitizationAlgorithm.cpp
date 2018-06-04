// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Digitization/ClusterDigitizationAlgorithm.hpp"

#include <iostream>
#include <stdexcept>

#include "ACTFW/Barcode/BarcodeSvc.hpp"
#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"
#include "ACTS/Detector/DetectorElementBase.hpp"
#include "ACTS/Digitization/PlanarModuleStepper.hpp"
#include "ACTS/Digitization/Segmentation.hpp"
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

std::vector<FW::ClusterDigitizationAlgorithm::SingleParticleCluster>
FW::ClusterDigitizationAlgorithm::divideSingleParticleClusters(
    const std::vector<SingleParticleCluster>& clusters) const
{
  std::vector<SingleParticleCluster> resolvedClusters;

  bool         commonEdge;
  unsigned int clusterSize = clusters.size();

  for (unsigned int i = 0; i < clusterSize; i++)
    for (unsigned int j = 0; j < clusters[i].usedCells.size(); j++) {
      commonEdge = false;
      for (unsigned int k = 0; k < clusters[i].usedCells.size(); k++) {
        if (j == k) continue;
        if ((clusters[i].usedCells[k].channel0
                 == clusters[i].usedCells[j].channel0 - 1
             || clusters[i].usedCells[k].channel0
                 == clusters[i].usedCells[j].channel0
             || clusters[i].usedCells[k].channel0
                 == clusters[i].usedCells[j].channel0 + 1)
            && (clusters[i].usedCells[k].channel1
                    >= clusters[i].usedCells[j].channel1 - 1
                || clusters[i].usedCells[k].channel1
                    >= clusters[i].usedCells[j].channel1
                || clusters[i].usedCells[k].channel1
                    >= clusters[i].usedCells[j].channel1 + 1)) {
          commonEdge = true;
          break;
        }
      }
      if (!commonEdge) {
        SingleParticleCluster tmpSc = clusters[i];
        tmpSc.usedCells.clear();
        tmpSc.usedCells.push_back(clusters[i].usedCells[j]);
        resolvedClusters.push_back(tmpSc);
        resolvedClusters.push_back(clusters[i]);
        resolvedClusters.back().usedCells.erase(
            resolvedClusters.back().usedCells.begin() + j);
      }
    }
  return std::move(resolvedClusters);
}

bool
FW::ClusterDigitizationAlgorithm::commonEdge(
    const SingleParticleCluster& cluster1,
    const SingleParticleCluster& cluster2) const
{
  for (auto& uC1 : cluster1.usedCells)
    for (auto& uC2 : cluster2.usedCells)
      if ((uC1.channel0 == uC2.channel0 - 1 || uC1.channel0 == uC2.channel0
           || uC1.channel0 == uC2.channel0 + 1)
          && (uC1.channel1 == uC2.channel1 - 1 || uC1.channel1 == uC2.channel1
              || uC1.channel1 == uC2.channel1 + 1))
        return true;
  return false;
}

std::
    vector<std::vector<FW::ClusterDigitizationAlgorithm::SingleParticleCluster>>
    FW::ClusterDigitizationAlgorithm::mergeSingleParticleClusters(
        const std::vector<SingleParticleCluster>& clusters) const
{
  std::vector<SingleParticleCluster> resolvedClusters
      = divideSingleParticleClusters(clusters);

  std::vector<std::vector<SingleParticleCluster>> mergedClusters;
  std::vector<unsigned int>                       mergedIndices;

  for (unsigned int i = 0; i < resolvedClusters.size(); i++) {
    mergedClusters.push_back({{resolvedClusters[i]}});
    for (unsigned int j = i + 1; j < resolvedClusters.size(); j++)
      if (commonEdge(resolvedClusters[i], resolvedClusters[j])) {
        mergedClusters.back().push_back(resolvedClusters[j]);
        mergedIndices.push_back(j);
      }
    std::vector<unsigned int>::reverse_iterator rit = mergedIndices.rbegin();
    for (; rit != mergedIndices.rend(); ++rit)
      resolvedClusters.erase(resolvedClusters.begin() + *rit);
    mergedIndices.clear();
  }

  return mergedClusters;
}

double
FW::ClusterDigitizationAlgorithm::meanLorentzShift(
    const std::vector<SingleParticleCluster>& clusters) const
{
  double meanLorentzShift = 0.;
  for (const SingleParticleCluster& spc : clusters)
    meanLorentzShift += spc.lorentzShift;

  return meanLorentzShift / clusters.size();
}

Acts::Vector2D
FW::ClusterDigitizationAlgorithm::localPosition(
    const std::vector<SingleParticleCluster>& clusters) const
{
  Acts::Vector2D locPos    = Acts::Vector2D::Zero(2);
  double         totalPath = 0;

  for (const SingleParticleCluster& spc : clusters) {
    locPos(0) += spc.localX;
    locPos(1) += spc.localY;
    totalPath += spc.totalPath;
  }
  locPos /= totalPath;

  return locPos;
}

Acts::ActsSymMatrixD<2>
FW::ClusterDigitizationAlgorithm::covariance(
    const std::vector<SingleParticleCluster>& clusters,
    const Acts::Vector2D                      mean) const
{
  Acts::ActsSymMatrixD<2> cov;

  cov << 0., 0., 0., 0.;

  double       totalPath = 0;
  unsigned int nSteps    = 0;

  for (const SingleParticleCluster& spc : clusters) {
    for (auto& dStep : spc.dSteps) {
      cov(0, 0) += (dStep.stepCellCenter.x() - mean(0))
          * (dStep.stepCellCenter.x() - mean(0)) * dStep.stepLength;
      cov(1, 0) += (dStep.stepCellCenter.x() - mean(0))
          * (dStep.stepCellCenter.y() - mean(1)) * dStep.stepLength;
      cov(1, 1) += (dStep.stepCellCenter.y() - mean(1))
          * (dStep.stepCellCenter.y() - mean(1)) * dStep.stepLength;
      nSteps++;
    }
    totalPath += spc.totalPath;
  }
  cov(0, 1) = cov(1, 0);

  cov /= totalPath;
  if (nSteps > 1)
    cov /= (double)nSteps;
  else  // TODO: cell size needs to be used
    cov << 0., 0., 0., 0.;

  return cov;
}

std::vector<Acts::PlanarModuleCluster>
FW::ClusterDigitizationAlgorithm::formClusters(
    const std::vector<SingleParticleCluster>& clusters) const
{
  std::vector<std::vector<SingleParticleCluster>> mergedClusters
      = mergeSingleParticleClusters(clusters);
  std::vector<Acts::PlanarModuleCluster> pClusters;

  for (auto& mCluster : mergedClusters) {
    Acts::Vector2D locPos = localPosition(mCluster);

    Acts::ActsSymMatrixD<2> cov = covariance(clusters, locPos);

    locPos(0) += meanLorentzShift(clusters);

    std::vector<Acts::DigitizationCell> usedCells;
    std::vector<Acts::ProcessVertex>    pVertex;

    for (SingleParticleCluster& spc : mCluster) {
      usedCells.insert(
          usedCells.end(), spc.usedCells.begin(), spc.usedCells.end());
      pVertex.push_back(*(spc.pVertex));
    }

    // create the planar cluster
    pClusters.push_back(
        Acts::PlanarModuleCluster(*(mCluster[0].hitSurface),
                                  Identifier(mCluster[0].geoID.value()),
                                  std::move(cov),
                                  locPos(0),
                                  locPos(1),
                                  std::move(usedCells),
                                  pVertex));
  }

  return std::move(pClusters);
}

void
FW::ClusterDigitizationAlgorithm::clusterize(
    FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>& planarClusters,
    const std::vector<SingleParticleCluster>& clusters,
    const geo_id_value                        volumeKey,
    const geo_id_value                        layerKey,
    const geo_id_value                        moduleKey) const
{
  if (clusters.empty()) return;

  std::vector<Acts::PlanarModuleCluster> pClusters = formClusters(clusters);

  for (Acts::PlanarModuleCluster pCluster : pClusters)
    // insert into the cluster map
    FW::Data::insert(
        planarClusters, volumeKey, layerKey, moduleKey, std::move(pCluster));
}

FW::ProcessCode
FW::ClusterDigitizationAlgorithm::execute(FW::AlgorithmContext ctx) const
{
  // prepare the input data
  const FW::DetectorData<geo_id_value,
                         std::pair<std::unique_ptr<const Acts::TrackParameters>,
                                   barcode_type>>* hitData
      = nullptr;
  // read and go
  if (ctx.eventStore.get(m_cfg.simulatedHitsCollection, hitData)
      == FW::ProcessCode::ABORT)
    return FW::ProcessCode::ABORT;

  ACTS_DEBUG("Retrieved hit data '" << m_cfg.simulatedHitsCollection
                                    << "' from event store.");

  // the particle mass table
  Acts::ParticleMasses pMasses;

  // prepare the output data: Clusters
  FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster> planarClusters;
  // perpare the second output data : SpacePoints
  FW::DetectorData<geo_id_value, Acts::Vector3D> spacePoints;

  // now digitise
  for (auto& vData : (*hitData)) {
    auto volumeKey = vData.first;
    ACTS_DEBUG("- Processing Volume Data collection for volume with ID "
               << volumeKey);
    for (auto& lData : vData.second) {
      auto layerKey = lData.first;
      ACTS_DEBUG("-- Processing Layer Data collection for layer with ID "
                 << layerKey);
      for (auto& sData : lData.second) {
        auto moduleKey = sData.first;
        ACTS_DEBUG("-- Processing Module Data collection for module with ID "
                   << moduleKey);

        std::vector<SingleParticleCluster> clusters;

        // get the hit parameters
        for (auto& hit : sData.second) {
          auto hitParameters   = hit.first.get();
          auto particleBarcode = hit.second;
          // get the surface
          const Acts::Surface& hitSurface = hitParameters->referenceSurface();
          // get the DetectorElement
          auto hitDetElement = hitSurface.associatedDetectorElement();
          if (hitDetElement) {
            // get the digitization module
            auto hitDigitizationModule = hitDetElement->digitizationModule();
            if (hitDigitizationModule) {
              struct SingleParticleCluster spc;
              spc.hitSurface = &hitSurface;

              // get the lorentz angle
              double lorentzAngle = hitDigitizationModule->lorentzAngle();
              double thickness    = hitDetElement->thickness();
              spc.lorentzShift    = thickness * tan(lorentzAngle);
              spc.lorentzShift *= -(hitDigitizationModule->readoutDirection());
              // parameters
              auto           pars     = hitParameters->parameters();
              auto           position = hitParameters->position();
              auto           momentum = hitParameters->momentum();
              Acts::Vector2D localIntersection(pars[Acts::ParDef::eLOC_0],
                                               pars[Acts::ParDef::eLOC_1]);
              Acts::Vector3D localDirection(
                  hitSurface.transform().inverse().linear() * momentum);
              // position
              spc.dSteps
                  = m_cfg.planarModuleStepper->cellSteps(*hitDigitizationModule,
                                                         localIntersection,
                                                         localDirection.unit());
              // everything under threshold or edge effects
              if (!spc.dSteps.size()) continue;
              /// let' create a cluster - centroid method
              // the cells to be used
              spc.usedCells.reserve(spc.dSteps.size());
              // loop over the steps
              for (auto dStep : spc.dSteps) {
                // @todo implement smearing
                spc.localX += dStep.stepLength * dStep.stepCellCenter.x();
                spc.localY += dStep.stepLength * dStep.stepCellCenter.y();
                spc.totalPath += dStep.stepLength;
                spc.usedCells.push_back(
                    std::move(Acts::DigitizationCell(dStep.stepCell.channel0,
                                                     dStep.stepCell.channel1,
                                                     dStep.stepLength)));
              }

              // get the segmentation & find the corresponding cell id
              const Acts::Segmentation& segmentation
                  = hitDigitizationModule->segmentation();
              auto           binUtility = segmentation.binUtility();
              Acts::Vector2D localPosition(spc.localX / spc.totalPath
                                               + spc.lorentzShift,
                                           spc.localY / spc.totalPath);
              // @todo remove unneccesary conversion
              size_t bin0          = binUtility.bin(localPosition, 0);
              size_t bin1          = binUtility.bin(localPosition, 1);
              size_t binSerialized = binUtility.serialize({bin0, bin1, 0});

              // create the indetifier
              spc.geoID.add(volumeKey, Acts::GeometryID::volume_mask);
              spc.geoID.add(layerKey, Acts::GeometryID::layer_mask);
              spc.geoID.add(moduleKey, Acts::GeometryID::sensitive_mask);
              spc.geoID.add(binSerialized, Acts::GeometryID::channel_mask);
              // create the truth for this - assume here muons
              Acts::ParticleProperties pProperties(
                  momentum, pMasses.mass[Acts::muon], 1., 13, particleBarcode);
              // the associated process vertex
              spc.pVertex = std::make_shared<Acts::ProcessVertex>(
                  Acts::ProcessVertex(position, 0., 0, {pProperties}, {}));

              // insert into the space point map
              FW::Data::insert(spacePoints,
                               volumeKey,
                               layerKey,
                               moduleKey,
                               hitParameters->position());

              clusters.push_back(spc);

            }  // hit moulde proection
          }    // hit element protection
        }      // hit loop

        clusterize(planarClusters, clusters, volumeKey, layerKey, moduleKey);

      }  // moudle loop
    }    // layer loop
  }      // volume loop

  // write the SpacePoints to the EventStore
  if (ctx.eventStore.add(m_cfg.spacePointsCollection, std::move(spacePoints))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }
  // write the clusters to the EventStore
  if (ctx.eventStore.add(m_cfg.clustersCollection, std::move(planarClusters))
      == FW::ProcessCode::ABORT) {
    return FW::ProcessCode::ABORT;
  }

  return FW::ProcessCode::SUCCESS;
}
