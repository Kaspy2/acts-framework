// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "ACTFW/EventData/DataContainers.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTS/Utilities/GeometryID.hpp"
#include "DigitizationAlgorithm.hpp"

namespace FW {

class ClusterDigitizationAlgorithm : public FW::DigitizationAlgorithm
{
public:
  /// @brief Collects recorded hits on surfaces, forms clusters from them and
  /// stores them on the whiteboard
  /// @param ctx access to the whiteboard data
  /// @return Indicator if the algorithm worked properly
  ProcessCode
  execute(AlgorithmContext ctx) const final override;

private:
  /// @brief This structure serves as data container. It allows moving data that
  /// belongs to a single particle on a single surface through functions without
  /// using quite long arguments.
  struct SingleParticleCluster
  {
    // Storage of the surface of the hit record
    Acts::Surface const* hitSurface = nullptr;
    // Storage of the Identifier of the hit
    Acts::GeometryID geoID;
    // Mean local x-coordinate of the hit on several digitization cells
    double localX = 0.;
    // Mean local y-coordinate of the hit on several digitization cells
    double localY = 0.;
    // Total path of the particle through the digitization module
    double totalPath = 0.;
    // Lorentz shift
    double lorentzShift = 0.;
    // Steps through the digitization module3
    std::vector<Acts::DigitizationStep> dSteps;
    // Digitization cells that recorded a signal from that particle
    std::vector<Acts::DigitizationCell> usedCells;
    // Truth vertex information
    std::shared_ptr<Acts::ProcessVertex> pVertex;
  };
  // TODO: config and constructor need to be modified (e.g. for commonCorners
  // selection of clustering)

  /// @brief Merges the hits of several particles to clusters
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @return Merged collection from @p clusters
  /// @note The structure of the return type is [Cluster index][Signals related
  /// to the cluster]
  std::vector<std::vector<SingleParticleCluster>>
  mergeSingleParticleClusters(
      const std::vector<SingleParticleCluster>& clusters) const;

  /// @brief Calculates the mean of all lorentz shifts. This is used for the
  /// case that it differs along the surface.
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @return Mean of the lorentz shifts
  double
  meanLorentzShift(const std::vector<SingleParticleCluster>& clusters) const;

  /// @brief Calculates the mean of the local postions of the clusters.
  /// @note This function does not apply a lorentz shift.
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @return Mean of the local postions of the cluster
  Acts::Vector2D
  localPosition(const std::vector<SingleParticleCluster>& clusters) const;

  /// @brief Calculates the covariance matrix of the cluster
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @param mean mean of the local postions of the cluster without lorentz
  /// shift correction
  /// @return Covariance matrix of the cluster
  Acts::ActsSymMatrixD<2>
  covariance(const std::vector<SingleParticleCluster>& clusters,
             const Acts::Vector2D                      mean) const;

  /// @brief Forms clusters and creates Acts::PlanarModuleCluster-objects out of
  /// a collection of hits on a surface
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @return Collection of Acts::PlanarModuleClusters-objects
  std::vector<Acts::PlanarModuleCluster>
  formClusters(const std::vector<SingleParticleCluster>& clusters) const;

  /// @brief Forms clusters out of hits on a surface and prepares their storage
  /// on the whiteboard
  /// @param planarClusters data container that will be written on the
  /// whiteboard carrying the data
  /// @param clusters collection of data about particles interacting with the
  /// surface
  /// @param volumeKey identifier of the volume of the surface
  /// @param layerKey identifier of the layer of the surface
  /// @param moduleKey identifier of the module of the surface
  /// @note the *Key parameters are used to construct a unique identifier
  void
  clusterize(
      FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>& planarClusters,
      const std::vector<SingleParticleCluster>& clusters,
      const geo_id_value                        volumeKey,
      const geo_id_value                        layerKey,
      const geo_id_value                        moduleKey) const;
};

}  // namespace FW
