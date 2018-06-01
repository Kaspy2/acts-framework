// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "DigitizationAlgorithm.hpp"
#include "ACTS/Utilities/GeometryID.hpp"
#include "ACTS/Digitization/DigitizationModule.hpp"
#include "ACTS/Digitization/PlanarModuleCluster.hpp"
#include "ACTFW/EventData/DataContainers.hpp"

namespace FW {

class ClusterDigitizationAlgorithm : public FW::DigitizationAlgorithm
{
public:

	/// @brief Collects recorded hits on surfaces, forms clusters from them and stores them on the whiteboard
	/// @param ctx access to the whiteboard data
	/// @return Indicator if the algorithm worked properly
	ProcessCode
	execute(AlgorithmContext ctx) const final override;
	
private:

	/// @brief This structure serves as data container. It allows moving data that belongs to a single particle on a single surface through functions without using quite long arguments.
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
		double totalPath = 0.;
		double lorentzShift = 0.;
		
		std::vector<Acts::DigitizationStep> dSteps;
		std::vector<Acts::DigitizationCell> usedCells;
		
		std::shared_ptr<Acts::ProcessVertex> pVertex;
		
	};
	// TODO: config and constructor need to be modified (e.g. for commonCorners selection of clustering)
	
	std::vector<std::vector<SingleParticleCluster>>
	mergeSingleParticleClusters(const std::vector<SingleParticleCluster>& clusters) const;
	
	double
	meanLorentzShift(const std::vector<SingleParticleCluster>& clusters) const;
	
	Acts::Vector2D
	localPosition(const std::vector<SingleParticleCluster>& clusters) const;

	Acts::ActsSymMatrixD<2>
	covariance(const std::vector<SingleParticleCluster>& clusters, const Acts::Vector2D mean) const;

	std::vector<Acts::PlanarModuleCluster>
	formClusters(const std::vector<SingleParticleCluster>& clusters) const;

	void
	clusterize(FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>& planarClusters, const std::vector<SingleParticleCluster>& clusters, const geo_id_value volumeKey, const geo_id_value layerKey, const geo_id_value moduleKey) const;
	
};

}  // namespace FW
