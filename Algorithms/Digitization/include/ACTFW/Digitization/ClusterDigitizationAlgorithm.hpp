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
	ProcessCode
	execute(AlgorithmContext ctx) const final override;
	
private:

	struct SingleParticleCluster
	{
		Acts::Surface const* hitSurface = nullptr;
		
		Acts::GeometryID geoID;
		
		double localX = 0.;
		double localY = 0.;
		double totalPath = 0.;
		double lorentzShift = 0.;
		
		std::vector<Acts::DigitizationCell> usedCells;
		
		std::shared_ptr<Acts::ProcessVertex> pVertex;
		
	};
	// TODO: config und constructor, damit commoncorners rein kann
	
	std::vector<std::vector<SingleParticleCluster>>
	mergeSingleParticleClusters(const std::vector<SingleParticleCluster>& clusters) const;
	
	double
	meanLorentzShift(const std::vector<SingleParticleCluster>& clusters) const;
	
	Acts::Vector2D
	localPosition(const std::vector<SingleParticleCluster>& clusters) const;

	Acts::ActsSymMatrixD<2>
	covariance(const std::vector<SingleParticleCluster>& cluster, Acts::Vector2D mean) const;

	std::vector<Acts::PlanarModuleCluster>
	formClusters(const std::vector<SingleParticleCluster>& clusters) const;

	void
	clusterize(FW::DetectorData<geo_id_value, Acts::PlanarModuleCluster>& planarClusters, const std::vector<SingleParticleCluster>& clusters, const geo_id_value volumeKey, const geo_id_value layerKey, const geo_id_value moduleKey) const;
	
};

}  // namespace FW
