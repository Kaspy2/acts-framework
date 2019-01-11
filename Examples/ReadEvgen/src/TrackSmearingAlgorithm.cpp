// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TrackSmearingAlgorithm.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "ACTFW/Random/RandomNumberDistributions.hpp"
#include "ACTFW/Random/RandomNumbersSvc.hpp"

#include <tbb/tbb.h>

#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
#include "Acts/Vertexing/FullVertexFitter.hpp"

struct Config
  {
    std::shared_ptr<FW::RandomNumbersSvc> randomNumbers = nullptr;

    std::array<double, 2> gaussParameters   = {{0., 1.}};
    std::array<double, 2> uniformParameters = {{0., 1.}};
    std::array<double, 2> landauParameters  = {{0., 1.}};
    std::array<double, 2> gammaParameters   = {{0., 1.}};
    int poissonParameter = 40;

    size_t drawsPerEvent = 0;
  };


FWE::TrackSmearingAlgorithm::TrackSmearingAlgorithm(const Config& cfg, Acts::Logging::Level level)
  : FW::BareAlgorithm("TrackSmearing", level), m_cfg(cfg)
{
}





FW::ProcessCode
FWE::TrackSmearingAlgorithm::execute(FW::AlgorithmContext context) const
{

	const double eta_cut = 3.0;

	// Define parameter for pt-dependent IP resolution
	// of the form sigma_d/z(p_t[GeV]) = A*exp(-B*p_t[GeV]) + C
	const double ipResA = 100.7439 * Acts::units::_um;
	const double ipResB = 0.23055;
	const double ipResC = 20. * Acts::units::_um;

	// Create and fill input event
	const std::vector<FW::Data::SimVertex<>>* inputEvent = nullptr;
	if (context.eventStore.get(m_cfg.collection, inputEvent) == FW::ProcessCode::ABORT)
	{
    	return FW::ProcessCode::ABORT;
	}

	// Define perigee surface center coordinates
	const Acts::Vector3D surfaceCenter(0.,0.,0.);

	std::shared_ptr<Acts::PerigeeSurface> perigeeSurface = std::make_shared<Acts::PerigeeSurface>(surfaceCenter);

	// Set up b-field and stepper
	Acts::ConstantBField bField(Acts::Vector3D(0.,0.,1.)*Acts::units::_T);
	Acts::EigenStepper<Acts::ConstantBField> stepper(bField);
	
	// Set up propagator with void navigator
	Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>> propagator(stepper);

	// Set up propagator options
	Acts::PropagatorOptions<> options;

	// Create random number generator and spawn gaussian distribution
	FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

	// typedef for simplicity
	using BoundParamsVector = std::vector<Acts::BoundParameters>;

	// Vector to store smrdTracksAtVtx for all vertices of event
	std::vector<BoundParamsVector> smrdTrackCollection;

	// Vector to store true vertices positions
	std::vector<Acts::Vector3D> trueVertices;

	int count_v = 0; // vertex count
	for (auto& vtx: (*inputEvent)){
		count_v++;
		// Take only first vertex now
		//if (count_v != 1) continue;

		// Vector to store smeared tracks at current vertex
		BoundParamsVector smrdTracksAtVtx;

		// Iterate over all particle emerging from current vertex
		for (auto const& particle : vtx.out){

			const Acts::Vector3D& ptclMom = particle.momentum();
			// Calculate pseudo-rapidity
			const double eta = Acts::VectorHelpers::eta(ptclMom);
			// Only charged particles for |eta| < 2.5
			if (particle.q() !=0 && std::abs(eta) < eta_cut) 
			{
				// Define start track params
				Acts::CurvilinearParameters 
					start(nullptr, particle.position(), ptclMom, particle.q());

				// Run propagator
				const auto result = propagator.propagate(start, *perigeeSurface, options);

				if (result.status == Acts::Status::SUCCESS){

					const auto& perigeeParameters = result.endParameters->parameters(); // (d0, z0, phi, theta,q/p)

					if (std::abs(perigeeParameters[0]) > 30 || std::abs(perigeeParameters[1]) > 200){
						continue;
					}

					// Calculate pt-dependent IP resolution
					const double pclPt = 
							Acts::VectorHelpers::perp(ptclMom)/Acts::units::_GeV;
					const double ipRes = ipResA * std::exp(-ipResB*pclPt) + ipResC;

					// except for IP resolution, following variances are rough guesses
					// Gaussian distribution for IP resolution
					FW::GaussDist gaussDist_IP(0., ipRes);
					// Gaussian distribution for angular resolution
					FW::GaussDist gaussDist_angular(0., 0.1);
					// Gaussian distribution for q/p (momentum) resolution
					FW::GaussDist gaussDist_qp(0., 0.1*perigeeParameters[4]);

					double rn_d0 = gaussDist_IP(rng);
					double rn_z0 = gaussDist_IP(rng);
					double rn_ph = gaussDist_angular(rng);
					double rn_th = gaussDist_angular(rng);
					double rn_qp = gaussDist_qp(rng);

					double smrd_d0 	= perigeeParameters[0] + rn_d0;
					double smrd_z0	= perigeeParameters[1] + rn_z0;
					double smrd_phi 	= perigeeParameters[2] + rn_ph;
					double smrd_theta	= perigeeParameters[3] + rn_th;
					double srmd_qp	= perigeeParameters[4] + rn_qp;

					// smearing can bring theta out of range ->close to beam line -> discard
					if(smrd_theta < 0 || smrd_theta > M_PI){
						continue;
					}

					double new_eta = -log(tan(smrd_theta/2));
					if(std::abs(new_eta) > eta_cut) continue;
					
					Acts::TrackParametersBase::ParVector_t paramVec;
					paramVec << smrd_d0, smrd_z0, smrd_phi, smrd_theta, srmd_qp;

					// Fill vector of smeared tracks
					std::unique_ptr<Acts::ActsSymMatrixD<5>> covMat = std::make_unique<Acts::ActsSymMatrixD<5>>();
					covMat->setZero();
					(*covMat)(0,0) = rn_d0*rn_d0;
					(*covMat)(1,1) = rn_z0*rn_z0;
					(*covMat)(2,2) = rn_ph*rn_ph;
					(*covMat)(3,3) = rn_th*rn_th;
					(*covMat)(4,4) = rn_qp*rn_qp;

					smrdTracksAtVtx.push_back(Acts::BoundParameters(std::move(covMat), paramVec, perigeeSurface));
				}
			}
		}

		if(!smrdTracksAtVtx.empty()){
			smrdTrackCollection.push_back(smrdTracksAtVtx);

			// Store true vertex position for debugging only
			trueVertices.push_back(vtx.position); 
		}
	}

	assert(smrdTrackCollection.size() == trueVertices.size());
	std::cout << "Total amount of vertices in event: " << trueVertices.size() << std::endl;

	// Set up vertex fitter
	Acts::FullVertexFitter<Acts::ConstantBField>::Config vertexFitterCfg(bField);
	Acts::FullVertexFitter<Acts::ConstantBField> vertexFitter(vertexFitterCfg);


	tbb::parallel_for(
      tbb::blocked_range<size_t>(0, smrdTrackCollection.size()),
      [&](const tbb::blocked_range<size_t>& r) {
        for (size_t event_idx = r.begin(); event_idx != r.end(); ++event_idx) {

			BoundParamsVector currentParamVectorAtVtx = smrdTrackCollection[event_idx];
			if (currentParamVectorAtVtx.size() > 1){

				Acts::Vector3D currentTrueVtx = trueVertices[event_idx];

				std::cout << "Vertex " << event_idx << " with " << currentParamVectorAtVtx.size() << " tracks" << std::endl;

				Acts::Vertex fittedVertex = vertexFitter.fit(currentParamVectorAtVtx);

				Acts::Vector3D diffVtx = currentTrueVtx - fittedVertex.position();

				std::cout << "true vertex:   "
					 << "(" << currentTrueVtx[0] << "," <<  currentTrueVtx[1] << ","<<   currentTrueVtx[2] << ")" << std::endl;
				std::cout << "fitted vertex: " 
					<< "(" << fittedVertex.position()[0] << "," <<  fittedVertex.position()[1] << ","<<   fittedVertex.position()[2] << ")\n" << std::endl;
			}
        }
      });

	// TODO: store fitted vertices and write to eventstore

	if(context.eventStore.add(m_cfg.collectionOut, std::move(smrdTrackCollection))
		!= FW::ProcessCode::SUCCESS)
	{
		return FW::ProcessCode::ABORT;
	}

	
	return FW::ProcessCode::SUCCESS;
}











