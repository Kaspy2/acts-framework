
#include "FullVertexFitter.hpp"
#include "LinearizedTrackFactory.hpp"
#include "TrackAtVertex.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace {

	struct BilloirTrack
	{
		BilloirTrack(const Acts::BoundParameters& params, LinearizedTrack* lTrack) :
			originalParams(params),
			linTrack(lTrack) 
		{}
                
        BilloirTrack(const BilloirTrack& arg) = default;
       
		const Acts::BoundParameters originalParams;
		LinearizedTrack* linTrack;
		double chi2;
		Acts::ActsMatrixD<5,3> Di_mat;
		Acts::ActsMatrixD<5,3> Ei_mat;
		Acts::ActsSymMatrixD<3> Gi_mat;
		Acts::ActsSymMatrixD<3> Bi_mat;     // Bi = Di.T * Wi * Ei
		Acts::ActsSymMatrixD<3> Ci_inv;     // Ci = (Ei.T * Wi * Ei)^-1
		Acts::Vector3D Ui_vec;     // Ui = Ei.T * Wi * dqi
		Acts::ActsSymMatrixD<3> BCi_mat;       // BCi = Bi * Ci^-1
		Acts::ActsVectorD<5> delta_q;
	};

	struct BilloirVertex
	{
                BilloirVertex()
                {
                  A_mat.setZero();
                  T_vec.setZero();
                  BCB_mat.setZero();
                  BCU_vec.setZero();
                };

		Acts::ActsSymMatrixD<3> A_mat;              // T  = sum{Di.T * Wi * Di}
		Acts::Vector3D T_vec;              // A  = sum{Di.T * Wi * dqi}
		Acts::ActsSymMatrixD<3> BCB_mat;       // BCB = sum{Bi * Ci^-1 * Bi.T}
		Acts::Vector3D BCU_vec;       // BCU = sum{Bi * Ci^-1 * Ui}
	};



} // end anonymous namespace

template <typename BField>
Vertex FullVertexFitter<BField>::fit(const std::vector<Acts::BoundParameters>& paramVector) const
{
	bool debug = false;
	
	double chi2 = 1E10;
	double newChi2 = 0;
	unsigned int nTracks = paramVector.size();

	// Factory for linearizing tracks
	typename LinearizedTrackFactory<BField>::Config lt_config(m_cfg.bField);
	LinearizedTrackFactory<BField> linFactory(lt_config);

	std::vector<BilloirTrack> billoirTracks;

	std::vector<Acts::Vector3D> trackMomenta;

	Acts::Vector3D linPoint(m_cfg.startingPoint);

	std::unique_ptr<Vertex> fittedVertex = std::make_unique<Vertex>();

	for(int n_iter = 0; n_iter < m_cfg.maxIterations; ++n_iter)
	{
		billoirTracks.clear();

		newChi2 = 0;

		BilloirVertex billoirVertex;

		int i_track = 0;
		// iterate over all tracks
		for(const auto& trackParams : paramVector)
		{
			if (n_iter == 0){

				double phi 		= trackParams.parameters()[Acts::ParID_t::ePHI];
				double theta 	= trackParams.parameters()[Acts::ParID_t::eTHETA];
				double qop 		= trackParams.parameters()[Acts::ParID_t::eQOP];
				trackMomenta.push_back(Acts::Vector3D(phi, theta, qop));
				}

			LinearizedTrack* linTrack = linFactory.linearizeTrack(&trackParams, linPoint);

			double d0 = linTrack->parametersAtPCA()[Acts::ParID_t::eLOC_D0];
			double z0 = linTrack->parametersAtPCA()[Acts::ParID_t::eLOC_Z0];
			double phi = linTrack->parametersAtPCA()[Acts::ParID_t::ePHI];
			double theta = linTrack->parametersAtPCA()[Acts::ParID_t::eTHETA];
			double qOverP = linTrack->parametersAtPCA()[Acts::ParID_t::eQOP];

			// calculate f(V_0,p_0)  f_d0 = f_z0 = 0
			double f_phi = trackMomenta[i_track][0];
			double f_theta =trackMomenta[i_track][1];
			double f_qOverP = trackMomenta[i_track][2];

			BilloirTrack currentBilloirTrack(trackParams, linTrack);

			//calculate delta_q[i]
			currentBilloirTrack.delta_q[0] = d0;
			currentBilloirTrack.delta_q[1] = z0;
			currentBilloirTrack.delta_q[2] = phi - f_phi;
			currentBilloirTrack.delta_q[3] = theta - f_theta;
			currentBilloirTrack.delta_q[4] = qOverP -f_qOverP;
			
			// position jacobian (D matrix)
			Acts::ActsMatrixD<5,3> D_mat;
			D_mat = linTrack->positionJacobian();

			// momentum jacobian (E matrix)
			Acts::ActsMatrixD<5,3> E_mat;
			E_mat = linTrack->momentumJacobian();

			// cache some matrix multiplications
			Acts::ActsMatrixD<3,5> Dt_W_mat;
			Dt_W_mat.setZero();
			Acts::ActsMatrixD<3,5> Et_W_mat;
			Et_W_mat.setZero();
			Dt_W_mat = D_mat.transpose() * (linTrack->covarianceAtPCA().inverse().eval());
			Et_W_mat = E_mat.transpose() * (linTrack->covarianceAtPCA().inverse().eval());

			// compute billoir tracks
			currentBilloirTrack.Di_mat = D_mat;
			currentBilloirTrack.Ei_mat = E_mat;
			currentBilloirTrack.Gi_mat = Et_W_mat*E_mat;
			currentBilloirTrack.Bi_mat = Dt_W_mat * E_mat; // Di.T * Wi * Ei
			currentBilloirTrack.Ui_vec = Et_W_mat * currentBilloirTrack.delta_q; // Ei.T * Wi * dqi
			currentBilloirTrack.Ci_inv = (Et_W_mat * E_mat).inverse().eval(); // (Ei.T * Wi * Ei)^-1

			// sum up over all tracks
			billoirVertex.T_vec += Dt_W_mat * currentBilloirTrack.delta_q; // sum{Di.T * Wi * dqi}
			billoirVertex.A_mat += Dt_W_mat * D_mat ; // sum{Di.T * Wi * Di}

			// remember those results for all tracks
			currentBilloirTrack.BCi_mat = currentBilloirTrack.Bi_mat * currentBilloirTrack.Ci_inv; // BCi = Bi * Ci^-1

			// and some summed results
			billoirVertex.BCU_vec += currentBilloirTrack.BCi_mat * currentBilloirTrack.Ui_vec; // sum{Bi * Ci^-1 * Ui}
			billoirVertex.BCB_mat += currentBilloirTrack.BCi_mat * currentBilloirTrack.Bi_mat.transpose();// sum{Bi * Ci^-1 * Bi.T}
			
			billoirTracks.push_back(currentBilloirTrack);
			++i_track;

		} // end loop tracks

		// calculate delta (billoirFrameOrigin-position), might be changed by the beam-const
		Acts::Vector3D    V_del     = billoirVertex.T_vec - billoirVertex.BCU_vec;      // V_del = T-sum{Bi*Ci^-1*Ui}
		Acts::ActsSymMatrixD<3> V_wgt_mat = billoirVertex.A_mat - billoirVertex.BCB_mat;      // V_wgt = A-sum{Bi*Ci^-1*Bi.T}

		// TODO: if constraint: change accordingly

		// cov(delta_V) = V_wgt^-1
		Acts::ActsSymMatrixD<3> cov_delta_V_mat = V_wgt_mat.inverse().eval() ;

		// delta_V = cov_(delta_V) * V_del;
		Acts::Vector3D delta_V = cov_delta_V_mat * V_del;

		//--------------------------------------------------------------------------------------
		// start momentum related calculations

		std::vector<std::unique_ptr<Acts::ActsSymMatrixD<5>>> cov_delta_P_mat ( nTracks );
		std::vector<double> chi2PerTrack;
		chi2PerTrack.clear();


		i_track = 0;
		for(auto& bTrack : billoirTracks){
	
			Acts::Vector3D deltaP = ( bTrack.Ci_inv ) * ( bTrack.Ui_vec - bTrack.Bi_mat.transpose() * delta_V );
	
			// update track momenta
			trackMomenta[i_track][0] += deltaP[0];

			trackMomenta[i_track][1] += deltaP[1];

			trackMomenta[i_track][2] += deltaP[2];

			// TODO: correct for 2PI / PI periodicity

			//calculate 5x5 cov_delta_P matrix
			// d(d0,z0,phi,theta,qOverP)/d(x,y,z,phi,theta,qOverP)-transformation matrix
			Acts::ActsMatrixD<5,6> trans_mat;
			trans_mat.setZero();
			trans_mat ( 0,0 ) = bTrack.Di_mat ( 0,0 ); trans_mat ( 0,1 ) = bTrack.Di_mat ( 0,1 );
			trans_mat ( 1,0 ) = bTrack.Di_mat ( 1,0 ); trans_mat ( 1,1 ) = bTrack.Di_mat ( 1,1 ); trans_mat ( 1,2 ) = 1.;
			trans_mat ( 2,3 ) = 1.; trans_mat ( 3,4 ) = 1.; trans_mat ( 4,5 ) = 1.;
			
			//some intermediate calculations to get 5x5 matrix
			//cov(V,V)
			Acts::ActsSymMatrixD<3> V_V_mat; 
			V_V_mat.setZero();  
			V_V_mat = cov_delta_V_mat;
			//std::cout<<"V_V_mat = "<<V_V_mat<<std::endl; 

			//cov(V,P)
			Acts::ActsSymMatrixD<3> 
			V_P_mat; V_P_mat.setZero(); 
			V_P_mat =  -cov_delta_V_mat*bTrack.Gi_mat*bTrack.Ci_inv;
			//std::cout<<"V_P_mat = "<<V_P_mat<<std::endl;

			//cov(P,P)
			Acts::ActsSymMatrixD<3> P_P_mat; 
			P_P_mat.setZero(); 
			P_P_mat =  bTrack.Ci_inv + bTrack.BCi_mat.transpose() * cov_delta_V_mat*bTrack.BCi_mat;
			//std::cout<<"P_P_mat = "<<P_P_mat<<std::endl;

			Acts::ActsSymMatrixD<6> cov_mat;
			cov_mat.setZero();
			cov_mat.block<3,3>(0,3) = V_P_mat;
			cov_mat.block<3,3>(3,0) = V_P_mat.transpose() ;
			cov_mat.block<3,3>(0,0) = V_V_mat ;
			cov_mat.block<3,3>(3,3) = P_P_mat ;

			// cov_delta_P calculation
			//          std::cout<<"trans_mat = "<<trans_mat<<std::endl;
			cov_delta_P_mat[i_track] = std::make_unique<Acts::ActsSymMatrixD<5>>(trans_mat*cov_mat*trans_mat.transpose());

			// Calculate chi2 per track. 
			bTrack.chi2= ( ( bTrack.delta_q - bTrack.Di_mat*delta_V - bTrack.Ei_mat*deltaP ).transpose() 
				* bTrack.linTrack->covarianceAtPCA().inverse().eval() 
				* ( bTrack.delta_q - bTrack.Di_mat*delta_V - bTrack.Ei_mat*deltaP ) )[0];
			newChi2 += bTrack.chi2;
			//std::cout << "track " << i_track << " with chi2 " << bTrack.chi2 << std::endl;

			++i_track;
		}

		// assign new linearization point (= new vertex position in global frame)
		linPoint += delta_V;

		if (newChi2 < chi2){
			chi2 = newChi2;

			Acts::Vector3D vertex(linPoint);

			if(debug) std::cout << "\tset new vertex " << vertex << " with chi2: " << chi2 << std::endl;


			fittedVertex->setPosition(vertex);
			fittedVertex->setCovariance(cov_delta_V_mat);

			std::vector<std::unique_ptr<TrackAtVertex>> tracksAtVertex;
			
			Acts::PerigeeSurface perigee(vertex);
			int i_track = 0;
			for(auto& bTrack : billoirTracks){

				// new refitted trackparameters
				Acts::TrackParametersBase::ParVector_t paramVec;
				paramVec << 0.,0.,trackMomenta[i_track](0),trackMomenta[i_track](1),trackMomenta[i_track](2);

				Acts::BoundParameters refittedParams(std::move(cov_delta_P_mat[i_track]), 
									  paramVec, perigee);

				std::unique_ptr<TrackAtVertex> trackVx
						= std::make_unique<TrackAtVertex>(bTrack.chi2, refittedParams, bTrack.originalParams); 
				tracksAtVertex.push_back(std::move(trackVx));
				++i_track;
			}

			fittedVertex->setTracksAtVertex(tracksAtVertex);

		}

	} // end loop iterations

	return *fittedVertex;
	

}