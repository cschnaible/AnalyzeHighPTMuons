#include "HighPTMuons/Analysis/include/FillHighPTInfo.h"

void FillEventInfo::fill(
		const edm::Event& iEvent, const reco::BeamSpot& beamSpot, const reco::Vertex& PV){
  reset();
  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.eventAuxiliary().luminosityBlock();
	beamSpot_pos.push_back(beamSpot.x0());
	beamSpot_pos.push_back(beamSpot.y0());
	beamSpot_pos.push_back(beamSpot.z0());
	primaryVertex_pos.push_back(PV.x());
	primaryVertex_pos.push_back(PV.y());
	primaryVertex_pos.push_back(PV.z());
}

void FillTrackInfo::fill(
		const reco::Track& track, const reco::BeamSpot& beamSpot, const reco::Vertex& PV,const float& scaleCov) {

	reset();

	good = true;
	K = track.qoverp();
	K_err = track.qoverpError();
	chi2 = track.chi2();
	pt = track.pt();
	nValidHits = track.numberOfValidHits();
	nValidMuonHits = track.hitPattern().numberOfValidMuonHits();
	nValidCSCHits = track.hitPattern().numberOfValidMuonCSCHits();
	nValidDTHits = track.hitPattern().numberOfValidMuonDTHits();
	nValidRPCHits = track.hitPattern().numberOfValidMuonRPCHits();
	nValidMuonStations = track.hitPattern().muonStationsWithValidHits();
	nValidCSCStations = track.hitPattern().cscStationsWithValidHits();
	nValidDTStations = track.hitPattern().dtStationsWithValidHits();
	nValidRPCStations = track.hitPattern().rpcStationsWithValidHits();
	dxyBS = track.dxy(beamSpot.position());
	dxyPV = track.dxy(PV.position());
	dszBS = track.dsz(beamSpot.position());
	dszPV = track.dsz(PV.position());
	phi = track.phi();
	eta = track.eta();
	theta = track.theta();
	ndof = track.ndof();
	pos.push_back(track.vx());
	pos.push_back(track.vy());
	pos.push_back(track.vz());
	mom.push_back(track.px());
	mom.push_back(track.py());
	mom.push_back(track.pz());
	AlgebraicSymMatrix55 Wmatrix = track.covariance();
	AlgebraicSymMatrix55 CovMatrix = track.covariance();
	CovMatrix*=scaleCov;
	Wmatrix*=scaleCov;
	if (!Wmatrix.Invert()) {
		std::cout << "Cannot invert cov matrix" << std::endl;
		std::cout << track.parameters() << std::endl;
		std::cout << CovMatrix << std::endl;
	}
	AlgebraicMatrix55 diagCov;
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (i==j) diagCov[i][j] = 1./std::sqrt(CovMatrix[i][j]);
			else diagCov[i][j] = 0.;
		}
	}
	AlgebraicSymMatrix55 CorrMatrix = ROOT::Math::Similarity(diagCov,CovMatrix);
	for (int i=0; i<5; i++) {
		par.push_back(track.parameters()[i]);
		std::vector<double> thisCovRow;
		std::vector<double> thisWRow;
		std::vector<double> thisCorrRow;
		for (int j=0; j<5; j++) {
			thisCovRow.push_back(CovMatrix[i][j]);
			thisWRow.push_back(Wmatrix[i][j]);
			thisCorrRow.push_back(CorrMatrix[i][j]);
		}
		cov.push_back(thisCovRow);
		w.push_back(thisWRow);
		corr.push_back(thisCorrRow);
	}
}

void FillCombinationInfo::fill(const reco::Track& tracker, const reco::Track& refit,const float& scaleCov, const float& correlation) {
	// Combination of two correlated 5D estimates
	// Using notation from 
	// https://www.sciencedirect.com/science/article/pii/S0168900203003292

	reset();

	AlgebraicSymMatrix55 covTrk = tracker.covariance();
	AlgebraicVector5 parTrk = tracker.parameters();

	AlgebraicSymMatrix55 covRef = refit.covariance();
	AlgebraicVector5 parRef = refit.parameters();

	typedef ROOT::Math::SMatrix<double,10,10,ROOT::Math::MatRepSym<double,10> > SymMatrix1010;
	typedef ROOT::Math::SMatrix<double,10,5,ROOT::Math::MatRepStd<double,10,5> > Matrix105;
	typedef ROOT::Math::SMatrix<double,5,10,ROOT::Math::MatRepStd<double,5,10> > Matrix510;
	typedef ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepStd<double,5,5> > Matrix55;
	typedef ROOT::Math::SVector<double,10> Vector10;

	// Set 10x1 measurement vector
	Vector10 y;
	for (int i=0; i<5; i++) {
		y(i) = parRef(i);
		y(i+5) = parTrk(i);
	}
	/*
	std::cout << "y" << std::endl;
	std::cout << y << std::endl;
	*/

	// Set 10x10 covariance matrix
	SymMatrix1010 M;
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (j>i) continue;
			// Upper left block is refit covariance
			M(i,j) = covRef(i,j)*scaleCov;
			// Bottom right block is tracker covariance
			M(i+5,j+5) = covTrk(i,j);
		}
	}
	// Off diagional blocks are correlation between tracker and refit
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (i==0 || j==0) {
				M(i,j+5) = 0.; // zero curvature correlation
				continue;
			}
			// correlation = 1 for all other parameters
			M(i,j+5) = correlation * std::sqrt(covRef(i,i)*scaleCov) * std::sqrt(covTrk(j,j));
		}
	}
	SymMatrix1010 Minv = M;
	if (!Minv.Invert()) std::cout << "Cannot invert 10x10 trk+ref covariance matrix" << std::endl;

	/*
	std::cout << "M" << std::endl;
	std::cout << M << std::endl;
	std::cout << "Minv" << std::endl;
	std::cout << Minv << std::endl;
	*/

	// Define U projection matrix
	Matrix105 U;
	U(0,0) = 1.; U(0,1) = 0.; U(0,2) = 0.; U(0,3) = 0.; U(0,4) = 0.;
	U(1,0) = 0.; U(1,1) = 1.; U(1,2) = 0.; U(1,3) = 0.; U(1,4) = 0.;
	U(2,0) = 0.; U(2,1) = 0.; U(2,2) = 1.; U(2,3) = 0.; U(2,4) = 0.;
	U(3,0) = 0.; U(3,1) = 0.; U(3,2) = 0.; U(3,3) = 1.; U(3,4) = 0.;
	U(4,0) = 0.; U(4,1) = 0.; U(4,2) = 0.; U(4,3) = 0.; U(4,4) = 1.;
	U(5,0) = 1.; U(5,1) = 0.; U(5,2) = 0.; U(5,3) = 0.; U(5,4) = 0.;
	U(6,0) = 0.; U(6,1) = 1.; U(6,2) = 0.; U(6,3) = 0.; U(6,4) = 0.;
	U(7,0) = 0.; U(7,1) = 0.; U(7,2) = 1.; U(7,3) = 0.; U(7,4) = 0.;
	U(8,0) = 0.; U(8,1) = 0.; U(8,2) = 0.; U(8,3) = 1.; U(8,4) = 0.;
	U(9,0) = 0.; U(9,1) = 0.; U(9,2) = 0.; U(9,3) = 0.; U(9,4) = 1.;
	Matrix510 UT = ROOT::Math::Transpose(U);

	/*
	std::cout << "U" << std::endl;
	std::cout << U << std::endl;
	std::cout << "UT" << std::endl;
	std::cout << UT << std::endl;
	*/

	Matrix55 UTMinvU = UT*(Minv*U);
	Matrix55 covComb = UTMinvU;
	if (!covComb.Invert()) 
		std::cout << "Cannot invert UT*Minv*U to get combined covariance matrix" << std::endl;
	Matrix510 UTMinv = UT*Minv;
	
	/*
	std::cout << "UTMinvU" << std::endl;
	std::cout << UTMinvU << std::endl;
	std::cout << "UTMinv" << std::endl;
	std::cout << UTMinv << std::endl;
	std::cout << "C*UTMinv" << std::endl;
	std::cout << C*UTMinv << std::endl;
	std::cout << "C*Cinv" << std::endl;
	std::cout << (C*UTMinv)*U << std::endl;
	*/

	AlgebraicVector5 parComb = (covComb*UTMinv)*y;

	/*
	std::cout << "Combined parameters" << std::endl;
	std::cout << parComb << std::endl;
	std::cout << "Combined covariance" << std::endl;
	std::cout << covComb << std::endl;
	*/

	/*
	// Calculate correlation matrix
	AlgebraicMatrix55 diagCovComb;
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (i==j) diagCovComb[i][j] = 1./std::sqrt(covComb[i][j]);
			else diagCovComb[i][j] = 0.;
		}
	}
	AlgebraicSymMatrix55 corrComb = ROOT::Math::Similarity(diagCovComb,covComb);
	*/

	// save stuff
	for (int i=0; i<5; i++) {
		par.push_back(parComb[i]);
		std::vector<double> thisCovRow;
		std::vector<double> thisWRow;
		//std::vector<double> thisCorrRow;
		for (int j=0; j<5; j++) {
			thisCovRow.push_back(covComb(i,j));
			thisWRow.push_back(UTMinvU(i,j));
			//thisCorrRow.push_back(corrComb(i,j));
		}
		cov.push_back(thisCovRow);
		w.push_back(thisWRow);
		//corr.push_back(thisCorrRow);
	}

	// Calculate match chi-square
	AlgebraicSymMatrix55 sumCov = covRef + covTrk;
	AlgebraicSymMatrix55 sumCovInv = sumCov;
	if (!sumCovInv.Invert()) std::cout << "Cannot invert trk+ref cov matrix" << std::endl;
	chi2 = ROOT::Math::Similarity(parRef-parTrk,sumCovInv);

	good = true;
	
	// double check with no correlation combination
	/*
	AlgebraicSymMatrix55 wRef = refit.covariance();
	AlgebraicSymMatrix55 wTrk = tracker.covariance();
	if (!wTrk.Invert()) std::cout << "Cannot invert covTrk" << std::endl;
	if (!wRef.Invert()) std::cout << "Cannot invert covRef" << std::endl;
	AlgebraicSymMatrix55 wComb = wTrk + wRef;
	AlgebraicSymMatrix55 covCombOld = wComb;
	if (!covCombOld.Invert()) std::cout << "Cannot invert wComb" << std::endl;
	AlgebraicVector5 parCombOld = covComb * (wTrk*parTrk + wRef*parRef);
	std::cout << "old combined parameters" << std::endl;
	std::cout << parCombOld << std::endl;
	std::cout << "old combined covariance" << std::endl;
	std::cout << covCombOld << std::endl;
	*/
}

void FillGenInfo::fill(const reco::GenParticle& muon) {
	reset();
	K = muon.charge()/muon.p();
	eta = muon.eta();
	phi = muon.phi();
	pt = muon.pt();
}

