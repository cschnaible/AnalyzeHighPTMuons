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
		const reco::Track& track, const reco::BeamSpot& beamSpot, const reco::Vertex& PV) {

	reset();

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
			thisCovRow.push_back(track.covariance()[i][j]);
			thisWRow.push_back(Wmatrix[i][j]);
			thisCorrRow.push_back(CorrMatrix[i][j]);
		}
		cov.push_back(thisCovRow);
		w.push_back(thisWRow);
		corr.push_back(thisCorrRow);
	}
}

void FillCombinationInfo::fill(const reco::Track& tracker, const reco::Track& refit) {

	reset();

	AlgebraicSymMatrix55 wTrk = tracker.covariance();
	AlgebraicSymMatrix55 covTrk = tracker.covariance();
	AlgebraicVector5 parTrk = tracker.parameters();
	if (!wTrk.Invert()) {
		std::cout << "Cannot invert tracker cov matrix" << std::endl;
		std::cout << covTrk << std::endl;
	}

	AlgebraicSymMatrix55 wRef = refit.covariance();
	AlgebraicSymMatrix55 covRef = refit.covariance();
	AlgebraicVector5 parRef = refit.parameters();
	if (!wRef.Invert()) {
		std::cout << "Cannot invert refit cov matrix" << std::endl;
		std::cout << covRef << std::endl;
	}

	// Calculate cov, cov^{-1} matrices
	AlgebraicSymMatrix55 wComb = wTrk + wRef;
	AlgebraicSymMatrix55 covComb = wComb;
	if (!covComb.Invert()) {
		std::cout << "Cannot invert combined W matrix" << std::endl;
	}
	// Calculate correlation matrix
	AlgebraicMatrix55 diagCovComb;
	for (int i=0; i<5; i++) {
		for (int j=0; j<5; j++) {
			if (i==j) diagCovComb[i][j] = 1./std::sqrt(covComb[i][j]);
			else diagCovComb[i][j] = 0.;
		}
	}
	AlgebraicSymMatrix55 corrComb = ROOT::Math::Similarity(diagCovComb,covComb);

	K = (parRef[0]/covRef[0][0] + parTrk[0]/covTrk[0][0]) / (1./covRef[0][0] + 1./covTrk[0][0]);
	K_err = 1. / (1./covRef[0][0] + 1./covTrk[0][0]);

	// save stuff
	AlgebraicVector5 parComb = covComb * (wTrk*parTrk + wRef*parRef);
	for (int i=0; i<5; i++) {
		par.push_back(parComb[i]);
		std::vector<double> thisCovRow;
		std::vector<double> thisWRow;
		std::vector<double> thisCorrRow;
		for (int j=0; j<5; j++) {
			thisWRow.push_back(wComb[i][j]);
			thisCovRow.push_back(covComb[i][j]);
			thisCorrRow.push_back(corrComb[i][j]);
		}
		w.push_back(thisWRow);
		cov.push_back(thisCovRow);
		corr.push_back(thisCorrRow);
	}
	// Calculate match chi-square
	AlgebraicSymMatrix55 sumCov = covRef + covTrk;
	AlgebraicSymMatrix55 sumCovInv = sumCov;
	if (!sumCovInv.Invert()) std::cout << "Cannot invert trk+ref cov matrix" << std::endl;
	chi2 = ROOT::Math::Similarity(parRef-parTrk,sumCovInv);
}

void FillGenInfo::fill(const reco::GenParticle& muon) {
	reset();
	K = muon.charge()/muon.p();
	eta = muon.eta();
	phi = muon.phi();
	pt = muon.pt();
}

