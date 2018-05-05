#ifndef HIGHPT_MUONTRACKTUPLES_FILLTRACKINFO_H
#define HIGHPT_MUONTRACKTUPLES_FILLTRACKINFO_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>

#include "TFile.h"
#include "TTree.h"
//#include "TString.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class HighPTTreeContainer {
	public:
		HighPTTreeContainer(std::string treeName, std::string treeTitle) {
			edm::Service<TFileService> fs;
			tree = fs->make<TTree>(TString(treeName),TString(treeTitle));
		}
		void write() {tree->Write();}
		void fill() {tree->Fill();}
		TTree * tree;
};

class FillHighPTInfo {
	public:
		FillHighPTInfo(HighPTTreeContainer& tree) : fTree(&tree) {reset();};
		virtual ~FillHighPTInfo() {};
	protected:
		virtual void reset() {};

		// Book single variable
		template<class T>
			void book(std::string name, T& var, std::string type) {
				fTree->tree->Branch(name.c_str(), &var, (name+"/"+type).c_str());
			}
		// Book vector
		template<class T>
			void book(std::string name, std::vector<T>& varv) {
				fTree->tree->Branch(name.c_str(), &varv);
			}

		HighPTTreeContainer * fTree;
};


class FillEventInfo : public FillHighPTInfo {
	public:

		FillEventInfo(HighPTTreeContainer& tree) :FillHighPTInfo(tree) {
			book("event",event,"I");
			book("run"  ,run,"I");
			book("lumi", lumi,"I");
			book("beamSpot_pos", beamSpot_pos);
			book("primaryVertex_pos", primaryVertex_pos);
		}
		virtual ~FillEventInfo() {};

	private:
		int event;
		int run;
		int lumi;
		std::vector<double> beamSpot_pos, primaryVertex_pos;

		virtual void reset(){
			event = -1;
			run = -1;
			lumi = -1;
			beamSpot_pos.clear();
			primaryVertex_pos.clear();
		}

		public:

		void fill(const edm::Event& iEvent, const reco::BeamSpot& beamSpot, const reco::Vertex& PV);

};


class FillTrackInfo : public FillHighPTInfo {
	public:
		FillTrackInfo(HighPTTreeContainer& tree, std::string trackType) : 
			FillHighPTInfo(tree),
			trackType(trackType)	
	{
			book(trackType+"_nValidHits", nValidHits, "I");
			book(trackType+"_nValidMuonHits", nValidMuonHits, "I");
			book(trackType+"_nValidCSCHits", nValidCSCHits, "I");
			book(trackType+"_nValidDTHits", nValidDTHits, "I");
			book(trackType+"_nValidRPCHits", nValidRPCHits, "I");
			book(trackType+"_nValidMuonStations", nValidMuonStations, "I");
			book(trackType+"_nValidCSCStations", nValidCSCStations, "I");
			book(trackType+"_nValidDTStations", nValidDTStations, "I");
			book(trackType+"_nValidRPCStations", nValidRPCStations, "I");
			book(trackType+"_ndof",ndof,"I");
			book(trackType+"_chi2",chi2,"D");
			book(trackType+"_K",K,"D");
			book(trackType+"_K_err",K_err,"D");
			book(trackType+"_phi", phi, "D");
			book(trackType+"_eta", eta, "D");
			book(trackType+"_theta", theta, "D");
			book(trackType+"_dxyBS", dxyBS, "D");
			book(trackType+"_dxyPV", dxyPV, "D");
			book(trackType+"_dszBS", dszBS, "D");
			book(trackType+"_dszPV", dszPV, "D");
			book(trackType+"_pt",pt,"D");
			book(trackType+"_par",par);
			book(trackType+"_cov",cov);
			book(trackType+"_pos",pos);
			book(trackType+"_mom",mom);
			book(trackType+"_w",w);
			book(trackType+"_corr",corr);
		}
		TString trackType;
		virtual ~FillTrackInfo() {};

	private:
		int nValidHits, nValidMuonHits, nValidCSCHits, nValidDTHits, nValidRPCHits;
		int nValidMuonStations, nValidCSCStations, nValidDTStations, nValidRPCStations;
		int ndof;
		double chi2,K,K_err,phi,eta,theta,dxyBS,dxyPV,dszBS,dszPV,pt;
		std::vector<double> par, mom, pos;
		std::vector< std::vector<double> > cov,w,corr;
		
		virtual void reset() {
			nValidHits = -1;
			nValidMuonHits = -1;
			nValidCSCHits = -1;
			nValidDTHits = -1;
			nValidRPCHits = -1;
			nValidMuonStations = -1;
			nValidCSCStations = -1;
			nValidDTStations = -1;
			nValidRPCStations = -1;
			ndof = -1;
			chi2 = -1;
			K = -999.;
			K_err = -999.;
			phi = -999.;
			eta = -999.;
			theta = -999.;
			dxyBS = -999.;
			dxyPV = -999.;
			dszBS = -999.;
			dszPV = -999.;
			pt = -999.;
			par.clear();
			mom.clear();
			pos.clear();
			cov.clear();
			w.clear();
			corr.clear();
		}
		
		public: 

		void fill(const reco::Track& track, const reco::BeamSpot& beamSpot, const reco::Vertex& PV);

};
	
class FillCombinationInfo : public FillHighPTInfo {
	public:
		FillCombinationInfo(HighPTTreeContainer& tree, std::string trackType) : 
			FillHighPTInfo(tree),trackType(trackType) 
		{
			book(trackType+"_K", K, "D");
			book(trackType+"_K_err", K_err, "D");
			book(trackType+"_chi2", chi2, "D");
			book(trackType+"_par", par);
			book(trackType+"_w", w);
			book(trackType+"_cov", cov);
			book(trackType+"_corr", corr);
		}
		TString trackType;
		virtual ~FillCombinationInfo() {};
	private:
		double K,K_err,chi2;
		std::vector<double> par;
		std::vector< std::vector<double> > w,cov,corr;

		virtual void reset() {
			K = -999.;
			K_err = -999.;
			par.clear();
			w.clear();
			cov.clear();
			corr.clear();
		}

	public:
		void fill(const std::vector<double>& par, const std::vector< std::vector<double> >& cov);
		void fill(const reco::Track& tracker, const reco::Track& refit);
};

class FillGenInfo : public FillHighPTInfo {
	public:
		FillGenInfo(HighPTTreeContainer& tree) : FillHighPTInfo(tree) {
			book("gen_K", K, "D");
			book("gen_eta", eta, "D");
			book("gen_phi", phi, "D");
			book("gen_pt", pt, "D");
		}
	private:
		double K,eta,phi,pt;

		virtual void reset() {
			K = -999.;
			eta = -999.;
			phi = -999.;
			pt = -999;
		}
	public:
		void fill(const reco::GenParticle& muon);
};

#endif
