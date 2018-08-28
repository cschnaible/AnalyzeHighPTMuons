// -*- C++ -*-
//
// Package:    AnalyzeTracks/AnalyzeTracks
// Class:      AnalyzeTracks
// 
/**\class AnalyzeTracks AnalyzeTracks.cc AnalyzeTracks/AnalyzeTracks/plugins/AnalyzeTracks.cc

 Description: Make ROOT TTree with RECO track parameters

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christian Schnaible (UCLA)
//         Created:  Fri, 09 Feb 2018 19:41:34 GMT
//
//


#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/OneToValue.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "HighPTMuons/Analysis/include/FillHighPTInfo.h"

//#include "TTree.h"
#include "TMath.h"

//
// class declaration
//

class AnalyzeTracks : public edm::EDAnalyzer {
   public:
      explicit AnalyzeTracks(const edm::ParameterSet&);
      ~AnalyzeTracks();

   typedef edm::AssociationMap<edm::OneToOne<std::vector<reco::Track>,std::vector<reco::Track> > > AssocTrackToTrack;

   private:
      virtual void beginJob() {};
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() {};
			template<typename M, typename V> int 
			matchMuonByDR(const M& mu, const V &v, const std::vector<int> &used) const;
			//int matchMuonByDR(
					//const reco::GenParticle&, const std::vector<reco::Track>&, const std::vector<int>&) const;
			int matchRecoToGenMuon(
					const reco::GenParticle&, const std::vector<reco::Muon>&, const std::vector<int>&) const;

      // ----------member data ---------------------------
	  
	  std::string filename;
		double correlation;
		double covScale;
	  
		edm::EDGetTokenT<std::vector<reco::Muon>> recoMuonToken_;

	  edm::EDGetTokenT<std::vector<reco::Track>> globalTrackToken_;
	  edm::EDGetTokenT<std::vector<reco::Track>> standAloneTrackToken_;
	  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> trackerOnlyMapToken_;

		edm::EDGetTokenT<AssocTrackToTrack> globalMuonOnlyMapToken_;
		edm::EDGetTokenT<AssocTrackToTrack> globalMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> pickyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> dytMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> globalTrackRankMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> globalTrackRankMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyTrackRankMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyTrackRankMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytTrackRankMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytTrackRankMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tunePTrackRankMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tunePTrackRankMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
	  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

		HighPTTreeContainer tree;
		FillEventInfo eventInfo;
		FillGenInfo genInfo;
		FillTrackInfo tunePInfo,   tunePMuonOnlyInfo,  tunePMuonOnlyUpdateInfo;
		FillTrackInfo globalInfo, globalMuonOnlyInfo, globalMuonOnlyUpdateInfo;
		FillTrackInfo pickyInfo, pickyMuonOnlyInfo, pickyMuonOnlyUpdateInfo;
		FillTrackInfo dytInfo, dytMuonOnlyInfo, dytMuonOnlyUpdateInfo;
		FillTrackInfo tpfmsInfo,  tpfmsMuonOnlyInfo,  tpfmsMuonOnlyUpdateInfo;
		FillTrackInfo globalTrackRankMuonOnlyInfo, globalTrackRankMuonOnlyUpdateInfo;
		FillTrackInfo pickyTrackRankMuonOnlyInfo, pickyTrackRankMuonOnlyUpdateInfo;
		FillTrackInfo dytTrackRankMuonOnlyInfo, dytTrackRankMuonOnlyUpdateInfo;
		FillTrackInfo tunePTrackRankMuonOnlyInfo, tunePTrackRankMuonOnlyUpdateInfo;
		FillTrackInfo standAloneInfo,   trackerInfo;
		FillCombinationInfo tunePCombInfo;
		FillCombinationInfo globalCombInfo;
		FillCombinationInfo pickyCombInfo;
		FillCombinationInfo dytCombInfo;
		FillCombinationInfo tpfmsCombInfo;
		FillCombinationInfo globalTrackRankCombInfo;
		FillCombinationInfo pickyTrackRankCombInfo;
		FillCombinationInfo dytTrackRankCombInfo;
		FillCombinationInfo tunePTrackRankCombInfo;
};

//
// constructors and destructor
//
AnalyzeTracks::AnalyzeTracks(const edm::ParameterSet& iConfig) : 
	tree("HighPTTrackTree","Tree holding HighPT Refit Tracks")
	// Turn on branch filling
	, eventInfo(tree)
	, genInfo(tree)
	, tunePInfo(tree,"tuneP")
	, tunePMuonOnlyInfo(tree,"tunePMuonOnly")
	, tunePMuonOnlyUpdateInfo(tree,"tunePMuonOnlyUpdate")
	, globalInfo(tree,"global")
	, globalMuonOnlyInfo(tree,"globalMuonOnly")
	, globalMuonOnlyUpdateInfo(tree,"globalMuonOnlyUpdate")
	, pickyInfo(tree,"picky")
	, pickyMuonOnlyInfo(tree,"pickyMuonOnly")
	, pickyMuonOnlyUpdateInfo(tree,"pickyMuonOnlyUpdate")
	, dytInfo(tree,"dyt")
	, dytMuonOnlyInfo(tree,"dytMuonOnly")
	, dytMuonOnlyUpdateInfo(tree,"dytMuonOnlyUpdate")
	, tpfmsInfo(tree,"tpfms")
	, tpfmsMuonOnlyInfo(tree,"tpfmsMuonOnly")
	, tpfmsMuonOnlyUpdateInfo(tree,"tpfmsMuonOnlyUpdate")
	, globalTrackRankMuonOnlyInfo(tree,"globalTrackRankMuonOnly")
	, globalTrackRankMuonOnlyUpdateInfo(tree,"globalTrackRankMuonOnlyUpdate")
	, pickyTrackRankMuonOnlyInfo(tree,"pickyTrackRankMuonOnly")
	, pickyTrackRankMuonOnlyUpdateInfo(tree,"pickyTrackRankMuonOnlyUpdate")
	, dytTrackRankMuonOnlyInfo(tree,"dytTrackRankMuonOnly")
	, dytTrackRankMuonOnlyUpdateInfo(tree,"dytTrackRankMuonOnlyUpdate")
	, tunePTrackRankMuonOnlyInfo(tree,"tunePTrackRankMuonOnly")
	, tunePTrackRankMuonOnlyUpdateInfo(tree,"tunePTrackRankMuonOnlyUpdate")
	, standAloneInfo(tree,"standAlone")
	, trackerInfo(tree,"tracker")
	, tunePCombInfo(tree,"tunePComb")
	, globalCombInfo(tree,"globalComb")
	, pickyCombInfo(tree,"pickyComb")
	, dytCombInfo(tree,"dytComb")
	, tpfmsCombInfo(tree,"tpfmsComb")
	, globalTrackRankCombInfo(tree,"globalTrackRankComb")
	, pickyTrackRankCombInfo(tree,"pickyTrackRankComb")
	, dytTrackRankCombInfo(tree,"dytTrackRankComb")
	, tunePTrackRankCombInfo(tree,"tunePTrackRankComb")
{
	//fileName = iConfig.getUntrackedParameter<std::string>("fileName");

	globalTrackToken_ = 
		consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("globalTracks"));
	standAloneTrackToken_ = 
		consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("standAloneTracks"));
	genParticlesToken_ = 
		consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));

	recoMuonToken_ = 
		consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("recoMuons"));

	trackerOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("trackerMapTag"));

	globalMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalMuonOnlyMapTag"));
	globalMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalMuonOnlyUpdateMapTag"));

	pickyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMapTag"));
	pickyMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMuonOnlyMapTag"));
	pickyMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMuonOnlyUpdateMapTag"));

	dytMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMapTag"));
	dytMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMuonOnlyMapTag"));
	dytMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMuonOnlyUpdateMapTag"));

	tpfmsMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMapTag"));
	tpfmsMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMuonOnlyMapTag"));
	tpfmsMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMuonOnlyUpdateMapTag"));

	globalTrackRankMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalTrackRankMuonOnlyMapTag"));
	globalTrackRankMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalTrackRankMuonOnlyUpdateMapTag"));
	pickyTrackRankMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyTrackRankMuonOnlyMapTag"));
	pickyTrackRankMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyTrackRankMuonOnlyUpdateMapTag"));
	dytTrackRankMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytTrackRankMuonOnlyMapTag"));
	dytTrackRankMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytTrackRankMuonOnlyUpdateMapTag"));
	tunePTrackRankMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tunePTrackRankMuonOnlyMapTag"));
	tunePTrackRankMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tunePTrackRankMuonOnlyUpdateMapTag"));

	beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"));
	vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

	correlation = iConfig.getParameter<double>("correlation");
	covScale = iConfig.getParameter<double>("covScale");

}


AnalyzeTracks::~AnalyzeTracks() {tree.write();}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnalyzeTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// Get the BeamSpot
	edm::Handle<reco::BeamSpot> beamspot;
	iEvent.getByToken(beamSpotToken_,beamspot);
	const reco::BeamSpot beamSpot = (*beamspot);

	// Get the Primary Vertex
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return; // skip the event if no PV is found
	const reco::Vertex &PV = vertices->front();

	// Fill the event info
	eventInfo.fill(iEvent,beamSpot,PV);

	// Reco muons
	edm::Handle< reco::MuonCollection > recoMuons;
	iEvent.getByToken(recoMuonToken_,recoMuons);

	// Global track
	edm::Handle< std::vector<reco::Track> > globalTracks;
	iEvent.getByToken(globalTrackToken_,globalTracks);

	// StandAlone track
	edm::Handle< std::vector<reco::Track> > standAloneTracks;
	iEvent.getByToken(standAloneTrackToken_,standAloneTracks);

	// Gen Particles
	edm::Handle< std::vector<reco::GenParticle> > genParticles;
	iEvent.getByToken(genParticlesToken_,genParticles);
	
	// Tracker-only
	edm::Handle<AssocTrackToTrack> trackerOnlyTracks;
	iEvent.getByToken(trackerOnlyMapToken_,trackerOnlyTracks);

	// Global No RPC track
	edm::Handle<AssocTrackToTrack> globalMuonOnlyTracks;
	iEvent.getByToken(globalMuonOnlyMapToken_,globalMuonOnlyTracks);
	// Global muon-only with update track
	edm::Handle<AssocTrackToTrack> globalMuonOnlyUpdateTracks;
	iEvent.getByToken(globalMuonOnlyUpdateMapToken_,globalMuonOnlyUpdateTracks);
	
	// picky track
	edm::Handle<AssocTrackToTrack> pickyTracks;
	iEvent.getByToken(pickyMapToken_,pickyTracks);
	// picky muon-only track
	edm::Handle<AssocTrackToTrack> pickyMuonOnlyTracks;
	iEvent.getByToken(pickyMuonOnlyMapToken_,pickyMuonOnlyTracks);
	// picky muon-only with update track
	edm::Handle<AssocTrackToTrack> pickyMuonOnlyUpdateTracks;
	iEvent.getByToken(pickyMuonOnlyUpdateMapToken_,pickyMuonOnlyUpdateTracks);

	// dyt track
	edm::Handle<AssocTrackToTrack> dytTracks;
	iEvent.getByToken(dytMapToken_,dytTracks);
	// dyt muon-only track
	edm::Handle<AssocTrackToTrack> dytMuonOnlyTracks;
	iEvent.getByToken(dytMuonOnlyMapToken_,dytMuonOnlyTracks);
	// dyt muon-only with update track
	edm::Handle<AssocTrackToTrack> dytMuonOnlyUpdateTracks;
	iEvent.getByToken(dytMuonOnlyUpdateMapToken_,dytMuonOnlyUpdateTracks);

	// tpfms track
	edm::Handle<AssocTrackToTrack> tpfmsTracks;
	iEvent.getByToken(tpfmsMapToken_,tpfmsTracks);
	// tpfms muon-only track
	edm::Handle<AssocTrackToTrack> tpfmsMuonOnlyTracks;
	iEvent.getByToken(tpfmsMuonOnlyMapToken_,tpfmsMuonOnlyTracks);
	// tpfms muon-only with update track
	edm::Handle<AssocTrackToTrack> tpfmsMuonOnlyUpdateTracks;
	iEvent.getByToken(tpfmsMuonOnlyUpdateMapToken_,tpfmsMuonOnlyUpdateTracks);

	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> globalTrackRankMuonOnlyTracks;
	iEvent.getByToken(globalTrackRankMuonOnlyMapToken_,globalTrackRankMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> globalTrackRankMuonOnlyUpdateTracks;
	iEvent.getByToken(globalTrackRankMuonOnlyUpdateMapToken_,globalTrackRankMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> pickyTrackRankMuonOnlyTracks;
	iEvent.getByToken(pickyTrackRankMuonOnlyMapToken_,pickyTrackRankMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> pickyTrackRankMuonOnlyUpdateTracks;
	iEvent.getByToken(pickyTrackRankMuonOnlyUpdateMapToken_,pickyTrackRankMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> dytTrackRankMuonOnlyTracks;
	iEvent.getByToken(dytTrackRankMuonOnlyMapToken_,dytTrackRankMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> dytTrackRankMuonOnlyUpdateTracks;
	iEvent.getByToken(dytTrackRankMuonOnlyUpdateMapToken_,dytTrackRankMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> tunePTrackRankMuonOnlyTracks;
	iEvent.getByToken(tunePTrackRankMuonOnlyMapToken_,tunePTrackRankMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> tunePTrackRankMuonOnlyUpdateTracks;
	iEvent.getByToken(tunePTrackRankMuonOnlyUpdateMapToken_,tunePTrackRankMuonOnlyUpdateTracks);


	// list of matched global/SA track indices
	std::vector<int> globalUsed = {-1};
	std::vector<int> recoMuonUsed= {-1};
	std::vector<int> standAloneUsed = {-1};



	// Loop on gen muons and match a global track
	// From the global track each additional refit
	// can be found via an AssociationMap
	//
	// Each entry in output tree is a single muon
	for (unsigned int g=0; g < (*genParticles).size(); g++) {

		// muon and status 1
		if (fabs((*genParticles)[g].pdgId())!=13) continue;
		if ((*genParticles)[g].status()!=1) continue;
		auto genMu = (*genParticles)[g];
		// Fill Gen Muon info
		genInfo.fill(genMu);
		
		// Match global track to gen muon
		int glb = matchMuonByDR(genMu,*globalTracks,globalUsed);
		bool doComb = false;
		bool doCombGlobal = false;
		bool doCombPicky = false;
		bool doCombDYT = false;
		bool doCombTPFMS = false;
		bool doCombGlobalTrackRank= false;
		bool doCombPickyTrackRank= false;
		bool doCombDYTTrackRank= false;
		bool doCombTunePTrackRank= false;
		std::string tunePchoice;
		// If global track is found then fill additional global track refits
		if (glb>=0) {

		  globalUsed.push_back(glb);
			globalInfo.fill((*globalTracks)[glb],beamSpot,PV);
			// tuneP
			int reco = matchMuonByDR((*globalTracks)[glb],*recoMuons,recoMuonUsed);
			if (reco>=0){
				recoMuonUsed.push_back(reco);
				if ((*recoMuons)[reco].muonBestTrackType()==1) {
					std::cout << "skip when tuneP is tracker track" << std::endl;
				}
				else {
					tunePInfo.fill(*(*recoMuons)[reco].bestTrack(),beamSpot,PV);
				}
			}

			// ------------------------------------------------------------------------ //
		  // Match tracker track to global track
		  //AssocTrackToTrack const & trackerMap = *trackerOnlyTracks;
		  try {
				trackerInfo.fill(*(*recoMuons)[reco].innerTrack(),beamSpot,PV);
				doComb = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tracker track with update to tracker track" << std::endl;
		  }

		  // Match Mu-only global track to global track
		  AssocTrackToTrack const & globalMuonOnlyMap = *globalMuonOnlyTracks;
		  try {
		    auto globalMuonOnlyTrack = 
					*globalMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				globalMuonOnlyInfo.fill(globalMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only global track glb track?" << std::endl;
		  }

		  // Match Mu-only global track with update to global track
		  AssocTrackToTrack const & globalMuonOnlyUpdateMap = *globalMuonOnlyUpdateTracks;
		  try {
			  auto globalMuonOnlyUpdateTrack = 
					*globalMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				globalMuonOnlyUpdateInfo.fill(globalMuonOnlyUpdateTrack,beamSpot,PV,covScale);
		    doCombGlobal = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only global track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match picky track to global track
		  AssocTrackToTrack const & pickyMap = *pickyTracks;
		  try {
		    auto pickyTrack = 
					*pickyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyInfo.fill(pickyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match picky track to glb track?" << std::endl;
		  }
		  
		  // Match Mu-only picky track to global track
		  AssocTrackToTrack const & pickyMuonOnlyMap = *pickyMuonOnlyTracks;
		  try {
		    auto pickyMuonOnlyTrack = 
					*pickyMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyMuonOnlyInfo.fill(pickyMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only picky track to glb track?" << std::endl;
		  }

		  // Match Mu-only picky track with update to picky track
		  AssocTrackToTrack const & pickyMuonOnlyUpdateMap = *pickyMuonOnlyUpdateTracks;
		  try {
			  auto pickyMuonOnlyUpdateTrack = 
					*pickyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyMuonOnlyUpdateInfo.fill(pickyMuonOnlyUpdateTrack,beamSpot,PV,covScale);
				doCombPicky = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only picky track with update to picky track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match dyt track to global track
		  AssocTrackToTrack const & dytMap = *dytTracks;
		  try {
		    auto dytTrack = 
					*dytMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytInfo.fill(dytTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match dyt track to glb track?" << std::endl;
		  }

		  // Match Mu-only dyt track to global track
		  AssocTrackToTrack const & dytMuonOnlyMap = *dytMuonOnlyTracks;
		  try {
		    auto dytMuonOnlyTrack = 
					*dytMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytMuonOnlyInfo.fill(dytMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only dyt track to glb track?" << std::endl;
		  }

		  // Match Mu-only dyt track with update to dyt track
		  AssocTrackToTrack const & dytMuonOnlyUpdateMap = *dytMuonOnlyUpdateTracks;
		  try {
			  auto dytMuonOnlyUpdateTrack = 
					*dytMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytMuonOnlyUpdateInfo.fill(dytMuonOnlyUpdateTrack,beamSpot,PV,covScale);
				doCombDYT = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only dyt track with update to dyt track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match tpfms track to global track
		  AssocTrackToTrack const & tpfmsMap = *tpfmsTracks;
		  try {
		    auto tpfmsTrack = 
					*tpfmsMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tpfmsInfo.fill(tpfmsTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match tpfms track to glb track?" << std::endl;
		  }
			
		  // Match Mu-only tpfms track to global track
		  AssocTrackToTrack const & tpfmsMuonOnlyMap = *tpfmsMuonOnlyTracks;
		  try {
		    auto tpfmsMuonOnlyTrack = 
					*tpfmsMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tpfmsMuonOnlyInfo.fill(tpfmsMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only tpfms track to glb track?" << std::endl;
		  }

		  // Match Mu-only tpfms track with update to tpfms track
		  AssocTrackToTrack const & tpfmsMuonOnlyUpdateMap = *tpfmsMuonOnlyUpdateTracks;
		  try {
			  auto tpfmsMuonOnlyUpdateTrack = 
					*tpfmsMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tpfmsMuonOnlyUpdateInfo.fill(tpfmsMuonOnlyUpdateTrack,beamSpot,PV,covScale);
				doCombTPFMS = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tpfms track with update to tpfms track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match Mu-only globalTrackRank track to global track
		  AssocTrackToTrack const & globalTrackRankMuonOnlyMap = *globalTrackRankMuonOnlyTracks;
		  try {
		    auto globalTrackRankMuonOnlyTrack = 
					*globalTrackRankMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				globalTrackRankMuonOnlyInfo.fill(globalTrackRankMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only globalTrackRank track to glb track?" << std::endl;
		  }

		  // Match Mu-only globalTrackRank track with update to globalTrackRank track
		  AssocTrackToTrack const & globalTrackRankMuonOnlyUpdateMap = *globalTrackRankMuonOnlyUpdateTracks;
		  try {
			  auto globalTrackRankMuonOnlyUpdateTrack = 
					*globalTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				globalTrackRankMuonOnlyUpdateInfo.fill(globalTrackRankMuonOnlyUpdateTrack,beamSpot,PV,covScale);
		    doCombGlobalTrackRank= true;
		  } catch(...) {
			  std::cout << "Could not match muon-only globalTrackRank track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match Mu-only pickyTrackRank track to global track
		  AssocTrackToTrack const & pickyTrackRankMuonOnlyMap = *pickyTrackRankMuonOnlyTracks;
		  try {
		    auto pickyTrackRankMuonOnlyTrack = 
					*pickyTrackRankMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyTrackRankMuonOnlyInfo.fill(pickyTrackRankMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only pickyTrackRank track to glb track?" << std::endl;
		  }

		  // Match Mu-only pickyTrackRank track with update to pickyTrackRank track
		  AssocTrackToTrack const & pickyTrackRankMuonOnlyUpdateMap = *pickyTrackRankMuonOnlyUpdateTracks;
		  try {
			  auto pickyTrackRankMuonOnlyUpdateTrack = 
					*pickyTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyTrackRankMuonOnlyUpdateInfo.fill(pickyTrackRankMuonOnlyUpdateTrack,beamSpot,PV,covScale);
		    doCombPickyTrackRank= true;
		  } catch(...) {
			  std::cout << "Could not match muon-only pickyTrackRank track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match Mu-only dytTrackRank track to global track
		  AssocTrackToTrack const & dytTrackRankMuonOnlyMap = *dytTrackRankMuonOnlyTracks;
		  try {
		    auto dytTrackRankMuonOnlyTrack = 
					*dytTrackRankMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytTrackRankMuonOnlyInfo.fill(dytTrackRankMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only dytTrackRank track to glb track?" << std::endl;
		  }

		  // Match Mu-only dytTrackRank track with update to dytTrackRank track
		  AssocTrackToTrack const & dytTrackRankMuonOnlyUpdateMap = *dytTrackRankMuonOnlyUpdateTracks;
		  try {
			  auto dytTrackRankMuonOnlyUpdateTrack = 
					*dytTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytTrackRankMuonOnlyUpdateInfo.fill(dytTrackRankMuonOnlyUpdateTrack,beamSpot,PV,covScale);
		    doCombDYTTrackRank= true;
		  } catch(...) {
			  std::cout << "Could not match muon-only dytTrackRank track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match Mu-only tunePTrackRank track to global track
		  AssocTrackToTrack const & tunePTrackRankMuonOnlyMap = *tunePTrackRankMuonOnlyTracks;
		  try {
		    auto tunePTrackRankMuonOnlyTrack = 
					*tunePTrackRankMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tunePTrackRankMuonOnlyInfo.fill(tunePTrackRankMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only tunePTrackRank track to glb track?" << std::endl;
		  }

		  // Match Mu-only tunePTrackRank track with update to tunePTrackRank track
		  AssocTrackToTrack const & tunePTrackRankMuonOnlyUpdateMap = *tunePTrackRankMuonOnlyUpdateTracks;
		  try {
			  auto tunePTrackRankMuonOnlyUpdateTrack = 
					*tunePTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tunePTrackRankMuonOnlyUpdateInfo.fill(tunePTrackRankMuonOnlyUpdateTrack,beamSpot,PV,covScale);
		    doCombTunePTrackRank= true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tunePTrackRank track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //


			// If a muon-only updated track refit is found then store combination with tracker track
			if (doComb){
				auto trackerTrack = *(*recoMuons)[reco].innerTrack();
				if (doCombGlobal) {
					auto globalMuonOnlyUpdateTrack = 
						*globalMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					globalCombInfo.fill(trackerTrack, globalMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombPicky) {
					auto pickyMuonOnlyUpdateTrack = 
						*pickyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					pickyCombInfo.fill(trackerTrack, pickyMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombDYT) {
					auto dytMuonOnlyUpdateTrack = 
						*dytMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					dytCombInfo.fill(trackerTrack, dytMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombTPFMS) {
					auto tpfmsMuonOnlyUpdateTrack = 
						*tpfmsMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tpfmsCombInfo.fill(trackerTrack, tpfmsMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombGlobalTrackRank) {
					auto globalTrackRankMuonOnlyUpdateTrack = 
						*globalTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					globalTrackRankCombInfo.fill(trackerTrack, globalTrackRankMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombPickyTrackRank) {
					auto pickyTrackRankMuonOnlyUpdateTrack = 
						*pickyTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					pickyTrackRankCombInfo.fill(trackerTrack, pickyTrackRankMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombDYTTrackRank) {
					auto dytTrackRankMuonOnlyUpdateTrack = 
						*dytTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					dytTrackRankCombInfo.fill(trackerTrack, dytTrackRankMuonOnlyUpdateTrack,covScale,correlation);
				}
				if (doCombTunePTrackRank) {
					auto tunePTrackRankMuonOnlyUpdateTrack = 
						*tunePTrackRankMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePTrackRankCombInfo.fill(trackerTrack, tunePTrackRankMuonOnlyUpdateTrack,covScale,correlation);
				}

				// 
				// tuneP
				//
				if ((*recoMuons)[reco].muonBestTrackType()==0)
					std::cout << "Why is there no best track choice?" << std::endl;
				else if ((*recoMuons)[reco].muonBestTrackType()==1) {
					//tunePCombInfo.fill(trackerTrack,trackerTrack);
					std::cout << " skipping tracker track tunePComb... " << std::endl;
				}
				else if ((*recoMuons)[reco].muonBestTrackType()==2)
					std::cout << "StandAlone track not a valid tuneP choice..." << std::endl;
				else if ((*recoMuons)[reco].muonBestTrackType()==3 && doCombGlobal) {
					auto globalMuonOnlyUpdateTrack = 
						*globalMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					auto globalMuonOnlyTrack = 
						*globalMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePMuonOnlyInfo.fill(globalMuonOnlyTrack,beamSpot,PV);
					tunePMuonOnlyUpdateInfo.fill(globalMuonOnlyUpdateTrack,beamSpot,PV,covScale);
					tunePCombInfo.fill(trackerTrack,globalMuonOnlyUpdateTrack,covScale,correlation);
				}
				else if ((*recoMuons)[reco].muonBestTrackType()==4 && doCombTPFMS) {
					auto tpfmsMuonOnlyUpdateTrack = 
						*tpfmsMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					auto tpfmsMuonOnlyTrack = 
						*tpfmsMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePMuonOnlyInfo.fill(tpfmsMuonOnlyTrack,beamSpot,PV);
					tunePMuonOnlyUpdateInfo.fill(tpfmsMuonOnlyUpdateTrack,beamSpot,PV,covScale);
					tunePCombInfo.fill(trackerTrack,tpfmsMuonOnlyUpdateTrack,covScale,correlation);
				}
				else if ((*recoMuons)[reco].muonBestTrackType()==5 && doCombPicky) {
					auto pickyMuonOnlyUpdateTrack = 
						*pickyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					auto pickyMuonOnlyTrack = 
						*pickyMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePMuonOnlyInfo.fill(pickyMuonOnlyTrack,beamSpot,PV);
					tunePMuonOnlyUpdateInfo.fill(pickyMuonOnlyUpdateTrack,beamSpot,PV,covScale);
					tunePCombInfo.fill(trackerTrack,pickyMuonOnlyUpdateTrack,covScale,correlation);
				}
				else if ((*recoMuons)[reco].muonBestTrackType()==6 && doCombDYT) {
					auto dytMuonOnlyUpdateTrack = 
						*dytMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					auto dytMuonOnlyTrack = 
						*dytMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePMuonOnlyInfo.fill(dytMuonOnlyTrack,beamSpot,PV);
					tunePMuonOnlyUpdateInfo.fill(dytMuonOnlyUpdateTrack,beamSpot,PV,covScale);
					tunePCombInfo.fill(trackerTrack,dytMuonOnlyUpdateTrack,covScale,correlation);
				}
				else {
					std::cout << "tuneP choice " 
						<< (*recoMuons)[reco].muonBestTrackType() 
						<< " did not have a valid refit" << std::endl;
				}

			}

		} // end if glb>=0
		else {
		  std::cout << "Could not match Global track to gen muon?" << std::endl;
		} // end global muon

		// Stand alone muons
		int sa = matchMuonByDR(genMu,*standAloneTracks,standAloneUsed);
		if (sa>=0) {
		  standAloneUsed.push_back(sa);
		  auto standAloneTrack = (*standAloneTracks)[sa];
		  standAloneInfo.fill(standAloneTrack,beamSpot,PV);
		} else {
		  std::cout << "Could not match to SA track to gen muon?" << std::endl;
		}

		tree.fill();

	 //std::cout << " ------ " << std::endl;
   } // end loop on gen muons
   globalUsed.clear();
   standAloneUsed.clear();
	 recoMuonUsed.clear();
	 //std::cout << " ****** " << std::endl;
	 //std::cout << " ****** " << std::endl;

   
}

template<class M, class V> int
AnalyzeTracks::matchMuonByDR(const M &mu, const V &tracks, const std::vector<int> &used) const
{
	int i = 0;
	for (auto track : tracks) {

		double this_dR = deltaR(mu,track);
		//std::cout << i << " " << this_dR << std::endl;
		/*
		if (std::find(used.begin(),used.end(),i) != used.end()) {
			i++;
			continue;
		}
		else*/ if (this_dR<0.5) {
			return i;
		}
		else {
			i++;
		}
	}
	return -1;
}

int
AnalyzeTracks::matchRecoToGenMuon(const reco::GenParticle &genMu, const std::vector<reco::Muon> &muons, const std::vector<int> &used) const
{
	//for (auto u : used) std::cout << u << " ";
	//std::cout << std::endl;
	int key = -1;
	for (const auto &muon : muons) {

		double this_dR = deltaR(genMu,muon);
		key = muon.combinedMuon().key();
		//std::cout << i << " " << this_dR << std::endl;
		if (std::find(used.begin(),used.end(),key) != used.end()) {
			continue;
		}
		if (this_dR<0.5) {
			return key;
		}
	}
	return -1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeTracks);
			// ------------------------------------------------------------------------ //
	/*
	, bestCombInfo(tree,"bestComb")
	, trkCurvCombInfo(tree,"trkCurvComb")
	, tunePCurvCombInfo(tree,"tunePCurvComb")
	, dxyCombInfo(tree,"dxyComb")
	, curvRelErrCombInfo(tree,"curvRelErrComb")
	*/
	/*
	, bestMuonOnlyInfo(tree,"bestMuonOnly")
	, bestMuonOnlyUpdateInfo(tree,"bestMuonOnlyUpdate")
	, trkCurvMuonOnlyInfo(tree,"trkCurvMuonOnly")
	, trkCurvMuonOnlyUpdateInfo(tree,"trkCurvMuonOnlyUpdate")
	, tunePCurvMuonOnlyInfo(tree,"tunePCurvMuonOnly")
	, tunePCurvMuonOnlyUpdateInfo(tree,"tunePCurvMuonOnlyUpdate")
	, dxyMuonOnlyInfo(tree,"dxyMuonOnly")
	, dxyMuonOnlyUpdateInfo(tree,"dxyMuonOnlyUpdate")
	, curvRelErrMuonOnlyInfo(tree,"curvRelErrMuonOnly")
	, curvRelErrMuonOnlyUpdateInfo(tree,"curvRelErrMuonOnlyUpdate")
	*/
		/*
		//FillCombinationInfo bestCombInfo;
		FillCombinationInfo trkCurvCombInfo;
		FillCombinationInfo tunePCurvCombInfo;
		FillCombinationInfo dxyCombInfo;
		FillCombinationInfo curvRelErrCombInfo;
		*/
		/*
		//FillTrackInfo bestMuonOnlyInfo, bestMuonOnlyUpdateInfo;
		FillTrackInfo trkCurvMuonOnlyInfo, trkCurvMuonOnlyUpdateInfo;
		FillTrackInfo tunePCurvMuonOnlyInfo, tunePCurvMuonOnlyUpdateInfo;
		FillTrackInfo dxyMuonOnlyInfo, dxyMuonOnlyUpdateInfo;
		FillTrackInfo curvRelErrMuonOnlyInfo, curvRelErrMuonOnlyUpdateInfo;
		*/
		/*
	  edm::EDGetTokenT<AssocTrackToTrack> bestMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> bestMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> trkCurvMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> trkCurvMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tunePCurvMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tunePCurvMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dxyMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dxyMuonOnlyUpdateMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> curvRelErrMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> curvRelErrMuonOnlyUpdateMapToken_;
		*/
	/*
	bestMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("bestMuonOnlyMapTag"));
	bestMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("bestMuonOnlyUpdateMapTag"));
	trkCurvMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("trkCurvMuonOnlyMapTag"));
	trkCurvMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("trkCurvMuonOnlyUpdateMapTag"));
	tunePCurvMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tunePCurvMuonOnlyMapTag"));
	tunePCurvMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tunePCurvMuonOnlyUpdateMapTag"));
	dxyMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dxyMuonOnlyMapTag"));
	dxyMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dxyMuonOnlyUpdateMapTag"));
	curvRelErrMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("curvRelErrMuonOnlyMapTag"));
	curvRelErrMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("curvRelErrMuonOnlyUpdateMapTag"));
	*/
	/*
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> bestMuonOnlyTracks;
	iEvent.getByToken(bestMuonOnlyMapToken_,bestMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> bestMuonOnlyUpdateTracks;
	iEvent.getByToken(bestMuonOnlyUpdateMapToken_,bestMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> trkCurvMuonOnlyTracks;
	iEvent.getByToken(trkCurvMuonOnlyMapToken_,trkCurvMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> trkCurvMuonOnlyUpdateTracks;
	iEvent.getByToken(trkCurvMuonOnlyUpdateMapToken_,trkCurvMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> tunePCurvMuonOnlyTracks;
	iEvent.getByToken(tunePCurvMuonOnlyMapToken_,tunePCurvMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> tunePCurvMuonOnlyUpdateTracks;
	iEvent.getByToken(tunePCurvMuonOnlyUpdateMapToken_,tunePCurvMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> dxyMuonOnlyTracks;
	iEvent.getByToken(dxyMuonOnlyMapToken_,dxyMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> dxyMuonOnlyUpdateTracks;
	iEvent.getByToken(dxyMuonOnlyUpdateMapToken_,dxyMuonOnlyUpdateTracks);
	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> curvRelErrMuonOnlyTracks;
	iEvent.getByToken(curvRelErrMuonOnlyMapToken_,curvRelErrMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> curvRelErrMuonOnlyUpdateTracks;
	iEvent.getByToken(curvRelErrMuonOnlyUpdateMapToken_,curvRelErrMuonOnlyUpdateTracks);
	*/

			/*
		  // Match Mu-only best track to global track
		  AssocTrackToTrack const & bestMuonOnlyMap = *bestMuonOnlyTracks;
		  try {
		    auto bestMuonOnlyTrack = 
					*bestMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				bestMuonOnlyInfo.fill(bestMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only best track to glb track?" << std::endl;
		  }

		  // Match Mu-only best track with update to best track
		  AssocTrackToTrack const & bestMuonOnlyUpdateMap = *bestMuonOnlyUpdateTracks;
		  try {
			  auto bestMuonOnlyUpdateTrack = 
					*bestMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				bestMuonOnlyUpdateInfo.fill(bestMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombBest = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only best track with update to global track" << std::endl;
		  }

			// ------------------------------------------------------------------------ //
		  // Match Mu-only trkCurv track to global track
		  AssocTrackToTrack const & trkCurvMuonOnlyMap = *trkCurvMuonOnlyTracks;
		  try {
		    auto trkCurvMuonOnlyTrack = 
					*trkCurvMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				trkCurvMuonOnlyInfo.fill(trkCurvMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only trkCurv track to glb track?" << std::endl;
		  }

		  // Match Mu-only trkCurv track with update to trkCurv track
		  AssocTrackToTrack const & trkCurvMuonOnlyUpdateMap = *trkCurvMuonOnlyUpdateTracks;
		  try {
			  auto trkCurvMuonOnlyUpdateTrack = 
					*trkCurvMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				trkCurvMuonOnlyUpdateInfo.fill(trkCurvMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombTrkCurv = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only trkCurv track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //

			// ------------------------------------------------------------------------ //
		  // Match Mu-only tunePCurv track to global track
		  AssocTrackToTrack const & tunePCurvMuonOnlyMap = *tunePCurvMuonOnlyTracks;
		  try {
		    auto tunePCurvMuonOnlyTrack = 
					*tunePCurvMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tunePCurvMuonOnlyInfo.fill(tunePCurvMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only tunePCurv track to glb track?" << std::endl;
		  }

		  // Match Mu-only tunePCurv track with update to tunePCurv track
		  AssocTrackToTrack const & tunePCurvMuonOnlyUpdateMap = *tunePCurvMuonOnlyUpdateTracks;
		  try {
			  auto tunePCurvMuonOnlyUpdateTrack = 
					*tunePCurvMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tunePCurvMuonOnlyUpdateInfo.fill(tunePCurvMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombTunePCurv = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tunePCurv track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //


			// ------------------------------------------------------------------------ //
		  // Match Mu-only dxy track to global track
		  AssocTrackToTrack const & dxyMuonOnlyMap = *dxyMuonOnlyTracks;
		  try {
		    auto dxyMuonOnlyTrack = 
					*dxyMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dxyMuonOnlyInfo.fill(dxyMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only dxy track to glb track?" << std::endl;
		  }

		  // Match Mu-only dxy track with update to dxy track
		  AssocTrackToTrack const & dxyMuonOnlyUpdateMap = *dxyMuonOnlyUpdateTracks;
		  try {
			  auto dxyMuonOnlyUpdateTrack = 
					*dxyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dxyMuonOnlyUpdateInfo.fill(dxyMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombDxy = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only dxy track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //


			// ------------------------------------------------------------------------ //
		  // Match Mu-only curvRelErr track to global track
		  AssocTrackToTrack const & curvRelErrMuonOnlyMap = *curvRelErrMuonOnlyTracks;
		  try {
		    auto curvRelErrMuonOnlyTrack = 
					*curvRelErrMuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				curvRelErrMuonOnlyInfo.fill(curvRelErrMuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only curvRelErr track to glb track?" << std::endl;
		  }

		  // Match Mu-only curvRelErr track with update to curvRelErr track
		  AssocTrackToTrack const & curvRelErrMuonOnlyUpdateMap = *curvRelErrMuonOnlyUpdateTracks;
		  try {
			  auto curvRelErrMuonOnlyUpdateTrack = 
					*curvRelErrMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				curvRelErrMuonOnlyUpdateInfo.fill(curvRelErrMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombCurvRelErr = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only curvRelErr track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //
			*/
				/*
				if (doCombBest) {
					auto bestMuonOnlyUpdateTrack = 
						*bestMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					bestCombInfo.fill(trackerTrack, bestMuonOnlyUpdateTrack);
				}
				if (doCombTrkCurv) {
					auto trkCurvMuonOnlyUpdateTrack = 
						*trkCurvMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					trkCurvCombInfo.fill(trackerTrack, trkCurvMuonOnlyUpdateTrack);
				}
				if (doCombTunePCurv) {
					auto tunePCurvMuonOnlyUpdateTrack = 
						*tunePCurvMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tunePCurvCombInfo.fill(trackerTrack, tunePCurvMuonOnlyUpdateTrack);
				}
				if (doCombDxy) {
					auto dxyMuonOnlyUpdateTrack = 
						*dxyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					dxyCombInfo.fill(trackerTrack, dxyMuonOnlyUpdateTrack);
				}
				if (doCombCurvRelErr) {
					auto curvRelErrMuonOnlyUpdateTrack = 
						*curvRelErrMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					curvRelErrCombInfo.fill(trackerTrack, curvRelErrMuonOnlyUpdateTrack);
				}
				*/
			/*
			// ------------------------------------------------------------------------ //
		  // Match Mu-only trackRank track to global track
			// f = 3
		  AssocTrackToTrack const & trackRank3MuonOnlyMap = *trackRank3MuonOnlyTracks;
		  try {
		    auto trackRank3MuonOnlyTrack = 
					*trackRank3MuonOnlyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				trackRank3MuonOnlyInfo.fill(trackRank3MuonOnlyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match muon-only trackRank3 track to glb track?" << std::endl;
		  }

		  // Match Mu-only trackRank3 track with update to trackRank3 track
		  AssocTrackToTrack const & trackRank3MuonOnlyUpdateMap = *trackRank3MuonOnlyUpdateTracks;
		  try {
			  auto trackRank3MuonOnlyUpdateTrack = 
					*trackRank3MuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				trackRank3MuonOnlyUpdateInfo.fill(trackRank3MuonOnlyUpdateTrack,beamSpot,PV);
		    doCombTrackRank3= true;
		  } catch(...) {
			  std::cout << "Could not match muon-only trackRank3 track with update to global track" << std::endl;
		  }
			// ------------------------------------------------------------------------ //
			*/
				/*
				if (doCombTrackRank3) {
					auto trackRank3MuonOnlyUpdateTrack = 
						*trackRank3MuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					trackRank3CombInfo.fill(trackerTrack, trackRank3MuonOnlyUpdateTrack);
				}
				*/
