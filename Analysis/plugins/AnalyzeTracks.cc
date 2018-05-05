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
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "TTree.h"
#include "TMath.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/OneToValue.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "HighPTMuons/Analysis/include/FillHighPTInfo.h"

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
			int matchToGenMuon(
					const reco::GenParticle&, const std::vector<reco::Track>&, const std::vector<int>&) const;

      // ----------member data ---------------------------
	  
	  std::string filename;
	  
	  edm::EDGetTokenT<std::vector<reco::Track>> globalTrackToken_;
	  edm::EDGetTokenT<std::vector<reco::Track>> standAloneTrackToken_;
	  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> trackerOnlyMapToken_;

		edm::EDGetTokenT<AssocTrackToTrack> globalNoRPCToken_;
		edm::EDGetTokenT<AssocTrackToTrack> globalMuonOnlyMapToken_;
		edm::EDGetTokenT<AssocTrackToTrack> globalMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> pickyMapToken_;
		edm::EDGetTokenT<AssocTrackToTrack> pickyNoRPCToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> pickyMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> dytMapToken_;
		edm::EDGetTokenT<AssocTrackToTrack> dytNoRPCToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> dytMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMapToken_;
		edm::EDGetTokenT<AssocTrackToTrack> tpfmsNoRPCToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> tpfmsMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<AssocTrackToTrack> bestMuonOnlyMapToken_;
	  edm::EDGetTokenT<AssocTrackToTrack> bestMuonOnlyUpdateMapToken_;

	  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
	  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;

		HighPTTreeContainer tree;
		FillEventInfo eventInfo;
		FillGenInfo genInfo;
		FillTrackInfo globalInfo, globalMuonOnlyInfo, globalMuonOnlyUpdateInfo;
		FillTrackInfo pickyInfo,  pickyMuonOnlyInfo,  pickyMuonOnlyUpdateInfo;
		FillTrackInfo dytInfo,    dytMuonOnlyInfo,    dytMuonOnlyUpdateInfo;
		FillTrackInfo tpfmsInfo,  tpfmsMuonOnlyInfo,  tpfmsMuonOnlyUpdateInfo;
		FillTrackInfo globalNoRPCInfo, pickyNoRPCInfo, dytNoRPCInfo, tpfmsNoRPCInfo;
		FillTrackInfo bestMuonOnlyInfo, bestMuonOnlyUpdateInfo;
		FillTrackInfo standAloneInfo,   trackerInfo;
		FillCombinationInfo globalMuonOnlyUpdateTrackerCombInfo;
		FillCombinationInfo pickyMuonOnlyUpdateTrackerCombInfo;
		FillCombinationInfo dytMuonOnlyUpdateTrackerCombInfo;
		FillCombinationInfo tpfmsMuonOnlyUpdateTrackerCombInfo;
		FillCombinationInfo bestMuonOnlyUpdateTrackerCombInfo;
};

//
// constructors and destructor
//
AnalyzeTracks::AnalyzeTracks(const edm::ParameterSet& iConfig) : 
	tree("HighPTTrackTree","Tree holding HighPT Refit Tracks")
	// Turn on branch filling
	, eventInfo(tree)
	, genInfo(tree)
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
	, globalNoRPCInfo(tree,"globalNoRPC")
	, pickyNoRPCInfo(tree,"pickyNoRPC")
	, dytNoRPCInfo(tree,"dytNoRPC")
	, tpfmsNoRPCInfo(tree,"tpfmsNoRPC")
	, bestMuonOnlyInfo(tree,"bestMuonOnly")
	, bestMuonOnlyUpdateInfo(tree,"bestMuonOnlyUpdate")
	, standAloneInfo(tree,"standAlone")
	, trackerInfo(tree,"tracker")
	, globalMuonOnlyUpdateTrackerCombInfo(tree,"globalComb")
	, pickyMuonOnlyUpdateTrackerCombInfo(tree,"pickyComb")
	, dytMuonOnlyUpdateTrackerCombInfo(tree,"dytComb")
	, tpfmsMuonOnlyUpdateTrackerCombInfo(tree,"tpfmsComb")
	, bestMuonOnlyUpdateTrackerCombInfo(tree,"bestComb")
{
	//fileName = iConfig.getUntrackedParameter<std::string>("fileName");

	globalTrackToken_ = 
		consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("globalTracks"));
	standAloneTrackToken_ = 
		consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("standAloneTracks"));
	genParticlesToken_ = 
		consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));

	trackerOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("trackerMapTag"));

	globalNoRPCToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalNoRPCMapTag"));
	globalMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalMuonOnlyMapTag"));
	globalMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("globalMuonOnlyUpdateMapTag"));

	pickyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMapTag"));
	pickyNoRPCToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyNoRPCMapTag"));
	pickyMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMuonOnlyMapTag"));
	pickyMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("pickyMuonOnlyUpdateMapTag"));

	dytMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMapTag"));
	dytNoRPCToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytNoRPCMapTag"));
	dytMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMuonOnlyMapTag"));
	dytMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("dytMuonOnlyUpdateMapTag"));

	tpfmsMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMapTag"));
	tpfmsNoRPCToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsNoRPCMapTag"));
	tpfmsMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMuonOnlyMapTag"));
	tpfmsMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("tpfmsMuonOnlyUpdateMapTag"));

	bestMuonOnlyMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("bestMuonOnlyMapTag"));
	bestMuonOnlyUpdateMapToken_ = 
		consumes<AssocTrackToTrack>(iConfig.getParameter<edm::InputTag>("bestMuonOnlyUpdateMapTag"));

	beamSpotToken_ = consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot"));
	vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

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
	edm::Handle<AssocTrackToTrack> globalNoRPCTracks;
	iEvent.getByToken(globalNoRPCToken_,globalNoRPCTracks);
	// Global muon-only track
	edm::Handle<AssocTrackToTrack> globalMuonOnlyTracks;
	iEvent.getByToken(globalMuonOnlyMapToken_,globalMuonOnlyTracks);
	// Global muon-only with update track
	edm::Handle<AssocTrackToTrack> globalMuonOnlyUpdateTracks;
	iEvent.getByToken(globalMuonOnlyUpdateMapToken_,globalMuonOnlyUpdateTracks);
	
	// picky track
	edm::Handle<AssocTrackToTrack> pickyTracks;
	iEvent.getByToken(pickyMapToken_,pickyTracks);
	// picky No RPC track
	edm::Handle<AssocTrackToTrack> pickyNoRPCTracks;
	iEvent.getByToken(pickyNoRPCToken_,pickyNoRPCTracks);
	// picky muon-only track
	edm::Handle<AssocTrackToTrack> pickyMuonOnlyTracks;
	iEvent.getByToken(pickyMuonOnlyMapToken_,pickyMuonOnlyTracks);
	// picky muon-only with update track
	edm::Handle<AssocTrackToTrack> pickyMuonOnlyUpdateTracks;
	iEvent.getByToken(pickyMuonOnlyUpdateMapToken_,pickyMuonOnlyUpdateTracks);

	// dyt track
	edm::Handle<AssocTrackToTrack> dytTracks;
	iEvent.getByToken(dytMapToken_,dytTracks);
	// dyt No RPC track
	edm::Handle<AssocTrackToTrack> dytNoRPCTracks;
	iEvent.getByToken(dytNoRPCToken_,dytNoRPCTracks);
	// dyt muon-only track
	edm::Handle<AssocTrackToTrack> dytMuonOnlyTracks;
	iEvent.getByToken(dytMuonOnlyMapToken_,dytMuonOnlyTracks);
	// dyt muon-only with update track
	edm::Handle<AssocTrackToTrack> dytMuonOnlyUpdateTracks;
	iEvent.getByToken(dytMuonOnlyUpdateMapToken_,dytMuonOnlyUpdateTracks);

	// tpfms track
	edm::Handle<AssocTrackToTrack> tpfmsTracks;
	iEvent.getByToken(tpfmsMapToken_,tpfmsTracks);
	// tpfms No RPC track
	edm::Handle<AssocTrackToTrack> tpfmsNoRPCTracks;
	iEvent.getByToken(tpfmsNoRPCToken_,tpfmsNoRPCTracks);
	// tpfms muon-only track
	edm::Handle<AssocTrackToTrack> tpfmsMuonOnlyTracks;
	iEvent.getByToken(tpfmsMuonOnlyMapToken_,tpfmsMuonOnlyTracks);
	// tpfms muon-only with update track
	edm::Handle<AssocTrackToTrack> tpfmsMuonOnlyUpdateTracks;
	iEvent.getByToken(tpfmsMuonOnlyUpdateMapToken_,tpfmsMuonOnlyUpdateTracks);

	// combinatoric muon-only track
	edm::Handle<AssocTrackToTrack> bestMuonOnlyTracks;
	iEvent.getByToken(bestMuonOnlyMapToken_,bestMuonOnlyTracks);
	// combinatoric muon-only with update track
	edm::Handle<AssocTrackToTrack> bestMuonOnlyUpdateTracks;
	iEvent.getByToken(bestMuonOnlyUpdateMapToken_,bestMuonOnlyUpdateTracks);


	// list of matched global/SA track indices
	std::vector<int> globalUsed = {-1};
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
		int glb = matchToGenMuon(genMu,*globalTracks,globalUsed);
		bool doComb = false;
		bool doCombGlobal = false;
		bool doCombPicky = false;
		bool doCombDYT = false;
		bool doCombTPFMS = false;
		bool doCombBest = false;
		// If global track is found then fill additional global track refits
		if (glb>=0) {

		  globalUsed.push_back(glb);
			globalInfo.fill((*globalTracks)[glb],beamSpot,PV);

		  // Match tracker track to global track
		  AssocTrackToTrack const & trackerMap = *trackerOnlyTracks;
		  try {
			  auto trackerTrack = 
					*trackerMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				trackerInfo.fill(trackerTrack,beamSpot,PV);
				doComb = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tracker track with update to tracker track" << std::endl;
		  }
		  
		  // Match no rpc global track to global track
		  AssocTrackToTrack const & globalNoRPCMap = *globalNoRPCTracks;
		  try {
		    auto globalNoRPCTrack = 
					*globalNoRPCMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				globalNoRPCInfo.fill(globalNoRPCTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match no rpc global track glb track?" << std::endl;
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
				globalMuonOnlyUpdateInfo.fill(globalMuonOnlyUpdateTrack,beamSpot,PV);
		    doCombGlobal = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only global track with update to global track" << std::endl;
		  }

		  // Match picky track to global track
		  AssocTrackToTrack const & pickyMap = *pickyTracks;
		  try {
		    auto pickyTrack = 
					*pickyMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyInfo.fill(pickyTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match picky track to glb track?" << std::endl;
		  }
		  
		  // Match no rpc picky track to global track
		  AssocTrackToTrack const & pickyNoRPCMap = *pickyNoRPCTracks;
		  try {
		    auto pickyNoRPCTrack = 
					*pickyNoRPCMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				pickyNoRPCInfo.fill(pickyNoRPCTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match no rpc picky track glb track?" << std::endl;
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
				pickyMuonOnlyUpdateInfo.fill(pickyMuonOnlyUpdateTrack,beamSpot,PV);
				doCombPicky = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only picky track with update to picky track" << std::endl;
		  }

		  // Match dyt track to global track
		  AssocTrackToTrack const & dytMap = *dytTracks;
		  try {
		    auto dytTrack = 
					*dytMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytInfo.fill(dytTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match dyt track to glb track?" << std::endl;
		  }
			
		  // Match no rpc dyt track to global track
		  AssocTrackToTrack const & dytNoRPCMap = *dytNoRPCTracks;
		  try {
		    auto dytNoRPCTrack = 
					*dytNoRPCMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				dytNoRPCInfo.fill(dytNoRPCTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match no rpc dyt track glb track?" << std::endl;
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
				dytMuonOnlyUpdateInfo.fill(dytMuonOnlyUpdateTrack,beamSpot,PV);
				doCombDYT = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only dyt track with update to dyt track" << std::endl;
		  }

		  // Match tpfms track to global track
		  AssocTrackToTrack const & tpfmsMap = *tpfmsTracks;
		  try {
		    auto tpfmsTrack = 
					*tpfmsMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tpfmsInfo.fill(tpfmsTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match tpfms track to glb track?" << std::endl;
		  }
			
		  // Match no rpc tpfms track to global track
		  AssocTrackToTrack const & tpfmsNoRPCMap = *tpfmsNoRPCTracks;
		  try {
		    auto tpfmsNoRPCTrack = 
					*tpfmsNoRPCMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				tpfmsNoRPCInfo.fill(tpfmsNoRPCTrack,beamSpot,PV);
		  } catch(...) {
				std::cout << "Could not match no rpc tpfms track glb track?" << std::endl;
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
				tpfmsMuonOnlyUpdateInfo.fill(tpfmsMuonOnlyUpdateTrack,beamSpot,PV);
				doCombTPFMS = true;
		  } catch(...) {
			  std::cout << "Could not match muon-only tpfms track with update to tpfms track" << std::endl;
		  }

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

			// If a muon-only updated track refit is found then store combination with tracker track
			if (doComb){
			  auto trackerTrack = 
					*trackerMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
				if (doCombGlobal) {
					auto globalMuonOnlyUpdateTrack = 
						*globalMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					globalMuonOnlyUpdateTrackerCombInfo.fill(trackerTrack, globalMuonOnlyUpdateTrack);
				}
				if (doCombPicky) {
					auto pickyMuonOnlyUpdateTrack = 
						*pickyMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					pickyMuonOnlyUpdateTrackerCombInfo.fill(trackerTrack, pickyMuonOnlyUpdateTrack);
				}
				if (doCombDYT) {
					auto dytMuonOnlyUpdateTrack = 
						*dytMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					dytMuonOnlyUpdateTrackerCombInfo.fill(trackerTrack, dytMuonOnlyUpdateTrack);
				}
				if (doCombTPFMS) {
					auto tpfmsMuonOnlyUpdateTrack = 
						*tpfmsMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					tpfmsMuonOnlyUpdateTrackerCombInfo.fill(trackerTrack, tpfmsMuonOnlyUpdateTrack);
				}
				if (doCombBest) {
					auto bestMuonOnlyUpdateTrack = 
						*bestMuonOnlyUpdateMap[edm::Ref<std::vector<reco::Track>>(globalTracks,glb)];
					bestMuonOnlyUpdateTrackerCombInfo.fill(trackerTrack, bestMuonOnlyUpdateTrack);
				}
			}

		} else {
		  std::cout << "Could not match Global track to gen muon?" << std::endl;
		} // end global muon

		// Stand alone muons
		int sa = matchToGenMuon(genMu,*standAloneTracks,standAloneUsed);
		if (sa>=0) {
		  standAloneUsed.push_back(sa);
		  auto standAloneTrack = (*standAloneTracks)[sa];
		  standAloneInfo.fill(standAloneTrack,beamSpot,PV);
		} else {
		  std::cout << "Could not match to SA track to gen muon?" << std::endl;
		}

		tree.fill();

   }
   globalUsed.clear();
   standAloneUsed.clear();

   
}

int
AnalyzeTracks::matchToGenMuon(const reco::GenParticle &genMu, const std::vector<reco::Track> &tracks, const std::vector<int> &used) const
{
	//for (auto u : used) std::cout << u << " ";
	//std::cout << std::endl;
	int i = 0;
	for (auto track : tracks) {

		double this_dR = deltaR(genMu,track);
		//std::cout << i << " " << this_dR << std::endl;
		if (std::find(used.begin(),used.end(),i) != used.end()) {
			//std::cout << "already seen " << i;
			i++;
			//std::cout << " " << i << std::endl;
			continue;
		}
		if (this_dR<0.3) {
			return i;
		}
		i++;
	}
	return -1;

	/*
	int i = 0;
	int closest_i = 0;
	double closest_dR = 9999.;
	for (auto u : used) std::cout << u << " ";
	std::cout << std::endl;
	for (auto track : tracks) {
		// If I've already used this track, skip it
		double this_dR = deltaR(genMu,track);
		//if (this_dR>0.3) continue;
		if (std::find(used.begin(),used.end(),i) != used.end()) continue;
		std::cout << this_dR << std::endl;
		if (i==0) {
			// First one is closest by default
			closest_i = i;
			closest_dR = this_dR;
		} else {
			if (this_dR < closest_dR) {
				// This one is now the closest
				closest_i = i;
				closest_dR = this_dR;
			} else {
				// Skip if it isn't closest
				continue;
			}
		}
		i++;
	}
	std::cout << closest_dR << std::endl;
	if (closest_dR>9998.) 
		return -1;
	else 
		return closest_i;
		*/
}


//define this as a plug-in
DEFINE_FWK_MODULE(AnalyzeTracks);
