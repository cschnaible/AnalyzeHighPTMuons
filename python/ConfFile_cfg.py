import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
		)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
		'file:DOESNOTEXIST.root'
    )
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.analyzer = cms.EDAnalyzer('AnalyzeTracks',

	vertices = cms.InputTag('offlinePrimaryVertices'),
	beamSpot = cms.InputTag('offlineBeamSpot'),

	globalTracks = cms.InputTag('globalMuons'),
	genParticles = cms.InputTag('genParticles'),
	standAloneTracks = cms.InputTag('standAloneMuons','UpdatedAtVtx'),

	trackerMapTag = cms.InputTag('tevTrackerRefit','trackerRefit'),
	#trackerMapTag = cms.InputTag('highPTMuonsRefit','tracker'),

	globalMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','default'),
	globalMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','default'),

	pickyMapTag = cms.InputTag('tevMuons','picky'),
	pickyMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','picky'),
	pickyMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','pickyVtxUpdate'),

	dytMapTag = cms.InputTag('tevMuons','dyt'),
	dytMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','dyt'),
	dytMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','dytVtxUpdate'),

	tpfmsMapTag = cms.InputTag('tevMuons','firstHit'),
	tpfmsMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','firstHit'),
	tpfmsMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','firstHitVtxUpdate'),

	bestMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','combinatoric'),
	bestMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','combinatoricVtxUpdate'),

	#filename = cms.untracked.string('DOESNOTEXIST.root'),
	#doPicky = cms.bool(False)
)

# Zpmm
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_refits_ZpMM_test.root')
#process.analyzer.filename = cms.untracked.string('AnalyzeTracks_ZpMM.root')
# ZMM
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_refits_ZMM_test.root')
#process.analyzer.filename = cms.untracked.string('AnalyzeTracks_ZMM.root')
# Muon gun pT = 10 - 1500 GeV
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_refits_muonMC_10_1500.root')
#process.analyzer.filename = cms.untracked.string('AnalyzeTracks_muonPT.root')
# Muon gun pT = 10 - 1500 GeV with Picky track refits
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/test.root')
#process.analyzer.filename = cms.untracked.string('AnalyzeTracks_muonPT_picky.root')
# test
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_refits_muonMC_1hit_outIn.root')
#process.analyzer.filename = cms.untracked.string('AnalyzeTracks_muonPT_1hit_outIn.root')
# No RPC hits, at least 3 DT+CSC
#process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_refits_muonMC_10_1500_noRPChits.root')
#process.analyzer.filename = cms.untracked.string('analyze_no_RPC_hits.root')

process.source.fileNames = cms.untracked.vstring('file:/home/cschnaib/931/highPT_test_50k.root')
process.TFileService = cms.Service('TFileService',fileName=cms.string('AnalyzeTracks_test_50k_04262018.root'))

process.p = cms.Path(process.analyzer)
