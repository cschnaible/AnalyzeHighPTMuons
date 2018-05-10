# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --datatier GEN-SIM-RECO,DQMIO --conditions auto:run1_mc -s RAW2DIGI,L1Reco,RECO,EI,VALIDATION:@standardValidationNoHLT,DQM:@standardDQMFakeHLT --eventcontent RECOSIM,DQM -n 100 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
		input = cms.untracked.int32(-1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/scratch3/HighPT/data/muonMC_10_1500_full.root'),
    secondaryFileNames = cms.untracked.vstring(),
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('/scratch3/HighPT/data/muonMC_10_1500_pruned.root'),
		# This forces output to be only 'official' recipe modules
    #outputCommands = process.RECOSIMEventContent.outputCommands,
		# Select branches to keep
    outputCommands = cms.untracked.vstring(
			'drop *_*_*_*',

			'keep *_muons_*_*',
			'keep *_genParticles_*_*',
			'keep *_generator_*_*',

			'keep *_offlinePrimaryVertices_*_*',
			'keep *_OfflinePrimaryVerticesBS_*_*',
			'keep *_offlineBeamSpot_*_*',

			'keep *_csc2DRecHits_*_*',
			'keep *_cscSegments_*_*',
			'keep *_dt4DSegments_*_*',
			'keep *_dt1DRecHits_*_*',
			'keep *_rpcRecHits_*_*',
			'keep *_dedxHitInfo_*_*',
			'keep *_siPixelClusters_*_*',
			'keep *_siStripClusters_*_*',

			'keep *_generalTracks_*_*',
			'keep *_standAloneMuons_*_*',
			'keep *_refittedStandAloneMuons_*_*',
			'keep *_globalMuons_*_*',

			'keep *_trackExtrapolator_*_*',
			'keep *_tevMuons_*_*',
			'keep *_tevMuonsNoRPC_*_*',
			'keep *_highPTMuonsRefit_*_*',
			'keep *_highPTMuonsTrackerRefit_*_*',
		),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load("RecoMuon.GlobalMuonProducer.tevMuons_cff")
process.load("RecoMuon.GlobalMuonProducer.highPT_cff")

# Reduce min hits in KF Fitter to 1
process.KFFitterForRefitInsideOut.minHits = cms.int32(1)
process.KFFitterForRefitOutsideIn.minHits = cms.int32(1)

# Tracker Refits
process.highPTMuonsTrackerRefit = process.highPTMuons.clone()
process.highPTMuonsTrackerRefit.Refits = cms.vstring('tracker')

# TeV Refits w/o RPC hits
process.tevMuonsNoRPC = process.tevMuons.clone()
process.tevMuonsNoRPC.RefitterParameters.RefitRPCHits = cms.bool(False)

# Path and EndPath definitions
process.p = cms.Path(process.highPTMuonsTrackerRefit * process.tevMuonsNoRPC)

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


# Schedule definition
process.schedule = cms.Schedule(process.p,process.RECOSIMoutput_step)

