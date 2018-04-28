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
		input = cms.untracked.int32(500)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:RelValZpMM.root'),
    #fileNames = cms.untracked.vstring('file:RelValZMM.root'),
    fileNames = cms.untracked.vstring('file:muonMC_10_1500_full.root'),
    secondaryFileNames = cms.untracked.vstring(),
	#eventsToProcess = cms.untracked.VEventRange('1:7:623','1:32:3141','1:33:3257')
	#eventsToProcess = cms.untracked.VEventRange('1:7:623')
	#eventsToProcess = cms.untracked.VEventRange('1:32:3141')
	#eventsToProcess = cms.untracked.VEventRange('1:33:3257')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #fileName = cms.untracked.string('highPT_refits_ZpMM_test.root'),
    #fileName = cms.untracked.string('highPT_refits_ZMM_test.root'),
	#fileName = cms.untracked.string('highPT_refits_muonMC_10_1500.root'),
	#fileName = cms.untracked.string('badevent_r1_l7_e623_outsideIn.root'),
	#fileName = cms.untracked.string('badevent_r1_l32_e3141.root'),
	#fileName = cms.untracked.string('highPT_refits_muonMC_10_1500_noRPChits.root'),
	#fileName = cms.untracked.string('bleh.root'),
	#fileName = cms.untracked.string('badevent_r1_l7_e623_outsideIn_test.root'),
	fileName = cms.untracked.string('highPT_test.root'),
	#fileName = cms.untracked.string('highPT_refits_muonMC_1hit_outIn.root'),
	# This forces output to be only 'official' recipe modules
    #outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_mcRun2_asymptotic_v2', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


#process.load("RecoMuon.GlobalMuonProducer.tevMuons_cff")
process.load("RecoMuon.GlobalMuonProducer.highPT_cff")

# Reduce min hits in KF Fitter to 1
process.KFFitterForRefitInsideOut.minHits = cms.int32(1)
process.KFFitterForRefitOutsideIn.minHits = cms.int32(1)

process.highPTMuonsRefit = process.highPTMuons.clone()
process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string('trackRank')
process.highPTMuonsRefit.UtilitiesParameters.trackRankFactor = cms.double(1.)

# do tracker only refit
#process.tevTrackerRefit = process.tevMuons.clone()
#process.tevTrackerRefit.RefitIndex = cms.vint32(6)
#process.tevTrackerRefit.Refits = cms.vstring('trackerRefit')
#process.tevTrackerRefit.RefitterParameters.printStuff = cms.bool(False)

# Path and EndPath definitions
process.p = cms.Path(process.highPTMuonsRefit)# * process.tevTrackerRefit)


process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


# Schedule definition
process.schedule = cms.Schedule(process.p,process.RECOSIMoutput_step)

