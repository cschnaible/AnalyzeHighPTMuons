# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --datatier GEN-SIM-RECO,DQMIO --conditions auto:run1_mc -s RAW2DIGI,L1Reco,RECO,EI,VALIDATION:@standardValidationNoHLT,DQM:@standardDQMFakeHLT --eventcontent RECOSIM,DQM -n 100 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process('RECO3')

options = VarParsing.VarParsing()
options.register('selector',\
		'curvPull',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Selector for best combinatoric refit (curvPull is default)')
options.register('MC',\
		'muonMC',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'MC to run on : ZpMM, ZMM, or muonMC (default)')
options.register('extra',\
		'',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Extra name to append to output file')
options.register('rankfactor',\
		1.,\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.float,\
		'Track rank = NHits * rankfactor - chi2')
options.register('curvPullCut',\
		1.,\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.float,\
		'Pre-select tracks with curvature pull < curvPullCut')
options.register('maxEvents',\
		# somewhere above 94000 i got a seg fault in muonMC?
		25000,\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.int,\
		'Maximum number of events to run')
options.parseArguments()

# some checks
selectors = ['dxy','curvPull','TEST','trackRank']
if options.selector not in selectors:
	raise ValueError(options.selector+' is not a valid selector')

# Set input file name
inputDir = 'file:/scratch3/HighPT/data/'
filename = ''
if options.MC=='ZpMM':
	filename = 'RelValZpMM.root'
elif options.MC == 'ZMM': 
	filename = 'RelValZMM.root'
elif options.MC == 'muonMC':
	#filename = 'muonMC_10_1500_pruned.root'
	filename = 'muonMC_10_1500_full.root'
else:
	raise ValueError(options.MC+' is not a valid MC')
if options.selector=='trackRank':
	extra = 'f'+str(int(options.rankfactor)) + ('_' + options.extra if options.extra else '')
else:
	extra = options.extra
inputFile = inputDir+filename
# Set output file name
outputFile = 'highPT_refits_'+options.MC+\
		('_'+options.selector)+\
		('_'+extra if extra else '')+\
		'.root'

print 'Input file :',inputFile
print 'Output file ',outputFile


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
		input = cms.untracked.int32(options.maxEvents)
		# Somewhere above 94000 there is a seg fault...
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
		SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(inputFile),
    secondaryFileNames = cms.untracked.vstring(),
	#eventsToProcess = cms.untracked.VEventRange('1:1:2')
	#eventsToProcess = cms.untracked.VEventRange('1:7:623','1:32:3141','1:33:3257')
	#eventsToProcess = cms.untracked.VEventRange('1:7:623')
	#eventsToProcess = cms.untracked.VEventRange('1:32:3141')
	#eventsToProcess = cms.untracked.VEventRange('1:33:3257')
	#eventsToProcess = cms.untracked.VEventRange('1:38:3770') # this one sucks
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
	#fileName = cms.untracked.string('highPT_refits_muonMC_curvPull.root'),
	fileName = cms.untracked.string(outputFile),
	# This forces output to be only 'official' recipe modules
    #outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
if options.MC=='ZpMM':
	process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')
elif options.MC=='ZMM':
	process.GlobalTag = GlobalTag(process.GlobalTag, '92X_mcRun2_asymptotic_v2', '')
elif options.MC=='muonMC':
	process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
else:
	raise ValueError(options.MC+' is not a valid MC')

process.load("RecoMuon.GlobalMuonProducer.tevMuons_cff")
process.load("RecoMuon.GlobalMuonProducer.highPT_cff")

# Reduce min hits in KF Fitter to 1
process.KFFitterForRefitInsideOut.minHits = cms.int32(1)
process.KFFitterForRefitOutsideIn.minHits = cms.int32(1)

# Muon-only refits
process.highPTMuonsRefit = process.highPTMuons.clone()
process.highPTMuonsRefit.Refits = cms.vstring('default','firstHit','picky','dyt','combinatoric')
#process.highPTMuonsRefit.Refits = cms.vstring('combinatoric')
process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string(options.selector)
#process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string('TEST')
#process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string('dxy')
#process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string('curvPull')
#process.highPTMuonsRefit.UtilitiesParameters.Selector = cms.string('trackRank')
# always set but not always used
process.highPTMuonsRefit.UtilitiesParameters.trackRankFactor = cms.double(options.rankfactor)
process.highPTMuonsRefit.UtilitiesParameters.curvPullCut = cms.double(options.curvPullCut)


# Path and EndPath definitions
process.p = cms.Path(process.highPTMuonsRefit)

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)


# Schedule definition
process.schedule = cms.Schedule(process.p,process.RECOSIMoutput_step)

