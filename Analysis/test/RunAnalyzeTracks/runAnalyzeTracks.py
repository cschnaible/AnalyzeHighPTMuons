import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()
options.register('MC',\
		'muonMC',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'MC to run on : ZpMM, ZMM, or muonMC (default)')
options.register('selector',\
		'asdf',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Selector for best combinatoric refit')
options.register('inputDir',\
		'/scratch3/HighPT/CMSSW_9_4_6_patch1/src/HighPTMuons/Production/RunHighPTRefitter/',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Directory for input file')
options.register('outputDir',\
		'',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Directory for output file')
options.register('extra',\
		'',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Extra label to put at end of output file name')
options.parseArguments()

inputFile = \
		options.inputDir+\
		'highPT_refits_'+\
		options.MC+\
		('_'+options.selector if options.selector else '')+\
		('_'+options.extra if options.extra else '')+\
		'.root'
outputFile = \
		options.outputDir+\
		'AnalyzerTracks_'+\
		options.MC+\
		('_'+options.selector if options.selector else '')+\
		('_'+options.extra if options.extra else '')+\
		'.root'

print 'Input file :',inputFile
print 'Output file :',outputFile

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

	trackerMapTag = cms.InputTag('highPTMuonsTrackerRefit','tracker'),

	globalNoRPCMapTag = cms.InputTag('tevMuonsNoRPC','default'),
	globalMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','default'),
	globalMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','default'),

	pickyMapTag = cms.InputTag('tevMuons','picky'),
	pickyNoRPCMapTag = cms.InputTag('tevMuonsNoRPC','picky'),
	pickyMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','picky'),
	pickyMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','pickyVtxUpdate'),

	dytMapTag = cms.InputTag('tevMuons','dyt'),
	dytNoRPCMapTag = cms.InputTag('tevMuonsNoRPC','dyt'),
	dytMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','dyt'),
	dytMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','dytVtxUpdate'),

	tpfmsMapTag = cms.InputTag('tevMuons','firstHit'),
	tpfmsNoRPCMapTag = cms.InputTag('tevMuonsNoRPC','firstHit'),
	tpfmsMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','firstHit'),
	tpfmsMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','firstHitVtxUpdate'),

	bestMuonOnlyMapTag = cms.InputTag('highPTMuonsRefit','combinatoric'),
	bestMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefit','combinatoricVtxUpdate'),

)

process.source.fileNames = cms.untracked.vstring('file:'+inputFile)
process.TFileService = cms.Service('TFileService',fileName=cms.string(outputFile))

process.p = cms.Path(process.analyzer)
