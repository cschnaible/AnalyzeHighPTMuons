import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()
options.register('MC',\
		'muonMC',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'MC to run on : ZpMM, ZMM, or muonMC (default)')
options.register('input',\
		'',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Extra name for input ROOT file')
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
options.register('output',\
		'',\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.string,\
		'Extra name for output ROOT file')
options.register('correlation',\
		1.0,\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.float,\
		'Correlation between lambda, phi, dxy, dsz of mu-only+vtx track and tracker track')
options.register('covScale',\
		1.0,\
		VarParsing.VarParsing.multiplicity.singleton,\
		VarParsing.VarParsing.varType.float,\
		'Scale factor applied to mu-only+vtx track covariance')
options.parseArguments()

inputFile = \
		options.inputDir+\
		'highPT_refits_'+\
		options.MC+\
		('_'+options.input if options.input else '')+\
		'.root'
outputFile = \
		options.outputDir+\
		'AnalyzeTracks_'+\
		options.MC+\
		('_'+options.output if options.output else '')+\
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

	correlation = cms.double(options.correlation),
	covScale = cms.double(options.covScale),

	vertices = cms.InputTag('offlinePrimaryVertices'),
	beamSpot = cms.InputTag('offlineBeamSpot'),

	recoMuons = cms.InputTag('muons'),
	globalTracks = cms.InputTag('globalMuons'),
	genParticles = cms.InputTag('genParticles'),
	standAloneTracks = cms.InputTag('standAloneMuons','UpdatedAtVtx'),

	trackerMapTag = cms.InputTag('highPTMuonsTrackerRefit','tracker'),

	globalMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitStd','default'),
	globalMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitStd','defaultVtxUpdate'),

	pickyMapTag = cms.InputTag('tevMuons','picky'),
	pickyMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitStd','picky'),
	pickyMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitStd','pickyVtxUpdate'),

	dytMapTag = cms.InputTag('tevMuons','dyt'),
	dytMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitStd','dyt'),
	dytMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitStd','dytVtxUpdate'),

	tpfmsMapTag = cms.InputTag('tevMuons','firstHit'),
	tpfmsMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitStd','firstHit'),
	tpfmsMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitStd','firstHitVtxUpdate'),

	globalTrackRankMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitGlobalTrackRank','combinatoric'),
	globalTrackRankMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitGlobalTrackRank','combinatoricVtxUpdate'),
	pickyTrackRankMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitPickyTrackRank','combinatoric'),
	pickyTrackRankMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitPickyTrackRank','combinatoricVtxUpdate'),
	dytTrackRankMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitDYTTrackRank','combinatoric'),
	dytTrackRankMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitDYTTrackRank','combinatoricVtxUpdate'),
	tunePTrackRankMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitTunePTrackRank','combinatoric'),
	tunePTrackRankMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitTunePTrackRank','combinatoricVtxUpdate'),

)

process.source.fileNames = cms.untracked.vstring('file:'+inputFile)
process.TFileService = cms.Service('TFileService',fileName=cms.string(outputFile))

process.p = cms.Path(process.analyzer)












#	trkCurvMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitCurvPullTrk','combinatoric'),
#	trkCurvMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitCurvPullTrk','combinatoricVtxUpdate'),
#	tunePCurvMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitCurvPullTuneP','combinatoric'),
#	tunePCurvMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitCurvPullTuneP','combinatoricVtxUpdate'),
#	dxyMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitDxy','combinatoric'),
#	dxyMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitDxy','combinatoricVtxUpdate'),
#	curvRelErrMuonOnlyMapTag = cms.InputTag('highPTMuonsRefitRelCurvErr','combinatoric'),
#	curvRelErrMuonOnlyUpdateMapTag = cms.InputTag('highPTMuonsRefitRelCurvErr','combinatoricVtxUpdate'),
