import sys
import ROOT as R
import math as math
import numpy as np
import argparse
R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)
parser = argparse.ArgumentParser()
parser.add_argument('-in','--infile',dest='infile',default='ZpMM',help='Which MC file to use')
parser.add_argument('-out','--outfile',dest='outfile',default='ZpMM',help='Name of output file')
parser.add_argument('-ex','--extra',default='',help='Extra name to output histogram file')
args = parser.parse_args()

class HistogramConfigurations(object):
	def __init__(self,trackName):

		self.stdfits = ['global','picky','dyt','tpfms','tuneP']
		self.trackName = trackName

		self.attr = ('name',\
				'xTitle','yTitle',\
				'nBins','binLow','binHigh')
		self.attr2D = ('name',\
				'xTitle','yTitle',\
				'nBinsX','binLowX','binHighX','nBinsY','binLowY','binHighY')

		# 1D
		self.data1D = {}
		self.data1D['genCurvRes'] = self.makeAttrDict((self.HName('genCurvRes'),\
				'#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}','Counts',\
				500,-5,5))
		self.data1D['genCurvPull'] = self.makeAttrDict((self.HName('genCurvPull'),\
				'#frac{#kappa - #kappa_{gen}}{#sigma(#kappa)}','Counts',\
				250,-5,5))
#		self.data1D['genCurvDiff'] = self.makeAttrDict((self.HName('genCurvDiff'),\
#				'#kappa - #kappa_{gen}','Counts',\
#				100,-0.0025,0.0025))
#		self.data1D['genCurvAbsDiff'] = self.makeAttrDict((self.HName('genCurvAbsDiff'),\
#				'|#kappa - #kappa_{gen}|','Counts',\
#				200,0,0.0025))
		# 2D
		self.data2D = {}
		self.data2D['genCurvResVsPt'] = self.makeAttrDict((self.HName('genCurvResVsPt'),\
				'gen p_{T} [GeV]','#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',\
				75,0,1500,500,-5,5))
		self.data2D['genCurvPullVsPt'] = self.makeAttrDict((self.HName('genCurvPullVsPt'),\
				'gen p_{T} [GeV]','#frac{#kappa - #kappa_{gen}}{#sigma(#kappa)}',\
				75,0,1500,250,-5,5))
		self.data2D['genCurvResVsEta'] = self.makeAttrDict((self.HName('genCurvResVsEta'),\
				'gen #eta','#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',\
				48,-2.4,2.4,500,-5,5))
		self.data2D['genCurvPullVsEta'] = self.makeAttrDict((self.HName('genCurvPullVsEta'),\
				'gen #eta','#frac{#kappa - #kappa_{gen}}{#sigma(#kappa)}',\
				48,-2.4,2.4,250,-5,5))
#		self.data2D['genCurvDiffVsPt'] = self.makeAttrDict((self.HName('genCurvDiffVsPt'),\
#				'gen p_{T} [GeV]','#kappa - #kappa_{gen}',\
#				75,0,1500,100,-0.0025,0.0025))
#		self.data2D['genCurvAbsDiffVsPt'] = self.makeAttrDict((self.HName('genCurvAbsDiffVsPt'),\
#				'gen p_{T} [GeV]','|#kappa - #kappa_{gen}|',\
#				75,0,1500,200,0,0.0025))
#		self.data2D['genCurvDiffVsEta'] = self.makeAttrDict((self.HName('genCurvDiffVsEta'),\
#				'gen #eta','#kappa - #kappa_{gen}',
#				48,-2.4,2.4,100,-0.0025,0.0025))
#		self.data2D['genCurvAbsDiffVsEta'] = self.makeAttrDict((self.HName('genCurvAbsDiffVsEta'),\
#				'gen #eta','|#kappa - #kappa_{gen}|',\
#				48,-2.4,2.4,200,0,0.0025))
		#
		for stdfit in self.stdfits:
			if stdfit+'Comb' == trackName:
				self.data1D['fullCurvRes'] = self.makeAttrDict((self.HName('fullCurvRes'),\
						'#frac{#kappa_{comb} - #kappa_{full}}{#kappa_{full}}','Counts',\
						500,-5,5))
				self.data1D['fullCurvPull'] = self.makeAttrDict((self.HName('fullCurvPull'),\
						'#frac{#kappa_{comb} - #kappa_{full}}{#sigma(#kappa_{comb})}','Counts',\
						250,-5,5))
				self.data1D['fullCurvDiff'] = self.makeAttrDict((self.HName('fullCurvDiff'),\
						'#kappa_{comb} - #kappa_{full}','Counts',\
						100,-0.0025,0.0025))
				self.data1D['fullCurvAbsDiff'] = self.makeAttrDict((self.HName('fullCurvAbsDiff'),\
						'#kappa_{comb} - #kappa_{full}','Counts',\
						200,0,0.0025))
				# 2D
				self.data2D['fullCurvResVsPt'] = self.makeAttrDict((self.HName('fullCurvResVsPt'),\
						'gen p_{T} [GeV]','#frac{#kappa_{comb} - #kappa_{full}}{#kappa_{full}}',\
						75,0,1500,500,-5,5))
				self.data2D['fullCurvPullVsPt'] = self.makeAttrDict((self.HName('fullCurvPullVsPt'),\
						'gen p_{T} [GeV]','#frac{#kappa_{comb} - #kappa_{full}}{#sigma(#kappa_{comb})}',
						75,0,1500,250,-5,5))
				self.data2D['fullCurvResVsEta'] = self.makeAttrDict((self.HName('fullCurvResVsEta'),\
						'gen #eta','#frac{#kappa_{comb} - #kappa_{full}}{#kappa_{full}}',\
						48,-2.4,2.4,500,-5,5))
				self.data2D['fullCurvPullVsEta'] = self.makeAttrDict((self.HName('fullCurvPullVsEta'),\
						'gen #eta','#frac{#kappa_{comb} - #kappa_{full}}{#sigma(#kappa_{comb})}',\
						48,-2.4,2.4,250,-5,5))
#				self.data2D['fullCurvDiffVsPt'] = self.makeAttrDict((self.HName('fullCurvDiffVsPt'),\
#						'gen p_{T} [GeV]','#kappa - #kappa_{gen}',\
#						75,0,1500,100,-0.0025,0.0025))
#				self.data2D['fullCurvAbsDiffVsPt'] = self.makeAttrDict((self.HName('fullCurvAbsDiffVsPt'),\
#						'gen p_{T} [GeV]','|#kappa - #kappa_{gen}|',\
#						75,0,1500,200,0,0.0025))
#				self.data2D['fullCurvDiffVsEta'] = self.makeAttrDict((self.HName('fullCurvDiffVsEta'),\
#						'gen #eta','#kappa - #kappa_{gen}',\
#						48,-2.4,2.4,100,-0.0025,0.0025))
#				self.data2D['fullCurvAbsDiffVsEta'] = self.makeAttrDict((self.HName('fullCurvAbsDiffVsEta'),\
#						'gen #eta','|#kappa - #kappa_{gen}|',\
#						48,-2.4,2.4,200,0,0.0025))

		# cuts
		self.histcuts = {
				'eta':{
					'all':'fabs(gen_eta)>=0',
					'barrel':'fabs(gen_eta)<1.2',
					'endcap':'fabs(gen_eta)>1.2',
					'fwdendcap':'fabs(gen_eta)>2.1',
					},
				'pt':{
					'all':'gen_pt>0',
					'200':'gen_pt>200',
					#'400':gen_pt>400',
					'600':'gen_pt>600',
					#'800':'gen_pt>800',
					'1000':'gen_pt>1000',
					}
				}

	def HName(self,key):
		return self.trackName+'_'+key

	def makeAttrDict(self,tup):
		if len(tup)==len(self.attr):
			return dict(zip(self.attr,tup))
		elif len(tup)==len(self.attr2D):
			return dict(zip(self.attr2D,tup))
		else:
			print len(tup) , len(self.attr), len(self.attr2D)
			raise ValueError('Cannot make attribute dictionary from tuple')

	def cuts(self,key):
		if key in self.data1D.keys():
			cutdict = {ptkey+'_'+etakey:{} for ptkey in self.histcuts['pt'].keys() for etakey in self.histcuts['eta'].keys()}
			for ptkey in self.histcuts['pt'].keys():
				for etakey in self.histcuts['eta'].keys():
					cutdict[ptkey+'_'+etakey] = self.histcuts['pt'][ptkey]+' && '+self.histcuts['eta'][etakey]
			return cutdict
		elif key in self.data2D.keys():
			y,x = key.split('Vs')
			if x=='Pt':
				return self.histcuts['eta']
			elif x=='Eta':
				return self.histcuts['pt']
			else: raise ValueError(key+' is not a valid histogram type')

	def getstuff(self,trackName,key,cutkey):
		if key in self.data1D.keys():
			AD = self.data1D[key]
			return AD['name']+'_'+cutkey, ';{};{}'.format(AD['xTitle'], AD['yTitle']), AD['nBins'], AD['binLow'], AD['binHigh']
		elif key in self.data2D.keys():
			AD = self.data2D[key]
			return AD['name']+'_'+cutkey, ';{};{}'.format(AD['xTitle'], AD['yTitle']), AD['nBinsX'], AD['binLowX'], AD['binHighX'], AD['nBinsY'], AD['binLowY'], AD['binHighY']


def Draw(t,HConfig,trackName,key,cutkey,expression):
	if 'full' in key: trackNameFull=trackName.strip('Comb')
	toPlot = expression.format(**locals())
	hName = trackName+'_'+key+'_'+cutkey
	cuts = HConfig.cuts(key)
	drawCMD = '{toPlot}>>{hName}'.format(toPlot=toPlot,hName=hName)
	drawCUTS = '{trackName}_good'.format(trackName=trackName) +' && '+cuts[cutkey].format(trackName=trackName)
	print hName
	print toPlot,drawCUTS
	t.Draw(drawCMD,drawCUTS)


def fillPlots(HISTS):
	# get file and tree
	inputFile = R.TFile.Open(args.infile)
	t = inputFile.Get('analyzer/HighPTTrackTree')
	for trackName in TrackNames:
		print trackName
		HConfig = HistogramConfigurations(trackName)
		for key in HConfig.data1D.keys():
			for cutkey in HConfig.cuts(key):
				if 'Comb' not in trackName and 'full' in key: continue
				HISTS[trackName][key+'_'+cutkey] = R.TH1D(*HConfig.getstuff(trackName,key,cutkey))
				print HConfig.getstuff(trackName,key,cutkey)
				Draw(t,HConfig,trackName,key,cutkey,HExpressions[key])
				HISTS[trackName][key+'_'+cutkey].SetDirectory(0)
		for key in HConfig.data2D.keys():
			for cutkey in HConfig.cuts(key):
				if 'Comb' not in trackName and 'full' in key: continue
				HISTS[trackName][key+'_'+cutkey] = R.TH2D(*HConfig.getstuff(trackName,key,cutkey))
				print HConfig.getstuff(trackName,key,cutkey)
				Draw(t,HConfig,trackName,key,cutkey,HExpressions[key])
				HISTS[trackName][key+'_'+cutkey].SetDirectory(0)
#
# Do stuff
#

HExpressions = {
		'genCurvRes'   :'({trackName}_par[0]-gen_K)/gen_K',
		'genCurvPull'  :'({trackName}_par[0]-gen_K)/sqrt({trackName}_cov[0,0])',
		'genCurvDiff'  :'{trackName}_par[0]-gen_K',
		'genCurvAbsDiff':'fabs({trackName}_par[0] - gen_K)',
		'genCurvResVsPt'  :'({trackName}_par[0] - gen_K) / gen_K:gen_pt',
		'genCurvResVsEta' :'({trackName}_par[0] - gen_K) / gen_K:gen_eta',
		'genCurvPullVsPt' :'({trackName}_par[0] - gen_K) / sqrt({trackName}_cov[0,0]):gen_pt',
		'genCurvPullVsEta':'({trackName}_par[0] - gen_K) / sqrt({trackName}_cov[0,0]):gen_eta',
		'genCurvDiffVsPt' :'({trackName}_par[0] - gen_K):gen_pt',
		'genCurvDiffVsEta':'({trackName}_par[0] - gen_K):gen_eta',
		'genCurvAbsDiffVsPt' :'fabs({trackName}_par[0] - gen_K):gen_pt',
		'genCurvAbsDiffVsEta':'fabs({trackName}_par[0] - gen_K):gen_eta',
		'fullCurvRes'   :'({trackName}_par[0] - {trackNameFull}_par[0]) / {trackNameFull}_par[0]',
		'fullCurvPull'  :'({trackName}_par[0] - {trackNameFull}_par[0]) / sqrt({trackName}_cov[0,0])',
		'fullCurvDiff'  :'{trackName}_par[0] - {trackNameFull}_par[0]',
		'fullCurvAbsDiff':'fabs({trackName}_par[0] - {trackNameFull}_par[0])',
		'fullCurvResVsPt'  :'({trackName}_par[0] - {trackNameFull}_par[0]) / {trackNameFull}_par[0]:gen_pt',
		'fullCurvResVsEta' :'({trackName}_par[0] - {trackNameFull}_par[0]) / {trackNameFull}_par[0]:gen_eta',
		'fullCurvPullVsPt' :'({trackName}_par[0] - {trackNameFull}_par[0]) / sqrt({trackName}_cov[0,0]):gen_pt',
		'fullCurvPullVsEta':'({trackName}_par[0] - {trackNameFull}_par[0]) / sqrt({trackName}_cov[0,0]):gen_eta',
		'fullCurvDiffVsPt' :'({trackName}_par[0] - {trackNameFull}_par[0]):gen_pt',
		'fullCurvDiffVsEta':'({trackName}_par[0] - {trackNameFull}_par[0]):gen_eta',
		'fullCurvAbsDiffVsPt' :'fabs({trackName}_par[0] - {trackNameFull}_par[0]):gen_pt',
		'fullCurvAbsDiffVsEta':'fabs({trackName}_par[0] - {trackNameFull}_par[0]):gen_eta',
		}

#BaseTrackNames = ['global','picky','dyt','tpfms','dxy','trkCurv','tunePCurv','curvRelErr','standAlone','tracker','tuneP','trackRank']
#BaseTrackNames = ['global','tuneP','trackRank','tracker','tpfms','dyt','picky']
BaseTrackNames = ['tunePTrackRank','globalTrackRank','pickyTrackRank','dytTrackRank',\
		'picky','dyt','tuneP','global','tracker','tpfms']
TrackTypes = ['','MuonOnly','MuonOnlyUpdate','Comb']
TrackNames = []
for baseName in BaseTrackNames:
	for trackType in TrackTypes:
		if baseName in ['dxy','trkCurv','tunePCurv','curvRelErr','trackRank'] and trackType=='': continue
		if 'TrackRank' in baseName and trackType=='': continue
		if baseName in ['standAlone','tracker'] and trackType!='': continue
		TrackNames.append(baseName+trackType)

HISTS = {trackName:{} for trackName in TrackNames}
fillPlots(HISTS)

outputFile = R.TFile(args.outfile,'RECREATE')
outputFile.cd()
for trackName in HISTS.keys():
	for key in HISTS[trackName].keys():
		HISTS[trackName][key].Write()

