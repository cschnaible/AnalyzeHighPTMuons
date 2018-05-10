import sys
import ROOT as R
import HighPTMuons.Analysis.tools as t
from HighPTMuons.Analysis.plot_styles import styles

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

MC = 'muonMC'

#selector = 'curvPull'
#selector = 'dxy'
f = str(15)
selector = 'trackRank_f'+f
extra = ''

fileDir = '/scratch3/HighPT/CMSSW_9_4_6_patch1/src/HighPTMuons/Analysis/test/RunAnalyzeTracks/'
#fileName = 'AnalyzeTracks_'+MC+'_'+selector+('_'+extra if extra else '')+'.root'
fileName = 'AnalyzeTracks_'+MC+'_'+selector+'.root'
inputFile = fileDir+fileName
outputFile = 'hists_curvature_'+MC+'_'+selector+('_'+extra if extra else '')+'.root'

print 'Input file :',inputFile
print 'Output file :',outputFile

inFile = R.TFile(inputFile)
tree = inFile.Get('analyzer/HighPTTrackTree')

outFile = R.TFile(outputFile,'recreate')

tracktypes = ['global','picky','tpfms','dyt','standAlone','best','tracker',\
		'globalNoRPC','pickyNoRPC','tpfmsNoRPC','dytNoRPC']
fits = ['','MuonOnly','MuonOnlyUpdate','Comb']

varList = {
		'gen_eta':{
			'bins':(60,-3,3),
			'label':'#eta',
			'log':False,
			},
		'gen_pt':{
			'bins':(150,0,1500),
			'label':'gen p_{T} [GeV]',
			'log':False,
			},
		'K_gen_res':{
			'bins':(100,-1,1),
			'label':'#kappa(track)-#kappa(gen)/#kappa(gen)',
			'log':False,
			},
		}

combTypes = ['simple','full']
ptcuts = [0,200,400,600,800,1000]
HISTS = {tracktype:{fit:{ptcut:{} for ptcut in ptcuts} for fit in fits} for tracktype in tracktypes}

for tracktype in tracktypes:
	for fit in fits:
		for combType in combTypes:
			for ptcut in ptcuts:
				htitle = tracktype+fit

				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType] =\
						R.TH1D(htitle+'_K_gen_res_'+combType+'_'+str(ptcut),'',*varList['K_gen_res']['bins'])

				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType+'_vs_gen_pt'] =\
						R.TH2D(htitle+'_K_gen_res_'+combType+'_vs_gen_pt_'+str(ptcut),'',\
								*(varList['gen_pt']['bins']+varList['K_gen_res']['bins']))

				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType+'_vs_gen_eta'] =\
						R.TH2D(htitle+'_K_gen_res_'+combType+'_vs_gen_eta_'+str(ptcut),'',\
								*(varList['gen_eta']['bins']+varList['K_gen_res']['bins']))
		

for l,entry in enumerate(tree):
	#if l>2000: break
	if l%1000==0: print 'Entry {l:>5}\r'.format(**locals())
	p = t.pyTree(tree)
	cuts = {
			0:p.gen_pt>0,
			200:p.gen_pt>200.,
			400:p.gen_pt>400.,
			600:p.gen_pt>600.,
			800:p.gen_pt>800.,
			1000:p.gen_pt>1000.
			}

	# Plots for all tracks
	for trackType in tracktypes:
		for fit in fits:
			if (trackType=='tracker' or trackType=='standAlone') and\
					(fit=='MuonOnly' or fit=='MuonOnlyUpdate' or fit=='Comb'):
				continue
			if (trackType=='best' and fit==''): continue
			if ('NoRPC' in trackType and fit!=''): continue
			trackname = trackType+fit
			track = t.Track(trackname,p)
			if track._bad==True: continue

			for ptcut in ptcuts:
				if not cuts[ptcut]: continue
				# Simple 1D combination
				if fit=='MuonOnlyUpdate':
					tracker = t.Track('tracker',p)
					k_simple_err = 1./(1./track.cov(0,0) + 1./tracker.cov(0,0))
					k_simple = (track.K()/track.cov(0,0) + tracker.K()/tracker.cov(0,0)) * k_simple_err
					res_simple = (k_simple-p.gen_K)/p.gen_K
					HISTS[trackType][fit][ptcut]['K_gen_res_simple'].Fill(res_simple)
					HISTS[trackType][fit][ptcut]['K_gen_res_simple_vs_gen_pt'].Fill(p.gen_pt, res_simple)
					HISTS[trackType][fit][ptcut]['K_gen_res_simple_vs_gen_eta'].Fill(p.gen_eta, res_simple)
				# Full 5D combination
				res_full = (track.K() - p.gen_K) / p.gen_K
				HISTS[trackType][fit][ptcut]['K_gen_res_full'].Fill(res_full)
				HISTS[trackType][fit][ptcut]['K_gen_res_full_vs_gen_pt'].Fill(p.gen_pt, res_full)
				HISTS[trackType][fit][ptcut]['K_gen_res_full_vs_gen_eta'].Fill(p.gen_eta, res_full)

outFile.cd()
for tracktype in tracktypes:
	for fit in fits:
		for combType in combTypes:
			for ptcut in ptcuts:
				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType].Write()
				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType+'_vs_gen_pt'].Write()
				HISTS[tracktype][fit][ptcut]['K_gen_res_'+combType+'_vs_gen_eta'].Write()

outFile.Close()
