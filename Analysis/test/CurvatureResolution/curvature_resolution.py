import sys
import ROOT as R
import HighPTMuons.Analysis.tools as t
from HighPTMuons.Analysis.plot_styles import styles

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

#f = str(1)

#fileName = '/scratch3/HighPT/CMSSW_9_4_6_patch1/src/HighPTMuons/Analysis/test/RunAnalyzeTracks/analyze_muonMC_f'+f+'_test.root'
fileName = '/scratch3/HighPT/CMSSW_9_4_6_patch1/src/HighPTMuons/Analysis/test/RunAnalyzeTracks/analyze_muonMC_bestPull_test.root'
#fileName = '/scratch3/HighPT/CMSSW_9_4_6_patch1/src/HighPTMuons/Analysis/test/RunAnalyzeTracks/analyze_muonMC_dxy.root'
inFile = R.TFile(fileName)
tree = inFile.Get('analyzer/HighPTTrackTree')
#outputFilename = 'hists_curvature_muonMC_f'+f+'.root'
outputFilename = 'hists_curvature_muonMC_bestPull_test_pt_gt600.root'
#outputFilename = 'hists_curvature_dxy.root'
outFile = R.TFile(outputFilename,'recreate')

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

HISTS = {tracktype:{fit:{} for fit in fits} for tracktype in tracktypes}

for tracktype in tracktypes:
	for fit in fits:
		htitle = tracktype+fit
		HISTS[tracktype][fit]['K_gen_res'] = R.TH1D(htitle+'_K_gen_res','',*varList['K_gen_res']['bins'])
		HISTS[tracktype][fit]['K_gen_res_vs_gen_pt'] = R.TH2D(htitle+'_K_gen_res_vs_gen_pt','',\
				*(varList['gen_pt']['bins']+varList['K_gen_res']['bins']))
		HISTS[tracktype][fit]['K_gen_res_vs_gen_eta'] = R.TH2D(htitle+'_K_gen_res_vs_gen_eta','',\
				*(varList['gen_eta']['bins']+varList['K_gen_res']['bins']))
		

for l,entry in enumerate(tree):
	#if l>2000: break
	sys.stdout.write('Entry {l:>5}\r'.format(**locals()))
	p = t.pyTree(tree)
	if p.gen_pt < 600: continue

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
			#print p.gen_pt, (track.K()-p.gen_K)/p.gen_K, trackType+fit
			HISTS[trackType][fit]['K_gen_res'].Fill((track.K()-p.gen_K)/p.gen_K)
			HISTS[trackType][fit]['K_gen_res_vs_gen_pt'].Fill(p.gen_pt, (track.K()-p.gen_K)/p.gen_K)
			HISTS[trackType][fit]['K_gen_res_vs_gen_eta'].Fill(p.gen_eta, (track.K()-p.gen_K)/p.gen_K)

outFile.cd()
for tracktype in tracktypes:
	for fit in fits:
		HISTS[tracktype][fit]['K_gen_res'].Write()
		HISTS[tracktype][fit]['K_gen_res_vs_gen_pt'].Write()
		HISTS[tracktype][fit]['K_gen_res_vs_gen_eta'].Write()
outFile.Close()
