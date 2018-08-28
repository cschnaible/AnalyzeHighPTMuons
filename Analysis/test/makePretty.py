import sys
import ROOT as R
import math as math
import numpy as np
import argparse
from HighPTMuons.Analysis.ResolutionAndPullFunctions import *
R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)
parser = argparse.ArgumentParser()
parser.add_argument('-in','--infile',dest='infile',default='ZpMM',help='Which MC file to use')
parser.add_argument('-out','--outfile',dest='outfile',default='ZpMM',help='Which MC file to use')
args = parser.parse_args()

infile = R.TFile(args.infile)
outfile = R.TFile.Open(args.outfile,'RECREATE')

gen1Dplots  = ['genCurvRes', 'genCurvPull']#, 'genCurvDiff', 'genCurvAbsDiff']
gen2Dplots  = ['genCurvResVsPt','genCurvPullVsPt','genCurvResVsEta','genCurvPullVsEta']#,
#				'genCurvDiffVsPt','genCurvAbsDiffVsPt','genCurvDiffVsEta','genCurvAbsDiffVsEta']
full1Dplots = ['fullCurvRes', 'fullCurvPull', 'fullCurvDiff', 'fullCurvAbsDiff']
full2Dplots = ['fullCurvResVsPt','fullCurvPullVsPt','fullCurvResVsEta','fullCurvPullVsEta',
				'fullCurvDiffVsPt','fullCurvAbsDiffVsPt','fullCurvDiffVsEta','fullCurvAbsDiffVsEta']
cutdict = {
		'eta':['all','barrel','endcap','fwdendcap'],
		'pt':['all','200','600','1000'],# 400, 800,
		}
cuts1D = [ptcut+'_'+etacut for etacut in cutdict['eta'] for ptcut in cutdict['pt']]

allTrackTypes = ['tracker',\
		'global','globalComb','globalMuonOnlyUpdate','globalTrackRankComb','globalTrackRankMuonOnlyUpdate',\
		'tuneP','tunePComb','tunePMuonOnlyUpdate','tunePTrackRankComb','tunePTrackRankMuonOnlyUpdate',\
		'picky','pickyComb','pickyMuonOnlyUpdate','pickyTrackRankComb','pickyTrackRankMuonOnlyUpdate',\
		'dyt','dytComb','dytMuonOnlyUpdate','dytTrackRankComb','dytTrackRankMuonOnlyUpdate',\
		'tpfms','tpfmsComb','tpfmsMuonOnlyUpdate']

hists = {}
for trackName in allTrackTypes:
	hists[trackName] = {}
	for gen1Dplot in gen1Dplots:
		hists[trackName][gen1Dplot] = {}
		for cut in cuts1D:
			hists[trackName][gen1Dplot][cut] = {}
			hname = trackName+'_'+gen1Dplot+'_'+cut
			print hname
			hists[trackName][gen1Dplot][cut]['plot'] = infile.Get(hname).Clone()
	for gen2Dplot in gen2Dplots:
		hists[trackName][gen2Dplot] = {}

		if 'Pt' in gen2Dplot:
			for cut in cutdict['eta']:
				hname = trackName+'_'+gen2Dplot+'_'+cut
				print hname
				hists[trackName][gen2Dplot][cut] = {}
				hists[trackName][gen2Dplot][cut]['plot'] = infile.Get(hname).Clone()
				hists[trackName][gen2Dplot][cut]['fitresults'] = R.TObjArray()
				if 'genCurvPull' in gen2Dplot:
					func = cruijff('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-5,5)
				elif 'genCurvRes' in gen2Dplot:
					func = cruijff('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-1,1)
				elif 'genCurvDiff' in gen2Dplot:
					func = gaus('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-2,2)
					func.labels = [('Normalization','norm'),('diff #mu','mu'),('diff #sigma','sigma'),('fit #chi^{2}','chi2')]
				else: continue
				hists[trackName][gen2Dplot][cut]['plot'].FitSlicesY(\
						func.func,0,hists[trackName][gen2Dplot][cut]['plot'].GetNbinsX()+1,0,'QNR',\
						hists[trackName][gen2Dplot][cut]['fitresults'])
				hists[trackName][gen2Dplot][cut]['ffunc'] = func

		if 'Eta' in gen2Dplot:
			for cut in cutdict['pt']:
				hname = trackName+'_'+gen2Dplot+'_'+cut
				print hname
				hists[trackName][gen2Dplot][cut] = {}
				hists[trackName][gen2Dplot][cut]['plot'] = infile.Get(hname).Clone()
				hists[trackName][gen2Dplot][cut]['fitresults'] = R.TObjArray()
				if 'genCurvPull' in gen2Dplot:
					func = cruijff('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-5,5)
				elif 'genCurvRes' in gen2Dplot:
					func = cruijff('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-1,1)
				elif 'genCurvDiff' in gen2Dplot:
					func = gaus('ffunc_'+trackName+'_'+gen2Dplot+'_'+cut,-2,2)
					func.labels = [('Normalization','norm'),('diff #mu','mu'),('diff #sigma','sigma'),('fit #chi^{2}','chi2')]
				else: continue
				hists[trackName][gen2Dplot][cut]['plot'].FitSlicesY(\
						func.func,0,hists[trackName][gen2Dplot][cut]['plot'].GetNbinsX()+1,0,'QNR',\
						hists[trackName][gen2Dplot][cut]['fitresults'])
				hists[trackName][gen2Dplot][cut]['ffunc'] = func


toPlot = {
#		'trkCurv':   {'std':'tuneP','comb':'trkCurvComb','vtxmu':'trkCurvMuonOnlyUpdate','tracker':'tracker'},
#		'tunePCurv': {'std':'tuneP','comb':'tunePCurvComb','vtxmu':'tunePCurvMuonOnlyUpdate','tracker':'tracker'},
#		'dxy':       {'std':'tuneP','comb':'dxyComb','vtxmu':'dxyMuonOnlyUpdate','tracker':'tracker'},
#		'curvRelErr':{'std':'tuneP','comb':'curvRelErrComb','vtxmu':'curvRelErrMuonOnlyUpdate','tracker':'tracker'},
#		'trackRank': {'std':'tuneP','comb':'trackRankComb','vtxmu':'trackRankMuonOnlyUpdate','tracker':'tracker'},
		'tuneP': {'std':'tuneP','comb':'tunePComb','vtxmu':'tunePMuonOnlyUpdate','tracker':'tracker'},
		'global': {'std':'global','comb':'globalComb','vtxmu':'globalMuonOnlyUpdate','tracker':'tracker'},
		'picky': {'std':'picky', 'comb':'pickyComb','vtxmu':'pickyMuonOnlyUpdate','tracker':'tracker'},
		'dyt': {'std':'dyt','comb':'dytComb','vtxmu':'dytMuonOnlyUpdate','tracker':'tracker'},
		'tpfms': {'std':'tpfms','comb':'tpfmsComb','vtxmu':'tpfmsMuonOnlyUpdate','tracker':'tracker'},
		'tunePTrackRank': {'std':'tuneP','comb':'tunePTrackRankComb','vtxmu':'tunePMuonOnlyUpdate','tracker':'tracker'},
		'globalTrackRank': {'std':'global','comb':'globalTrackRankComb','vtxmu':'globalTrackRankMuonOnlyUpdate','tracker':'tracker'},
		'pickyTrackRank': {'std':'picky', 'comb':'pickyTrackRankComb','vtxmu':'pickyTrackRankMuonOnlyUpdate','tracker':'tracker'},
		'dytTrackRank': {'std':'dyt','comb':'dytTrackRankComb','vtxmu':'dytTrackRankMuonOnlyUpdate','tracker':'tracker'},
		#'refits':     {'dxy':'dxyComb','trkCurv':'trkCurvComb','trackRank':'trackRankComb','std':'tuneP'},
		}
prettyLeg = {
	'tuneP':'tuneP track',
	'global':'Global track',
	'picky':'Picky track',
	'dyt':'DYT track',
	'tpfms':'TPFMS track',
	'tracker':'Tracker track',
	'trkCurvComb':'Tracker + mu-only vtx track',
	'tunePCurvComb':'Tracker + mu-only vtx track',
	'dxyComb':'Tracker + mu-only vtx track',
	'trackRankComb':'Tracker + mu-only vtx track',
	'curvRelErrComb':'Tracker + mu-only vtx track',
	'globalComb':'Tracker + mu-only vtx track',
	'pickyComb':'Tracker + mu-only vtx track',
	'dytComb':'Tracker + mu-only vtx track',
	'tpfmsComb':'Tracker + mu-only vtx track',
	'trkCurvMuonOnlyUpdate':'mu-only vtx track',
	'tunePCurvMuonOnlyUpdate':'mu-only vtx track',
	'dxyMuonOnlyUpdate':'mu-only vtx track',
	'trackRankMuonOnlyUpdate':'mu-only vtx track',
	'curvRelErrMuonOnlyUpdate':'mu-only vtx track',
	'globalMuonOnlyUpdate':'mu-only vtx track',
	'pickyMuonOnlyUpdate':'mu-only vtx track',
	'dytMuonOnlyUpdate':'mu-only vtx track',
	'tpfmsMuonOnlyUpdate':'mu-only vtx track',
	'tunePComb':'Tracker + mu-only vtx track',
	'tunePMuonOnlyUpdate':'mu-only vtx track',
	'globalTrackRankMuonOnlyUpdate':'mu-only vtx track',
	'pickyTrackRankMuonOnlyUpdate':'mu-only vtx track',
	'dytTrackRankMuonOnlyUpdate':'mu-only vtx track',
	'tunePTrackRankMuonOnlyUpdate':'mu-only vtx track',
	'globalTrackRankComb':'Tracker + mu-only vtx track',
	'pickyTrackRankComb':'Tracker + mu-only vtx track',
	'dytTrackRankComb':'Tracker + mu-only vtx track',
	'tunePTrackRankComb':'Tracker + mu-only vtx track',
	}
prettyTitle = {
	'trkCurv':'Tracker #kappa pull',   
	'tunePCurv':'tuneP #kappa pull',
	'tuneP':'tuneP Track',
	'dxy':'mu-only d_{xy} pull',
	'trackRank':'Track Rank',
	'curvRelErr':'|#sigma(#kappa)/#kappa|',
	'global':'Global Track',
	'picky':'Picky Track',
	'dyt':'DYT Track',
	'tpfms':'TPFMS Track',
	'globalTrackRank':'Track Rank using only muon hits from Global track',
	'pickyTrackRank':'Track Rank using only muon hits from Picky track',
	'dytTrackRank':'Track Rank using only muon hits from DYT track',
	'tunePTrackRank':'Track Rank using only muon hits tuneP track',
	}

# tracks to draw together
allfits = ['tracker','vtxmu','comb','std']
colors = {
		'std':R.kBlue,
		'comb':R.kViolet,
		'vtxmu':R.kGreen+3,
		'tracker':R.kRed,
		}

def prettyCut(cut):
	cuts = cut.split('_')
	pt = cuts[0]
	eta = cuts[1]
	ptcut = ''
	etacut = ''
	if pt!='all': ptcut += 'p_{T} > '+str(pt)+' GeV'
	if eta=='barrel': etacut += '|#eta| < 1.2'
	elif eta=='endcap': etacut += '|#eta| > 1.2'
	elif eta=='fwdendcap': etacut += '|#eta| > 2.1'
	cuttitle = ''
	if   eta=='all' and pt!='all': cuttitle+= ' : '+ptcut
	elif pt=='all'  and eta!='all': cuttitle+=' : '+etacut
	elif pt!='all'  and eta!='all': cuttitle += ' : '+ptcut+' && '+etacut
	return cuttitle

def draw1Dplot(plotKey,plot,cut):
	outfile.cd()
	c = R.TCanvas()
	leg = R.TLegend(0.15,0.6,0.44,0.8)
	fits = list(allfits)
	maxlist = []
	for f,fit in enumerate(fits):
		trackType = toPlot[plotKey][fit]
		hist = hists[trackType][plot][cut]['plot']
		hist.SetMinimum(0)
		if plot=='genCurvRes' or plot=='fullCurvRes':
			hist.GetXaxis().SetRangeUser(-0.5,0.5)
		hist.SetLineColor(colors[fit])
		hist.SetLineWidth(2)
		hist.SetMarkerColor(colors[fit])
		draw = '' if f==0 else 'sames'
		hist.Draw(draw)
		leg.AddEntry(hist,prettyLeg[toPlot[plotKey][fit]],'l')
		maxlist.append(hist.GetMaximum())
		R.gPad.Update()
	leg.Draw()
	for f,fit in enumerate(fits):
		trackType = toPlot[plotKey][fit]
		hist = hists[trackType][plot][cut]['plot']
		hist.SetMaximum(max(maxlist)*1.1)
		statbox = hist.GetListOfFunctions().FindObject('stats')
		statbox.SetOptStat(1110)
		statbox.SetTextColor(colors[fit])
		statbox.SetX1NDC(statbox.GetX1NDC())
		statbox.SetX2NDC(statbox.GetX2NDC())
		y = 0.7 - f*0.15  - (0.025*f if f!=0 else 0.0)
		statbox.SetY1NDC(y)
		statbox.SetY2NDC(y + 0.15)
		statbox.Draw()
	c.Write(plotKey+'_'+plot+'_'+cut)

def draw2Dplot(plotKey,plot,cut):
	if 'full' in plot: fits = [fit for fit in allfits if fit not in ('vtxmu','std','tracker')]
	else: fits = list(allfits)
	outfile.cd()
	tmpTrackName = toPlot[plotKey][fits[0]]
	fitresults = len(hists[tmpTrackName][plot][cut]['fitresults'])
	for fitresult in range(fitresults):
		maxlist = []
		leg = R.TLegend(0.2,0.6,0.5,0.8)
		c = R.TCanvas()
		print plotKey, plot, cut, fitresult
		for f,fit in enumerate(fits):
			trackType = toPlot[plotKey][fit]
			hist = hists[trackType][plot][cut]['fitresults'][fitresult]
			hist.SetMinimum(0)
			if plot=='genCurvRes' or plot=='fullCurvRes':
				hist.GetXaxis().SetRangeUser(-0.5,0.5)
			hist.SetLineColor(colors[fit])
			hist.SetLineWidth(2)
			hist.SetMarkerColor(colors[fit])
			hist.GetXaxis().SetTitle(hists[trackType][plot][cut]['plot'].GetXaxis().GetTitle())
			func = hists[trackType][plot][cut]['ffunc']
			hist.GetYaxis().SetTitle(func.labels[fitresult][0])
			hist.SetTitle(prettyTitle[plotKey])
			hist.SetStats(R.kFALSE)
			draw = '' if f==0 else 'sames'
			hist.Draw(draw)
			leg.AddEntry(hist,prettyLeg[toPlot[plotKey][fit]],'lep')
			maxlist.append(hist.GetMaximum())
			R.gPad.Update()
		for f,fit in enumerate(fits):
			trackType = toPlot[plotKey][fit]
			hist = hists[trackType][plot][cut]['fitresults'][fitresult]
			hist.SetMaximum(max(maxlist)*1.1)
		leg.Draw()
		c.Write(plotKey+'_'+plot+'_'+cut+'_'+hists[trackType][plot][cut]['ffunc'].labels[fitresult][1])
			
for trackType in toPlot.keys():
	for plot in gen1Dplots:  
		for cut in cuts1D:
			print trackType,plot,cut
			draw1Dplot(trackType,plot,cut)
	for plot in gen2Dplots:  
		if 'Pt' in plot:
			for cut in cutdict['eta']:
				print trackType,plot,cut
				draw2Dplot(trackType,plot,cut)
		if 'Eta' in plot:
			for cut in cutdict['pt']:
				print trackType,plot,cut
				draw2Dplot(trackType,plot,cut)
#	if trackType in ['global','picky','dyt','tpfms']:
#		# plot only mu-only + vtx combination?
#		for plot in full1Dplots: 
#			for cut in cuts1D:
#				draw1Dplot(trackType,plot,cut)
#		for plot in full2Dplots:
#			if 'Pt' in plot:
#				for cut in cutdict['eta']:
#					draw2Dplot(trackType,plot,cut)
#			if 'Eta' in plot:
#				for cut in cutdict['pt']:
#					draw2Dplot(trackType,plot,cut)
#
