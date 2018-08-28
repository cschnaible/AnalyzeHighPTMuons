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



gen1Dplots  = ['genCurvRes', 'genCurvPull', 'genCurvDiff', 'genCurvAbsDiff']
gen2Dplots  = ['genCurvResVsPt','genCurvPullVsPt','genCurvResVsEta','genCurvPullVsEta',
				'genCurvDiffVsPt','genCurvAbsDiffVsPt','genCurvDiffVsEta','genCurvAbsDiffVsEta']
full1Dplots = ['fullCurvRes', 'fullCurvPull', 'fullCurvDiff', 'fullCurvAbsDiff']
full2Dplots = ['fullCurvResVsPt','fullCurvPullVsPt','fullCurvResVsEta','fullCurvPullVsEta',
				'fullCurvDiffVsPt','fullCurvAbsDiffVsPt','fullCurvDiffVsEta','fullCurvAbsDiffVsEta']

# tracks to draw together
allfits = ['tracker','vtxmu','comb','std']
colors = {
		'std':R.kBlue,
		'comb':R.kViolet,
		'vtxmu':R.kGreen+3,
		'tracker':R.kRed,
		}

cutdict = {
		'eta':['all','barrel','endcap','fwdendcap'],
		'pt':['all','200','600','1000'],# 400, 800,
		}

cuts1D = [ptcut+'_'+etacut for etacut in cutdict['eta'] for ptcut in cutdict['pt']]

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

def draw1Dplot(trackType,plot,cut):
	hists = {}
	fits = list(allfits)
	for fit in fits:
		print toPlot[trackType][fit]+'_'+plot+'_'+cut
		hists[fit] = infile.Get(toPlot[trackType][fit]+'_'+plot+'_'+cut).Clone()
		# shrink x-axis limits from [-5, 5] to [-0.25,0.25] 
		if plot=='genCurvRes' or plot=='fullCurvRes':
			hists[fit].GetXaxis().SetRangeUser(-0.5,0.5)
	c = R.TCanvas()
	for f,fit in enumerate(fits):
		hists[fit].SetLineWidth(2)
		hists[fit].SetLineColor(colors[fit])
		hists[fit].SetMarkerColor(colors[fit])
		hists[fit].SetTitle(prettyTitle[trackType]+prettyCut(cut))
		draw = '' if f==0 else 'sames'
		hists[fit].Draw(draw)
		R.gPad.Update()
	for f,fit in enumerate(fits):
		statbox = hists[fit].GetListOfFunctions().FindObject('stats')
		statbox.SetOptStat(1110)
		statbox.SetTextColor(colors[fit])
		statbox.SetX1NDC(statbox.GetX1NDC())
		statbox.SetX2NDC(statbox.GetX2NDC())
		y = 0.7 - f*0.15  - (0.025*f if f!=0 else 0.0)
		statbox.SetY1NDC(y)
		statbox.SetY2NDC(y + 0.15)
		statbox.Draw()
		
	hmax = max([hists[fit].GetMaximum() for fit in fits])
	hists[fits[0]].SetMaximum(hmax*1.1)
	hists[fits[0]].SetMinimum(0)
	leg = R.TLegend(0.15,0.6,0.44,0.8)
	for f,fit in enumerate(fits): leg.AddEntry(hists[fit],prettyLeg[toPlot[trackType][fit]],'l')
	leg.Draw()
	c.Write(trackType+'_'+plot+'_'+cut)

def draw2Dplot(trackType,plot,cut):
	if 'full' in plot: fits = [fit for fit in allfits if fit not in ('vtxmu','std','tracker')]
	else: fits = list(allfits)
	histSlices = {fit:R.TObjArray() for fit in fits}
	hists = {}
	for fit in fits:
		hists[fit] = infile.Get(toPlot[trackType][fit]+'_'+plot+'_'+cut).Clone()
		if 'genCurvPull' in plot:
			#func = gaus('ffunc_'+trackType+'_'+plot+'_'+cut,-2,2)
			#func.labels = [('Normalization','norm'),('pull #mu','mu'),('pull #sigma','sigma')]
			func = cruijff('ffunc_'+trackType+'_'+plot+'_'+cut,-5,5)
		elif 'genCurvRes' in plot:
			#func = gausSum('ffunc_'+trackType+'_'+plot+'_'+cut,-1,1)
			#func = gaus('ffunc_'+trackType+'_'+plot+'_'+cut,hists[fit],GetMean()-hists[fit].GetRMS(),GetMean()+hists[fit].GetRMS())
			func = cruijff('ffunc_'+trackType+'_'+plot+'_'+cut,-1,1)
		elif 'genCurvDiff' in plot:
			func = gaus('ffunc_'+trackType+'_'+plot+'_'+cut,-2,2)
			func.labels = [('Normalization','norm'),('diff #mu','mu'),('diff #sigma','sigma')]
		else: return
		print toPlot[trackType][fit]+'_'+plot+'_'+cut
		hists[fit].FitSlicesY(func.func,0,hists[fit].GetNbinsX()+1,0,'QNR',histSlices[fit])


	outfile.cd()
	for p in range(func.GetNpar()):
		c = R.TCanvas()
		for f,fit in enumerate(fits):
			histSlices[fit][p].SetLineWidth(2)
			histSlices[fit][p].SetLineColor(colors[fit])
			histSlices[fit][p].SetMarkerColor(colors[fit])
			histSlices[fit][p].GetXaxis().SetTitle(hists[fit].GetXaxis().GetTitle())
			histSlices[fit][p].GetYaxis().SetTitle(func.labels[p][0])
			histSlices[fit][p].SetTitle(prettyTitle[trackType])
			histSlices[fit][p].SetStats(R.kFALSE)
			#histSlices[fit][p].Write(hists[fit].GetName()+'_'+labels[p][1])
			draw = '' if f==0 else 'same'
			histSlices[fit][p].Draw(draw)
		hmax = max([histSlices[fit][p].GetMaximum() for fit in fits])
		histSlices[fits[0]][p].SetMaximum(hmax*1.1)
		if p!=1:
			histSlices[fits[0]][p].SetMinimum(0)
		leg = R.TLegend(0.2,0.6,0.5,0.8)
		for f,fit in enumerate(fits): leg.AddEntry(histSlices[fit][p],prettyLeg[toPlot[trackType][fit]],'lep')
		leg.Draw()
		c.Write(trackType+'_'+plot+'_'+cut+'_'+func.labels[p][1])

			
for trackType in toPlot.keys():
	for plot in gen1Dplots:  
		for cut in cuts1D:
			draw1Dplot(trackType,plot,cut)
	for plot in gen2Dplots:  
		if 'Pt' in plot:
			for cut in cutdict['eta']:
				draw2Dplot(trackType,plot,cut)
		if 'Eta' in plot:
			for cut in cutdict['pt']:
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
