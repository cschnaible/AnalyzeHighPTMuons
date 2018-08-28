import sys
import ROOT as R
import math as math
import numpy as np
import argparse
from HighPTMuons.Analysis.ResolutionAndPullFunctions import *
R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)
R.gStyle.SetOptFit(1)
parser = argparse.ArgumentParser()
parser.add_argument('-in','--infile',dest='infile',default='ZpMM',help='Which MC file to use')
parser.add_argument('-out','--outfile',dest='outfile',default='ZpMM',help='Which MC file to use')
args = parser.parse_args()

infile = R.TFile(args.infile)
outfile = R.TFile.Open(args.outfile,'RECREATE')

trackNameBase = ['tracker','tuneP','picky','dyt','tpfms',\
		'globalTrackRank','pickyTrackRank','dytTrackRank','tunePTrackRank']
trackNameExtra = ['','MuonOnly','MuonOnlyUpdate','Comb']
cutdict = {
		'eta':['all','barrel','endcap','fwdendcap'],
		'pt':['all','200','600','1000'],# 400, 800,
		}
gen2Dplots  = ['genCurvResVsPt','genCurvResVsEta','genCurvPullVsPt','genCurvPullVsEta']
#				'genCurvDiffVsPt','genCurvDiffVsEta','genCurvAbsDiffVsPt','genCurvAbsDiffVsEta']

def plotHistSliceFits(trackName,plot,cut):
	histname = trackName+'_'+plot+'_'+cut
	print histname
	infile.cd()
	hist2d = infile.Get(histname).Clone() 
	outfile.cd()
	outfile.mkdir(histname)
	outfile.cd(histname)
	for ibin in range(1,hist2d.GetNbinsX()+1):
		c = R.TCanvas()
		name = histname+'_'+str(hist2d.GetXaxis().GetBinLowEdge(ibin))+'_'+str(hist2d.GetXaxis().GetBinLowEdge(ibin+1))
		name = name.replace('-','n')
		name = name.replace('.','')
		#ffunc = gaus('ffunc_'+name,-0.1,0.1)
		#ffunc = gausSum('ffunc_'+name,-0.5,0.5)
		ffunc = cruijff('ffunc_'+name,-1,1) 
		hist = hist2d.ProjectionY(name,ibin,ibin+1)
		hist.SetStats(1)
		hist.Fit(ffunc.GetName(),'QNR')
		hist.Draw()
		ffunc.Draw('same')
		hist.SetMinimum(0)
		c.Write(name+'_proj')

for base in trackNameBase:
	for extra in trackNameExtra:
		trackName = base+extra
		if 'TrackRank' in base and extra=='': continue
		if base=='tracker' and extra!='': continue
		for plot in gen2Dplots:
			for cuttype in ['eta','pt']:
				for cut in cutdict[cuttype]:
					if 'Pt' in plot and cuttype=='pt': continue
					if 'Eta' in plot and cuttype=='eta': continue
					plotHistSliceFits(trackName,plot,cut)
