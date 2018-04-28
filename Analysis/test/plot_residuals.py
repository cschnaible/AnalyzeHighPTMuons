import sys
import ROOT as R
import math as math
import numpy as np
from scipy.stats import norm
import argparse
import tools as t

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',default='ZpMM',help='Which MC file to use')
parser.add_argument('-ex','--extra',default='',help='Extra name to output histogram file')
parser.add_argument('-r','--recreate',default=False,action='store_true',help='Recreate histogram file')
args = parser.parse_args()

fileName = 'AnalyzeTracks_'+str(args.file)+'.root'
nametmp = 'hists_residuals_'+str(args.file)+('' if args.extra=='' else '_'+str(args.extra))
outfilebasename = nametmp+'_base.root'
outfilename     = nametmp+'.root'
#outFile = R.TFile(outfilename,'recreate')
RECREATE = args.recreate

#dxy = R.TH1D('dxy','',100,-10,10)

if RECREATE:
	rFile = R.TFile(fileName)

	tree = rFile.Get('t')
	#outFileBase = R.TFile(outfilebasename,'recreate')

	for e, entry in enumerate(tree):

		p = t.pyTree(tree)
		if not p.allOkay: continue
		muon = t.Track('refit_noUpdate',p)
		full = t.Track('global',p)
		tracker = t.Track('tracker',p)
		comb = t.Track('comb',p)
		refit = t.Track('refit',p)

		#res_par,res_cov = t.kalman_prefit_residuals(muon._par,muon._cov,full._par,full._cov)
		#print res_par[2], full.dxy()-muon.dxy()
		#dxy.Fill(res_par[2])
		refit_par,refit_cov = t.kalman_filter_update(muon._par,muon._cov,full._par,full._cov)
		print 'refit'
		print refit_par[0],refit.K()
		print refit_par[1],refit.Lambda()
		print refit_par[2],refit.phi()
		print refit_par[3],refit.dxy()
		print refit_par[4],refit.dsz()
		comb_par,comb_cov = t.chisq_comb(refit_par,refit_cov,tracker._par,tracker._cov)
		print 'comb'
		print comb_par[0],comb.K()
		print comb_par[1],comb.Lambda()
		print comb_par[2],comb.phi()
		print comb_par[3],comb.dxy()
		print comb_par[4],comb.dsz()

		print
		break

#	outFile.cd()
#	dxy.Write()


