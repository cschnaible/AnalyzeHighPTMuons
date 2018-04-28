import ROOT as R
import numpy as np
import math as math
import tools as t

fin = R.TFile('AnalyzeTracks_muonPT.root')
tree = fin.Get('t')
step = 0.1
uplim = 2.0
slist = np.arange(0,uplim+step,step)
title = 'Comb = %(s)s*ref + trk;#frac{#kappa-#kappa_{gen}}{#kappa_{gen}};Counts'
reshists = {s:{R.TH1D('comb_scale_'+str(s),title%(str(s)),200,-1,1)} for s in slist}

for e,entry in enumerate(tree):
	p = t.PyTree(tree)
	glb = t.Track('global',p)
	trk = t.Track('tracker',p)
	glb_ref = t.Track('global_refit',p)
	glb_comb = t.Track('global_comb',p)

	for s in slist:
		par,cov = t.chisq_comb(glb_ref.par(),s*glb_ref._cov,trk.par(),trk._cov)



