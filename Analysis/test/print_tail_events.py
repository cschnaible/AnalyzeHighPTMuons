import ROOT as R
import numpy as np
import tools as t

fin = R.TFile('testing_again.root','read')
data = R.TFile('simple_res_plot.root','read')

glbreshist = data.Get('glbreshist_2')
combreshist = data.Get('combreshist_2')


tree = fin.Get('t')
ncombbad_1 = 0
ncombbad_2 = 0
ncombbad_3 = 0
ncombbad_4 = 0
ncombbad_5 = 0
nglbbad_1 = 0
nglbbad_2 = 0
nglbbad_3 = 0
nglbbad_4 = 0
nglbbad_5 = 0

for e,entry in enumerate(tree):
	p = t.pyTree(tree)
	glb = t.Track('global',p)
	glb_comb = t.Track('global_comb',p)
	glb_muOnly = t.Track('global_refit_noUpdate',p)
	glb_refit = t.Track('global_refit',p)
	tracker = t.Track('tracker',p)
	if len(glb_comb.par())<5: continue

	glbres = (glb.K()-p.gen_K)/p.gen_K
	combres = (glb_comb.K()-p.gen_K)/p.gen_K

	ptres = combreshist.GetBinContent(combreshist.FindBin(p.gen_pt))

	#print p.gen_pt,ptres,glbres/ptres,combres/ptres
	
	if abs(combres/ptres) > 1: 
		ncombbad_1+=1
		if abs(combres/ptres) > 2: 
			ncombbad_2+=1
			if abs(combres/ptres) > 3:
				ncombbad_3+=1
				if abs(combres/ptres) > 4:
					ncombbad_4+=1
					if abs(combres/ptres) > 5:
						ncombbad_5+=1
	if abs(glbres/ptres) > 1: 
		nglbbad_1+=1
		if abs(glbres/ptres) > 2: 
			nglbbad_2+=1
			if abs(glbres/ptres) > 3:
				nglbbad_3+=1
				if abs(glbres/ptres) > 4:
					nglbbad_4+=1
					if abs(glbres/ptres) > 5:
						nglbbad_5+=1

n = float(tree.GetEntries())
print 'ncombbad',ncombbad_1/n,ncombbad_2/n,ncombbad_3/n,ncombbad_4/n,ncombbad_5/n
print 'nglbbad',nglbbad_1/n,nglbbad_2/n,nglbbad_3/n,nglbbad_4/n,nglbbad_5/n
