import ROOT as R
import numpy as np

fin = open('muonRefitLog_chi2andprob','read')
fout = R.TFile('chi2ProbPlots.root','recreate')

hists = {}
cutvals = range(1,12)
plots = ['12','23','34','45','56','67','78','89']

Pbins = np.logspace(-3,1.75,101)
Phists = {plot:R.TH2D('prob_h'+plot,'-ln(Prob('+plot[1]+' hit)) vs. -ln(Prob('+plot[0]+' hit));-ln(Prob('+plot[0]+' hit));-ln(Prob('+plot[1]+' hit))',100,Pbins,100,Pbins) for plot in plots}

chi2bins = np.logspace(-2,2,101)
Chi2histslog = {plot:R.TH2D('chi2_h'+plot+'_log','#chi^{2}_{'+plot[1]+'} vs. #chi^{2}_{'+plot[0]+'};#chi^{2}_{'+plot[0]+'};#chi^{2}_{'+plot[1]+'}',100,chi2bins,100,chi2bins) for plot in plots}
Chi2hists = {plot:R.TH2D('chi2_h'+plot,'#chi^{2}_{'+plot[1]+'} vs. #chi^{2}_{'+plot[0]+'};#chi^{2}_{'+plot[0]+'};#chi^{2}_{'+plot[1]+'}',100,0,50,100,0,50) for plot in plots}
dChi2hists = {plot:R.TH2D('delta_chi2_h'+plot,'#chi^{2}_{'+plot[0]+'} vs. #Delta#chi^{2}_{'+plot[0]+','+plot[1]+'};#Delta#chi^{2}_{'+plot[0]+','+plot[1]+'};#chi^{2}_{'+plot[0]+' hit}',100,0,50,100,0,50) for plot in plots}
for l,line in enumerate(fin):
	cols = line.strip('\n').split()
	name = cols[0]
	vals = cols[1:]
	nHits = len(vals)
	if nHits==1: continue
	if name=='-ln(P)':
		for v,val in enumerate(vals):
			if v>=nHits-1: continue
			Pl = float(val)
			Pu = float(vals[v+1])
			if Pl==-1 or Pu==-1: continue
			Phists[str(v+1)+str(v+2)].Fill(Pl,Pu)
	if name=='chi2':
		for v,val in enumerate(vals):
			if v>=nHits-1: continue
			chi2l = float(val)
			chi2u = float(vals[v+1])
			if chi2l==-1 or chi2u==-1: continue
			Chi2histslog[str(v+1)+str(v+2)].Fill(chi2l,chi2u)
			Chi2hists[str(v+1)+str(v+2)].Fill(chi2l,chi2u)
			dChi2hists[str(v+1)+str(v+2)].Fill(chi2u-chi2l,chi2l)

fout.cd()
for plot in plots:
	Phists[plot].Write()
	Chi2hists[plot].Write()
	dChi2hists[plot].Write()
	Chi2histslog[plot].Write()
fout.Close()
