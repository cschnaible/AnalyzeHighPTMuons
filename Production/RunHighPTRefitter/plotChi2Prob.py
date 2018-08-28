import ROOT as R
import numpy as np

fin = open('muonRefitLog_50K','read')
fout = R.TFile('chi2ProbPlots.root','recreate')

hists = {}
cutvals = range(1,12)
bins = np.logspace(-3,1.75,101)
allhist = R.TH2D('hAll',';N(hits);-ln(P)',11,0,11,100,bins)
hists = {}
for val in cutvals:
	hists[val] = R.TH1D('h_'+str(val),';-ln(P);',100,bins)

for l,line in enumerate(fin):
	cols = line.strip('\n').split()
	nHits=0
	for i in str(cols[2]):
		try:
			nHits+=int(i)
		except:
			print l
			exit()
	chi2Prob = float(cols[6])
	allhist.Fill(nHits,chi2Prob)
	hists[nHits].Fill(chi2Prob)

fout.cd()
allhist.Write()
for val in cutvals:
	hists[val].Write()
fout.Close()
