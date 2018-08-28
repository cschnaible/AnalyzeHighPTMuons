import ROOT as R

fin = open('muonRefitLog_30k','read')
fout = R.TFile('chi2plots.root','recreate')

hists = {nHit:R.TH1D('h'+str(nHit),'',100,0,25) for nHit in [1,2,3,4,5,6,7,8,9]}

for line in fin:
	cols = line.strip('\n').split()
	nHits=0
	for i in cols[2]:
		nHits+=int(i)
	hists[nHits].Fill(float(cols[5]))

fout.cd()
for nHit in hists.keys():
	hists[nHit].Write()
fout.Close()
