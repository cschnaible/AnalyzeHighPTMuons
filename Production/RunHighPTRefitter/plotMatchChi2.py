import ROOT as R

fin = open('log_res_matchchi2','read')
fout = R.TFile('matchChi2plots.root','recreate')

hists = {}
cuttypes = ['gt','lt']
cutvals = [1,2,3,4,5]
allhist = R.TH2D('hAll',';#kappa-#kappa_{gen}/#sigma_{#kappa};match #chi^{2}',100,0,10,100,0,10)
for cut in cuttypes:
	hists[cut] = {}
	for val in cutvals:
		hists[cut][val] = R.TH1D('h_'+cut+str(val),'',100,0,10)

for line in fin:
	cols = line.strip('\n').split()
	if cols[0]!='FORPLOT':continue
	print cols
	nHits=0
	allhist.Fill(float(cols[1]),float(cols[2]))
	
	cuts = {val:'gt' if float(cols[1])>val else 'lt' for val in cutvals}
	for val in cuts.keys():
		hists[cuts[val]][val].Fill(float(cols[2]))

fout.cd()
for cut in cuttypes:
	for val in cutvals:
		hists[cut][val].Write()
fout.Close()
