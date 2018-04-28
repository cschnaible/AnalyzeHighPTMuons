import ROOT as R
import tools as t

f = open('globalHitChi2','r')
rout = R.TFile('hit_chi2.root','recreate')

DT2 = R.TH1D('DT2','DT 2-dim RH;#chi^{2};N(hits)',200,0,20)
DT2prob = R.TH1D('DT2prob','DT 2-dim RH;#chi^{2} Probability;N(hits)',100,0,1)
DT4 = R.TH1D('DT4','DT 4-dim RH;#chi^{2};N(hits)',200,0,20)
DT4prob = R.TH1D('DT4prob','DT 4-dim RH;#chi^{2} Probability;N(hits)',100,0,1)
CSC4 = R.TH1D('CSC4','CSC 4-dim RH;#chi^{2};N(hits)',200,0,20)
CSC4prob = R.TH1D('CSC4prob','CSC 4-dim RH;#chi^{2} Probability;N(hits)',100,0,1)
RPCb = R.TH1D('RPCb','RPC barrel 2-dim RH;#chi^{2};N(hits)',200,0,20)
RPCbprob = R.TH1D('RPCbprob','RPC barrel 2-dim RH;#chi^{2} Probability;N(hits)',100,0,1)
RPCe = R.TH1D('RPCe','RPC endcap 2-dim RH;#chi^{2};N(hits)',200,0,20)
RPCeprob = R.TH1D('RPCeprob','RPC endcap 2-dim RH;#chi^{2} Probability;N(hits)',100,0,1)
allhitprob = R.TH1D('allhitprob','All hits;#chi^{2} Probability;N(hits)',100,0,1)

for line in f:
	cols = line.strip('\n').split()
	if len(cols)<1: continue
	det = str(cols[0])
	chi2 = float(cols[1])
	dim = int(cols[2])
	if det=='DT' and dim==2:
		DT2.Fill(chi2)
		DT2prob.Fill(R.TMath.Prob(chi2,1))
	elif det=='DT' and dim==4:
		DT4.Fill(chi2)
		DT4prob.Fill(R.TMath.Prob(chi2,3))
	elif det=='CSC' and dim==4:
		CSC4.Fill(chi2)
		CSC4prob.Fill(R.TMath.Prob(chi2,3))
	elif det=='RPCBarrel' and dim==2:
		RPCb.Fill(chi2)
		RPCbprob.Fill(R.TMath.Prob(chi2,1))
	elif det=='RPCEndcap' and dim==2:
		RPCe.Fill(chi2)
		RPCeprob.Fill(R.TMath.Prob(chi2,1))
	else:
		print det,chi2,dim
	allhitprob.Fill(R.TMath.Prob(chi2,dim-1))

rout.cd()
DT2.Write()
DT4.Write()
CSC4.Write()
RPCb.Write()
RPCe.Write()
DT2prob.Write()
DT4prob.Write()
CSC4prob.Write()
RPCbprob.Write()
RPCeprob.Write()
allhitprob.Write()
