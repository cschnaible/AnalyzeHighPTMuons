import ROOT as R
import tools as t
import math as math

#f = R.TFile('testing_again.root')
f = R.TFile('AnalyzeTracks_muonPT.root')
tree = f.Get('t')

fout = R.TFile('simple_res_plot.root','recreate')

combreshist = R.TH2F('combreshist',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
combreshist2 = R.TH2F('combreshist2',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
combreshist3 = R.TH2F('combreshist3',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
combreshist4 = R.TH2F('combreshist4',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
glbreshist = R.TH2F('glbreshist',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
glbreshist3 = R.TH2F('glbreshist3',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)
glbreshist4 = R.TH2F('glbreshist4',';gen p_{T} [GeV];#frac{#kappa - #kappa_{gen}}{#kappa_{gen}}',150,0,1500,100,-0.5,0.5)

n = 0
for e,entry in enumerate(tree):
	p = t.pyTree(tree)
	glb_comb = t.Track('global_comb',p)
	glb = t.Track('global',p)
	tracker = t.Track('tracker',p)
	glb_refit_muOnly = t.Track('global_refit_noUpdate',p)
	glb_refit_update = t.Track('global_refit',p)
	
	if len(glb_refit_muOnly._par)<5: continue
	#if abs(glb_refit_muOnly.dxy())<2.0:
	if math.sqrt(glb_refit_muOnly.cov(0,0))/abs(glb_refit_muOnly.K())<0.5:
		combreshist2.Fill(p.gen_pt,(glb_comb.K()-p.gen_K)/p.gen_K)
		combreshist3.Fill(p.gen_pt,(glb_comb.K()-p.gen_K)/p.gen_K)
		glbreshist3.Fill(p.gen_pt,(glb.K()-p.gen_K)/p.gen_K)
	#elif abs(glb_refit_muOnly.dxy())>=2.0:
	elif math.sqrt(glb_refit_muOnly.cov(0,0))/abs(glb_refit_muOnly.K())>=0.5:
		combreshist2.Fill(p.gen_pt,(tracker.K()-p.gen_K)/p.gen_K)
		combreshist4.Fill(p.gen_pt,(tracker.K()-p.gen_K)/p.gen_K)
		glbreshist4.Fill(p.gen_pt,(glb.K()-p.gen_K)/p.gen_K)
	else:
		print 'why',abs(glb_refit_muOnly.dxy())
	glbreshist.Fill(p.gen_pt,(glb.K()-p.gen_K)/p.gen_K)
	combreshist.Fill(p.gen_pt,(glb_comb.K()-p.gen_K)/p.gen_K)
	#if abs(p.gen_pt-1200)<100 and abs(glb_refit_muOnly.dxyBS())>7:
	#	print entry.run, entry.lumi, entry.event

fout.cd()

for hist in [combreshist,combreshist2,combreshist3,combreshist4,glbreshist,glbreshist3,glbreshist4]:
	hist.Write()
	harr = R.TObjArray()
	func = R.TF1('name','gaus',-0.5,0.5)
	hist.FitSlicesY(func,0,-1,0,'QNR',harr)
	for h in harr:
		h.Write()
