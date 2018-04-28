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
nametmp = 'hists_'+str(args.file)+('' if args.extra=='' else '_'+str(args.extra))
outfilebasename = nametmp+'_base.root'
outfilename     = nametmp+'.root'
outFile = R.TFile(outfilename,'recreate')
RECREATE = args.recreate


styles = {
		'refit':{
			'color':R.kBlue,
			'name':'ref',
			'pretty':'#mu-only track + vtx pos & dir',
			},
		'global':{
			'color':R.kGreen+1,
			'name':'glb',
			'pretty':'Global track',
			},
		'tracker':{
			'color':R.kRed,
			'name':'trk',
			'pretty':'Tracker track',
			},
		'standAlone':{
			'color':R.kOrange,
			'name':'sa',
			'pretty':'Stand-alone track',
			},
		'refit_noUpdate':{
			'color':R.kMagenta+2,
			'name':'refNoUp',
			'pretty':'#mu-only track',
			},
		'comb':{
			'color':R.kMagenta+1,
			'name':'comb',
			'pretty':'Tracker track + #mu-only + vtx pos & dir track',
			},
		}

histDict = {
		'K':{
			'bins':(100,-0.01,0.01),
			'X':'#Kappa [Gev^{-1}]',
			'Y':'Counts',
			},
		'K_res':{
			'bins':(100,-0.5,0.5),
			'X':'#Kappa(track)-#Kappa(gen)/#Kappa(gen)',
			'Y':'Counts',
			},
		'K_b':{
			'bins':(100,-0.001,0.001),
			'X':'#Kappa(track)-#Kappa(gen)/#Kappa(gen)',
			'Y':'Counts',
			},
		'K_rel_err':{
			'bins':(100,0,3),
			'X':'#delta_{#Kappa}/|#Kappa|',
			'Y':'Counts',
			},
		'K_rel_err_vs_pt':{
			'bins':(150,0,1500,50,0,0.05),
			'X':'p_{T} [GeV]',
			'Y':'#delta_{#Kappa}/|#Kappa|',
			},
		'K_res_vs_pt':{
			'bins':(150,0,1500,50,-0.5,0.5),
			'X':'p_{T} [GeV]',
			'Y':'#Kappa(track)-#Kappa(gen)/#Kappa(gen)',
			},
		'K_rel_err_vs_pt2':{
			'bins':(50,0,10000000,50,0,0.05),
			'X':'p_{T}^{2} [GeV^{2}]',
			'Y':'#delta_{#Kappa}/|#Kappa|',
			},
		'K_res_vs_pt2':{
			'bins':(50,0,10000000,50,-0.5,0.5),
			'X':'p_{T}^{2} [GeV^{2}]',
			'Y':'#Kappa(track)-#Kappa(gen)/#Kappa(gen)',
			},
		'P_diff_vs_gen_P':{
			'bins':(75,0,2500,30,-1500,1500),
			'X':'Gen p_{T} [GeV]',
			'Y':'#Delta P [GeV]',
			},
		'K_res_vs_comb_chi2':{
			'bins':(25,0,25,50,-0.5,0.5),
			'X':'Combination #chi^{2}',
			'Y':'#Kappa(track)-#Kappa(gen)/#Kappa(gen)',
			},
		}

tracklist = ['global','comb','refit','tracker','standAlone','refit_noUpdate']

# K_res = K(track) - K(gen) / K(gen)
# K = q/P
# K_rel_err = sigma(K) / K
# K_b = K(track) - K(gen)
baselist = ['K_res','K','K_rel_err','K_b']
Kb_eta_list = ['minus','plus','mid_plus','mid_minus','fwd_plus','fwd_minus']
Kb_phi_list = ['all','1','2','3','4']

basehists = {track:{hist:{} for hist in baselist} for track in tracklist}
Kbhists = {track:{eta:{phi:{} for phi in Kb_phi_list} for eta in Kb_eta_list} for track in tracklist}

# 1/abs(K) = P
# K_rel_err_vs_pt  = sigma(K)/K vs. gen pt
# K_rel_err_vs_pt2 = sigma(K)/K vs. gen pt^2
# K_res_vs_pt  = K(track)-K(gen)/K(gen) vs. gen pt
# K_res_vs_pt2 = K(track)-K(gen)/K(gen) vs. gen pt^2
histlist2D = ['K_rel_err_vs_pt2','K_res_vs_pt2', 'K_rel_err_vs_pt','K_res_vs_pt','K_res_vs_comb_chi2']
hists2D = {track:{hist:{} for hist in histlist2D} for track in tracklist}

# frachist
# K(comb,f) = f*K(ref) + (1-f)*K(trk), sigma(K(comb,f)-K(gen) / K(gen)) vs. f

# P_diff_vs_gen_P
# P(1) - P(2) vs. P(gen)

# f1_f2_res
# K(comb,f1,f2) = (C(ref)^-1 + C(trk)^-1)^-1 * (f1 * C(ref)^-1 * P(ref) + f2 * C(trk)^-1 * P(trk))
# f1 vs. f2; sigma(K(comb,f1,f2)-K(gen)/K(gen))
difflist = ['ref_trk','trk_gen','trk_glb','ref_gen','ref_glb']
other = {
	'frachist':{},
	'P_diff_vs_gen_P':{diff:{} for diff in difflist},
	'michalis':{},
	'michalis2':{},
	}
fRvsfT = {}

if RECREATE:
	rFile = R.TFile(fileName)

	tree = rFile.Get('t')
	outFileBase = R.TFile(outfilebasename,'recreate')

	# Initialize histograms
	for track in tracklist:
		for hist in baselist:
			hname = hist+'_'+styles[track]['name']
			basehists[track][hist] = R.TH1D(hname,'',*histDict[hist]['bins'])
		for eta in Kb_eta_list:
			for phi in Kb_phi_list:
				hname = 'K_b_'+eta+'_'+phi+'_'+styles[track]['name']
				Kbhists[track][eta][phi] = R.TH1D(hname,'',*histDict[hist]['bins'])
		for hist in histlist2D:
			hname = hist+'_'+styles[track]['name']
			hists2D[track][hist] = R.TH2D(hname,'',*histDict[hist]['bins'])
	fracs = [i*0.05 for i in range(21)]
	other['frachist'] = R.TH2D('frachist','#Kappa_{comb} = f #times #Kappa_{refit} + (1-f) #times #Kappa_{tracker}',21,0,1.05,100,-0.5,0.5)
	for diff in difflist:
		other['P_diff_vs_gen_P'][diff] = R.TH2D('P_diff_vs_gen_P_'+diff,'',*histDict['P_diff_vs_gen_P']['bins'])
	other['michalis'] = R.TH2D('michalis',';1/p_{T,gen} [GeV^{-1}];#frac{1/p_{T,tracker}}{1/p_{T,gen}}',50,0,0.1,50,0,2)
	other['michalis2'] = R.TH2D('michalis2',';#kappa(gen) [GeV^{-1}];#kappa(tracker) / #kappa(gen)',50,0,0.1,50,0,2)

	'''
	# fR vs fT vs mu,sigma(K-Kgen/Kgen)
	step = 0.0002
	steps = np.arange(0.999,1.001+step,step)
	resdict = {fR:{fT:np.array([]) for fT in steps} for fR in steps}
	'''

	# loop on tree
	for l,entry in enumerate(tree):
		#if l>2000: continue
		sys.stdout.write('Entry {l:>5}\r'.format(**locals()))
		p = t.pyTree(tree)
		if not p.allOkay: continue
		for trackType in tracklist:
			track = t.Track(trackType,p)
			if track.cov(0,0)<0.:
				print "Negative covariance? WTF?"
				continue
			if track.eta()==-999.:
				print "track lambda is strange?"
				continue

			# Fill base histograms
			basehists[trackType]['K'].Fill(track.K())
			basehists[trackType]['K_res'].Fill((track.K()-p.gen_K)/p.gen_K)
			basehists[trackType]['K_b'].Fill(track.K()-p.gen_K)
			basehists[trackType]['K_rel_err'].Fill(math.sqrt(track.cov(0,0))/abs(track.K()))

			# Fill K_b plots
			etacuts = {
					'plus'     :track.eta()>   0.,
					'minus'    :track.eta()<   0.,
					'mid_plus' :track.eta()< 1.4 and track.eta()>0.,
					'mid_minus':track.eta()>-1.4 and track.eta()<0.,
					'fwd_plus' :track.eta()> 1.4,
					'fwd_minus':track.eta()<-1.4,
					}
			phicuts = {
					'all':True,
					'1':track.phi()<-math.pi/2.,
					'2':track.phi()>-math.pi/2 and track.phi()<0.,
					'3':track.phi()>0.         and track.phi()<math.pi/2,
					'4':track.phi()>math.pi/2,
					}
			for etacut in etacuts:
				for phicut in phicuts:
					if etacuts[etacut] and phicuts[phicut]:
						Kbhists[trackType][etacut][phicut].Fill(track.K()-p.gen_K)

			# Fill 2D histograms
			hists2D[trackType]['K_rel_err_vs_pt2'].Fill(p.gen_pt**2, track.cov(0,0)/(abs(track.K())**2))
			hists2D[trackType]['K_res_vs_pt2'].Fill(p.gen_pt**2, (track.K()-p.gen_K)/p.gen_K)
			hists2D[trackType]['K_rel_err_vs_pt'].Fill(p.gen_pt, track.cov(0,0)/(abs(track.K())**2))
			hists2D[trackType]['K_res_vs_pt'].Fill(p.gen_pt, (track.K()-p.gen_K)/p.gen_K)
			hists2D[trackType]['K_res_vs_comb_chi2'].Fill(p.comb_chi2,(track.K()-p.gen_K)/p.gen_K)

		refit = t.Track('refit',p)
		tracker = t.Track('tracker',p)
		glb = t.Track('global',p)
		# Fill res vs f histogram
		other['P_diff_vs_gen_P']['ref_trk'].Fill((1./abs(p.gen_K)),abs(1./refit.K())-abs(1./tracker.K()))
		other['P_diff_vs_gen_P']['ref_glb'].Fill((1./abs(p.gen_K)),abs(1./refit.K())-abs(1./glb.K()))
		other['P_diff_vs_gen_P']['ref_gen'].Fill((1./abs(p.gen_K)),abs(1./refit.K())-abs(1./p.gen_K))
		other['P_diff_vs_gen_P']['trk_gen'].Fill((1./abs(p.gen_K)),abs(1./tracker.K())-abs(1./p.gen_K))
		other['P_diff_vs_gen_P']['trk_glb'].Fill((1./abs(p.gen_K)),abs(1./tracker.K())-abs(1./glb.K()))
		for f in fracs:
			#Kcomb = f*p.refit_par[0] + (1.-f)*p.tracker_par[0]
			Kcomb = f*refit.K() + (1.-f)*tracker.K()
			other['frachist'].Fill(f,(Kcomb-p.gen_K)/p.gen_K)
		#print 1./p.gen_pt, (1./tracker.pt())/(1./p.gen_pt)
		if abs(p.gen_eta)>1.8:
			other['michalis'].Fill(1./p.gen_pt, (1./tracker.pt())/(1./p.gen_pt))
			other['michalis2'].Fill(p.gen_K , tracker.K()/p.gen_K)

		'''
		# do fR vs fT
		for fR in steps:
			for fT in steps:
				Ccov = (refit._w+tracker._w).getI()
				Cpar =  ((refit.par()*(fR*refit._w)) + (tracker.par()*(fT*tracker._w)))*Ccov
				# Cpar = [[K,lambda,phi,dxy,dsz]]
				# Cpar[0,0] = K
				res = (Cpar[0,0]-p.gen_K)/p.gen_K
				resdict[fR][fT] = np.append(resdict[fR][fT],res)
		'''
	

	# Write base histograms
	outFileBase.cd()
	outFileBase.mkdir('base')
	outFileBase.cd('base')
	for trackType in tracklist:
		for hist in baselist:
			basehists[trackType][hist].Write()
	# Write K_b histograms
	outFileBase.cd()
	outFileBase.mkdir('K_b')
	for trackType in tracklist:
		outFileBase.mkdir('K_b/'+styles[trackType]['name'])
		outFileBase.cd('K_b/'+styles[trackType]['name'])
		for eta in Kb_eta_list:
			for phi in Kb_phi_list:
				Kbhists[trackType][eta][phi].Write()
	for hist in histlist2D:
		outFileBase.cd()
		outFileBase.mkdir(hist)
		outFileBase.cd(hist)
		for track in tracklist:
			hists2D[track][hist].Write(hist+'_'+styles[track]['name'])
	# Write res vs f histogram
	outFileBase.cd()
	outFileBase.mkdir('res_vs_f')
	outFileBase.cd('res_vs_f')
	other['frachist'].Write()
	# Write delta P vs gen P histogram
	outFileBase.cd()
	outFileBase.mkdir('P_diff')
	outFileBase.cd('P_diff')
	for diff in difflist:
		other['P_diff_vs_gen_P'][diff].Write()
	'''
	# Write fR vs. fT
	fRvsfT['res'] = R.TH2D('fRvsfT_res','',11,0.999,1.0012,11,0.999,1.0012)
	fRvsfT['mu'] = R.TH2D('fRvsfT_mu','',11,0.999,1.0012,11,0.999,1.0012)
	for fR in steps:
		for fT in steps:
			mu,res = norm.fit(resdict[fR][fT])
			print fR,fT,mu,res
			fRvsfT['mu'].Fill(fR,fT,mu)
			fRvsfT['res'].Fill(fR,fT,res)
	outFileBase.cd()
	outFileBase.mkdir('fR_vs_fT')
	outFileBase.cd('fR_vs_fT')
	fRvsfT['res'].Write()
	fRvsfT['mu'].Write()
	'''
	outFileBase.cd()
	outFileBase.mkdir('test')
	outFileBase.cd('test')
	other['michalis'].Write()
	other['michalis2'].Write()

else:
	# Get file containing base histograms
	outFileBase = R.TFile(outfilebasename)
	# Write base histograms to file
	outFile.cd()
	outFile.mkdir('base')
	outFile.cd('base')
	for track in tracklist:
		for hist in baselist:
			basehists[track][hist] = outFileBase.Get('base/'+hist+'_'+styles[track]['name'])
			basehists[track][hist].Write(hist+'_'+styles[track]['name'])
	for track in tracklist:
		outFile.cd()
		outFile.mkdir('K_b/'+styles[track]['name'])
		outFile.cd('K_b/'+styles[track]['name'])
		for eta in Kb_eta_list:
			for phi in Kb_phi_list:
				hname = 'K_b_'+eta+'_'+phi+'_'+styles[track]['name']
				Kbhists[track][eta][phi] = outFileBase.Get('K_b/'+styles[track]['name']+'/'+hname)
				Kbhists[track][eta][phi].Write(hname)
	for hist in histlist2D:
		outFile.cd()
		outFile.mkdir(hist)
		outFile.cd(hist)
		for track in tracklist:
			hists2D[track][hist] = outFileBase.Get(hist+'/'+hist+'_'+styles[track]['name'])
			hists2D[track][hist].Write(hist+'_'+styles[track]['name'])
	other['frachist'] = outFileBase.Get('res_vs_f/frachist')
	outFile.cd()
	outFile.mkdir('res_vs_f')
	outFile.cd('res_vs_f')
	other['frachist'].Write('frachist')
	outFile.cd()
	outFile.mkdir('P_diff')
	outFile.cd('P_diff')
	for diff in difflist:
		other['P_diff_vs_gen_P'][diff] = outFileBase.Get('P_diff/P_diff_vs_gen_P_'+diff)
		other['P_diff_vs_gen_P'][diff].Write('P_diff_vs_gen_P_'+diff)

	'''
	fRvsfT['res'] = outFileBase.Get('fR_vs_fT/fRvsfT_res')
	fRvsfT['mu'] = outFileBase.Get('fR_vs_fT/fRvsfT_mu')
	outFile.cd()
	outFile.mkdir('fR_vs_fT')
	outFile.cd('fR_vs_fT')
	fRvsfT['res'].Write()
	fRvsfT['mu'].Write()
	'''
	other['michalis'] = outFileBase.Get('test/michalis')
	other['michalis2'] = outFileBase.Get('test/michalis2')
	outFile.cd()
	outFile.mkdir('test')
	outFile.cd('test')
	other['michalis'].Write()
	other['michalis2'].Write()
	
# ************************** #
# Plot Resolution Histograms #
# ************************** #

def plot_res(tracklist):
	c = R.TCanvas()
	name = 'res'
	names = [styles[track]['name'] for track in tracklist]
	for n in names: name += '_'+n
	plots = {track:basehists[track]['K_res'].Clone() for track in tracklist}
	for track in tracklist:
		# styling
		plots[track].SetLineColor(styles[track]['color'])
		plots[track].SetLineWidth(2)
		plots[track].Scale(1./plots[track].Integral())
	hmax = max([plots[track].GetMaximum() for track in tracklist])
	plots[tracklist[0]].SetMaximum(1.1*hmax)
	plots[tracklist[0]].GetXaxis().SetTitle('track(#Kappa)-gen(#Kappa)/gen(#Kappa)')
	plots[tracklist[0]].GetYaxis().SetTitle('A.U.')
	for t,track in enumerate(tracklist):
		optionstring = 'SAME' if t!=0 else ''
		plots[track].Draw('hist'+optionstring)
	c.Write(name)

outFile.cd()
outFile.mkdir('res')
outFile.cd('res')
plot_res(['global','refit','tracker','standAlone'])
plot_res(['refit','standAlone'])
plot_res(['global','comb'])
plot_res(['refit','tracker','comb'])
plot_res(['refit','tracker'])
for track in ['global','refit','tracker','standAlone','refit_noUpdate','comb']:
	plot_res([track])

# ***************************************** #
# Plot Resolution eta,phi for tracker/refit #
# ***************************************** #

def plot_K_b_grid(tracklist):
	c = R.TCanvas()
	plots = {eta:{phi:{track:Kbhists[track][eta][phi].Clone() for track in tracklist} for phi in Kb_phi_list} for eta in Kb_eta_list}
	gridplots = {track:R.TH2D('K_b_grid_'+track,'',4,-math.pi,math.pi,4,np.array([-2.4,-1.4,0,1.4,2.4])) for track in tracklist}
	diffname = 'K_b_grid_diff'
	names = [styles[track]['name'] for track in tracklist]
	for n in names:
		diffname += '_'+n
	griddiff = R.TH2D(diffname,'',4,-math.pi,math.pi,4,np.array([-2.4,-1.4,0,1.4,2.4]))
	gridtxt = {track:{eta:{phi:{} for phi in Kb_phi_list} for eta in Kb_eta_list} for track in tracklist}
	sdir = 'K_b/'
	for i,n in enumerate(names):
		if i==0: sdir += n
		else: sdir += '_'+n
	outFile.cd()
	outFile.mkdir(sdir)
	outFile.cd(sdir)
	for eta in Kb_eta_list:
		for phi in Kb_phi_list:
			name = 'K_b_'+eta+'_'+phi
			for n in names: name += '_'+n
			plot_K_b(tracklist,plots[eta][phi],sdir,name)
			for track in tracklist:
				gridtxt[track][eta][phi] = fit_K_b(plots[eta][phi][track],eta,phi,track)
	for e,eta in enumerate(Kb_eta_list):
		for p,phi in enumerate(Kb_phi_list):
			for track in tracklist:
				etabin = e+1
				phibin = p+1
				gridplots[track].SetBinContent(p,e,gridtxt[track][eta][phi][0])
				gridplots[track].SetBinError(p,e,gridtxt[track][eta][phi][1])
			griddiff.SetBinContent(p,e,gridtxt['tracker'][eta][phi][0]-gridtxt['refit'][eta][phi][0])
			griddiff.SetBinError(p,e,math.sqrt(gridtxt['tracker'][eta][phi][1]**2+gridtxt['refit'][eta][phi][1]**2))
	outFile.cd()
	outFile.cd('K_b')
	for track in tracklist:
		gridplots[track].Write()
	griddiff.Write()
	

def plot_K_b(tracklist,plots,sdir,name):
	c = R.TCanvas()
	for track in tracklist:
		plots[track].SetLineWidth(2)
		plots[track].SetLineColor(styles[track]['color'])
		plots[track].Scale(1./plots[track].Integral())
	hmax = max([plots[track].GetMaximum() for track in tracklist])
	for t,track in enumerate(tracklist):
		plots[track].Draw('hist'+('same' if t>0 else ''))
	plots[tracklist[0]].SetMaximum(1.1*hmax)
	c.Write(name)

def fit_K_b(plot,eta,phi,track):
	ffunc = R.TF1('name','gaus',-0.01,0.01)
	plot.Fit('name','QNR')
	mean = ffunc.GetParameter(1)
	sigma = ffunc.GetParameter(2)
	return [mean,sigma]
			
outFile.cd()
plot_K_b_grid(['tracker','refit'])

# ********************* #
# Plot Resolution vs. f #
# ********************* #

def plot_res_vs_f(frachist):
	outFile.cd()
	outFile.mkdir('res_vs_f/results')
	outFile.cd('res_vs_f/results')
	histSlices = R.TObjArray()
	func = R.TF1('name','gaus',-0.15,0.15)
	frachist.FitSlicesY(func,0,21,0,'QNR',histSlices)
	labels = ['Normalization','Gaussian Mean #mu','Gaussian Core #sigma','Gaussian Fit #chi^{2}']
	for h,hist in enumerate(histSlices):
		c = R.TCanvas()
		hist.SetLineWidth(2)
		hist.GetXaxis().SetTitle('f')
		hist.GetYaxis().SetTitle(labels[h])
		hist.Write()

def plot_hists_for_f(frachist):
	fracs = [i*0.05 for i in range(21)]
	outFile.cd()
	outFile.mkdir('res_vs_f/fits')
	outFile.cd('res_vs_f/fits')
	for i in range(frachist.GetNbinsX()):
		func = R.TF1('ff_'+str(fracs[i]),'gaus',-0.15,0.15)
		hist = frachist.ProjectionY('hres_'+str(fracs[i]),i+1,i+2)
		hist.Fit('ff_'+str(fracs[i]),'QR')
		hist.GetXaxis().SetTitle('#Kappa_{comb}-#Kappa_{gen}/#Kappa_{gen}')
		hist.GetYaxis().SetTitle('Counts')
		func.Draw('same')
		hist.Write()

plot_res_vs_f(other['frachist'])
plot_hists_for_f(other['frachist'])

# ************ #
# Plot Y vs. X #
# ************ #
# Does a FitSlicesY on the TH2 
# Plots the fits for each slice 
# Plots the mean, sigma, and sigma^2 for a few tracks together

def plot_fit_slices_Y_vs_X(hists2D,hist):
	outFile.cd()
	outFile.mkdir(hist+'/results')
	outFile.cd(hist+'/results')
	for track in tracklist:
		histSlices = R.TObjArray()
		func = R.TF1('name','gaus',-0.5,0.5)
		hists2D[track][hist].FitSlicesY(func,0,hists2D[track]['K_res_vs_pt'].GetNbinsX()+1,0,'QNR',histSlices)
		labels = ['Normalization','Gaussian Mean #mu','Gaussian Core #sigma','Gaussian Fit #chi^{2}']
		for h,fhist in enumerate(histSlices):
			c = R.TCanvas()
			fhist.SetLineWidth(2)
			fhist.GetXaxis().SetTitle(histDict[hist]['X'])
			fhist.GetYaxis().SetTitle(labels[h])
			fhist.Write()
			# plot sigma^2 vs pt
			if fhist.GetName()==hist+'_'+styles[track]['name']+'_1':
				fitplots['mean'][track][hist] = fhist.Clone()
			if fhist.GetName()==hist+'_'+styles[track]['name']+'_2':
				fitplots['sig'][track][hist] = fhist.Clone()
				sigsqfhist = fhist.Clone()
				fitplots['sigsq'][track][hist] = make_Y2_vs_X(sigsqfhist,hist,track,'sigsq')

def make_Y2_vs_X(hist,hname,track,what):
	for i in range(hist.GetNbinsX()):
		b = i+1
		sig2 = hist.GetBinContent(b)*hist.GetBinContent(b)
		esig2 = hist.GetBinError(b)*hist.GetBinContent(b)
		hist.SetBinContent(b,sig2)
		hist.SetBinError(b,esig2)
	hist.GetYaxis().SetTitle('#sigma_{fit}^{2}')
	hist.Write(track+'_'+hname+'_'+what+'_'+styles[track]['name']+'_2')
	return hist

def plot_hists_for_Y_vs_X(hists2D,hist):
	for track in tracklist:
		outFile.cd()
		outFile.mkdir(hist+'/'+track+'_fits')
		outFile.cd(hist+'/'+track+'_fits')
		dhist = hists2D[track][hist]
		for i in range(dhist.GetNbinsX()):
			func = R.TF1('ff_'+hist+'_'+str(i)+'_'+track,'gaus',-0.5,0.5)
			phist = dhist.ProjectionY(hist+'_'+track+'_'+str(i),i+1,i+2)
			phist.Fit('ff_'+hist+'_'+str(i)+'_'+track,'QR')
			phist.GetXaxis().SetTitle(histDict[hist]['Y'])
			phist.GetYaxis().SetTitle('Counts')
			func.Draw('same')
			phist.Write()

def plot_Y_vs_X_together(tracklist,fitplots,together,name):
	outFile.cd(name+'/results')
	c = R.TCanvas()
	hmaxes = []
	for t,track in enumerate(tracklist):
		hist = fitplots[track][name]
		hist.SetLineColor(styles[track]['color'])
		hist.Draw('same' if t>0 else '')
		hmaxes.append(hist.GetMaximum())
	hmax = max(hmaxes)
	fitplots[tracklist[0]][name].SetMaximum(hmax*1.1)
	for track in tracklist:
		name += '_'+styles[track]['name']
	name += '_'+together
	c.Write(name)

togetherlist = ['sigsq','sig','mean']
fitplots = {together:{track:{hist:{} for hist in histlist2D} for track in tracklist} for together in togetherlist}
for hist in histlist2D:
	plot_fit_slices_Y_vs_X(hists2D,hist)
	plot_hists_for_Y_vs_X(hists2D,hist)
	for together in togetherlist:
		plot_Y_vs_X_together(['tracker','refit'],fitplots[together],together,hist)
		plot_Y_vs_X_together(['global','comb'],fitplots[together],together,hist)
		plot_Y_vs_X_together(['global','comb','tracker','refit'],fitplots[together],together,hist)

# ******************************* #
# Plot Profile of P diff vs gen P #
# ******************************* #

def plot_profileY(hist):
	profhist = hist.ProfileY()
	profhist.Write()
def plot_profileX(hist):
	profhist = hist.ProfileX()
	profhist.Write()

outFile.cd()
outFile.cd('P_diff')
for diff in difflist:
	plot_profileY(other['P_diff_vs_gen_P'][diff])
	plot_profileX(other['P_diff_vs_gen_P'][diff])
