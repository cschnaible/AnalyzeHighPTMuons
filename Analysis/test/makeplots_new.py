import sys
import ROOT as R
import math as math
import numpy as np
from scipy.stats import norm
import argparse
import tools as t
from HighPT.Analysis.plot_styles import styles

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file',default='ZpMM',help='Which MC file to use')
parser.add_argument('-ex','--extra',default='',help='Extra name to output histogram file')
parser.add_argument('-r','--recreate',default=False,action='store_true',help='Recreate histogram file')
args = parser.parse_args()

#fileName = 'AnalyzeTracks_'+str(args.file)+'.root'
fileName = 'testing_again.root'
nametmp = 'hists_'+str(args.file)+('' if args.extra=='' else '_'+str(args.extra))
outfilebasename = nametmp+'_base.root'
outfilename     = nametmp+'.root'
outFile = R.TFile(outfilename,'recreate')
RECREATE = args.recreate

varList = {
		# K = q/P
		'K':{
			'bins':(100,-0.01,0.01),
			'label':'#kappa [Gev^{-1}]',
			'log':False,
			},
		# lambda
		'lambda':{
			'bins':(30,-1.57,1.57),
			'label':'#lambda',
			'log':False,
			},
		# phi
		'phi':{
			'bins':(30,-3.15,3.15),
			'label':'#phi',
			'log':False,
			},
		# d_xy
		'dxy':{
			'bins':(100,-10,10),# binning is for refit_noUpdate tracks, others should be -0.2 to +0.2
			'label':'d_{xy} [cm]',
			'log':False,
			},
		# d_sz
		'dsz':{
			'bins':(100,-20,20),
			'label':'d_{sz} [cm]',
			'log':False,
			},
		# eta = -ln(tan((pi/2-lambda)/2))
		'full_eta':{
			'bins':(60,-3,3),
			'label':'#eta',
			'log':False,
			},
		'gen_eta':{
			'bins':(60,-3,3),
			'label':'#eta',
			'log':False,
			},
		'eta':{
			'bins':(60,-3,3),
			'label':'#eta',
			'log':False,
			},
		# pT
		'gen_pt':{
			'bins':(150,0,1500),
			'label':'gen p_{T} [GeV]',
			'log':False,
			},
		# pT
		'full_pt':{
			'bins':(200,0,2000),
			'label':'p_{T} [GeV]',
			'log':False,
			},
		# pT
		'pt':{
			'bins':(200,0,2000),
			'label':'p_{T} [GeV]',
			'log':False,
			},
#		# pT^2
#		'pt2':{
#			'bins':(50,0,10000000),
#			'label':'p_{T}^{2} [GeV^{2}]',
#			'log':False,
#			},
		# K_gen_res = K(track) - K(gen) / K(gen)
		'K_gen_res':{
			'bins':(100,-1,1),
			'label':'#kappa(track)-#kappa(gen)/#kappa(gen)',
			'log':False,
			},
		# K_gen_res = K(track) - K(gen) / K(gen)
		'K_full_res':{
			'bins':(100,-5,5),
			'label':'#kappa(refit)-#kappa(full)/#kappa(full)',
			'log':False,
			},
#		# K_b = K(track) - K(gen)
#		'K_b':{
#			'bins':(100,-0.001,0.001),
#			'label':'#kappa(track)-#kappa(gen)/#kappa(gen)',
#			'log':False,
#			},
		# K_rel_err = sigma(K) / K
		'K_rel_err':{
			'bins':(100,0,10),
			'label':'#delta_{#kappa}/|#kappa|',
			'log':False,
			},
#		# P_diff = P(track1) - P(track2)
#		'P_diff':{
#			'bins':(30,-1500,1500),
#			'label':'Gen p_{T} [GeV]',
#			'log':False,
#			},
		# P
#		'gen_P':{
#			'bins':(75,0,2500),
#			'label':'Gen p_{T} [GeV]',
#			'log':False,
#			},
		# combination chi^2
		'chi2':{
			'bins':(100,0,100),
			'label':'Track #chi^{2}',
			'log':False,
			},
		# lambda(track1) - lambda(track2)
		'lambda_diff':{
			'bins':(30,-1.6,1.6),
			'label':'#Delta #lambda',
			'log':False,
			},
		# phi(track1) - phi(track2)
		'phi_diff':{
			'bins':(30,-6.3,6.3),
			'label':'#Delta #phi',
			'log':False,
			},
		# d_xy(track1) - d_xy(track2)
		'dxy_diff':{
			'bins':(30,-10,10),
			'label':'#Delta d_{xy} [cm]',
			'log':False,
			},
		'dxyBS_diff':{
			'bins':(30,-1,1),
			'label':'#Delta d_{xy} [cm]',
			'log':False,
			},
		'dxyPV_diff':{
			'bins':(30,-1,1),
			'label':'#Delta d_{xy} [cm]',
			'log':False,
			},
		# d_sz(track1) - d_sz(track2)
		'dsz_diff':{
			'bins':(30,-10,10),
			'label':'#Delta d_{sz} [cm]',
			'log':False,
			},
		'dszBS_diff':{
			'bins':(30,-1,1),
			'label':'#Delta d_{sz} [cm]',
			'log':False,
			},
		'dszPV_diff':{
			'bins':(30,-1,1),
			'label':'#Delta d_{sz} [cm]',
			'log':False,
			},
		# 3D distance between track reference points
		'pca_diff':{
			'bins':(50,np.logspace(-2,3,51)),
			'label':'|#Delta v| [cm]',
			'log':True,
			},
		# Number of valid hits used in track
		'nValidHits':{
			'bins':(80,0,80),
			'label':'N(valid hits)',
			'log':False,
			},
		# Number of lost hits used in track
		'nLostHits':{
			'bins':(80,0,80),
			'label':'N(lost hits)',
			'log':False,
			},
		# Number of degrees of freedom in track
		'ndof':{
			'bins':(80,0,80),
			'label':'N(dof)',
			'log':False,
			},
		}

tracklist = ['global','global_comb','global_refit','global_refit_noUpdate',\
			 'picky', 'picky_comb', 'picky_refit', 'picky_refit_noUpdate',\
			 'tracker','standAlone']

plotlist1D = ['K','lambda','phi','dxy','dsz',\
		'K_gen_res','K_full_res','K_rel_err',\
		'eta','pt','chi2',\
		'phi_diff','dxy_diff','dsz_diff',\
		'nValidHits','ndof','pca_diff',\
		'dxyBS_diff','dxyPV_diff','dszBS_diff','dszPV_diff']

plotlist2d = [
	'K_gen_res_vs_gen_pt', 'K_gen_res_vs_gen_eta',
	'K_gen_res_vs_dxy_diff', 'K_gen_res_vs_dsz_diff',
	'K_gen_res_vs_dxyBS_diff', 'K_gen_res_vs_dszBS_diff',
	'K_gen_res_vs_dxyPV_diff', 'K_gen_res_vs_dszPV_diff',
	'K_gen_res_vs_phi_diff', 
	'K_gen_res_vs_nValidHits',
	'K_gen_res_vs_pca_diff',

	'K_full_res_vs_full_pt', 'K_full_res_vs_full_eta',
	'K_full_res_vs_dxy_diff', 'K_full_res_vs_dsz_diff', 
	'K_full_res_vs_dxyBS_diff', 'K_full_res_vs_dszBS_diff', 
	'K_full_res_vs_dxyPV_diff', 'K_full_res_vs_dszPV_diff', 
	'K_full_res_vs_phi_diff', 
	'K_full_res_vs_nValidHits',
	'K_full_res_vs_pca_diff',

	'K_rel_err_vs_gen_pt', 'K_rel_err_vs_gen_eta',
	'K_rel_err_vs_dxy_diff', 'K_rel_err_vs_dsz_diff',
	'K_rel_err_vs_dxyBS_diff', 'K_rel_err_vs_dszBS_diff',
	'K_rel_err_vs_dxyPV_diff', 'K_rel_err_vs_dszPV_diff',
	'K_rel_err_vs_phi_diff', 
	'K_rel_err_vs_nValidHits',
	'K_rel_err_vs_pca_diff',

	'nValidHits_vs_gen_eta', 'nValidHits_vs_gen_pt',
	'phi_diff_vs_gen_pt', 'phi_diff_vs_gen_eta', 'phi_diff_vs_nValidHits',
	'pca_diff_vs_gen_pt', 'pca_diff_vs_gen_eta', 'pca_diff_vs_nValidHits',
	'dxy_diff_vs_gen_pt', 'dxy_diff_vs_gen_eta', 'dxy_diff_vs_nValidHits',
	'dsz_diff_vs_gen_pt', 'dsz_diff_vs_gen_eta', 'dsz_diff_vs_nValidHits',
	'dxyBS_diff_vs_gen_pt', 'dxyBS_diff_vs_gen_eta', 'dxyBS_diff_vs_nValidHits',
	'dszBS_diff_vs_gen_pt', 'dszBS_diff_vs_gen_eta', 'dszBS_diff_vs_nValidHits',
	'dxyPV_diff_vs_gen_pt', 'dxyPV_diff_vs_gen_eta', 'dxyPV_diff_vs_nValidHits',
	'dszPV_diff_vs_gen_pt', 'dszPV_diff_vs_gen_eta', 'dszPV_diff_vs_nValidHits'
	]

hists1D = {track:{} for track in tracklist}
hists2D = {track:{} for track in tracklist}



# frachist
# K(comb,f) = f*K(ref) + (1-f)*K(trk), sigma(K(comb,f)-K(gen) / K(gen)) vs. f
other = {
		'global':{'frachist':{}},
		'picky':{'frachist':{}},
	}

if RECREATE:
	rFile = R.TFile(fileName)

	tree = rFile.Get('t')
	outFileBase = R.TFile(outfilebasename,'recreate')

	for track in tracklist:
		for plot1D in plotlist1D:
			if ((track in ['global','picky','standAlone'])\
					and ('full' in plot1D or 'diff' in plot1D)):
				continue
			if (('comb' in track or 'refit'==track[-5:])\
					and ('pca' in plot1D or 'nValidHits' in plot1D)):
				continue
			name = track+'_'+plot1D
			title = ';'+varList[plot1D]['label']+';Counts'
			hists1D[track][plot1D] = R.TH1D(name,title,*varList[plot1D]['bins'])
		for plot2D in plotlist2d:
			if ((track in ['global','picky','standAlone'])\
					and ('full' in plot2D or 'diff' in plot2D)):
				continue
			if (('comb' in track or 'refit'==track[-5:])\
					and ('pca' in plot2D or 'nValidHits' in plot2D)):
				continue
			name = track+'_'+plot2D
			x,y = plot2D.split('_vs_')[1],plot2D.split('_vs_')[0]
			title = ';'+varList[x]['label']+';'+varList[y]['label'].strip('_gen')
			hists2D[track][plot2D] = R.TH2D(name,title,*varList[x]['bins']+varList[y]['bins'])

	fracs = [i*0.05 for i in range(21)]
	other['global']['frachist'] = R.TH2D('global_frachist','#kappa_{comb} = f #times #kappa_{refit} + (1-f) #times #kappa_{tracker}',21,0,1.05,100,-0.5,0.5)
	other['picky']['frachist'] = R.TH2D('picky_frachist','#kappa_{comb} = f #times #kappa_{refit} + (1-f) #times #kappa_{tracker}',21,0,1.05,100,-0.5,0.5)

	# loop on tree
	for l,entry in enumerate(tree):
		#if l>2000: continue
		sys.stdout.write('Entry {l:>5}\r'.format(**locals()))
		p = t.pyTree(tree)
		if not p.allOkay: continue

		# Plots for all tracks
		for trackType in tracklist:
			track = t.Track(trackType,p)
			if track.cov(0,0) < 0: continue
			hists1D[trackType]['K'].Fill(track.K())
			hists1D[trackType]['lambda'].Fill(track.Lambda())
			hists1D[trackType]['phi'].Fill(track.phi())
			hists1D[trackType]['dxy'].Fill(track.dxy())
			hists1D[trackType]['dsz'].Fill(track.dsz())
			hists1D[trackType]['eta'].Fill(track.eta())
			hists1D[trackType]['pt'].Fill(track.pt())
			hists1D[trackType]['K_gen_res'].Fill((track.K()-p.gen_K)/p.gen_K)
			hists1D[trackType]['K_rel_err'].Fill(math.sqrt(track.cov(0,0))/abs(track.K()))
			hists1D[trackType]['chi2'].Fill(track.chi2())
			hists2D[trackType]['K_gen_res_vs_gen_pt'].Fill(p.gen_pt, (track.K()-p.gen_K)/p.gen_K)
			hists2D[trackType]['K_gen_res_vs_gen_eta'].Fill(p.gen_eta, (track.K()-p.gen_K)/p.gen_K)
			hists2D[trackType]['K_rel_err_vs_gen_pt'].Fill(p.gen_pt, math.sqrt(track.cov(0,0))/abs(track.K()))
			hists2D[trackType]['K_rel_err_vs_gen_eta'].Fill(p.gen_eta, math.sqrt(track.cov(0,0))/abs(track.K()))

		# Plots for some tracks
		for full in ['picky','global']:
			fullTrack = t.Track(full,p)
			if fullTrack.cov(0,0) < 0: continue
			hists1D[full]['nValidHits'].Fill(fullTrack.nValidHits())
			hists1D[full]['ndof'].Fill(fullTrack.ndof())
			hists2D[full]['nValidHits_vs_gen_pt'].Fill(p.gen_pt, fullTrack.nValidHits())
			hists2D[full]['nValidHits_vs_gen_eta'].Fill(p.gen_eta, fullTrack.nValidHits())
			hists2D[full]['K_gen_res_vs_nValidHits'].Fill(fullTrack.nValidHits(), (fullTrack.K()-p.gen_K)/p.gen_K)
			hists2D[full]['K_rel_err_vs_nValidHits'].Fill(fullTrack.nValidHits(), math.sqrt(fullTrack.cov(0,0)))
			for trackType in [full+'_refit',full+'_refit_noUpdate',full+'_comb','tracker']:
				track = t.Track(trackType,p)
				if track.cov(0,0) < 0: continue
				#hists1D[trackType]['lambda_diff'].Fill(track.Lambda()-full.Lambda())
				hists1D[trackType]['phi_diff'].Fill(track.phi()-fullTrack.phi())
				hists1D[trackType]['dxy_diff'].Fill(track.dxy()-fullTrack.dxy())
				hists1D[trackType]['dxyBS_diff'].Fill(track.dxy()-fullTrack.dxyBS())
				hists1D[trackType]['dxyPV_diff'].Fill(track.dxy()-fullTrack.dxyPV())
				hists1D[trackType]['dsz_diff'].Fill(track.dsz()-fullTrack.dsz())
				hists1D[trackType]['dszBS_diff'].Fill(track.dsz()-fullTrack.dszBS())
				hists1D[trackType]['dszPV_diff'].Fill(track.dsz()-fullTrack.dszPV())
				hists1D[trackType]['K_full_res'].Fill((track.K()-fullTrack.K())/fullTrack.K())
				
				hists2D[trackType]['K_full_res_vs_full_pt'].Fill(fullTrack.pt(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_full_eta'].Fill(fullTrack.eta(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['phi_diff_vs_gen_pt'].Fill(p.gen_pt, track.phi()-fullTrack.phi())
				hists2D[trackType]['phi_diff_vs_gen_eta'].Fill(p.gen_eta, track.phi()-fullTrack.phi())
				hists2D[trackType]['dxy_diff_vs_gen_pt'].Fill(p.gen_pt, track.dxy()-fullTrack.dxy())
				hists2D[trackType]['dxy_diff_vs_gen_eta'].Fill(p.gen_eta, track.dxy()-fullTrack.dxy())
				hists2D[trackType]['dsz_diff_vs_gen_pt'].Fill(p.gen_pt, track.dsz()-fullTrack.dsz())
				hists2D[trackType]['dsz_diff_vs_gen_eta'].Fill(p.gen_eta, track.dsz()-fullTrack.dsz())

				hists2D[trackType]['K_gen_res_vs_phi_diff'].Fill(track.phi()-fullTrack.phi(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dxy_diff'].Fill(track.dxy()-fullTrack.dxy(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dxyBS_diff'].Fill(track.dxy()-fullTrack.dxyBS(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dxyPV_diff'].Fill(track.dxy()-fullTrack.dxyPV(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dsz_diff'].Fill(track.dsz()-fullTrack.dsz(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dszBS_diff'].Fill(track.dsz()-fullTrack.dszBS(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_dszPV_diff'].Fill(track.dsz()-fullTrack.dszPV(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_rel_err_vs_phi_diff'].Fill(track.phi()-fullTrack.phi(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dxy_diff'].Fill(track.dxy()-fullTrack.dxy(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dxyBS_diff'].Fill(track.dxy()-fullTrack.dxyBS(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dxyPV_diff'].Fill(track.dxy()-fullTrack.dxyPV(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dsz_diff'].Fill(track.dsz()-fullTrack.dsz(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dszBS_diff'].Fill(track.dsz()-fullTrack.dszBS(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_dszPV_diff'].Fill(track.dsz()-fullTrack.dszPV(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_full_res_vs_phi_diff'].Fill(track.phi()-fullTrack.phi(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dxy_diff'].Fill(track.dxy()-fullTrack.dxy(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dxy_diff'].Fill(track.dxyBS()-fullTrack.dxyBS(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dxy_diff'].Fill(track.dxyPV()-fullTrack.dxyPV(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dsz_diff'].Fill(track.dsz()-fullTrack.dsz(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dsz_diff'].Fill(track.dszBS()-fullTrack.dszBS(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_dsz_diff'].Fill(track.dszPV()-fullTrack.dszPV(), (track.K()-fullTrack.K())/fullTrack.K())

			for trackType in [full+'_refit_noUpdate','tracker']:
				track = t.Track(trackType,p)
				hists1D[trackType]['ndof'].Fill(track.ndof())
				hists1D[trackType]['nValidHits'].Fill(track.nValidHits())
				hists1D[trackType]['pca_diff'].Fill(t.pos_diff(track._pos,fullTrack._pos))

				hists2D[trackType]['pca_diff_vs_gen_pt'].Fill(p.gen_pt, t.pos_diff(track._pos,fullTrack._pos))
				hists2D[trackType]['pca_diff_vs_gen_eta'].Fill(p.gen_eta, t.pos_diff(track._pos,fullTrack._pos))
				hists2D[trackType]['pca_diff_vs_nValidHits'].Fill(track.nValidHits(), t.pos_diff(track._pos,fullTrack._pos))
				hists2D[trackType]['dxy_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dxy()-fullTrack.dxy())
				hists2D[trackType]['dxyBS_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dxy()-fullTrack.dxyBS())
				hists2D[trackType]['dxyPV_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dxy()-fullTrack.dxyPV())
				hists2D[trackType]['dsz_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dsz()-fullTrack.dsz())
				hists2D[trackType]['dszBS_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dsz()-fullTrack.dszBS())
				hists2D[trackType]['dszPV_diff_vs_nValidHits'].Fill(track.nValidHits(), track.dsz()-fullTrack.dszPV())
				hists2D[trackType]['K_gen_res_vs_nValidHits'].Fill(track.nValidHits(), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_gen_res_vs_pca_diff'].Fill(t.pos_diff(track._pos,fullTrack._pos), (track.K()-p.gen_K)/p.gen_K)
				hists2D[trackType]['K_rel_err_vs_pca_diff'].Fill(t.pos_diff(track._pos,fullTrack._pos), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_rel_err_vs_nValidHits'].Fill(track.nValidHits(), math.sqrt(track.cov(0,0))/abs(track.K()))
				hists2D[trackType]['K_full_res_vs_nValidHits'].Fill(track.nValidHits(), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['K_full_res_vs_pca_diff'].Fill(t.pos_diff(track._pos,fullTrack._pos), (track.K()-fullTrack.K())/fullTrack.K())
				hists2D[trackType]['phi_diff_vs_nValidHits'].Fill(track.nValidHits(), track.phi()-fullTrack.phi())
				hists2D[trackType]['nValidHits_vs_gen_pt'].Fill(p.gen_pt, track.nValidHits())
				hists2D[trackType]['nValidHits_vs_gen_eta'].Fill(p.gen_eta, track.nValidHits())


		global_refit = t.Track('global_refit',p)
		picky_refit = t.Track('picky_refit',p)
		tracker = t.Track('tracker',p)
		for f in fracs:
			#Kcomb = f*p.refit_par[0] + (1.-f)*p.tracker_par[0]
			KcombGlobal = f*global_refit.K() + (1.-f)*tracker.K()
			KcombPicky = f*picky_refit.K() + (1.-f)*tracker.K()
			other['global']['frachist'].Fill(f,(KcombGlobal-p.gen_K)/p.gen_K)
			other['picky']['frachist'].Fill(f,(KcombPicky-p.gen_K)/p.gen_K)

	outFileBase.cd()
	for track in tracklist:
		outFileBase.cd()
		outFileBase.mkdir(track)
		outFileBase.cd(track)
		for plot in hists1D[track].keys():
			hists1D[track][plot].Write()
		for plot in hists2D[track].keys():
			hists2D[track][plot].Write()

else:
	# Get file containing base histograms
	outFileBase = R.TFile(outfilebasename)
	# Write base histograms to file
	# Set dictionary with histograms
	for track in tracklist:
		outFile.cd()
		outFile.mkdir(track)
		outFile.cd(track)
		for plot1D in plotlist1D:
			if ((track in ['global','picky','standAlone'])\
					and ('full' in plot1D or 'diff' in plot1D)):
				continue
			if (('comb' in track or 'refit'==track[-5:])\
					and ('pca' in plot1D or 'nValidHits' in plot1D)):
				continue
			name = track+'_'+plot1D
			hists1D[track][plot1D] = outFileBase.Get(track+'/'+name)
			hists1D[track][plot1D].Write()
		for plot2D in plotlist2d:
			if ((track in ['global','picky','standAlone'])\
					and ('full' in plot2D or 'diff' in plot2D)):
				continue
			if (('comb' in track or 'refit'==track[-5:])\
					and ('pca' in plot2D or 'nValidHits' in plot2D)):
				continue
			name = track+'_'+plot2D
			hists2D[track][plot2D] = outFileBase.Get(track+'/'+name)
			hists2D[track][plot2D].Write()
	outFile.cd()


# Plots the gaussian fit normalization, mean, and sigma for y-axis slices
# Combines plots together for tracks in tracklist
def plot_gaussian_fits(tracklist,htype='K_gen_res_vs_gen_pt'):
	'''
	'''
	c = R.TCanvas()
	name = 'gaus_fit_res'
	names = [styles[track]['name'] for track in tracklist]
	for n in names: name += '_'+n
	slices = {track:R.TObjArray() for track in tracklist}
	# Make fit slices
	for track in tracklist:
		fitmin,fitmax = varList[htype.split('_vs_')[0]]['bins'][1],varList[htype.split('_vs_')[0]]['bins'][2]
		func = R.TF1('name','gaus',fitmin,fitmax)
		hists2D[track][htype].FitSlicesY(func,0,-1,0,'QNR',slices[track])
	legpos = {
			0:(0.45, 0.6, 0.85, 0.7),
			1:(0.15, 0.3, 0.45, 0.4),
			2:(0.15, 0.7, 0.45, 0.8),
			3:(0.45, 0.7, 0.85, 0.8),
			}
	for h in range(4):
		c = R.TCanvas()
		hmax = max([slices[track][h].GetMaximum() for track in tracklist])
		leg = R.TLegend(*legpos[h])
		for t,track in enumerate(tracklist):
			slices[track][h].SetLineColor(styles[track]['color'])
			slices[track][h].Draw('same' if t>0 else '')
			leg.AddEntry(slices[track][h],styles[track]['pretty'])
		leg.Draw()
		slices[tracklist[0]][h].SetMaximum(hmax*1.1)
		c.Write(name+'_'+str(h))

for htype in ['K_gen_res_vs_gen_pt']:#,'K_full_res_vs_full_pt','K_gen_res_vs_gen_eta']:
	outFile.cd()
	outFile.mkdir(htype)
	outFile.cd(htype)
	for track in tracklist:
		plot_gaussian_fits([track],htype)
	plot_gaussian_fits(['global','picky'],htype)
	plot_gaussian_fits(['global','global_comb'],htype)
	plot_gaussian_fits(['picky','picky_comb'],htype)
	plot_gaussian_fits(['global_comb','picky_comb'],htype)


# Plots simple distributions on top of each other
def plot_dist(tracklist, var):
	c = R.TCanvas()
	name = var
	names = [styles[track]['name'] for track in tracklist]
	for n in names: name += '_'+n
	for t,track in enumerate(tracklist):
		hists1D[track][var].SetLineColor(styles[track]['color'])
		hists1D[track][var].Draw('same' if t>0 else '')
	hmax = max([hists1D[track][var].GetMaximum() for track in tracklist])
	hists1D[tracklist[0]][var].SetMaximum(hmax*1.1)
	c.Write(name)

for htype in ['dxy','dsz','phi','lambda','eta','K_gen_res','K_rel_err']:
	outFile.cd()
	outFile.mkdir(htype)
	outFile.cd(htype)
	plot_dist(['global_refit_noUpdate','global'],htype)
	plot_dist(['picky_refit_noUpdate','picky'],htype)
	
exit()



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
	plots[tracklist[0]].GetXaxis().SetTitle('track(#kappa)-gen(#kappa)/gen(#kappa)')
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
		hist.GetXaxis().SetTitle('#kappa_{comb}-#kappa_{gen}/#kappa_{gen}')
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
