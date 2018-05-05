import sys
import ROOT as R
#import HighPT.Analysis.tools as t
#from HighPT.Analysis.plot_styles import styles

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

rankflist = [1,5,10,15,20,25,30,35,50]
tracktypes = ['tracker','standAlone','global','picky','tpfms','dyt']
fittypes = ['','NoRPC','Comb']

widths = ['rms','sigma','sigma_core']

data = {track:{width:{f:{} for f in rankflist} for width in widths} for track in ['bestComb','bestMuonOnlyUpdate']}

for trackName in ['bestComb','bestMuonOnlyUpdate']:
	#for f in rankflist:
	#inputFilename = 'hists_curvature_muonMC_f'+str(f)+'.root'
	#inputFilename = 'hists_curvature_muonMC_bestPull.root'
	#inputFilename = 'hists_curvature_dxy.root'
	#inputFilename = 'hists_curvature_muonMC_bestPull_test.root'
	inputFilename = 'hists_curvature_muonMC_bestPull_test_pt_gt200.root'
	inputFile = R.TFile(inputFilename)

	#trackName = 'bestComb'
	hname = trackName+'_K_gen_res'
	resplot = inputFile.Get(hname)

	#data[trackName]['rms'][f] = resplot.GetRMS()
	data[trackName]['rms'] = resplot.GetRMS()

	#fitname = trackName+'_'+str(f)
	fitname = trackName#+'_'+str(f)
	fit = R.TF1(fitname,'gaus',-0.5,0.5)
	resplot.Fit(fitname,'R')
	#data[trackName]['sigma'][f] = fit.GetParameter(2)
	data[trackName]['sigma'] = fit.GetParameter(2)

	fitcore = R.TF1(fitname+'_core','gaus',-0.1,0.1)
	resplot.Fit(fitname+'_core','R')
	#data[trackName]['sigma_core'][f] = fitcore.GetParameter(2)
	data[trackName]['sigma_core'] = fitcore.GetParameter(2)


trackdata = {track+ftype:{width:{} for width in widths} for track in tracktypes for ftype in fittypes}
for track in tracktypes:
	for ftype in fittypes:
		if track=='bestMuonOnlyUpdate' and ftype=='Comb': continue
		if (track=='tracker' or track=='standAlone') and (ftype=='NoRPC' or ftype=='Comb'): continue
		trackname = track+ftype
		plt = inputFile.Get(trackname+'_K_gen_res').Clone()
		plt2 = inputFile.Get(trackname+'_K_gen_res').Clone()
		trackdata[trackname]['rms'] = plt.GetRMS()

		thefit = R.TF1(trackname+'_fit','gaus',-0.5,0.5)
		plt.Fit(trackname+'_fit','R')
		trackdata[trackname]['sigma'] = thefit.GetParameter(2)

		core = R.TF1(trackname+'_fitcore','gaus',-0.1,0.1)
		plt2.Fit(trackname+'_fitcore','R')
		trackdata[trackname]['sigma_core'] = core.GetParameter(2)

print
print '*'*10
print
print '   {0:7} {1:7} {2:7}'.format('rms', 'sigma', 'sigcore')
for track in tracktypes:
	for ftype in fittypes:
		if track=='bestMuonOnlyUpdate' and ftype=='Comb': continue
		if (track=='tracker' or track=='standAlone') and (ftype=='NoRPC' or ftype=='Comb'): continue
		trackname = track+ftype
		print trackname
		print '   {0:7.4f} {1:7.4f} {2:7.4f}'.format(trackdata[trackname]['rms'], trackdata[trackname]['sigma'], trackdata[trackname]['sigma_core'])
	print
for trackName in ['bestMuonOnlyUpdate','bestComb']:
	print trackName
#	for f in rankflist:
#		print '{0:2} {1:7.4f} {2:7.4f} {3:7.4f}'.format(int(f), data[trackName]['rms'][f], data[trackName]['sigma'][f], data[trackName]['sigma_core'][f])
	print '   {0:7.4f} {1:7.4f} {2:7.4f}'.format(data[trackName]['rms'], data[trackName]['sigma'], data[trackName]['sigma_core'])
