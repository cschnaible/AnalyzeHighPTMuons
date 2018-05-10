import sys
import ROOT as R
#import HighPT.Analysis.tools as t
#from HighPT.Analysis.plot_styles import styles

R.gROOT.SetBatch(True)
R.gStyle.SetPadTickX(1)
R.gStyle.SetPadTickY(1)

rankflist = [1,5,10,15,20,25,30,35,50]
tracktypes = ['tracker','standAlone','global','picky','tpfms','dyt','bestMuonOnlyUpdate']
fittypes = ['','NoRPC','Comb']
#combine = 'simple'
combine = 'full'

#inputFilename = 'hists_curvature_muonMC_'+selector+('_'+extra if extra else '')+'.root'
selectors = ['dxy','curvPull','trackRank']
extras = {'dxy':[''],'curvPull':[''],\
		'trackRank':['f1','f3','f5','f10','f15']}
ptcuts = ['0','200','400','600','800','1000']

widths = ['rms','sigma','sigma_core']
bestdata = {track:{selector+extra:{ptcut:{width:{} for width in widths} for ptcut in ptcuts} for selector in selectors for extra in extras[selector]} for track in ['bestComb','bestMuonOnlyUpdate']}

for trackName in ['bestComb','bestMuonOnlyUpdate']:
	for selector in selectors:
		for extra in extras[selector]:
			for ptcut in ptcuts:
				inputFilename = 'hists_curvature_muonMC_'+selector+('_'+extra if extra else '')+'.root'
				inputFile = R.TFile(inputFilename)
				print inputFilename

				hname = trackName+'_K_gen_res_'+combine+'_'+ptcut
				resplot = inputFile.Get(hname)

				bestdata[trackName][selector+extra][ptcut]['rms'] = resplot.GetRMS()

				fitname = hname+'_fit'
				fit = R.TF1(fitname,'gaus',-0.5,0.5)
				resplot.Fit(fitname,'RQ')
				bestdata[trackName][selector+extra][ptcut]['sigma'] = fit.GetParameter(2)

				fitcore = R.TF1(fitname+'_core','gaus',-0.1,0.1)
				resplot.Fit(fitname+'_core','RQ')
				bestdata[trackName][selector+extra][ptcut]['sigma_core'] = fitcore.GetParameter(2)


print '   {0:7} {1:7} {2:7}'.format('rms', 'sigma', 'sigcore')

trackdata = {track+ftype:{ptcut:{} for ptcut in ptcuts} for track in tracktypes for ftype in fittypes}
for track in tracktypes:
	for ftype in fittypes:
		for ptcut in ptcuts:
			print track, ftype, ptcut
			if track=='bestMuonOnlyUpdate' and ftype!='': continue
			if (track=='tracker' or track=='standAlone') and (ftype=='NoRPC' or ftype=='Comb'): continue
			trackname = track+ftype
			inputFile = R.TFile('hists_curvature_muonMC_curvPull.root')
			plt = inputFile.Get(trackname+'_K_gen_res_'+combine+'_'+ptcut).Clone()
			plt2 = inputFile.Get(trackname+'_K_gen_res_'+combine+'_'+ptcut).Clone()
			trackdata[track+ftype][ptcut]['rms'] = plt.GetRMS()
			rms = plt.GetRMS()

			thefit = R.TF1(trackname+'_fit','gaus',-0.5,0.5)
			plt.Fit(trackname+'_fit','RQ')
			sigma = thefit.GetParameter(2)
			trackdata[track+ftype][ptcut]['sigma'] = thefit.GetParameter(2)

			core = R.TF1(trackname+'_fitcore','gaus',-0.1,0.1)
			plt2.Fit(trackname+'_fitcore','RQ')
			sig_core = core.GetParameter(2)
			trackdata[track+ftype][ptcut]['sigma_core'] = core.GetParameter(2)
			print '   {0:7.4f} {1:7.4f} {2:7.4f}'.format(rms, sigma, sig_core)
		print
	print
	print '*'*10
	print

for trackName in ['bestMuonOnlyUpdate','bestComb']:
	for ptcut in ptcuts:
		print
		print ptcut
		for selector in selectors:
			for extra in extras[selector]:
				print trackName, selector, extra
				print '   {0:7.4f} {1:7.4f} {2:7.4f}'.format(bestdata[trackName][selector+extra][ptcut]['rms'], bestdata[trackName][selector+extra][ptcut]['sigma'], bestdata[trackName][selector+extra][ptcut]['sigma_core'])
	print
	print '*'*10
	print



for width in widths:
	print
	print '{0:>12}{1:>7} {2:>7} {3:>7} {4:>7} {5:>7} {6:>7}'.format(width,'>10','>200','>400','>600','>800','>1000')
	for track in ['global','globalComb','bestMuonOnlyUpdate']:
		toprint = '{0:>12}'.format(track)
		for ptcut in ptcuts:
			res = trackdata[track][ptcut][width]
			toprint += '{0:7.4f} '.format(res)
		print toprint
	print
	for selector in selectors:
		for extra in extras[selector]:
			toprint = '{0:>12}'.format(selector+extra)
			for ptcut in ptcuts:
				res = bestdata['bestComb'][selector+extra][ptcut][width]
				toprint += '{0:7.4f} '.format(res)
			print toprint
