import numpy as np
import ROOT as R
import math as math
import array as array


def printMatrix(matrix):
	for i in range(5):
		for j in range(5):
			matrix_ij = matrix[i,j]
			print '{matrix_ij:>11.4e}'.format(**locals()),
		print
	print

def printVector(vector):
	for i in range(5):
		vector_i = vector[i]
		print '{vector_i:>11.4e}'.format(**locals()),
	print

class pyTree(object):
	''' 
	Class dedicated to converting all stl vectors into numpy objects
	'''
	def __init__(self,t,DecList=[\
				'global','global_refit','global_refit_noUpdate','global_comb',\
				'picky','picky_refit','picky_refit_noUpdate','picky_comb',\
				'tracker','standAlone','gen']
			):
		bad = False
		trackList = DecList[:]
		trackList.remove('gen')
		for track in trackList:
			if len(getattr(t,track+'_par')) < 5: bad = True
			#if getattr(t,track+'_par')[0] < -998: bad = True
			setattr(self,track+'_par',  np.array(getattr(t,track+'_par')))
			setattr(self,track+'_cov',  np.matrix(getattr(t,track+'_cov')))
			setattr(self,track+'_w',    np.matrix(getattr(t,track+'_w')))
			setattr(self,track+'_corr', np.matrix(getattr(t,track+'_corr')))
			setattr(self,track+'_chi2', getattr(t,track+'_chi2'))
			if hasattr(t,track+'_nValidHits'):
				setattr(self,track+'_nValidHits', getattr(t,track+'_nValidHits'))
			if hasattr(t,track+'_nLostHits'):
				setattr(self,track+'_nLostHits', getattr(t,track+'_nLostHits'))
			if hasattr(t,track+'_pos'):
				setattr(self,track+'_pos', np.array(getattr(t,track+'_pos')))
			if hasattr(t,track+'_mom'):
				setattr(self,track+'_mom', np.array(getattr(t,track+'_mom')))
			if hasattr(t,track+'_ndof'):
				setattr(self,track+'_ndof', getattr(t,track+'_ndof'))
			if hasattr(t,track+'_dxyBS'):
				setattr(self,track+'_dxyBS', getattr(t,track+'_dxyBS'))
			if hasattr(t,track+'_dxyPV'):
				setattr(self,track+'_dxyPV', getattr(t,track+'_dxyPV'))
			if hasattr(t,track+'_dszBS'):
				setattr(self,track+'_dszBS', getattr(t,track+'_dszBS'))
			if hasattr(t,track+'_dszPV'):
				setattr(self,track+'_dszPV', getattr(t,track+'_dszPV'))
		if 'gen' in DecList:
			self.gen_K = t.gen_K
			self.gen_eta = t.gen_eta
			self.gen_phi = t.gen_phi
			self.gen_pt = t.gen_pt
		self.allOkay = True if bad==False else False
		self.run = t.run
		self.lumi = t.lumi
		self.event = t.event

class Track(object):
	def __init__(self,track,pyTree):
		if track not in [\
				'global','global_refit','global_comb','global_refit_noUpdate',\
				'picky','picky_refit','picky_comb','picky_refit_noUpdate',\
				'tracker','standAlone']:
			raise NameError(track+' is not a valid track name')
		self._trackType = track
		self.__parList = ['K','lambda','phi','dxy','dsz']
		self._par  = getattr(pyTree,track+'_par')
		self._cov  = getattr(pyTree,track+'_cov')
		self._w    = getattr(pyTree,track+'_w')
		self._corr = getattr(pyTree,track+'_corr')
		self._chi2 = getattr(pyTree,track+'_chi2')
		if track in ['global','picky','global_refit_noUpdate','picky_refit_noUpdate',\
				'tracker','standAlone']:
			self._nValidHits = getattr(pyTree,track+'_nValidHits')
			self._nLostHits = getattr(pyTree,track+'_nLostHits')
			self._pos = getattr(pyTree,track+'_pos')
			self._mom = getattr(pyTree,track+'_mom')
			self._ndof = getattr(pyTree,track+'_ndof')
			self._dxyBS = getattr(pyTree,track+'_dxyBS')
			self._dxyPV = getattr(pyTree,track+'_dxyPV')
			self._dszBS = getattr(pyTree,track+'_dszBS')
			self._dszPV = getattr(pyTree,track+'_dszPV')
		else:
			self._nValidHits = -1
			self._nLostHits = -1
			self._pos = np.array([-999,-999,-999])
			self._mom = np.array([-999,-999,-999])
			self._ndof = -1
			self._dxyBS = -999
			self._dxyPV = -999
			self._dszBS = -999
			self._dszPV = -999

	def trackType(self):   return self._trackType
	def par(self):    return self._par
	def K(self):      return self._par[0]
	def Lambda(self): return self._par[1]
	def phi(self):    return self._par[2]
	def dxy(self):    return self._par[3]
	def dxyBS(self):    return self._dxyBS
	def dxyPV(self):    return self._dxyPV
	def dsz(self):    return self._par[4]
	def dszBS(self):    return self._dszBS
	def dszPV(self):    return self._dszPV
	def cov(self,i,j):  return self._cov[i,j]
	def w(self,i,j):    return self._w[i,j]
	def corr(self,i,j): return self._corr[i,j]
	def eta(self):    
		if abs(self._par[1])>math.pi/2:
			return -999.
		else:
			return -math.log(math.tan((math.pi/2 - self._par[1])/2))
	def pt(self): return abs(1./self._par[0])*math.cos(self._par[1])
	def pos(self,i): return self._pos[i]
	def mom(self,i): return self._mom[i]
	def nValidHits(self): return self._nValidHits
	def nLostHits(self): return self._nLostHits
	def chi2(self): return self._chi2
	def ndof(self): return self._ndof

	def __getitem__(self,par):
		if par not in self.__parList:
			raise NameError(par+' is not a valid parameter name')
		else:
			return self.par[self.__parList.index(par)]

	def __str__(self):
		out =  '\nTrack Type : {self._trackType}\n'.format(**locals())
		out += '\n'
		out += 'Track Parameters\n'
		for p,par in enumerate(self.__parList):
			par_p = self._par[p]
			out += '{par:<6} = {par_p:11.4e}\n'.format(**locals())
		out += '\n'
		out += 'Track Covariance\n'
		for i in range(5):
			for j in range(5):
				cov_ij = self._cov[i,j]
				out += '{cov_ij:11.4e} '.format(**locals())
			out += '\n'
		out += 'Track pT  = {:11.4e}\n'.format(self.pt())
		out += 'Track eta = {:4.2f}\n'.format(self.eta())
		out += 'Number of valid hits {:2}\n'.format(self._nValidHits)
		out += 'Number of lost hits  {:2}\n'.format(self._nLostHits)
		return out

def pos_diff(pos1,pos2):
	return math.sqrt(sum((pos1-pos2)**2))
def pos_diff_T(pos1,pos2):
	return math.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)

def draw_overflow(hist):
    # Function to paint the histogram with an extra bin for overflow
    nx = hist.GetNbinsX()+1
    xbins = array.array('d',[])
    for ibin in range(nx+1):
        if ibin==nx+1: xbins.append(xbins[nx-1]+hist.GetBinWidth(nx))
        else: xbins.append(hist.GetBinLowEdge(ibin+1))
    # Book a temporary histogram having extra bins for overflow
    htmp = R.TH1F(hist.GetName()+'_o', hist.GetTitle(), nx, xbins)
    htmp.Sumw2()
    # Fill the new histogram including the overflow
    for ibin in range(1,nx+1):
        htmp.SetBinContent(htmp.FindBin(htmp.GetBinCenter(ibin)),hist.GetBinContent(ibin))
        htmp.SetBinError(htmp.FindBin(htmp.GetBinCenter(ibin)),hist.GetBinError(ibin))
    htmp.SetBinContent(htmp.FindBin(htmp.GetBinCenter(1)-1),hist.GetBinContent(0))
    htmp.SetBinError(htmp.FindBin(htmp.GetBinCenter(1)-1),hist.GetBinError(0))
    # Restore the number of entries
    htmp.SetEntries(hist.GetEffectiveEntries())
    return htmp

def cumulative_histogram(h, type='ge'):
    """Construct the cumulative histogram in which the value of each
    bin is the tail integral of the given histogram.
    """
    
    nb = h.GetNbinsX()
    name = h.GetName()+'_cumulative_'+type
    hc = R.TH1F(name,'',nb,h.GetXaxis().GetXmin(),h.GetXaxis().GetXmax())
    hc.Sumw2()
    if type == 'ge':
        first, last, step = nb+1, 0, -1
    elif type == 'le':
        first, last, step = 0, nb+1, 1
    else:
        raise ValueError('type %s not recognized' % type)
    for i in xrange(first, last, step):
        prev = 0 if i == first else hc.GetBinContent(i-step)
        c = h.GetBinContent(i) + prev
        hc.SetBinContent(i, c)
        if c > 0:
            hc.SetBinError(i, c**0.5)
        else:
            hc.SetBinError(i, 0.)
    return hc

def poisson_interval(nobs, alpha=(1-0.6827)/2, beta=(1-0.6827)/2):
    lower = 0
    if nobs > 0:
        lower = 0.5 * R.Math.chisquared_quantile_c(1-alpha, 2*nobs)
    elif nobs == 0:
        beta *= 2
    upper = 0.5 * R.Math.chisquared_quantile_c(beta, 2*(nobs+1))
    return lower, upper

def clopper_pearson(n_on, n_tot, alpha=1-0.6827, equal_tailed=True):
    if equal_tailed:
        alpha_min = alpha/2
    else:
        alpha_min = alpha

    lower = 0
    upper = 1

    if n_on > 0:
        lower = R.Math.beta_quantile(alpha_min, n_on, n_tot - n_on + 1)
    if n_tot - n_on > 0:
        upper = R.Math.beta_quantile_c(alpha_min, n_on + 1, n_tot - n_on)

    if n_on == 0 and n_tot == 0:
        return 0, lower, upper
    else:
        return float(n_on)/n_tot, lower, upper

def clopper_pearson_poisson_means(x, y, alpha=1-0.6827):
	r, rl, rh = clopper_pearson(x, x+y, alpha)
	return r/(1 - r), rl/(1 - rl), rh/(1 - rh)

def kalman_filter_update(x,P,Z,R):
	'''
	By hand implementation of global track position and direction Kalman Filter
	update of muon only refits
	x = state of muon only refit
	P = Covariance of muon only refit
	Z = (full) measurement state of global track
	R = (full) covariance of global track
	'''
	H = np.matrix('0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1')
	HT = H.transpose()

	# y = H*Z - H*x # differs because Z is 5 par
	y = Z[1:].reshape(4,1) - x[1:].reshape(4,1)
	# S = H*R*HT + H*P*HT # differs because R is 5x5
	S = R[1:,1:] + P[1:,1:]
	SInv = np.linalg.inv(S)

	K = P*HT*SInv
	x_update = x + np.array(K*y).reshape(1,5)[0]
	P_update = P - K*H*P
	return x_update,P_update

def kalman_prefit_residuals(x,P,Z,R):
	'''
	Return only pre-fit residuals (exclude the q/P estimate!)
	only [lambda, phi, dxy, dsz]
	'''
	y = Z[1:] - x[1:]
	S = R[1:,1:] + P[1:,1:]
	return y.reshape(1,4)[0],S
	
def chisq_comb(x1,C1,x2,C2):
	W1 = np.linalg.inv(C1)
	W2 = np.linalg.inv(C2)
	C = np.linalg.inv(W1+W2)
	P = np.array(C*(W1*(x1.reshape(5,1))+W2*(x2.reshape(5,1)))).reshape(1,5)[0]
	return P,C


##########################################################################

if __name__=='__main__':
	import math as math
	# interactive testing
	f = R.TFile('testing.root')
	t = f.Get('t')
	for e,entry in enumerate(t):
		if e>0: break
		P = pyTree(t)
		globalTrack = Track('global',P)
		print globalTrack
		print '*'*30
		globalCombTrack = Track('global_comb',P)
		print globalCombTrack
		print '*'*30
		pickyTrack = Track('picky',P)
		print pickyTrack
		print '*'*30
		pickyCombTrack = Track('picky_comb',P)
		print pickyCombTrack
		print '*'*30


