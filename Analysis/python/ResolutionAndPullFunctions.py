import ROOT as R
import math as math

class gaus(object):
	def __init__(self,name,low,high):
		self.name = name
		self.func = R.TF1(name,self.__call__,low,high,3)
		self.func.SetParameters(1,0,0.05)
		self.func.SetParLimits(2,1E-5,10)
		self.labels = [('Normalization','norm'),('#mu','mu'),('#sigma','sigma'),('fit #chi^{2}','chi2')]
	
	def __getattr__(self,name):
		return getattr(self.func,name)

	def __call__(self,x,par):
		# par[0] = normalization
		# par[1] = mean
		# par[2] = sigma
		return par[0]*math.exp(-0.5*((x[0]-par[1])/par[2])**2)

class gaus_sum(object):
	def __init__(self,name,low,high):
		self.name = name
		self.func = R.TF1(self.name,self.__call__,low,high,5)
		self.func.SetParameters(1,0,0.05,1,3)
		self.func.SetParLimits(2,1E-5,10)
		self.func.SetParLimits(4,1E-5,10)
		self.labels = [('Normalization1','norm1'),('#mu','mu'),('#sigma1','sigma1'),\
				('Normalization2','norm2'),('#sigma2','sigma2'),('fit #chi^{2}','chi2')]

	def __getattr__(self,name):
		return getattr(self.func,name)

	def __call__(self,x,par):
		# par[0] = normalization gaussian 1
		# par[1] = mean gaussian 1 & 2
		# par[2] = sigma gaussian 1
		# par[3] = normalization gaussian 2
		# par[4] = sigma gaussian 2
		return par[0]*math.exp(-0.5*((x[0]-par[1])/par[2])**2) + \
				par[3]*math.exp(-0.5*((x[0]-par[1])/par[4])**2)

	# assumes the function has been fitted -- not sure if this works
	def sigma(self,x,par):
		return par[2]*par[2] + par[4]*par[4]

class dsCrystalBall(object):
	def __init__(self,name,low,high):
		self.name = name
		self.func = R.TF1(self.name,self.dsCrystalBall,low,high,5)

	def __getattr__(self,name):
		return getattr(self.func,name)
	
	def __call__(self,x,par):
		pass

class cruijff(object):
	def __init__(self,name,low,high):
		self.name = name
		self.func = R.TF1(self.name,self.__call__,low,high,5)
		self.func.SetParameters(1.,0.0,0.05,0.1,0.1)
		self.func.SetParLimits(2,1E-5,10)
		self.func.SetParLimits(3,0,10)
		self.func.SetParLimits(4,0,10)
		self.labels = [('Normalization','norm'),('#mu','mu'),('#sigma','sigma'),\
				('expR','expR'),('expL','expL'),('fit #chi^{2}','chi2')]

	def __getattr__(self,name):
		return getattr(self.func,name)

	def __call__(self,x,par):
		# par[0] = normalization
		# par[1] = gaussian mean
		# par[2] = gaussian sigma
		# par[3] = exp left tail
		# par[4] = exp right tail
		dx = x[0] - par[1]
		sigma = par[2]
		alpha = par[3] if dx<0 else par[4]
		f = 2*sigma*sigma + alpha*dx*dx
		return par[0] * math.exp(-dx*dx/f)
