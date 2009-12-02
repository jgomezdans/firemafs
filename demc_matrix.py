# -*- coding: utf-8 -*-
import numpy
import numpy.linalg
import scipy.stats
import sys, pdb
from pymc import inverse_wishart_like, rinverse_wishart

def full_gauss_den(x, mu, va, log):
	""" This function is the actual implementation
	of gaussian pdf in full matrix case.

	It assumes all args are conformant, so it should
	not be used directly Call gauss_den instead

	Does not check if va is definite positive (on inversible
	for that matter), so the inverse computation and/or determinant
	would throw an exception."""
	d       = mu.size
	inva    = numpy.linalg.inv(va)
	fac     = 1 / numpy.sqrt( (2*numpy.pi) ** d * numpy.fabs(numpy.linalg.det(va)))

	# we are using a trick with sum to "emulate"
	# the matrix multiplication inva * x without any explicit loop
	#y   = -0.5 * N.sum(N.dot((x-mu), inva) * (x-mu), 1)
	y   = -0.5 * numpy.dot(numpy.dot((x-mu), inva) * (x-mu),
			numpy.ones((mu.size, 1), x.dtype))[:, 0]

	if not log:
		y   = fac * numpy.exp(y)
	else:
		y   = y + numpy.log(fac)

	return y

class DEMC_sampler(object):
	"""
	=======================
	Python DEMC sampler
	=======================
	A python implementation of the Differential Evoluation MCMC sampler of ter Braak (2008) (ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain
	with snooker updater and fewer chains. Statistics and Computing http://dx.doi.org/10.1007/s11222-008-9104-9). It appears that this implementation (and its daughter, DREAM) provide a really nice and efficient MCMC sampling scheme. An additional plus is that the algorithm can be easily parallelised.

	The DEMC_sampler object *should not be used directly* (only for testing purposes), but it should be subclassed. Subclassing implies adding your own likelihood function, rather than the provided likelihoood function. The following python packages are required
		- numpy ( & linalg, usually comes as standard)
		- scipy (we need the stats objects from there)

	An example run


	Class usage is very simple. The following example defines a DEMC sampler

	.. python::
		DEMC = DEMC_sampler ( 10, n_generations=1000 )
		parameter_list=[['x1', 'scipy.stats.uniform(-20, 40)'], ['x2', 'scipy.stats.uniform(-20, 40)']]
		parameters = ['x1','x2']
		DEMC.prior_distributions ( parameter_list, parameters )
		Z = DEMC.ProposeStartingMatrix ( 150 )
		(Z_out, accept_rate) = DEMC.demc_zs ( Z )

	The first call instantiates the DEMC_sampler class. It is initialised with population size (10 in this example), and with the number of generations (or iterations).  The second line defines a list of pairs of parameter names, distributions. The distributions are distribution objects from scipy.stats, or code that implements the same methods that the user provides. Note that we do not provide a reference to these objects, but write the reference as a string (the reference is obtained internally). We also provide a list with the parameter names in the same order they are required in the M{\theta} parameter vector. These two lists are then included in the object by calling the prior_distributions() method. The sampler is initialised with a M{Z} matrix (num_params x [...]), sampled from the prior distributions. Or, you can provide your own. Finally, the sampler is started and it returns the Z_out matrix (the samples for the different parameters) and the accepctance rate (try to get it ~0.2-0.4)
	"""
	def __init__ ( self, num_population, CR=1.0, F=2.38, pSnooker=0.1, pGamma1=0.01, n_generations=1000, n_thin=5, n_burnin=200, eps_mult=0.1, eps_add=0) :
		"""
		The class creator. You should pass it the populations (say between 10 or 100), and the sampler parameters, such as n_generations, n_thing and n_burnin. The other parameters are DEMC-related parameters and shouldn't need changing accroding to Vrugt and ter Braak. CR is a crossover rate, and at present, is not used in this implementation of the code (but could be used in a future implementation that does DREAM)
		"""
		self.num_population = num_population
		self.CR = CR
		self.F = F
		self.pSnooker = pSnooker
		self.pGamma1 = pGamma1
		self.n_generations = n_generations
		self.n_thin = n_thin
		self.n_burnin =n_burnin
		self.eps_mult = eps_mult
		self.eps_add = eps_add
		self.CreateObservations ()
		
	def CreateObservations ( self, x=3, b=0.1):
		x=numpy.r_[-5:5:.2]
		a=3
		b=0.1
		x1=a*x+b
		x2=a*x+b
		mu=numpy.array([0,0])
		cov=numpy.array([[1.,0.1*numpy.sqrt(3.2)],[0.1*numpy.sqrt(3.2),3.2]])
		d=numpy.random.multivariate_normal (mu,cov,size=x1.shape)
		X1=x1+d[:,0] #+(numpy.random.rand(x1.shape[0])-0.5)*8#Actually x1==x2, but nevermind!
		X2=x2+d[:,1] #+(numpy.random.rand(x1.shape[0])-0.5)*8
		self.observations = numpy.c_[ X1, X2 ]
		print self.observations
		
	def choose_without_replacement(self, m,n,repeats=None):
		"""Choose n nonnegative integers less than m without replacement
		Returns an array of shape n, or (n,repeats).

		This code is from Anne Archiebald, as numpy/scipy don't seem to have this facility.... It looks like it might not be the most efficient code in the world
		"""
		if repeats is None:
			r = 1
		else:
			r = repeats
		if n>m:
			raise ValueError, "Cannot find %d nonnegative integers less than %d" %  (n,m)
		if n>m/2:
			res = numpy.sort(numpy.random.rand(m,r).argsort(axis=0)[:n,:],axis=0)
		else:
			res = numpy.random.random_integers(m,size=(n,r))
			while True:
				res = numpy.sort(res,axis=0)
				w = numpy.nonzero(numpy.diff(res,axis=0)==0)
				nr = len(w[0])
				if nr==0:
					break
				res[w] = numpy.random.random_integers(m,size=nr)
		if repeats is None:
			return res[:,0]
		else:
			return res

	def prior_distributions ( self, parameter_list, parameters ):
		"""This function defines the prior distributions from a (python) list. Useful to quickly add/update these parameters. The parameters list (2nd argument) is there to have the parameter names in the same order that vector M{\Theta} will have. So if theta required by the model is [ par1, par2, par3], then parameters=['par1','par2','par3']

		 Example parameter list:
		 parameter_list = [['fuel_1hr',"scipy.stats.uniform ( 0,100.)"],\
		 ['fuel_10hr', "scipy.stats.uniform ( 0,100.)"],\
		 ['fuel_100hr', "scipy.stats.uniform ( 0,100.)"] ]

		 The way this method works is by using setattr to add the prior methods to the self object. Note that in this particular implementation, the parameter distributions are independent (i.e., you can't have correlation between parameters). If this is an issue, you'll need to change the code.

		 In Bayesian statistics, the prior distributions encapsulate whatever knowledge we have about the distribution of a parameter prior to doing any experiment. On the one hand, they add a degree of subjectiveness, as different people may translate their knowledge in different distribution shapes. Also, it is sometimes useful to use non-informative priors to indicate ignorance or "objectivity" (yeah, right...). In practice, broad uniform distributions, truncated normals, lognormal distributions are often used. In sequential applications of a calibration exercise, the posterior of the previous run can be used as the prior for the current state.
		"""
		self.parameters = parameters # List is in the order of the parameter vector
		for [k,v] in parameter_list:
			setattr(self, k, eval( v ))
	def Wishart_Params ( self, DoF, Tau ):
		"""
		This method sets the parameters for the Wishart/Inverse
		Wishart prior: the degrees of freedom (n) and the
		(inverse) scale matrix (p). I think that n>p-1,
		and Tau has to be positive definite.
		"""
		assert ( DoF>(Tau.shape[0]-1))
		self.wishart_dof = DoF
		self.Tau = Tau
		
	def prior_probabilities ( self, theta, W ):
		""" The method that calculates the prior (log) probabilities. This is based on the prior distributions given in prior_distributions, and assumes independence, so we just add them up. 
		"""
		p = numpy.array([ numpy.log ( getattr ( self, self.parameters[i]).pdf ( theta[i])) for i in xrange(len(self.parameters)) ]).sum()
		p += inverse_wishart_like ( W, self.wishart_dof, self.Tau )
		#if numpy.isneginf(p):
			#p = numpy.log(1.0E-300)
		if p<1.0e-100:
			p = -numpy.infty
			
		return p

	def likelihood_function ( self, theta, W):

		x = numpy.r_[-5:5:.2]
		fwd = theta[0]*x + theta[1]
		D = fwd[:,numpy.newaxis]-self.observations
		
		p = - numpy.log(2.*numpy.pi)  - numpy.log(numpy.sqrt ( numpy.linalg.det (W))) -\
				0.5*numpy.array([ numpy.dot(numpy.dot(D[i,:],numpy.linalg.inv(W)),
				D[i,:]) for i in xrange(D.shape[0])]).sum()
		if numpy.isneginf(p):
			p = numpy.log(1.0E-300)
		return p

	def fitness ( self, theta, W ):
		"""
		The new posterior probability in log. Convenience, really
		"""
		return self.likelihood_function ( theta,W ) + self.prior_probabilities ( theta, W )
		
	###def MonitorChains ( self, x ):
		###alpha = 0.05                     # 95% intervals
		###(n, m) = x.shape
		###xdot = numpy.mean ( x, axis=0)
		###s2 = numpy.var ( x, axis=0)
		###W = numpy.mean(s2)
		###B = n*numpy.var(xdot)
		###muhat = numpy.mean(xdot)
		###varW = numpy.var ( s2 )/m
		###varB = (B**2)*2/(m-1.)
		###covWB = (n/m)*(numpy.cov ( s2, xdot**2)-2.*muhat*numpy.cov(s2,xdot))
		###sig2hat = ((n-1)*W+B)/n
		###quantiles =[ scipy.stats.scoreatpercentile ( x.flatten(1), percentile) for percentile in [2.5, 25., 50., 75., 97.5] ]
		###quantiles = numpy.array ( quantiles )
		####pdb.set_trace()
		###if (W>1.e-8): # Non-degenerate case
			###postvar = sig2hat + B/(m*n)
			###varpostvar = max ( 0, ((((n-1)**2)*varW + (1+1/m)**2*varB + 2.*(n-1)*(1+1/m)*covWB)/n**2)[0,0]) # CHECK!
			###post_df = min ( 2*(postvar**2/varpostvar), 1000)
			###post_range = muhat + numpy.sqrt ( postvar ) *scipy.stats.t.ppf ( 1-alpha/2, post_df)*numpy.array([-1,0,1])
			###verlo_df = 2*(W**2/varW)
			###confshrink_range = numpy.sqrt  ( numpy.array([ postvar/W, (n-1)/n + (1+1/m)*(1/n)*(B/W)*scipy.stats.f.ppf ( 97.5, m-1, verlo_df)])*(post_df+3)/(post_df+1))
			###n_eff = m*n*min(sig2hat/B,1)
			###print post_range
			###print quantiles
			###print confshrink_range
			###print n_eff
			###return ( post_range, quantiles, confshrink_range, n_eff)
		###else: # Degenerate case
			###return ( numpy.ones(3), quantiles, numpy.ones(2),1)
			
			#varlo.df <- chisqdf (W, varW)
			#confshrink.range <- sqrt (c(postvar/W,
			#(n-1)/n + (1+1/m)*(1/n)*(B/W) * qf(.975, m-1, varlo.df)) *
			#(post.df+3)/(post.df+1))
	def MonitorChains ( self, X ):
		(Npar, I, T2) = X.shape # (number of parameters, number of chains, number of iterations)
		T = T2/2 # Only use the latter half of the itarations. Deals with burn-in issues
		#We do it for each parameter, I guess?
		rhat = 0.
		for par in xrange(Npar):
			x = X[par, :, (T+1):] # For convenience, subset per parameter. x is now (chains, iterations)
			psi_i = numpy.mean ( x,axis=1) # Check axis!
			psi = numpy.mean ( psi_i )
			B = numpy.sum ( (psi_i-psi)*(psi_i-psi))*(1/(I-1.))
			S = numpy.mean( (x-psi_i[:,numpy.newaxis])**2,axis=1)
			W = numpy.mean(S)
			V_hat = ((T-.1)/T)*W + (1+1./I)*B
			#pdb.set_trace()
			if rhat<numpy.sqrt(V_hat/W):
				rhat = numpy.sqrt(V_hat/W)
			print "+--------------------------------------------------------------------------------------------------+"
			print "| %12.4g "%(numpy.sqrt(V_hat/W)),
			for percentile in [2.5, 25, 50, 75, 97.5]:
				print "|%8.2g "%(scipy.stats.scoreatpercentile ( x.flatten(1), percentile)),
			#print "%8g | %8g | %8g | %8g | %8g" % [scipy.stats.scoreatpercentile ( x.flatten(1), percentile) for percentile in [2.5, 25., 50., 75., 97.5] ]
			print "| %8.2g | %8.2g"%(x.flatten(1).mean(), x.flatten(1).std())
		print "==============================================================================="
		return rhat
			
			
			
		

	def ProposeStartingMatrix ( self, m0):
		"""
		The  proposed starting matrix. We initialise the chain with this matrix. Usually, a draw from the prior distribution is suitable
		"""
		if m0<=(self.num_population+len(self.parameters)):
			print "The size of the Z matrix needs to be larger than the population size + the number of parameters"
			sys.exit()
		Z = []
		for i in xrange(len(self.parameters)):
			Z.append(getattr ( self, self.parameters[i]).rvs ( size = m0))
		Zw = [ rinverse_wishart ( self.wishart_dof, self. Tau) for i in xrange(m0)]
		Zw = numpy.array(Zw)
		Zw=numpy.rollaxis(numpy.rollaxis(Zw,1),-1)
		return numpy.array(Z), Zw

	def demc_zs ( self,  Z, Zw):
		#CR=self.CR, F=self.F, pSnooker=0.1, pGamma1=0.1, n_generations=10000, n_thin=5, n_burnin=2000, eps_mult=0.1, eps_add=0
		#Z = numpy.zeros ( d, m0)
		#X = numpy.zeros (d, num_population)
		#pSnooker=0.
		d = 2
		X = Z[:,:self.num_population]
		Xw = Zw[:,:,:self.num_population]
		# We want to get the lower diagonal and main diagonal elements.
		# Not the cleanest implementation...
		diag_elements = numpy.nonzero (numpy.tril( numpy.ones(Xw[:,:,0].shape)))
		m0 = Z.shape[1]
		self.discard = int(m0+self.num_population*numpy.floor(self.n_burnin/self.n_thin))
		mZ = Z.shape[1]
		Npar = X.shape[0]
		Npar = 2
		Npar12 = (Npar-1)/2. # Factor for Metropolis ratio DE snooker update
		Npar12 = 2.
		F2 = self.F/numpy.sqrt ( 2.*Npar)
		F2 = self.F/numpy.sqrt(10.)
		F1 = 1.0
		accept = numpy.zeros ( self.n_generations )
		#iseq = numpy.arange(1, num_population)
		rr = 0.0 ; r_extra = 0#numpy.log(0.)
		#print int(m0+num_population*numpy.floor(n_burnin/n_thin))
		#Calculate the starting posteriors...
		logfitness_x = [ self.fitness(X[:,i], Xw[:,:,i]) \
								for i in xrange(self.num_population) ]
		# The diagnostic matrix
		#Augmented by the relevant elements of the covariance matrix
		Z_diagnostic = numpy.zeros ( (Npar+diag_elements[0].shape[0],\
				 self.num_population, self.n_generations))
		#Start of main loop
		iteration = -1
		#passer is our variable to add some hysterisis to convergence
		passer = 0
		while True:
			# We clear the acceptance counter
			accepti = 0
			for i in xrange(self.num_population):
				#Start of different chains loop
				#First decide whether this is a snooker update or not
				#pdb.set_trace()
				if (numpy.random.random()<self.pSnooker):
					#Snooker update
					#Select three chains
					while True:
						rr = self.choose_without_replacement(mZ, 3, repeats=None)
						rr-=1
						if not (numpy.any(rr==i)): break
					z = Z[:,rr[2]]
					zw = Zw[:,:,rr[2]]
					x_z = X[:,i] - z
					x_zw =Xw[:,:,i] - zw
					#Difference between the current point and one of the 3 chains
					D2 = max(numpy.sum(x_z*x_z), 1.0e-300) # This is the distance. Could do it with dot?
					D2w = max(numpy.sum(x_zw*x_zw), 1.0e-300)
					gamma_snooker = numpy.random.random()+0.2 # Snooker stochastic
					proj_diff = numpy.sum((Z[:,rr[0]] - Z[:,rr[1]])*x_z)/D2 # Project the difference onto x_z. Normalize by x_z's norm
					proj_diff_w = numpy.sum((Zw[:,:,rr[0]] - Zw[:,:,rr[1]])*x_zw)/D2w # Project the difference onto
					x_prop = X[:,i] + (gamma_snooker * proj_diff) * x_z # Proposed point
					xw_prop = Xw[:,:,i] + (gamma_snooker * proj_diff_w) * x_zw # Proposed matrix
					x_z = x_prop - z # update x_z
					x_zw = xw_prop - zw # update x_z
					#pdb.set_trace()
					D2prop = max( numpy.dot(x_z, x_z), 1.0e-30) # Calculate D2prop
#					D2propw = max( numpy.dot(x_zw, x_zw), 1.0e-30) # Calculate D2prop
					r_extra = Npar12*(numpy.log(D2prop) - numpy.log(D2))
					
				else:
					if (numpy.random.random()<self.pGamma1):
						gamma_par = F1
					else:
						gamma_par = F2 * numpy.random.uniform( low=1-self.eps_mult, high=1+self.eps_mult, size=Npar)
					
					while True:
						rr = self.choose_without_replacement(mZ, 2, repeats=None)-1
						if not ( numpy.any(rr==i)): break
					#pdb.set_trace()
					if (self.eps_add==0):
						x_prop = X[:,i] + gamma_par * ( Z[:,rr[0]] - Z[:,rr[1]])
						xw_prop = Xw[:,:,i] + gamma_par * ( Zw[:,:,rr[0]] - Zw[:,:,rr[1]])
					else:
						x_prop = X[:,i] + gamma_par * ( Z[:,rr[0]] - Z[:,rr[1]]) + self.eps_add*numpy.random.randn(Npar)
						
						xw_prop = Xw[:,:,i] + gamma_par * ( Zw[:,:,rr[0]] - Zw[:,:,rr[1]]) + self.eps_add*numpy.tril(numpy.random.randn(self.Tau.flatten().shape[0]).reshape(self.Tau.shape))
					#r_extra = 0
				#W_test=numpy.array([[1,-0.8*numpy.sqrt(3.2)],[-0.8*numpy.sqrt(3.2),3.2]])
				logfitness_x_prop = self.fitness ( x_prop, xw_prop )
				logr = logfitness_x_prop - logfitness_x[i]
				if (logr + r_extra)>numpy.log ( numpy.random.random()):
					accepti += 1
					X[:,i] = x_prop
					Xw[:,:,i] = xw_prop
					logfitness_x[i] = logfitness_x_prop
				#Store the samples in the diagnostic array
				#Note that we also add the covar matrix terms.
				Z_diagnostic [:, i, iteration] = numpy.r_[ X[:,i], \
						Xw[diag_elements[0], diag_elements[1],i] ]
				
			accept[iteration] = accepti
			iteration += 1
			#print logfitness_x[i], logfitness_x_prop ,logr,r_extra,accepti,x_prop
			#print iteration, accepti,X.mean(axis=1),Xw.mean(axis=2)
			
			if iteration%self.n_thin==0:
				Z = numpy.c_[X,Z]
				Zw = numpy.c_[Xw,Zw]
				mZ = Z.shape[1]
				#print "iteration ",iteration
				#print "\t",X[0,:]
				#print "\t",X[1,:]
			if (iteration% (self.n_generations)==0) and iteration>0:
				rhat = self.MonitorChains ( Z_diagnostic)
				print iteration, rhat
				new_pop = iteration
				if rhat<=1.3:
					passer +=1
				if (passer>=2) and (rhat>1.3):
					passer -=1
				if (rhat<=1.3) and (passer>3):
					break
				
				accept = numpy.append ( accept, numpy.zeros (self.n_generations+1), 0)
				Z_diagnostic = numpy.append(Z_diagnostic, \
					numpy.zeros ((Npar+diag_elements[0].shape[0], \
					self.num_population, self.n_generations+1)), 2)
		return (Z, Zw, accept/self.num_population)
		#return (Z[:,:int(m0+iteration*numpy.floor(self.n_burnin/self.n_thin))], accept/self.num_population)

if __name__=="__main__":
	DEMC = DEMC_sampler ( 7, n_generations=1000, n_burnin=1, n_thin=1)
	parameter_list=[['x1', 'scipy.stats.norm(3, 1)'], ['x2', 'scipy.stats.uniform(0, 1)']]
	parameters = ['x1','x2']
	DEMC.prior_distributions ( parameter_list, parameters )
	cov=numpy.array([[1,0.8*numpy.sqrt(3.2)],[0.8*numpy.sqrt(3.2),3.2]])
	#cov = numpy.array([[1.,0.1],[0.1,1]])
	DEMC.Wishart_Params ( 2, cov)
	(Z, Zw) = DEMC.ProposeStartingMatrix ( 35 )
	(Z_out, Zw_out, accept_rate) = DEMC.demc_zs ( Z, Zw )
