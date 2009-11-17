# -*- coding: utf-8 -*-
import numpy, pdb
def choose_without_replacement(m,n,repeats=None):
	"""Choose n nonnegative integers less than m without replacement
	Returns an array of shape n, or (n,repeats).
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
def fitness ( u ):
		x1 = u[0] ; x2 = u[1]
		sigma1=1.0
		sigma2=1.0
		mu1=-5.0
		mu2=5.0
		ret=1.0/numpy.sqrt(2*numpy.pi*sigma1**2) * numpy.exp(-1/(2*sigma1**2) * (mu1-x1)**2)\
			*1/numpy.sqrt(2*numpy.pi*sigma2**2) * numpy.exp(-1/(2*sigma2**2) * (mu2-x2)**2) ;
		return numpy.log(ret)
		
def demc_zs ( num_population, Z, CR=1.0, F=2.38, pSnooker=0.1, pGamma1=0.1, n_generations=10000, n_thin=5, n_burnin=2000, eps_mult=0.1, eps_add=0):
	#Z = numpy.zeros ( d, m0)
	#X = numpy.zeros (d, num_population)
	d = 2
	X = Z[:,:num_population]
	m0 = Z.shape[1]
	mZ = Z.shape[1]
	Npar = X.shape[0]
	Npar12 = (Npar-1)/2. # Factor for Metropolis ratio DE snooker update
	F2 = F/numpy.sqrt ( 2.*Npar)
	F1 = 1.0
	accept = numpy.zeros ( n_generations )
	#iseq = numpy.arange(1, num_population)
	rr = 0.0 ; r_extra = 0
	print int(m0+num_population*numpy.floor(n_burnin/n_thin))
	logfitness_x = [ fitness(X[:,i]) for i in xrange(num_population) ]
	for iteration in xrange(n_generations):
		accepti = 0
		for i in xrange(num_population):
			if (numpy.random.random()<pSnooker):
				rr = choose_without_replacement(mZ-1, 3, repeats=None)
				z = Z[:,rr[2]]
				x_z = X[:,i] - z
				D2 = max(numpy.sum(x_z*x_z), 1.0e-300)
				gamma_snooker = numpy.random.random()+1.2
				proj_diff = numpy.sum((Z[:,rr[0]] - Z[:,rr[1]])*x_z)/D2
				x_prop = X[:,i] + (gamma_snooker * proj_diff) * x_z
				x_z = x_prop - z
				#pdb.set_trace()
				D2prop = max( numpy.dot(x_z, x_z), 1.0e-30)
				r_extra = Npar12*(numpy.log(D2prop) - numpy.log(D2))
			else:
				if (numpy.random.random()<pGamma1):
					gamma_par = F1
				else:
					gamma_par = F2 * numpy.random.uniform( low=1-eps_mult, high=1+eps_mult, size=Npar)
				rr = choose_without_replacement(mZ-1, 2, repeats=None)
				if (eps_add==0):
					x_prop = X[:,i] + gamma_par * ( Z[:,rr[0]] - Z[:,rr[1]])
				else:
					x_prop = X[:,i] + gamma_par * ( Z[:,rr[0]] - Z[:,rr[1]]) + eps_add*numpy.random.randn(Npar)
				r_extra = 0
			logfitness_x_prop = fitness ( x_prop )
			logr = logfitness_x_prop - logfitness_x[i]
			if (logr + r_extra)>numpy.log ( numpy.random.random()):
				accepti += 1
				X[:,i] = x_prop
				logfitness_x[i] = logfitness_x_prop
		accept[iteration] = accepti
		#print logfitness_x[i], logfitness_x_prop ,logr,r_extra,accepti
		if iteration%n_thin==0:
			Z = numpy.c_[X,Z]
			mZ = Z.shape[1]
			
		#if iteration%thin==0:
			#list(Draws= Z[,-(1:(M0 + Npop* floor(n.burnin/n.thin)))] , accept.prob= accept/Npop, X.final = X, logfitness.X.final = logfitness_X)
	return (Z[:,:int(m0+num_population*numpy.floor(n_burnin/n_thin))], accept/num_population)

if __name__=="__main__":
	import pylab
	num_population = 100
	Z = numpy.c_[numpy.random.random(150)*20-10, numpy.random.random(150)*20-10].T
	(Zp, accept) = demc_zs ( num_population, Z )
	pylab.hist ( Zp[0,:], bins=30)
	pylab.hist ( Zp[1,:], bins=30)
	pylab.show()