# -*- coding: utf-8 -*-
from demc import DEMC_sampler
import numpy
import scipy.stats

class CalibrateSPITFIRE ( DEMC_sampler ):
	def __init__ ( self, num_population, grid_id, year, CR=1.0, F=2.38, pSnooker=0.1, pGamma1=0.1, n_generations=1000, n_thin=5, n_burnin=200, eps_mult=0.1, eps_add=0):
		self.grid_id = grid_id
		self.year = year
		self._ReadData ( )
		
	def _ReadData ( self ):
		"""
		This function reads the observations data from the database, and stores the results in the class.
		"""
		s = SpitfireDB()
		self.datos = s.GetSpinUpVars ( self.grid_id, self.year )
		self.modis_ba = s.GetMODISBA ( self.grid_id, self.year )
		self.modis_num_fires = s.GetMODIS_FireNum ( self.grid_id, self.year )
		self.modis_fmc = s.GetMODISFMC ( self.grid_id, self.year )

	def likelihood_function ( self, theta ):
		"""
		Simple (log)likelihood function. It assumes independent distributions for annual mismatch between model and data. The newer version ought to have some MVG version of this, but
		it probably implies assuming stuff on priors...
		"""
		error_ba = self.modis_ba.sum()/10.
		error_fmc = self.modis_fmc.sum()/20.
		error_num_fires = self.modis_num_fires.sum()/20.
		(fwd_model_ba, fwd_model_fmc, fwd_model_numfires ) = spitfire ( self.grid_id, self.year, theta, self.datos )
		modis_ba_diff = ( fwd_model_ba - self.modis_ba ).sum()
		modis_fmc_diff = ( fwd_model_fmc - self.modis_fmc ).sum()
		modis_num_fires = ( fwd_modis_num_fires - self.modis_num_fires ).sum()
		p_ba = numpy.log ( scipy.stats.norm.pdf ( modis_ba_diff, 0.0, error_ba ))
		p_fmc = numpy.log ( scipy.stats.norm.pdf ( modis_ba_fmc, 0.0, error_fmc ))
		p_num_fires = numpy.log ( scipy.stats.norm.pdf ( modis_num_fires, 0.0, error_num_fires ))
		return p_ba + p_fmc + p_num_fires

	#def fwd_model ( self ):
		