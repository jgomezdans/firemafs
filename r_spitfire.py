# -*- coding: utf-8 -*-
import numpy
import scipy
import spytfire

def RunSpytfire (  grid_id, year, theta, Datos, theta0 = numpy.array([0.39, 14.,  12.,  23.,  28.,  22.,  18.,  16.,   4.,   2.,\
                66.0, 3.58, 0.98]) ):

	#Obtain the parameters from subroutine...

	(pftpar, sla, tree, evergreen, summergreen, raingreen,\
		needle, boreal, lm_sapl, \
		sm_sapl, hm_sapl, rm_sapl, latosa, allom1, allom2,\
		allom3, allom4, wooddens, reinickerp) =\
                                spytfire.pftparameters()

	# theta is defined as [a_nd, fuel_bulk_density[pfts_present], sigma_1hr, sigma_10hr, sigma_100hr, alpha, moisture_extinction, moisture_conversion, duration_exponential]

	datos = deepcopy(Datos)
	a_nd = theta[0]
	cnt = 1
	for ipft in xrange(9):
		if datos['present'][ipft]==1:
			pftpar[ipft, 36] = theta[cnt]
			cnt+=1
		else:
			pftpar[ipft, 36] = theta0[ipft+1]
	sigma_1hr = theta[cnt] ; sigma_10hr = theta[cnt+1] ; sigma_100hr=theta[cnt+2]
	fbd_C3_livegrass = 4. ; fbd_C4_livegrass = 4.
	datos['a_nd'] = a_nd*numpy.ones((1,12),dtype=numpy.float32)

	( d_fdi, acflux_fire, mcflux_fire, afire_frac, num_fire, annum_fire, \
	area_burnt, an_areafires, mfdi,an_fdi,\
	an_fseason, mcflux_trace, acflux_trace, m_fc_crown, an_fc_crown, m_i_surface,\
	an_i_surface, dlm_lg, dlm_1hr, dlm_10hr, nesterov,  fire_durat, d_area_burnt, \
	d_numfire, ros_f, ros_b,lb,d_fuel_consumed, d_i_surface, cf) = \
	spytfire.fire ( year, pftpar, datos['dtemp'], datos['dtemp_min'], datos['dtemp_max'],\
	datos['dprec'], datos['dwindsp'], datos['lightn'], \
	datos['dphen'], datos['litter_ag'], datos['litter_bg'], \
	datos['fuel_1hr'], datos['fuel_10hr'], datos['fuel_100hr'], datos['fuel_1000hr'],\
	datos['lm_ind'], datos['rm_ind'], datos['sm_ind'], datos['hm_ind'], datos['nind'],\
	datos['dw1'],datos['present'],\
	tree, datos['lat'], datos['mw1'], datos['fpc_grid'], datos['popden'],\
	datos['a_nd'], datos['height'], datos['height_class'],\
	datos['dbh'],datos['tau_c'],datos['cl_t'], fbd_C3_livegrass, fbd_C4_livegrass,\
	sigma_1hr,sigma_10hr,sigma_100hr,\
 alpha, moisture_extinction, moisture_conversion, duration_exponential )

	return (area_burnt, dlm_lg, num_fire, ros_f, ros_b, dlm_10hr, fire_durat)