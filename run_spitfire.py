import spytfire
import numpy,pdb
from RetrieveDB import SpitfireDB
from copy import deepcopy

def RunSpytfire (  grid_id, year, theta, Datos ):
	(pftpar, sla, tree, evergreen, summergreen, raingreen,\
		needle, boreal, lm_sapl, \
		sm_sapl, hm_sapl, rm_sapl, latosa, allom1, allom2,\
		allom3, allom4, wooddens, reinickerp) =\
                                spytfire.pftparameters()
	
	#s = SpitfireDB()
	
	datos = deepcopy(Datos)
	a_nd = theta[0]
	#pdb.set_trace()
	cnt = 1
	for ipft in xrange(9):
		if (ipft+1) in datos['present']:
			pftpar[ ipft,36 ] = theta [ cnt ]
			cnt += 1	
	####pftpar[0,36] = theta[1] # Fuel Bulk Density
	####pftpar[1,36] = theta[2]
	####pftpar[7,36] = theta[3]
	####pftpar[8,36] = theta[4]
	sigma_1hr = theta[cnt] ; sigma_10hr = theta[cnt+1] ; sigma_100hr=theta[cnt+2]
	#moistfactor_livegrass = theta[13] ; moistfactor_1hr = theta[14]
	####moistfactor_10hr = theta[15]
	####moistfactor_100hr = theta[16]
	####moistfactor_1000hr = theta[17]
	moistfactor_livegrass = 0.2 ; moistfactor_1hr = 0.404
	moistfactor_10hr = 0.487
	moistfactor_100hr = 0.525
	moistfactor_1000hr = 0.544
	fbd_a = 1.2 ; fbd_b = 1.4
	fbd_C3_livegrass = 4. ; fbd_C4_livegrass = 4.
	#fbd_a = theta[1] ; fbd_b = theta[2] #1.2, 1.4
	#fbd_C3_livegrass = theta[3] ; fbd_C4_livegrass = theta[4] # 4.0, 4.0
	#datos = s.GetSpinUpVars ( grid_id, year )
	#modis_ba = s.GetMODISBA ( grid_id, year )
	datos['a_nd'] = a_nd*numpy.ones((1,12),dtype=numpy.float32)
	#datos['lightn'] =numpy.array( datos['lightn'])*kappa
	#pdb.set_trace()
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
	datos['dbh'],datos['tau_c'],datos['cl_t'], fbd_a, fbd_b, fbd_C3_livegrass,fbd_C4_livegrass,\
	sigma_1hr,sigma_10hr,sigma_100hr,\
	moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr,\
	moistfactor_100hr,moistfactor_1000hr ) 
	
	return (area_burnt)


def RunSpytfire_aNd (  grid_id, year, theta, Datos ):
	(pftpar, sla, tree, evergreen, summergreen, raingreen,\
		needle, boreal, lm_sapl, \
		sm_sapl, hm_sapl, rm_sapl, latosa, allom1, allom2,\
		allom3, allom4, wooddens, reinickerp) =\
                                spytfire.pftparameters()
	
	#s = SpitfireDB()
	theta0 = numpy.array([0.39, 14.,  12.,  23.,  28.,  22.,  18.,  16.,   4.,   2.,\
		66.0, 3.58, 0.98, 0.2, 0.404, 0.487, 0.525, 0.544])
	datos = deepcopy(Datos)
	a_nd = theta
	#pdb.set_trace()
	pftpar[:,36] = theta0[1:10] # Fuel Bulk Density
	sigma_1hr = theta0[10] ; sigma_10hr = theta0[11] ; sigma_100hr=theta0[12]
	moistfactor_livegrass = theta0[13] ; moistfactor_1hr = theta0[14]
	moistfactor_10hr = theta0[15]
	moistfactor_100hr = theta0[16]
	moistfactor_1000hr = theta0[17]
	
	fbd_a = 1.2 ; fbd_b = 1.4
	fbd_C3_livegrass = 4. ; fbd_C4_livegrass = 4.
	#fbd_a = theta[1] ; fbd_b = theta[2] #1.2, 1.4
	#fbd_C3_livegrass = theta[3] ; fbd_C4_livegrass = theta[4] # 4.0, 4.0
	#datos = s.GetSpinUpVars ( grid_id, year )
	#modis_ba = s.GetMODISBA ( grid_id, year )
	datos['a_nd'] = a_nd*numpy.ones((1,12),dtype=numpy.float32)
	#datos['lightn'] =numpy.array( datos['lightn'])*kappa
	#pdb.set_trace()
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
	datos['dbh'],datos['tau_c'],datos['cl_t'], fbd_a, fbd_b, fbd_C3_livegrass,fbd_C4_livegrass,\
	sigma_1hr,sigma_10hr,sigma_100hr,\
	moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr,\
	moistfactor_100hr,moistfactor_1000hr ) 
	
	return (area_burnt)

def RunSpytfire_OK (  grid_id, year, theta, Datos ):
	(pftpar, sla, tree, evergreen, summergreen, raingreen,\
		needle, boreal, lm_sapl, \
		sm_sapl, hm_sapl, rm_sapl, latosa, allom1, allom2,\
		allom3, allom4, wooddens, reinickerp) =\
                                spytfire.pftparameters()
	
	#s = SpitfireDB()
	theta0 = numpy.array([0.39, 14.,  12.,  23.,  28.,  22.,  18.,  16.,   4.,   2.,\
		66.0, 3.58, 0.98, 0.2, 0.404, 0.487, 0.525, 0.544])
	datos = deepcopy(Datos)
	a_nd = theta[0]
	#pdb.set_trace()
	pftpar[0,36] = theta[1] # Fuel Bulk Density
	pftpar[1,36] = theta[2]
	pftpar[7,36] = theta[3]
	pftpar[8,36] = theta[4]
	sigma_1hr = theta[5] ; sigma_10hr = theta[6] ; sigma_100hr=theta[7]
	moistfactor_livegrass = theta0[13] ; moistfactor_1hr = theta0[14]
	moistfactor_10hr = theta0[15]
	moistfactor_100hr = theta0[16]
	moistfactor_1000hr = theta0[17]
	
	fbd_a = 1.2 ; fbd_b = 1.4
	fbd_C3_livegrass = 4. ; fbd_C4_livegrass = 4.
	fbd_a = theta[1] ; fbd_b = theta[2] #1.2, 1.4
	fbd_C3_livegrass = theta[3] ; fbd_C4_livegrass = theta[4] # 4.0, 4.0
	#datos = s.GetSpinUpVars ( grid_id, year )
	#modis_ba = s.GetMODISBA ( grid_id, year )
	datos['a_nd'] = a_nd*numpy.ones((1,12),dtype=numpy.float32)
	#datos['lightn'] =numpy.array( datos['lightn'])*kappa
	#pdb.set_trace()
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
	datos['dbh'],datos['tau_c'],datos['cl_t'], fbd_a, fbd_b, fbd_C3_livegrass,fbd_C4_livegrass,\
	sigma_1hr,sigma_10hr,sigma_100hr,\
	moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr,\
	moistfactor_100hr,moistfactor_1000hr ) 
	
	return (area_burnt)


def RunSpytfire_selective (  grid_id, year, theta, Datos, theta0 ):
	(pftpar, sla, tree, evergreen, summergreen, raingreen,\
		needle, boreal, lm_sapl, \
		sm_sapl, hm_sapl, rm_sapl, latosa, allom1, allom2,\
		allom3, allom4, wooddens, reinickerp) =\
                                spytfire.pftparameters()
	
	#s = SpitfireDB()
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
	moistfactor_livegrass = 0.2 ; moistfactor_1hr = 0.404
	moistfactor_10hr = 0.487
	moistfactor_100hr = 0.525
	moistfactor_1000hr = 0.544
	fbd_a = 1.2 ; fbd_b = 1.4
	fbd_C3_livegrass = 4. ; fbd_C4_livegrass = 4.
	#fbd_a = theta[1] ; fbd_b = theta[2] #1.2, 1.4
	#fbd_C3_livegrass = theta[3] ; fbd_C4_livegrass = theta[4] # 4.0, 4.0
	#datos = s.GetSpinUpVars ( grid_id, year )
	#modis_ba = s.GetMODISBA ( grid_id, year )
	datos['a_nd'] = a_nd*numpy.ones((1,12),dtype=numpy.float32)
	#pdb.set_trace()
	#datos['lightn'] =numpy.array( datos['lightn'])*kappa
	#pdb.set_trace()
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
	datos['dbh'],datos['tau_c'],datos['cl_t'], fbd_a, fbd_b, fbd_C3_livegrass,fbd_C4_livegrass,\
	sigma_1hr,sigma_10hr,sigma_100hr,\
	moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr,\
	moistfactor_100hr,moistfactor_1000hr ) 
	
	return (area_burnt)