c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE FIRE
c     Biomass destruction through disturbance by fire

      subroutine fire(year,pftpar,dtemp,dtemp_min,dtemp_max,dprec,
     *  dwindsp,lightn,dphen,litter_ag,litter_bg,fuel_1hr,fuel_10hr,
     *  fuel_100hr, d_fdi,
     *  fuel_1000hr,acflux_fire,mcflux_fire,afire_frac,lm_ind,rm_ind,
     *  sm_ind,hm_ind,nind,dw1,present,tree,lat,mw1,fpc_grid, popden,
     *  a_nd_array,height,height_class,dbh,tau_c,cl_t,num_fire,
     *  annum_fire,
     *  area_burnt,an_areafires,mfdi,an_fdi,an_fseason,mcflux_trace,
     *  acflux_trace,m_fc_crown,an_fc_crown,m_i_surface,an_i_surface,
     *  dlm_lg, dlm_1hr,dlm_10hr,
     * nesterov, fire_durat, d_area_burnt, d_numfire, ros_f, ros_b,d_lb,
     * d_fuel_consumed, d_i_surface,cf, fbd_a, fbd_b, fbd_C3_livegrass,
     * fbd_C4_livegrass,sigma_1hr,sigma_10hr,sigma_100hr,
     * moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr,
     * moistfactor_100hr,moistfactor_1000hr, alpha,
     * moisture_extinction, moisture_conversion, duration_exponential )

cf2py intent(in) dw1, present, tree, lat, mw1
cf2py intent(in) popden, a_nd, height, height_class, dbh
cf2py intent(inout) lm_ind, rm_ind, sm_ind, hm_ind, nind
cf2py intent(in) tau_c, cl_t
cf2py intent(out) acflux_fire, afire_frac
cf2py intent(inout) fpc_grid
cf2py intent(in) year, pftpar
cf2py intent(in) dtemp, dtemp_min, dtemp_max, dprec, dwindsp
cf2py intent(in) lightn, dphen
cf2py intent(inout) litter_ag, litter_bg, fuel_1hr, fuel_10hr
cf2py intent(inout) fuel_100hr, fuel_1000hr
cf2py intent(out) num_fire,annum_fire
cf2py intent(out) area_burnt, mcflux_fire
cf2py intent(out) an_areafires,mfdi,an_fdi,an_fseason,mcflux_trace
cf2py intent(out) acflux_trace,m_fc_crown,an_fc_crown,m_i_surface,an_i_surface
cf2py intent(out) d_fdi, dlm_lg, dlm_1hr, dlm_10hr, nesterov
cf2py intent(out) fire_durat, d_area_burnt, d_numfire, ros_f, ros_b,d_lb
cf2py intent(out) d_fuel_consumed, d_i_surface, cf
cf2py intent(in) fbd_a, fbd_b, fbd_C3_livegrass, fbd_C4_livegrass
cf2py intent(in) sigma_1hr,sigma_10hr,sigma_100hr
cf2py intent(in) moistfactor_100hr,moistfactor_1000hr
cf2py intent(in) moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr
cf2py intent(in) alpha, moisture_extinction
cf2py intent(in) moisture_conversion, duration_exponential
      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar
         parameter (npft=9,npftpar=50,nsoilpar=7)
      real pi
         parameter (pi=3.14159265)
      real minfuel
         parameter (minfuel=100.0)  
c     SPITFIRE
      real wind_speed
*         parameter (wind_speed=76.8)    
c     total mineral content
      real MINER_TOT
         parameter(MINER_TOT=0.055)
c     calorific heat content (kJ/kg)
      real H
         parameter(H=18000.0)
c     surface-area-to-volume ratio (cm/cm��)
c     This comes from python
      real sigma_1hr,sigma_10hr,sigma_100hr
      real alpha
c         parameter(sigma_1hr=66.0,sigma_10hr=3.58,sigma_100hr=0.98)
c     ALLAN
       real sigma_1000hr
       parameter(sigma_1000hr = 0.5) 
       real sigma_livegrass
       parameter (sigma_livegrass=80.0)
c       Implement formula for Me after P&R(1986). Later, weight by dead fuel classes.

c      discussion 25/07/05: physical vs. chemical composition of dead fuel on flammability.
c      With this approach no PFT influence on fire ignition, only physical characteristics
c             of dead fuel
c      Check this approach vs. me as a PFT parameter.
       real moistfactor_livegrass,moistfactor_1hr, moistfactor_10hr
       real moistfactor_100hr,moistfactor_1000hr
       real moisture_conversion, duration_exponential
c       parameter (moistfactor_livegrass = 0.398)   
!        parameter (moistfactor_livegrass=0.2)
!        parameter (moistfactor_1hr = 0.404)         
!        parameter (moistfactor_10hr = 0.487)        
!        parameter (moistfactor_100hr = 0.525)      
!        parameter (moistfactor_1000hr = 0.544)    

c     fire duration (min)
      real fire_durat
c         parameter(fire_durat=480.0)   
c       parameter (fire_durat=240.0) 
c	LPJ limited to daily timestep.  This is assumed to be the  
c	maximum fire duration during the 
c	mid-high latitude summer or sub- /equatorial dry season.	
       real moisture_extinction ( 1:npft)
       real d_fuel_consumed(1:365)
       real fbd_a  
       real fbd_b  
c       parameter(fbd_a = 1.2, fbd_b = 1.4)

       real fbd_C3_livegrass,fbd_C4_livegrass
c       parameter (fbd_C3_livegrass=4.0,fbd_C4_livegrass=4.0) 

	real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12 
	parameter(p1=0.19,p2=0.17,p3=0.14,p4=0.04,
     *  p5=0,p6=0,p7=0,p8=0,
     *  p9=0.01,p10=0.08,p11=0.17,p12=0.21) 

c     ARGUMENTS
      real a_nd_array(1:12)
      real pftpar(1:npft,1:npftpar)
      real dtemp(1:365),dtemp_min(1:365),dtemp_max(1:365)
      real dprec(1:365),dwindsp(1:365)
      real dphen(1:365,1:npft)
      real litter_ag(1:npft),litter_bg(1:npft)
      real fuel_1hr(1:npft),fuel_10hr(1:npft),fuel_100hr(1:npft)
      real fuel_1000hr(1:npft)
      real acflux_fire,acflux_trace(1:6),mcflux_fire(1:12)
      real mcflux_fire_pft(1:12,1:npft),mcflux_trace(1:12,1:6)
      real mcflux_trace_pft(1:12,1:6,1:npft)
      real lm_ind(1:npft),rm_ind(1:npft)
      real sm_ind(1:npft),hm_ind(1:npft)
      real nind(1:npft)
      real dw1(1:365),mw1(1:12),dlm(1:365) 
      real afire_frac,fire_frac
      real fire_length,an_fseason
      real height(1:npft),dbh(1:npft)       
      real tau_c(0:4,1:npft)                
      real cl_t(0:4,1:npft)
      logical present(1:npft),tree(1:npft)
      integer year

c     SPITFIRE 
      real mtemp(1:12),mtemp_dmin(1:12),mtemp_dmax(1:12)
      real d_numfire(1:365),d_numf_old(1:365),num_fire(1:12),annum_fire
      real d_area_burnt(1:365),area_burnt(1:12),an_areafires
      real d_fdi(1:365),mfdi(1:12),an_fdi
      real fpc_grid(1:npft)
      real lat
      real popden,popden_yr(1:102),a_nd

c     LOCAL VARIABLES
      integer pft,d,m,n,i,month(1:12),month_length(1:12),x,count
      integer count_int,count_fdi,count_yr,n_pft
      real fuel,fire_prob
      real fire_index,disturb
      real fuel_1hr_total,fuel_10hr_total,fuel_100hr_total
      real fuel_1000hr_total
      real ef_trace(1:npft,1:6),acflux_fire_pft(1:npft)
      real dcflux_fire_pft(1:365,1:npft)
      real resist(1:npft)
      real moistfactor,me,litter_ag_total
      real fire_term,area, area_ha
c     SPITFIRE 
      real m_firelength(1:12),length(0:12)
c      real dtemp_min,dtemp_max
      real U_back,U_front,gamma

      real ros_f(1:365),ros_b(1:365)
      real df(1:365),db(1:365)
      real lightn(1:12)                      
      real lightn1  		
      real FDI                         
      real human_ign                   
      real popden1
      real fpc_tree_total,fpc_grass_total,lb_grass,lb
	 real d_lb(1:365)
      real dens_fuel(1:npft),dens_fuel_tot,dens_fuel_ave,sigma
      real net_fuel,dead_fuel,livegrass
      real wind_forward,wind_backward,back_ws
      real d_i_surface(1:365),m_i_surface(1:12),an_i_surface  
      real fuel_consum,an_fuel_consum
*      real pot_fc_1hr(1:npft),pot_fc_10hr(1:npft),pot_fc_100hr(1:npft)
      real pot_fc_lg(1:npft), ratio_lg_reduc
      real fc_1hr_total,fc_10hr_total,fc_100hr_total,fc_lg_total
      real fc_1hr(1:npft),fc_10hr(1:npft),fc_100hr(1:npft)
      real fc_lg(1:npft),fc_1000hr(1:npft)
      real sh(1:npft),f(1:npft),g(1:npft) 
      real crown(1:npft),ck(1:npft),cl(1:npft)
      real ck_t(1:5,1:npft)  
      real fc_crown(1:365),m_fc_crown(1:12),an_fc_crown
      real wc_live(1:npft),postf_mort(1:npft)
      real cf(1:npft),tau_l(1:npft)
      real pm_tau(1:npft),pm_tau_class(0:4,1:npft)
      real pm_ck(1:npft),r_ck(1:npft),p(1:npft)
      real nind_fuel(1:npft),nind_old(1:npft),nind_kill(1:npft)






       real ratio_fbd   
       real ratio_dead_fuel
       real ratio_live_fuel
       real dens_deadgrass_ave
       real dens_livegrass_ave
       real char_dens_fuel_ave
       real char_net_fuel
       real deadgrass
       real char_sigma
       real char_moistfactor
       real ratio_C3_livegrass
       real ratio_C4_livegrass
c       real moistfactor_C3_livegrass
c       real moistfactor_C4_livegrass
c       parameter (moistfactor_C3_livegrass = 0.1, 
c     *       moistfactor_C4_livegrass =0.1)

c       real moistfactor_1hr, moistfactor_10hr
c       real moistfactor_100hr, moistfactor_1000hr
       real   dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365) 
       real nesterov(1:365) 
       real   dlm_100hr (1:365), dlm_1000hr (1:365)
c       real fbd_C3_livegrass
       real nind_fa(1:npft) 
       real nind_nkill(1:npft) 
       real ag _non_leaf_ind(1:npft) 
       real fuel_update_1hr(1:npft),fuel_update_10hr(1:npft)
       real fuel_update_100hr(1:npft),fuel_update_1000hr(1:npft)  
       real prop_fa_nk(1:npft) 
       real afire_frac_temp,an_areafires_temp
       real class, height_class(0:4,1:npft)
       real pot_fc_lg_temp(1:npft)
       real ignitions 

c     Initialise
	  month(1)=31
	  month(2)=59
	  month(3)=90
	  month(4)=120
	  month(5)=151
       month(6)=181
       month(7)=212
	  month(8)=243
	  month(9)=273
	  month(10)=304
	  month(11)=334
	  month(12)=365



      fire_prob=0.0
      fire_index=0.0
      fire_length=0.0
      fuel=0.0
      disturb=0.0
      acflux_fire=0.0
      annum_fire=0.0
      an_areafires=0.0
      an_fdi=0.0
      an_fseason=0.0
      sigma=0.0
      lb_grass=0.0
      lb=0.0
      wind_speed=0.0  

      fpc_tree_total=0.0
      fpc_grass_total=0.0
      wind_forward=0.0
*      p=0.0
      back_ws=0.0
      area_ha=0.0
      net_fuel=0.0
      dead_fuel=0.0
      dens_fuel_tot=0.0
      dens_fuel_ave=0.0
      livegrass=0.0
      fire_frac=0.0
c      sh=0.0
      an_i_surface=0.0
      an_fc_crown=0.0
*      r_ck=0.0
*      p=0.0

      do m=1,12
        length(m)=0.0
        m_firelength(m)=0.0
        num_fire(m)=0.0
        area_burnt(m)=0.0
        mfdi(m)=0.0
        m_fc_crown(m)=0.0
        m_i_surface(m)=0.0
        mcflux_fire(m)=0.0
        do pft=1,npft
          mcflux_fire_pft(m,pft)=0.0
        enddo
        do x=1,6
          mcflux_trace(m,x)=0.0
          do pft=1,npft
             mcflux_trace_pft(m,x,pft)=0.0
          enddo
        enddo
      enddo

      do i=1,365
        d_fdi(i)=0.0
        d_numfire(i)=0.0
        dlm(i)=0.0
        d_numf_old(i)=0.0
        ros_f(i)=0.0
        ros_b(i)=0.0
        d_area_burnt(i)=0.0
        d_i_surface(i)=0.0
        df(i)=0.0
        db(i)=0.0
        fc_crown(i)=0.0
      enddo

      do x=1,6
       acflux_trace(x)=0.0
      enddo

      do pft=1,npft
        do d=1,365
          dcflux_fire_pft(d,pft)=0.0
        enddo
       acflux_fire_pft(pft)=0.0
       wc_live(pft)=0.0
       pot_fc_lg(pft)=0.0
       fc_1hr(pft)=0.0
       fc_10hr(pft)=0.0
       fc_100hr(pft)=0.0
       fc_1000hr(pft)=0.0
       sh(pft)=0.0
       nind_old(pft)=0.0
       nind_kill(pft)=0.0
       postf_mort(pft)=0.0
       pot_fc_lg(pft)=0.0
       fuel_update_1hr(pft)=0.0
       fuel_update_10hr(pft)=0.0
       fuel_update_100hr(pft)=0.0
       fuel_update_1000hr(pft)=0.0
      enddo

c     Assign a minimum fire fraction (for presentational purposes)

      afire_frac=0.001

c    ASSIGN PFT PARAMETER REQUIRED IN THE FIRE ROUTINE

c    ASSIGN PFT PARAMETER FOR Glob-FIRM AND SPITFIRE

c          Calculate total above-ground litter

      litter_ag_total=0.0
      do pft=1,npft
        litter_ag_total=litter_ag_total+litter_ag(pft)
      enddo

c        Calculate litter moisture weighting factor (moisture of extinction me)

      moistfactor=0.0

      do pft=1,npft
!        me=pftpar(pft,6)
         me = moisture_extinction ( pft )
        if (litter_ag_total.gt.0.0) then
          moistfactor=moistfactor+(litter_ag(pft)/litter_ag_total)*me
        else
          moistfactor=0.0
        endif
      enddo


c        Assign emission factors for trace gas emissions resulting from acflux_fire_pft
      do pft=1,npft
         ef_trace(pft,1)=pftpar(pft,38)/1000.0   
         ef_trace(pft,2)=pftpar(pft,39)/1000.0   
         ef_trace(pft,3)=pftpar(pft,40)/1000.0   
         ef_trace(pft,4)=pftpar(pft,41)/1000.0   
         ef_trace(pft,5)=pftpar(pft,42)/1000.0   
         ef_trace(pft,6)=pftpar(pft,43)/1000.0   
       enddo
c     ASSIGN PFT PARAMETER FOR Glob-FIRM ONLY 
c     Assign PFT resistance to fire

      do pft=1,npft
        resist(pft)=pftpar(pft,8)
      enddo

c     START OF REGIONAL FIRE MODEL SPITFIRE, when year of observed fire reached
c     otherwise Glob-FIRM
      if (year.lt.1001) then

c     Calculate the length of the fire season (units=days)
      i=1
      do d=1,365

c       Calculate today's fire probability, fire_prob
c       Assume fire is only possible when temperature is above zero

        if (dtemp(d).gt.0.0.and.moistfactor.gt.0.0) then
          fire_prob=EXP((-pi)*(dw1(d)/moistfactor)**2)
        else
          fire_prob=0.0
        endif

        fire_length=fire_length+fire_prob

      enddo
      an_fseason=fire_length

c     Calculate annual fire index

      fire_index=fire_length/365.0

c     Calculate the available fuel (above-ground litter) to carry the fire

      do pft=1,npft
        fuel=fuel+litter_ag(pft)
      enddo


      fire_term=fire_index-1.0

      afire_frac=fire_index*EXP(fire_term/(0.45*fire_term**3
     *    + 2.83*fire_term**2 + 2.96*fire_term + 1.04))

      call pixelarea(lat,area)
        an_areafires=(afire_frac*area)/10000.0

c     Reduce fraction of grid cell affected by fire when fuel
c     becomes limiting (reduced carrying capacity)

        if (fuel.lt.minfuel) then
           an_areafires=0.0
           afire_frac=0.0
         endif

      if (afire_frac.lt.0.001) afire_frac=0.001

c     Implement the effect of the fire on vegetation structure and litter
c     in the disturbed fraction.

c     Each PFT is assigned a resistance to fire, representing the fraction of
c     the PFT which survives a fire. Grasses assumed already to have completed
c     their life cycle and thus are not affected by fire, giving them
c     a competitive advantage against woody PFTs.

c     Calculate trace gas emission resulting from biomass burning

      do pft=1,npft

        if (present(pft).and.tree(pft)) then

c         Calculate the fraction of individuals in grid cell which die

          disturb=(1.0-resist(pft))*afire_frac

c         Calculate carbon flux to atmosphere (gC/m2) due to burnt biomass

          acflux_fire=acflux_fire+disturb*(nind(pft)*
     *      (lm_ind(pft)+sm_ind(pft)+hm_ind(pft)+rm_ind(pft)))

          acflux_fire_pft(pft)=disturb*(nind(pft)*
     *       (lm_ind(pft)+sm_ind(pft)+hm_ind(pft)+rm_ind(pft)))

c         Update the individual density

          nind(pft)=nind(pft)*(1.0-disturb)


        endif

c       Add combusted litter to carbon flux to atmosphere term

        acflux_fire=acflux_fire+(afire_frac*litter_ag(pft))


        acflux_fire_pft(pft)=acflux_fire_pft(pft)+
     *   (afire_frac*litter_ag(pft))

c         Calculate trace gas emissions (gSpecies/m��)

        do x=1,6
           if (acflux_fire.gt.0.0) then
               acflux_trace(x)=acflux_trace(x)+
     *                (acflux_fire_pft(pft)*ef_trace(pft,x)/0.45)
           else
               acflux_trace(x)=0.0
           endif

*        print*,'acfluxfire',year,pft,x,ef_trace(pft,x),
*     *    ef_trace(pft,x)/0.45,acflux_fire_pft(pft),acflux_trace(x)

        enddo

c       Update the above ground litter term

        litter_ag(pft)=(1.0-afire_frac)*litter_ag(pft)

c       fuel classes
          fuel_1hr(pft)=(1.0-afire_frac)*fuel_1hr(pft)

          fuel_10hr(pft)=(1.0-afire_frac)*fuel_10hr(pft)
          fuel_100hr(pft)=(1.0-afire_frac)*fuel_100hr(pft)
          fuel_1000hr(pft)=(1.0-afire_frac)*fuel_1000hr(pft)


      enddo

      else 


c  START OF SPITFIRE
c    ASSIGN PFT PARAMETER VALUES FOR SPITFIRE ONLY
     
c    fuel bulk density
*      dens_fuel_ave=0.0
      do pft=1,npft
        if (present(pft)) dens_fuel(pft)=pftpar(pft,37)
      enddo

c    for calculation of fire perimeter calculate tree and grass coverage
        do pft=1,npft
          if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
            else
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)

            endif
           endif
        enddo


c      f=0.0
      do pft=1,npft

c     crown length as a proportion of tree height
c         crown(pft)=pftpar(pft,44)

c     scorch height parameter for crown fire
c         f=f+(fpc_grid(pft)/fpc_tree_total)*pftpar(pft,45)
          f(pft)=pftpar(pft,45)
          g(pft)=pftpar(pft,46)
c     r(ck) and p parameter for postfire mortality as a result of crown damage
         r_ck(pft)=pftpar(pft,49)
         p(pft)=pftpar(pft,50)

      enddo

c     population density
c       popden1=popden
c       print*,year,popden 


       call pixelarea(lat,area)
       area_ha=(area/10000.0)

c       do pft=1,npft
c       if (present(pft)) then
c       print*,'fuel 1hr',pft,fuel_1hr(pft)
c       print*,'fuel 10hr',pft,fuel_10hr(pft)
c       print*,'fuel 100hr',pft,fuel_100hr(pft)
c       print*,'fuel 1000hr',pft,fuel_1000hr(pft)
c       endif
c       enddo
c---------------------------------------
c  SPITFIRE at daily time step
c---------------------------------------
        m=1
        count=1
        count_int=1
        count_yr=1
        count_fdi=1
        afire_frac=0.0

       do d=1,365

c	print*,'-------- '
c        print*,'day= ',d
c         if (year.ge.1050) pause

c         d_numf_old(d)=0.0
         d_i_surface(d)=0.0
         fc_crown(d)=0.0
         wind_speed=0.0
         fire_frac=0.0
         fuel_consum=0.0
        livegrass=0.0
c        deadgrass=0.0  
	ratio_lg_reduc=0.0
        ratio_dead_fuel=0.0
        ratio_live_fuel=0.0

	do pft=1,npft
	fc_lg(pft)=0.0
	fc_1hr(pft)=0.0
	fc_10hr(pft)=0.0
	fc_100hr(pft)=0.0
	fc_1000hr(pft)=0.0
        nind_kill(pft)=0.0
        nind_fa(pft)=0.0
        nind_nkill(pft)=0.0
        pot_fc_lg(pft)=0.0
c        do class=1,5
          ck(pft)=0.0
c        enddo
        pm_ck(pft)=0.0
        pm_tau(pft)=0.0
        postf_mort(pft)=0.0
        sh(pft)=0.0
	enddo

c influence of wind_speed on rate of forward spread
c wind speed from NCEP reanalysis data is in m/s, ROS is in m/min

c	if (year.ge. 1050) print*,'day= ',d,'wspeed= ',dwindsp(d) 

        wind_speed=dwindsp(d)*60  
c  introduce reduction factor for wind speed with respect to forest or grass fpc 
c       ref: Rothermel 1983, Pyne 1996

        wind_speed=(fpc_tree_total*dwindsp(d)*60.0*0.4)+
     *               (fpc_grass_total*dwindsp(d)*60.0*0.6)
c        wind_speed=(fpc_tree_total*1.28*60.0*0.4)+
c     *               (fpc_grass_total*1.28*60.0*0.6)
c        print*,'wind_speed',fpc_tree_total,fpc_grass_total,wind_speed 
c        pause
c	converts wind_speed (m/min) to ft/min
c 	for input into Rothermel's formula for phi_wind in the ROS S/R
        wind_forward=3.281*wind_speed       

c        back_ws=wind_speed*exp(-0.05039*wind_speed)     
c        wind_backward=3.281*back_ws

        
    
c       weight lb according to presence of grasses and forests
c       wind speed in Canadian formula is in km/h

c	CHECK UNITS
c 	check fpc- proportion?

       If (wind_speed.lt.16.67) then
                lb = 1
       else 
              lb=fpc_tree_total*
     *   (1.0+(8.729*((1.0-(exp(-0.03*0.06*wind_speed)))**2.155)))
     *    +(fpc_grass_total*(1.1+((0.06*wind_speed)**0.0464)))
       endif

       if (lb.gt.8)  lb = 8
	  d_lb(d) = lb
c Kirsten: check alternatively:
c         else 
c              lb=min(8,fpc_tree_total*
c     *   (1.0+(8.729*((1.0-(exp(-0.03*0.06*wind_speed)))**2.155)))
c     *    +(fpc_grass_total*(1.1+((0.06*wind_speed)**0.0464))))
c       endif
  

c    FUEL CHARACTERISTICS

c    net fuel load and total amount of dead fuel per fuel class
c    and amount of livegrass
c    ACHTUNG: take grass phenology into account (cf. grass curing) per time step
c             to calculate amount of live grass
        dead_fuel=0.0
        fuel_1hr_total=0.0
        fuel_10hr_total=0.0
        fuel_100hr_total=0.0
        dens_fuel_ave=0.0

        do pft=1,npft
          if (present(pft)) then

            fuel_1hr_total=fuel_1hr_total+fuel_1hr(pft)/0.45

            if (tree(pft)) then
              fuel_10hr_total=fuel_10hr_total+fuel_10hr(pft)/0.45
              fuel_100hr_total=fuel_100hr_total+fuel_100hr(pft)/0.45
            else 
c   KIRSTEN:  take proportion of grass leafmass, when green grass leaves are on
c         todays amount of green grass leaves: [gC/m2],influence on ROS only through moist_lg_1hr
              livegrass=livegrass+
     *                        (lm_ind(pft)/0.45)*nind(pft)*dphen(d,pft)  

c              print*,'livegrass',d,pft,livegrass,dphen(d,pft)
c              used in fire effects section only
               pot_fc_lg(pft)=(lm_ind(pft)/0.45*nind(pft))*dphen(d,pft)  

c  Allan
c              livegrass=livegrass+
c     *                    (lm_ind(pft)/0.45)*nind(pft)  
c               pot_fc_lg(pft)=(lm_ind(pft)/0.45)*nind(pft)   

c               deadgrass = deadgrass + (fuel_1hr(pft)/0.45) 	
            endif

          endif
        enddo


        dead_fuel = fuel_1hr_total + fuel_10hr_total 
     *     + fuel_100hr_total  

c       net fuel load
        net_fuel=0.0
        if (dead_fuel.gt.0.0)
     *    net_fuel=(1.0-MINER_TOT)*(dead_fuel/1000.0)  

c    fuel bulk density, weighted per fuel class and fuel load
c    ACHTUNG: WEIGHTING per fpc AND fuel load or per fuel load only? Reg-FIRM: per FPC       
        do pft=1,npft
          if (present(pft)) then
            if (dead_fuel.gt.0.0) then
              ratio_fbd = ((fuel_1hr(pft) + fbd_a * fuel_10hr(pft) + 
     *                     fbd_b * fuel_100hr(pft)) / 0.45) / dead_fuel 
              dens_fuel_ave = dens_fuel_ave + dens_fuel(pft) * ratio_fbd    
            else
              dens_fuel_ave=0.0
            endif
          endif
        enddo


c    livegrass

       if (livegrass.gt.0.0) then
c    KIRSTEN: calculate leaf moisture content of grasses from dw1
c    ACHTUNG: change with water scalar value for grasses
c      ORIGINAL representation
c       dlm_lg(d) = max(0.0,((10.0/9.0)*dw1(d)-(1.0/9.0)))
       dlm_lg(d) = max(0.0,(( moisture_conversion )*dw1(d)-
     * (moisture_conversion-1)))
       ratio_C3_livegrass = pot_fc_lg(8) / livegrass
       ratio_C4_livegrass = pot_fc_lg(9) / livegrass
       else
         dlm_lg(d)=0.0
       endif


c      influence of livegrass on FBD
c    Kirsten: presence of livegrass additionally lowers FBD, when trees mixed 
c    with grasses
       dens_livegrass_ave =  
     *       fbd_C3_livegrass *  ratio_C3_livegrass +
     *      fbd_C4_livegrass *  ratio_C4_livegrass


        if (dead_fuel.gt.0.0 .or. livegrass.gt.0.0) then 
         ratio_dead_fuel = dead_fuel  / (dead_fuel + livegrass)
         ratio_live_fuel = livegrass / (dead_fuel + livegrass)
         char_dens_fuel_ave = dens_fuel_ave* ratio_dead_fuel
     *    + dens_livegrass_ave * ratio_live_fuel	
        endif

c         print*,'chardensfuel',d,dead_fuel,char_dens_fuel_ave,
c     *    dens_fuel_ave,livegrass,ratio_live_fuel,ratio_dead_fuel

         char_net_fuel = net_fuel +(1.0-MINER_TOT)*livegrass/1000.0 

       

c      moistfactor
c     Kirsten: original method, where moistfactor is a PFT parameter 
c       	moistfactor = (moistfactor_1hr * fuel_1hr_total
c     *          + moistfactor_10hr * fuel_10hr_total
c     *          + moistfactor_100hr * fuel_100hr_total)/dead_fuel

c Kirsten: why should livegrass increase flammability? 
c  moistfactor for livegrass less than for dead fuel
        char_moistfactor = moistfactor *  ratio_dead_fuel
     *       + moistfactor_livegrass * ratio_live_fuel


c    surface-area-to-volume ratio weighted among fuel classes
c     In future, account for SAV diffs between pine needles and
c     broad-leaves.  Previous studies report order
c     of magnitude difference in SAV b/w the two leaf types.
	
       if (dead_fuel.gt.0.0) then
         sigma=(fuel_1hr_total * sigma_1hr +
     *         fuel_10hr_total * sigma_10hr +
     *         fuel_100hr_total * sigma_100hr) / dead_fuel	

c    higher sigma value for livegrass increases fire 
c       char_sigma =sigma*ratio_dead_fuel+sigma_livegrass*ratio_live_fuel

       else
         sigma=0.00001
c          char_sigma=0.00001
       endif
c        print*,'charsigma',char_sigma,sigma
c         if (d.eq.365) then
c        print*,'year=',year,'char_moist= ',char_moistfactor,moistfactor
c        pause
c         endif

c-----------------------------------------------------------
c      end preparing daily variables - start of simulation
c-----------------------------------------------------------

c       climatic fire danger
c	OLD: Assume moisture content of livegrass close to 1.0 ? NO
         call fire_danger_index
     *    (d_fdi,dlm,dlm_lg,dtemp_min,dtemp_max,dprec,
     *    d,moistfactor,fuel_1hr_total,fuel_10hr_total, nesterov,
     *    fuel_100hr_total,
     *    dead_fuel,char_moistfactor, ratio_dead_fuel,ratio_live_fuel,
     *    dlm_1hr,dlm_10hr, dlm_100hr, dlm_1000hr,year,
     *    alpha, sigma_1hr, sigma_10hr, sigma_100hr ) ! changed to account for alpha calibration

         if (d_fdi(d).gt.0.0) then
           mfdi(m)= mfdi(m)+d_fdi(d)
           an_fdi=an_fdi+d_fdi(d)
           count=count+1
         endif

c  number of fires. only ignitions, when enough dead fuel load
c        old: if(fuel.gt.minfuel)
         if (net_fuel.gt.0.001) then   
c	CHECK
c	CHECK

c	RUN4 Allan


c	lightn1=lightn*365.0

	SELECT CASE(d)
  	CASE(1:31)
     		lightn1 = lightn(1)
			a_nd = a_nd_array(1)
  	CASE(32:59)
     		lightn1 = lightn(2)
			a_nd = a_nd_array(2)

	CASE(60:90)
      		lightn1 = lightn(3)
			a_nd = a_nd_array(3)

  	CASE(91:120)
     		lightn1 = lightn(4)
			a_nd = a_nd_array(4)

  	CASE(121:151)
     		lightn1 = lightn(5)
			a_nd = a_nd_array(5)

  	CASE(152:181)
     		lightn1 = lightn(6)
			a_nd = a_nd_array(6)

  	CASE(182:212)
      		lightn1 = lightn(7)
			a_nd = a_nd_array(7)

  	CASE(213:243)
     		lightn1 = lightn(8)
			a_nd = a_nd_array(8)

  	CASE(244:273)
     		lightn1 = lightn(9)
			a_nd = a_nd_array(9)

  	CASE(274:304)
      		lightn1 = lightn(10)
			a_nd = a_nd_array(10)

  	CASE(305:334)
     		lightn1 = lightn(11)
			a_nd = a_nd_array(11)

  	CASE(335:365)
     		lightn1 = lightn(12)
			a_nd = a_nd_array(12)

  	END SELECT

c	end RUN4 changes

c399	continue

c        ignitions = human_ign(popden,a_nd,year)+lightn 



	ignitions = human_ign(popden,a_nd,year)+lightn1
c	print*,d,human_ign(popden,a_nd,year),lightn
c        ignitions = human_ign(popden,a_nd,year) RUN3 Allan
c	if (ignitions.lt.1) then
c	 	ignitions = 1
c	endif
c        if (year.eq.1002) then
c        print*,d,a_nd,human_ign(popden,a_nd,year),lightn,ignitions
c        endif
           d_numfire(d)=d_fdi(d)*
     *       ignitions*0.000001*area_ha 

c            continued burning of yesterdays ignitions, not accounting in the # fires statistics
*             if (d_i_surface(d-1).gt.1000.0.and.dprec(d).le.3.0)
*     *          d_numf_old(d)=d_numfire(d-1)
         else  
           d_numfire(d)=0.0
         endif

*	 if (year.eq.1070) then
c      print*,'d_fdi(d)',d,lightn,d_numfire(d)
c      pause
*        endif


c    area burnt

c rate of spread
       if (net_fuel.gt.0.0) then

c          ros_f(d)=U_front(wind_forward,dens_fuel_ave,sigma,dlm,d,
c     *             net_fuel,moistfactor,H, char_dens_fuel_ave, 
c     *             char_sigma, char_net_fuel, char_moistfactor)


c     calculation of backward ROS taken out, because of inconsistencies
c     Allan: calculated after Hirsch (1996)
c         call rate_of_spread(U_front,wind_backward,dens_fuel_ave,sigma,
c     *       dlm,d,net_fuel,moistfactor,H, char_dens_fuel_ave, 
c     *             char_sigma, char_net_fuel, char_moistfactor,gamma)
        
c          ros_b(d)=U_front

       call rate_of_spread(U_front,wind_forward,dens_fuel_ave,sigma,dlm,
     *            d, net_fuel,moistfactor,H, char_dens_fuel_ave, 
     *            char_sigma, char_net_fuel, char_moistfactor,gamma)

       ros_f(d)=U_front


 

c	CHECK FORMULAE: HIRSCH (1996)
c     this doubles ros_b(d) compare to the old method to estimate backward wind speed
       ros_b(d) = ros_f(d) * exp(-0.012 * wind_speed) 
           if (ros_b(d).lt. 0.05) ros_b(d) = 0.0 
     

c	check the parameter value
c     fire duration as a function of d_fdi
c          fire_durat=301.0/(1.0+(((301.0/1.)-1.)*exp(-11.4*d_fdi(d))))
c           fire_durat=361.0/(1.0+(((361.0/1.)-1.)*exp(-11.76*d_fdi(d))))
c        ORIGINAL
c          fire_durat=241.0/(1.0+(((241.0/1.)-1.)*exp(-10.95*d_fdi(d))))
          fire_durat = 241./( 1.0 +
     *    (240.* exp(-duration_exponential*d_fdi(d))))
          db(d)=ros_b(d)*fire_durat 
          df(d)=ros_f(d)*fire_durat
       else
          ros_b(d)=0.0
          db(d)=0.0
          ros_f(d)=0.0
          df(d)=0.0
       endif

c  area burnt from todays ignitions

         if(lb.gt.0.0) then
         d_area_burnt(d)=d_numfire(d)*
     *        ((pi/(4.0*lb))*((df(d)+db(d))**2.0))/10000.0 

c       print*,'ROS',d,ros_b(d),ros_f(d),db(d),df(d),
c     *    d_area_burnt(d)
c     *  'ros_f(d)= ',ros_f(d),
c     *  'df(d)= ', df(d),
c     *  'ros_b(d)= ',ros_b(d),
c     *  'db(d)= ', db(d)
c       pause

        fire_frac=(d_area_burnt(d)*10000.0)/area

c	ALLAN: Do not update annual burnt area stats yet. Need to assess FI threshold first.
c	Put info into temp vars for now. KIRSTEN: May not be necessary
c         afire_frac=afire_frac+fire_frac
          afire_frac_temp=afire_frac+fire_frac
          an_areafires_temp=an_areafires+d_area_burnt(d)

c--------------------------------------------------------
c    prevent to burn the entire cell AND ensure fire_frac is consistent with its later use
         if ((an_areafires_temp.gt.area_ha).or.
     *      (d_area_burnt(d).gt.area_ha)) then
               d_area_burnt(d)=area_ha-an_areafires_temp+d_area_burnt(d)
               fire_frac=1.0-afire_frac_temp+fire_frac
c               an_areafires_temp=area_ha
c               afire_frac_temp=1.0
         endif
c-------------------------------------------------------
c	if (year.ge.1069.and.year.lt.1071) then
c        print*,'prevent burnt > grid cell',year,d,
c     *   ' d_area_burnt(d)= ', d_area_burnt(d),
c     *   ' fire_frac= ',fire_frac
c     *   ' afire_frac= ', afire_frac
c             check references for numbers
c     *   'area_ha= ',area_ha 
c        endif

         endif 


c dead fuel consumption: in g biomass per m�� - not gC/m�� 

      call fuel_consumption(npft,present,fuel_consum,fire_frac,
     *  fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,pot_fc_lg,
     *  tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,fc_1hr,fc_lg,
     *  fc_10hr,fc_100hr,fc_1000hr,cf,char_moistfactor,
     *   dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *   moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *   moistfactor_100hr, moistfactor_1000hr)

      d_fuel_consumed(d) = fuel_consum
c    surface fire intensity
        if ((fire_frac.gt.0.0).and.(fuel_consum.gt.0.0)) then

         d_i_surface(d)=H*(fuel_consum/1000.0/fire_frac)*(ros_f(d)/60.0)
c        print*,'i_surface',fuel_consum,d_i_surface(d),ros_f(d)

c       scorch height
c          sh=f*(d_i_surface(d)**0.667)

c    update of no. fires, area burnt, fire_frac, sh
c    ACHTUNG: there might be two influences of surface fire intensity 
c             a) as a stopping rule and b) as criteria for multiple-day burning
c             check references for numbers
c    Fires may stop but continue to smoulder; or stop and extinguish.
c    Pyne et al (1986): 50 and 300 respectively. CHECK

        endif 

         if (d_i_surface(d).lt.50.0) then
           d_numfire(d)=0.0
           d_area_burnt(d)=0.0
           fire_frac=0.0
           d_i_surface(d)=0.0
c          sh=0.0
         else
           count_int=count_int+1
         endif

          num_fire(m)=num_fire(m)+d_numfire(d)
          annum_fire=annum_fire+d_numfire(d)

          an_areafires=an_areafires+d_area_burnt(d)
          afire_frac=afire_frac+fire_frac

          m_i_surface(m)=m_i_surface(m)+d_i_surface(d)
          an_i_surface=an_i_surface+d_i_surface(d)

c     prevent to burn the entire cell. DOUBLE CHECK NECESSARY??????
c         if (an_areafires.gt.area_ha.or.d_area_burnt(d).gt.area_ha) then
c               d_area_burnt(d)=area_ha-an_areafires+d_area_burnt(d)
c               fire_frac=1.0-afire_frac+fire_frac
c               an_areafires=area_ha
c               afire_frac=1.0
c         endif

          area_burnt(m)=area_burnt(m)+d_area_burnt(d)

c	goto 6601 

c   FIRE EFFECTS - THIS IS NEW cf REG-FIRM: CROWN FIRES
      do pft=1,npft
        
       if (present(pft)) then

c    Fire effects for trees, when surface fire intensity is greater than 50 W/m
c    At the moment it has the same threshold as for extinguishing fires

       if (d_i_surface(d).ge.50.0) then


c   surface fire: update fuel load per dead fuel class and litter_ag in general
         fuel_1hr(pft)=fuel_1hr(pft)-(fc_1hr(pft)*0.45) 
         fuel_10hr(pft)=fuel_10hr(pft)-(fc_10hr(pft)*0.45)
          fuel_100hr(pft)=fuel_100hr(pft)-(fc_100hr(pft)*0.45)
         fuel_1000hr(pft)=fuel_1000hr(pft)-fc_1000hr(pft)

         litter_ag(pft)=litter_ag(pft)-(fc_1hr(pft)*0.45)-
     *         (fc_10hr(pft)*0.45)-(fc_100hr(pft)*0.45)-fc_1000hr(pft)

c   start of crown fires
       if (tree(pft)) then

c   scorch height per PFT
          sh(pft)=f(pft)*(d_i_surface(d)**0.667)
c      print*,'sh',pft,g(pft),sh(pft)

c      post-fire mortality from cambial damage
c      tau_r=cf/Gamma, tau_l=tau_r, if not tau_l=2*tau_r
          tau_l(pft)=2.0*(cf(pft)/gamma)
c      print*,'tau_l',pft,tau_l(pft)
c      pause
c   uniform distribution of height (5 classes), not needed with real age classes
       do class=0,4

c   crown kill in [%] assuming the crown shape being a cylinder
c   crown height as a fraction of tree height definded per PFT
c   propn of canopy burnt = (SH - (height - cl))/cl = (SH - height + cl)/cl 


        if (sh(pft).lt.(height_class(class,pft)-cl_t(class,pft)))
     *            ck_t(class,pft)=0.0

        if (sh(pft).ge.(height_class(class,pft)-cl_t(class,pft))
     *       .and.sh(pft).lt.height_class(class,pft))
     *   ck_t(class,pft)=
     *   ((sh(pft)-height_class(class,pft)+cl_t(class,pft))
     *          /cl_t(class,pft))

          if (sh(pft).ge.height_class(class,pft)) ck_t(class,pft)=1.0

          ck(pft)=ck(pft)+ck_t(class,pft)
          
c       post-fire mortality from crown scorching
          pm_ck(pft)=pm_ck(pft)+(r_ck(pft)*(ck_t(class,pft)**p(pft)))

c      print*,'height',pft,class,height_class(class,pft),ck_t(class,pft)

c       post-fire mortality from cambial damage
c         Allan's version after Peterson&Ryan

          if ((tau_l(pft)/tau_c(class,pft)).ge.2.0) then
               pm_tau_class(class,pft)=1.0
          else
             if ((tau_l(pft)/tau_c(class,pft)).gt.0.22) then
           pm_tau_class(class,pft)=
     *              (0.563*(tau_l(pft)/tau_c(class,pft)))-0.125
             else
               pm_tau_class(class,pft)=0.0
             endif
          endif

          pm_tau(pft)=pm_tau(pft)+pm_tau_class(class,pft)

c        print*,'day ',d,pft,class,height_class(class,pft),
c     *    ' pm_tau(pft)',pm_tau_class(class,pft),pm_ck(pft)
c     *    ,ck(pft),pm_tau(pft),pm_ck(pft)
c        pause

        enddo 
        
          ck(pft)=ck(pft)/5.0
          pm_ck(pft)=pm_ck(pft)/5.0
          pm_tau(pft)=pm_tau(pft)/5.0
           
c       Calculate total post-fire mortality from crown scorching AND cambial kill

       postf_mort(pft)=pm_tau(pft)+pm_ck(pft)-(pm_tau(pft)*pm_ck(pft))
c       print*,'postf mort.',d,pft,d_i_surface(d)
c     *    ,sh(pft),postf_mort(pft),ck(pft),pm_tau(pft),
c     *  pm_ck(pft)
c       pause
c	number of indivs affected by fire in grid cell 
        nind_fa(pft)=fire_frac*nind(pft)
c        print*,'nind pot. affected by fire',pft,nind(pft),nind_fa(pft)
c     ACHTUNG: moisture content of live fuel can be calculated depending from
c     canopy conductance etc.
*              wc_live(pft)=1.0       

c          amount of LIVE FUEL combusted in a crown fire
c	ALLAN: Assume 100% of leaves, 100% of small twigs, 100% of large twigs
c	and 5% of small branches combusted (Stocks et al 2004 etc).
        disturb=0.0

        ag_non_leaf_ind(pft) = sm_ind(pft) + hm_ind(pft)

c	Kirsten: the fractions of hm_ind and sm_ind to fuel classes should always 
c       add up to 1.0.
c       no 1000hr fuel involved in crown kill
       disturb = nind_fa(pft) * ck(pft) * (lm_ind(pft) +
     *    0.045 * ag_non_leaf_ind(pft) +
     *    0.075 * ag_non_leaf_ind(pft) +
     *    0.21 * 0.05 * ag_non_leaf_ind(pft)) 

       dcflux_fire_pft(d,pft)=dcflux_fire_pft(d,pft)+disturb
       fc_crown(d)=fc_crown(d)+disturb


c	Calculate total number of indivs killed for fire affected area

       nind_kill(pft) = postf_mort(pft) * nind_fa(pft)
c       print*,'nind_kill',pft,nind_kill(pft)
       
c	Send a/g non-combusted biomass of dead trees to fuel cats & a/g litter 
c       But make them available for burning in the next year 

       fuel_update_1hr(pft) = fuel_update_1hr(pft) + 
     *                (nind_kill(pft) * (1-ck(pft)) *
     *                (lm_ind(pft) + 0.045 * ag_non_leaf_ind(pft)))

       fuel_update_10hr(pft) = fuel_update_10hr(pft) + 
     *               (nind_kill(pft) * (1-ck(pft)) *
     *               (0.075 * ag_non_leaf_ind(pft)))


c     cambial damage and left-overs from crown scorching
       fuel_update_100hr(pft)  = fuel_update_100hr(pft) + 
     *                (nind_kill(pft) * (1-ck(pft)) * 
     *               (0.21 * ag_non_leaf_ind(pft)) +  
     *		     nind_kill(pft) * ck(pft) *
     *               (0.95 * 0.21 * ag_non_leaf_ind(pft)))


c      cambial damage and crown scorching together
       fuel_update_1000hr(pft) = fuel_update_1000hr(pft) + 
     *                 nind_kill(pft) * (0.67 * ag_non_leaf_ind(pft)) 


c	Send roots of dead trees to b/g litter, Kirsten: this can be done daily, no effect on fire spread
       litter_bg(pft) = litter_bg(pft) + rm_ind(pft) * nind_kill(pft)



c       update individual density and send dead biomass as a result of postfire mort. to
c       litter and fuel classes, resp. 
        
c    add DEAD FUEL combustion to fire carbon flux per PFT

        dcflux_fire_pft(d,pft)=dcflux_fire_pft(d,pft)+
     *   ((fc_1hr(pft)+fc_10hr(pft)+fc_100hr(pft))*0.45)+fc_1000hr(pft)

c    total carbon flux from life fuel and dead fuel combustion
        mcflux_fire_pft(m,pft)=mcflux_fire_pft(m,pft)+
     *               dcflux_fire_pft(d,pft)
        acflux_fire_pft(pft)=acflux_fire_pft(pft)+dcflux_fire_pft(d,pft)
c       print*,'acflux',d,pft,acflux_fire_pft(pft),dcflux_fire_pft(d,pft)

c	Number indivs NOT killed by fire
       nind_nkill(pft) = nind(pft) - nind_kill(pft)
c        print*,'nind not killed',pft,nind_nkill(pft),
c     *                               nind(pft),nind_kill(pft)

c	ALLAN: If at least partial crown kill occurs, and some fire-affected
c       trees survive fire, and 
c	number of trees not killed by fire is non-zero, then redistribute
c	leaf, sapwood, and heartwood and rootmass among survivors

c   Kirsten: error in setting the brackets: the last one wasn't seen by the
c   compiler. Check: performance ok now?

      if ((ck(pft).gt.0.0).and.
     *  ((nind_fa(pft)-nind_kill(pft)).gt.0.0) 
     *  .and. (nind_nkill(pft) .gt. 0.0)) then

       prop_fa_nk(pft) = (nind_fa(pft)-nind_kill(pft))/nind_nkill(pft)

       lm_ind(pft)=lm_ind(pft)-prop_fa_nk(pft)*ck(pft)*lm_ind(pft)

       sm_ind(pft) = sm_ind(pft) - prop_fa_nk(pft) * ck(pft) *
     *   (0.045 + 0.075 + (0.21 * 0.05)) * sm_ind(pft)

       hm_ind(pft) = hm_ind(pft) - prop_fa_nk(pft) * ck(pft) *
     *   (0.045 + 0.075 + (0.21 * 0.05)) * hm_ind(pft)

c   KIRSTEN: where is the rootmass??? No rootmass treatment, otherwise carbon balance not closed
c    ALLAN: see my original code below. 'No rootmass treatment' still applies
c       rm_ind(pft) = rm_ind(pft) - prop_fa_nk(pft)*ck(pft)*rm_ind(pft) 

       endif 

c	ALLAN:  Assume the proportional decrease in individual root mass due 
c	to stress of surviving a crown-fire among fire-affected indivs = 0%. 
c	So, rm_ind(pft) stays unchanged.  
c	Otherwise, implement proportional reduction (z) = 0.05, say:
c	rm_ind(pft) = rm_ind(pft) -  (prop_fa_nk(pft) * ck(pft) * 
c     *                rm_ind(pft) * 0.05) 
c	litter_bg(pft) = litter_bg(pft) + (prop_fa_nk(pft) * ck(pft) *
c     *                   rm_ind(pft) * 0.05)


c	Update number of indivs surviving fire, ready for the next day. 
       nind(pft)=nind_nkill(pft)
c       print*,'nind',pft,nind(pft)
c      pause
c6603	continue

      else 

c       dead grass
        dcflux_fire_pft(d,pft)=dcflux_fire_pft(d,pft)+(fc_1hr(pft)*0.45)
c        fuel_1hr(pft)=fuel_1hr(pft)-(fc_1hr(pft)*0.45)
c        litter_ag(pft)=litter_ag(pft)-(fc_1hr(pft)*0.45)


c      live grass        
       dcflux_fire_pft(d,pft)=dcflux_fire_pft(d,pft)+fc_lg(pft)
c  KIRSTEN: no update of rootmass and litter_bg here, otherwise carbon balance not closed
c       lm_ind(pft)=lm_ind(pft)-fc_lg(pft)

c    Kirsten: update of livegrass only through lm_ind, due to combination of dphen and livegrass 
c       consumption
c       livegrass=livegrass-fc_lg(pft) 
        pot_fc_lg_temp(pft)=lm_ind(pft)-fc_lg(pft) 
c        lm_ind(pft)=lm_ind(pft)-fc_lg(pft) 
 
c     Kirsten: never forget to protect variables against division by zero
        ratio_lg_reduc=0.0 
c        if (pot_fc_lg_temp(pft).gt.0.0)
c     *      ratio_lg_reduc=lm_ind(pft)/pot_fc_lg_temp(pft)

	if (lm_ind(pft).gt.0) 
     *     ratio_lg_reduc=pot_fc_lg_temp(pft)/lm_ind(pft) 
           lm_ind(pft)=lm_ind(pft)-fc_lg(pft) 

c  KIRSTEN: no update of rootmass and litter_bg here, otherwise carbon balance not closed
c        rm_ind(pft)=rm_ind(pft)-(ratio_lg_reduc*rm_ind(pft))
c        litter_bg(pft)=litter_bg(pft)+(ratio_lg_reduc*rm_ind(pft))

c       Kirsten: this way we loose all the roots
c       lm_ind(pft)=pot_fc_lg_temp(pft)*0.45 
        litter_bg(pft)=litter_bg(pft) + (1-ratio_lg_reduc) * rm_ind(pft)
        rm_ind(pft)=ratio_lg_reduc * rm_ind(pft)



c    total carbon flux from life fuel and dead fuel combustion
        mcflux_fire_pft(m,pft)=mcflux_fire_pft(m,pft)+
     *               dcflux_fire_pft(d,pft)
        acflux_fire_pft(pft)=acflux_fire_pft(pft)+dcflux_fire_pft(d,pft)
C       print*,'acflux',d,pft,acflux_fire_pft(pft),dcflux_fire_pft(d,pft)
C       print*,'d,pft,acflux_fire_pft(pft),dcflux_fire_pft(d,pft)' 
C       print*,d,pft,acflux_fire_pft(pft),dcflux_fire_pft(d,pft) 

c6602	continue
      endif   

      endif 


      endif     

      enddo     

c6601	continue

         m_fc_crown(m)=m_fc_crown(m)+fc_crown(d)
         an_fc_crown=an_fc_crown+fc_crown(d)


c    Kirsten:  adding up of daily to monthly values
         if (d.eq.month(m)) then
            if (count.gt.1) then
               mfdi(m)=mfdi(m)/count
               count_fdi=count_fdi+count
            endif
            if (count_int.gt.1) then
               m_i_surface(m)=m_i_surface(m)/count_int
               count_yr=count_yr+count_int
            endif
            m=m+1
            count_int=0
            count=0
          elseif (d.le.month(m).and.afire_frac.eq.1.0) then
            mfdi(m)=mfdi(m)/count
            count_fdi=count_fdi+count
            m_i_surface(m)=m_i_surface(m)/count_int
            count_yr=count_yr+count_int
            count=0
            count_int=0
          endif

c     stop daily loop if entire grid cell burnt
         if (afire_frac.eq.1.0) then
c            print*,year,' day= ',d,afire_frac
            goto 200
         endif

c       pause
       enddo 

200   continue

c  send dead biomass from cambial damage to fuel classes at the end of the year
c  to avoid that the dead fuel influence fire spread of the actual fire season (standing dead biomass)
       do pft=1,npft
       if (present(pft)) then
       fuel_1hr(pft) = fuel_1hr(pft) + fuel_update_1hr(pft)
       fuel_10hr(pft) = fuel_10hr(pft) + fuel_update_10hr(pft)
       fuel_100hr(pft) = fuel_100hr(pft) + fuel_update_100hr(pft)
       fuel_1000hr(pft) = fuel_1000hr(pft) + fuel_update_1000hr(pft)

       litter_ag(pft) = litter_ag(pft) + fuel_update_1hr(pft) +
     *                  fuel_update_10hr(pft) + fuel_update_100hr(pft) +
     *                  fuel_update_1000hr(pft)    
       endif
       enddo

c    summation of monthly and annual variables

       do m=1,12 
        do pft=1,npft  
           mcflux_fire(m)=mcflux_fire(m)+mcflux_fire_pft(m,pft)  
        enddo  
       enddo  

      do pft=1,npft
         acflux_fire=acflux_fire+acflux_fire_pft(pft)
      enddo

        an_fdi=an_fdi/count_fdi
        an_i_surface=an_i_surface/count_yr

c        

c      ACHTUNG: Check this, is this doubling the effort?
c         an_areafires is in ha, area in m��
c          afire_frac=((an_areafires*10000.0)/area)

*	  if (year.gt.1075)
c       print*,'fdi,nof,areab,ffrac,int_surface',year,an_fdi,
c     * annum_fire,an_areafires,afire_frac,an_i_surface
c      pause
c     Calculate monthly trace gas emission resulting from biomass burning (gSpecies/m��)


       do m=1,12 
        do x=1,6 

         do pft=1,npft
         if (present(pft)) then

            if (mcflux_fire_pft(m,pft).gt.0.0) then 
c            print*,'mcflux_fire (m,pft) gt zero

            mcflux_trace_pft(m,x,pft)=(mcflux_fire_pft(m,pft)*
     *                                    ef_trace(pft,x)/0.45)
          mcflux_trace(m,x)=mcflux_trace(m,x)+mcflux_trace_pft(m,x,pft)
           else
               mcflux_trace(m,x)=0.0
           endif

         endif 

        enddo
           acflux_trace(x)=acflux_trace(x)+mcflux_trace(m,x)
       enddo
      enddo
*            print*,'acflux_trace(x)',(acflux_trace(x),x=1,6)
c            if (year.eq.1102) then
c               print*,'mcflux_trace(m,x)','year=',year
c               do x=1,6
c                  print*,(mcflux_trace(m,x),m=1,12)
c               enddo
c               print*,'acflux_trace(x),(acflux_trace(x),x=1,6)',
c     *                'year=',year 
c               print*,(acflux_trace(x),x=1,6) 
c            endif

*      print*,'litc, reg end',year,(litter_ag(pft),pft=1,npft)
*      if (year.ge.1000) print*,'firec', year,acflux_fire,afire_frac
*        print*,'fuel_1hr reg end',(fuel_1hr(pft),pft=1,npft)
*        print*,'fuel_10hr reg end',(fuel_10hr(pft),pft=1,npft)
*        print*,'fuel_100hr reg end',(fuel_100hr(pft),pft=1,npft)
*        print*,'fuel_1000hr reg end',(fuel_1000hr(pft),pft=1,npft)
*        print*,'lm reg end',(lm_ind(pft),pft=1,npft)
*        print*,'rm reg end',(rm_ind(pft),pft=1,npft)
*        print*,'hm reg end',(hm_ind(pft),pft=1,npft)
*        print*,'sm reg end',(sm_ind(pft),pft=1,npft)
*        print*,'litter_bg reg end',(litter_bg(pft),pft=1,npft)
*        print*,'nind',(nind(pft),pft=1,npft) 
*        print*,'nind_old',(nind_old(pft),pft=1,npft)
*      pause

      endif 

*      print*,'subr fire',year,acflux_fire

c	ALLAN: For the sake of neatness, account for rounding errors 
c	leading to slightly negative balances.

	do pft=1,npft
           if (present(pft)) then
	     	if (litter_ag(pft).lt.0.0) then
		litter_ag(pft)=0.0
	     	endif
	     	if (litter_bg(pft).lt.0.0) then
                litter_bg(pft)=0.0
             	endif
		if (fuel_1hr(pft).lt.0.0) then
		fuel_1hr(pft)=0.0
		endif
                if (fuel_10hr(pft).lt.0.0) then
                fuel_10hr(pft)=0.0
                endif
                if (fuel_100hr(pft).lt.0.0) then
                fuel_100hr(pft)=0.0
                endif
                if (fuel_1000hr(pft).lt.0.0) then
                fuel_1000hr(pft)=0.0
                endif
		if (lm_ind(pft).lt.0.0) then
		lm_ind(pft)=0.0
		endif
                if (sm_ind(pft).lt.0.0) then
                sm_ind(pft)=0.0
                endif
                if (hm_ind(pft).lt.0.0) then
                hm_ind(pft)=0.0
                endif
                if (rm_ind(pft).lt.0.0) then
                rm_ind(pft)=0.0
                endif
	   endif
	enddo


*      print*,'subr fire',year,acflux_fire
      return
      end


c----------------------------------------------------------------------
c     FUNCTION FIRE DANGER INDEX
c     Calculation of the fire danger index using the Nesterov Index
c


      subroutine fire_danger_index(d_fdi,dlm,dlm_lg,dtemp_min,
     *    dtemp_max,dprec,d,moistfactor,fuel_1hr_total,
     *    fuel_10hr_total, nesterov,
     *    fuel_100hr_total,dead_fuel,char_moistfactor,ratio_dead_fuel,
     *    ratio_live_fuel,dlm_1hr,dlm_10hr, dlm_100hr, dlm_1000hr,
     *    year, alpha, sigma_1hr, sigma_10hr, sigma_100hr )

      implicit none

      real dtemp_min(1:365),dtemp_max(1:365),dprec(1:365)
      real dw1(1:365),dlm(1:365),d_fdi(1:365),dlm_lg(1:365)
      real nesterov(1:365)
      real moistfactor
      real fuel_1hr_total,fuel_10hr_total,fuel_100hr_total
      real dead_fuel
      integer d,year
      real sigma_1hr, sigma_100hr, sigma_10hr ! Needed to calculate \alpha for dead fuel moisture
      real alpha
c        parameter (alpha=0.0015)
      real alpha_1hr, alpha_10hr,alpha_100hr
      real alpha_livegrass,alpha_1000hr

c        parameter (alpha_1hr=0.001,alpha_10hr=0.00005424) 
c        parameter (alpha_100hr=0.00001485, alpha_1000hr = 0.000001) 
c        parameter (alpha_livegrass=0.0005) 
c     Allan: Alpha values for livegrass and 1000hr dead fuels are particularly
c            subjective 
c     KIRSTEN:i.e. values made up, revisit alpha_livegrass with livegrass
c             moisture driven upper soil moisture
      real   dlm_1hr(1:365), dlm_10hr(1:365)  
      real   dlm_100hr(1:365), dlm_1000hr (1:365)
      real char_moistfactor, ratio_dead_fuel,ratio_live_fuel

c     LOCAL VARIABLES
      real ni_acc
      real d_NI,alpha_fuel,char_alpha_fuel
      real fdi_spread

c     initialise

      d_NI=0.0
      fdi_spread=0.0
      alpha_fuel=0.0
      alpha_1hr = alpha
      alpha_10hr = alpha/(sigma_1hr/sigma_10hr)
      alpha_100hr = alpha/(sigma_1hr/sigma_100hr)
c     calculate Nesterov Index   equ. (2)
        if (dprec(d).le.3.0.and.(dtemp_min(d)-4.0).ge.0.0) then

          d_NI=dtemp_max(d)*(dtemp_max(d)-dtemp_min(d)-4.0)
        else
          d_NI=0.0
        endif

        if (d_NI.gt.0.0) then
          ni_acc=ni_acc+d_NI
        else
          ni_acc=0.0
        endif
        nesterov(d) = ni_acc
c   litter moisture index, weighted per dead fuel class and
c   fuel load per dead fuel class; version with 3 alpha now confirmed

*        dlm(d)=exp(-alpha*ni_acc)
       if (dead_fuel.gt.0.0) then
        alpha_fuel=(alpha_1hr*fuel_1hr_total+
     *              alpha_10hr*fuel_10hr_total+
     *              alpha_100hr*fuel_100hr_total)/dead_fuel 


c   backcalculation of alpha_livegrass from livegrass moisture
        if (ni_acc.gt.0.0) then
          if (dlm_lg(d).gt.0.0) then
           alpha_livegrass =(log(dlm_lg(d))/ni_acc)*(-1.0)
         else
           alpha_livegrass = 0.0
          endif
        endif
c         print*,'alpha_livegrass',alpha_livegrass,ni_acc,dlm_lg(d)
c         pause

c    Allan
        char_alpha_fuel=alpha_fuel * ratio_dead_fuel
     *          + alpha_livegrass * ratio_live_fuel
       else
c         alpha_fuel=0.00001
         char_alpha_fuel=0.00001
       endif

c        dlm(d)=exp(-alpha_fuel*ni_acc)
        dlm(d)=exp(-char_alpha_fuel*ni_acc)  
       
       dlm_1hr(d) = exp(-alpha_1hr * ni_acc)
       dlm_10hr(d) = exp(-alpha_10hr * ni_acc)
       dlm_100hr(d) = exp(-alpha_100hr * ni_acc)
       dlm_1000hr(d) = exp(-alpha_1000hr * ni_acc)
      

c probability of fire spread
*        if (dlm(d).le.moistfactor) then
*           fdi_spread=1.0-(dlm(d)/moistfactor)
*            fdi_spread=max(0.0,(1.0-2.59*(dlm(d)/moistfactor)+
*     *             5.11*(dlm(d)/moistfactor)**2.0-
*     *             3.52*(dlm(d)/moistfactor)**3.0))
*        else
*           fdi_spread=0.0
*        endif

c     calculate Fire Danger Index equ. (9)
        if (d_NI.le.0.0) then
           d_fdi(d)=0.0
        else
c        d_fdi(d)=max(0.0,
c     *           (1.0-((1.0/moistfactor)*(exp(-alpha_fuel*ni_acc)))))
c
*     *           (1.0-((1.0/moistfactor)*dlm(d))))

        d_fdi(d)=max(0.0,(1.0-((1.0/char_moistfactor)*
     *                   (exp(-char_alpha_fuel*ni_acc)))))


        endif

*        print*,'ni_acc,d_fdi,dlm',d,ni_acc,moistfactor,d_fdi(d),dlm(d)
*        pause

      return
      end


c-------------------------------------------------------------------------------
c     FUNCTION HUMAN IGNITION
c     Calculation of the numbers of ignition caused by humans, which is
c     determined by population density, see equ. (13), where the value of
c     a(N_d) and k(pop_den) are equal to their mathematical expectation

      real function human_ign(popden,a_nd,year)
cf2py intent(in) popden,a_nd,year
cf2py intent(out) human_ign
      implicit none

      integer year
      real popden

      real a_nd_const,a_nd
      parameter (a_nd_const=0.251)

c       if (popden.lt.1.0) popden=1.0
c       if (a_nd.lt.0.11) a_nd=0.165
c      popden=7.5
       
c       human_ign=gamma*popden1**0.43
c       human_ign=6.8*(exp(-0.5*(popden1**0.5)))*a_nd_const*popden1
c       human_ign=20.4*(exp(-0.5*(popden**0.5)))*a_nd*popden
c        human_ign=6.8*a_nd*popden**0.43
	human_ign=30.0*(exp(-0.5*(popden**0.5)))*a_nd*popden 
c        human_ign=100
      
      end

c-------------------------------------------------------------------------------
c     FUNCTION LIGHTNING
c     Calculation of lightning strikes. Constant at the moment, in the future
c     depending on latitude and other factors
c

c      real function lightn_ign(lat)

c      implicit none

c      real lat

c      lightn=0.02
c     for case of Brandenburg
c       lightn_ign=0.2

c      end


c-------------------------------------------------------------------------------
c     subroutine RATE OF forward SPREAD

      subroutine rate_of_spread(U_front,base_wind,dens_fuel_ave,sigma,
     *            dlm,d,net_fuel,moistfactor,H, char_dens_fuel_ave, 
     *	          char_sigma, char_net_fuel, char_moistfactor,gamma)


      implicit none

      integer npft,npftpar
        parameter (npft=9,npftpar=50)
      integer d,m
      real U_front,base_wind
      real mw1(1:12),dw1(1:365),dlm(1:365)
      real dens_fuel_ave,sigma,H
      real net_fuel,moistfactor
      real dummy,dummy2
       real char_dens_fuel_ave,char_sigma
       real char_net_fuel, char_moistfactor

c     emissitivity of flames unitless
*      real emiss
*        parameter(emiss=0.3)

c     Rothermal fire spread
	
	real beta,ir,xi,eps,q_ig,phi_wind
      real gamma,gamma_aptr,moist_damp,mw_weight
	real gamma_max,bet,beta_op,a,c,b,e

c     mineral dampening coefficient
      real MINER_DAMP
        parameter(MINER_DAMP=0.41739)
      real pi
        parameter (pi=3.14159265)
      real part_dens
        parameter (part_dens=513.0)
c    initialise
      bet=0.0
      q_ig=0.0
      eps=0.0
      a=0.0
      b=0.0
      c=0.0
      e=0.0
      phi_wind=0.0
      xi=0.0
      gamma_max=0.0
      gamma_aptr=0.0
      mw_weight=0.0
      moist_damp=0.0
      ir=0.0
      dummy=0.0
      dummy2=0.0

c     start of function

c        bet=dens_fuel_ave*0.200395*(sigma**(-0.8189))/513.0     
       beta = char_dens_fuel_ave / part_dens
       beta_op = 0.200395*(sigma**(-0.8189))
c       beta_op = 0.200395*(char_sigma**(-0.8189))
       bet = beta / beta_op

c     heat of pre-ignition
c     ACHTUNG: check influence of litter moisture
       q_ig=581.0+2594.0*dlm(d)
c     effective heating number
       eps=exp(-4.528/sigma)
c        eps=exp(-4.528/char_sigma)
c     influence of wind speed
           b=0.15988*(sigma**0.54)
           c=7.47*(exp(-0.8711*(sigma**0.55)))
           e=0.715*(exp(-0.01094*sigma))
c           b = 0.15988 * (char_sigma**0.54)
c           c = 7.47 * (exp(-0.8711 * (char_sigma**0.55)))
c           e = 0.715 * (exp(-0.01094 * char_sigma))

       phi_wind=c*(base_wind**b)*(bet**(-e))
c     propagating flux
      if (sigma.le.0.00001) then
c       if (char_sigma.le.0.00001) then

        xi=0.0
       else
c        xi=(exp((0.792+3.7597*(sigma**0.5))*
c     *     ((dens_fuel_ave/513.0)+0.1)))/
c     *     (192.0+7.9095*sigma)
        xi = (exp((0.792 + 3.7597 * (sigma**0.5)) *
     *     (beta + 0.1))) / (192 + 7.9095 * sigma)

       endif

c    reaction intensity
           a=8.9033*(sigma**(-0.7913))
c           a = 8.9033 * (char_sigma**(-0.7913))
          if (sigma.le.0.00001) then
c           if (char_sigma.le.0.00001) then
            dummy=0.0
           else
            dummy=exp(a*(1.0-bet))
           endif
           gamma_max=1.0/(0.0591+2.926*(sigma**(-1.5)))
c           gamma_max = 1.0 / (0.0591 + 2.926 * (char_sigma**(-1.5)))

*        gamma=gamma_max*(bet**a)*exp(a*(1.0-bet))

         gamma_aptr=gamma_max*(bet**a)*dummy

c        if (moistfactor.gt.0.0) then
        if (char_moistfactor.gt.0.0) then

c          mw_weight=dlm(d)/moistfactor
         mw_weight=dlm(d)/char_moistfactor
        else
          mw_weight=0.0
        endif

        moist_damp=max(0.0,(1.0-(2.59*mw_weight)+
     *     (5.11*(mw_weight**2.0))-(3.52*(mw_weight**3.0))))

c        ir=gamma*net_fuel*H*moist_damp*MINER_DAMP
         ir=gamma_aptr * char_net_fuel * H * moist_damp * MINER_DAMP

c        for use in postfire mortality (tau_r)
         gamma=gamma_aptr*moist_damp*MINER_DAMP
      
*        print*,'ir',ir,gamma,net_fuel,H,moist_damp,MINER_DAMP
c    reaction intensity end

c        if(dens_fuel_ave.le.0.0.or.eps.le.0.0.or.q_ig.le.0.0) then
         if ((char_dens_fuel_ave.le.0.0).or.
     *      (eps.le.0.0).or.(q_ig.le.0.0)) then
          U_front=0.0
        else
c          U_front=(ir*xi*(1.0+phi_wind))/(dens_fuel_ave*eps*q_ig)
          U_front=(ir * xi * (1.0 + phi_wind)) /
     *            (char_dens_fuel_ave * eps * q_ig)

        endif

*       if (d.gt.220)
c       print*,'U_front',d,U_front,
c     *      xi,sigma,
c     *      dens_fuel_ave
*     *           c,base_wind,b,bet,e,sigma

*      pause
      return
      end

c-------------------------------------------------------------------------------
c    subroutine for calculation of fuel consumption in the area affected by fire

      subroutine fuel_consumption(npft,present,fuel_consum,fire_frac,
     *  fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,pot_fc_lg,
     *  tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,fc_1hr,fc_lg,
     *  fc_10hr,fc_100hr,fc_1000hr,cf,char_moistfactor,
     *   dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *   moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *   moistfactor_100hr, moistfactor_1000hr)


      implicit none

      integer d,pft,npft
      real fuel_1hr_total,livegrass
      real moistfactor,dw1(1:365),fire_frac,dlm(1:365)
      real pot_fc_lg(1:npft),fuel_1hr(1:npft),fuel_10hr(1:npft)
      real fuel_100hr(1:npft),fuel_1000hr(1:npft)
      real MINER_TOT
      real fuel_consum,fc_1hr(1:npft),fc_lg(1:npft),fc_10hr(1:npft)
      real fc_100hr(1:npft),fc_1000hr(1:npft)
      real tau_l(1:npft),cf(1:npft)
      real char_moistfactor
      real moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr
      real moistfactor_100hr, moistfactor_1000hr
      real  dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365)  
      real  dlm_100hr(1:365), dlm_1000hr (1:365)
      logical present(1:npft),tree(1:npft)

c     local variables
       real moist_lg_d1hr,moist_1hr,moist_10_100hr
      real net_moist_lg, net_moist_1hr, net_moist_10hr
      real net_moist_100hr, net_moist_1000hr
      real pot_fc_1hr_total,mw_weight
      real cf_lg, cf_1hr,cf_10hr,cf_100hr,cf_1000hr

c     initialise
      fuel_consum=0.0
      moist_lg_d1hr=0.0
      moist_1hr=0.0
      moist_10_100hr=0.0
      mw_weight=0.0
      cf_1hr=0.0
      cf_10hr=0.0
      cf_100hr=0.0

       cf_1000hr=0
       cf_lg=0

c       net_moist_1hr = dlm_1hr(d) / moistfactor_1hr
c       net_moist_10hr = dlm_10hr(d) / moistfactor_10hr
c       net_moist_100hr = dlm_100hr(d) / moistfactor_100hr
c       net_moist_1000hr = dlm_1000hr(d) / moistfactor_1000hr

c     fuel consumption depending on fuel moisture
c     influence of livefuel on 1hr fuel moisture content
c     ACHTUNG: change with dlm_lg as a function of upper soil moisture
c     CHECK method

       moist_lg_d1hr=dlm_1hr(d)+(dlm_lg(d)*(livegrass/fuel_1hr_total))
c       print*,'moist_lg',d,moist_lg_d1hr,dlm_1hr(d),dlm_lg(d),
c     *        (livegrass/fuel_1hr_total),livegrass
c       pause
c      moist_lg_d1hr=dlm(d)

c     ACHTUNG: check with litter moisture
      if(moistfactor.gt.0.0) then
        moist_1hr=moist_lg_d1hr/char_moistfactor
        moist_10_100hr=dlm(d)/char_moistfactor
      else
        moist_1hr=1.0
        mw_weight=0.0
        moist_10_100hr=1.0
      endif


      do pft=1,npft
        fc_1hr(pft)=0.0
        fc_lg(pft)=0.0
        fc_10hr(pft)=0.0
        fc_100hr(pft)=0.0
        fc_1000hr(pft)=0.0
c        tau_l(pft)=0.0 

       if (present(pft)) then

c     1hr fuel consumption
      if(moist_1hr.le.0.18) then
         fc_1hr(pft)=1.0*(1.0-MINER_TOT)*(fuel_1hr(pft)/0.45)*fire_frac
         cf_1hr=1.0  
      else
         if(moist_1hr.gt.0.18.and.moist_1hr.le.0.73) then
            fc_1hr(pft)=(1.10-0.62*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft)/0.45)*fire_frac
            cf_1hr=1.10-0.62*moist_1hr
         else
           if(moist_1hr.gt.0.73.and.moist_1hr.le.1.0) then
            fc_1hr(pft)=(2.45-2.45*moist_1hr)*(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft)/0.45)*fire_frac
            cf_1hr=2.45-2.45*moist_1hr
           else
            fc_1hr(pft)=0.0
            cf_1hr=0.0
           endif
         endif
      endif


c    time required for cambial kill, 1hr fuel [gBiomass/cm��]
c	/1e4 to get from 1m2 to cm2, the units used by P&R(1986) 
c         tau_l(pft)=tau_l(pft)+((fuel_1hr(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_1hr)**0.5)))

c     Cf only for tau_r
          cf(pft)=cf(pft)+cf_1hr

*      print*,'tau_l, 1hr',d,dlm(d),pft,(fuel_1hr(pft)/0.45),cf_1hr,
*     * tau_l(pft)


c       print*,'fuel_1hr',pft,fuel_1hr(pft),fc_1hr(pft)
c     livegrass consumption
      if (.not.tree(pft)) then
c      if(moist_1hr.le.0.18) then

c  ACHTUNG can we use dlm_lg directly? Thresholds don't match, we set 
c  dlm_lg to zero at dw1=0.3

      if(dlm_lg(d).le.0.18) then
         fc_lg(pft)=1.0*(1.0-MINER_TOT)*pot_fc_lg(pft)*fire_frac
      else
cc         if(moist_1hr.gt.0.18.and.moist_1hr.le.0.73) then
         if(dlm_lg(d).gt.0.18.and.dlm_lg(d).le.0.73) then
           fc_lg(pft)=(1.10-0.62*dlm_lg(d))*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft)*fire_frac
        else
cc           if(moist_1hr.gt.0.73.and.moist_1hr.le.1.0) then
           if(dlm_lg(d).gt.0.73.and.dlm_lg(d).le.1.0) then
            fc_lg(pft)=(2.45-2.45*dlm_lg(d))*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft)*fire_frac
           else
            fc_lg(pft)=0.0
           endif
        endif
      endif
      endif 

c     10hr fuel consumption
      if(moist_10_100hr.le.0.12) then
        fc_10hr(pft)=1.0*(1.0-MINER_TOT)*(fuel_10hr(pft)/0.45)*fire_frac
        cf_10hr=1.0
      else
       if(moist_10_100hr.gt.0.12.and.moist_10_100hr.le.0.51) then
          fc_10hr(pft)=(1.09-0.72*moist_10_100hr)*
     *               (1.0-MINER_TOT)*(fuel_10hr(pft)/0.45)*fire_frac
          cf_10hr=1.09-0.72*moist_10_100hr
         else
           if(moist_10_100hr.gt.0.51.and.moist_10_100hr.le.1.0) then
             fc_10hr(pft)=(1.47-1.47*moist_10_100hr)*
     *                (1.0-MINER_TOT)*(fuel_10hr(pft)/0.45)*fire_frac
             cf_10hr=1.47-1.47*moist_10_100hr
           else
             fc_10hr(pft)=0.0
             cf_10hr=0.0
           endif
         endif
      endif


c    time required for cambial kill, 10hr fuel
c         tau_l(pft)=tau_l(pft)+((fuel_10hr(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_10hr)**0.5)))

      

c     100hr fuel consumption
      if(moist_10_100hr.le.0.38) then
        fc_100hr(pft)=(0.98-0.85*moist_10_100hr)*(1.0-MINER_TOT)
     *                *(fuel_100hr(pft)/0.45)*fire_frac
        cf_100hr=0.98-0.85*moist_10_100hr
      else
        if(moist_10_100hr.gt.0.38.and.moist_10_100hr.le.1.0) then
          fc_100hr(pft)=(1.06-1.06*moist_10_100hr)*(1.0-MINER_TOT)
     *               *(fuel_100hr(pft)/0.45)*fire_frac
          cf_100hr=1.06-1.06*moist_10_100hr
        else
           fc_100hr(pft)=0.0
          cf_100hr=0.0
        endif
      endif


c    time required for cambial kill, 100hr fuel
c         tau_l(pft)=tau_l(pft)+((fuel_100hr(pft)/0.45/1e4)*
c     *                     (1.0-((1.0-cf_100hr)**0.5)))

c    calculating Cf from fuel classes, tau_r=Cf/Gamma tau_l=2*tau_r
        cf(pft)=(cf_1hr+cf_10hr+cf_100hr)/3.0


c    Peterson& Ryan: 39.4*tau_l, assuming particle density of 510, 
c    we have 513, thus 39.09*tau_l
c         tau_l(pft)=39.09*tau_l(pft)

c    1000hr fuel consumption, not influencing rate of spread or I_surface (Rothermel 1972)

        fc_1000hr(pft)=(-0.8*mw_weight+0.8)*(1.0-MINER_TOT)*
     *                    fuel_1000hr(pft)*fire_frac

c	Allan: Approximate form. No data.

c    total fuel consumption (without 1000hr fuel) in g biomass per m��
c    Used to calculate fire intensity in the FLAMING FRONT.
      fuel_consum=fuel_consum+fc_1hr(pft)+fc_10hr(pft)+fc_100hr(pft)
                 

*       print*,'fuel classes',pft,pot_fc_1hr(pft),pot_fc_10hr(pft),
*     *           pot_fc_100hr(pft)
c      if (.not.tree(pft).and.fire_frac.gt.0.0)
c     *         print*,'fuel_consum',d,pft,fc_lg(pft)
c     *      fc_1hr(pft),fc_10hr(pft),fc_100hr(pft)
c      pause
       endif  

      enddo 


      return
      end

c-------------------------------------------------------------------------------
c    subroutine for calculation of fuel consumption in the area affected by fire

      subroutine fuel_consumption_netmoist(npft,present,fuel_consum,
     *  fire_frac,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,
     *  pot_fc_lg,tree,fuel_1hr_total,moistfactor,dlm,d,MINER_TOT,fc_1hr
     *  ,fc_lg,fc_10hr,fc_100hr,fc_1000hr,tau_l,char_moistfactor,
     *   dlm_lg, dlm_1hr, dlm_10hr, dlm_100hr, dlm_1000hr,
     *   moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr,
     *   moistfactor_100hr, moistfactor_1000hr)


      implicit none

      integer d,pft,npft
      real fuel_1hr_total,livegrass
      real moistfactor,dw1(1:365),fire_frac,dlm(1:365)
      real pot_fc_lg(1:npft),fuel_1hr(1:npft),fuel_10hr(1:npft)
      real fuel_100hr(1:npft),fuel_1000hr(1:npft)
      real MINER_TOT
      real fuel_consum,fc_1hr(1:npft),fc_lg(1:npft),fc_10hr(1:npft)
      real fc_100hr(1:npft),fc_1000hr(1:npft)
      real tau_l(1:npft)
      real char_moistfactor
      real moistfactor_livegrass, moistfactor_1hr,moistfactor_10hr
      real moistfactor_100hr, moistfactor_1000hr
      real  dlm_lg (1:365), dlm_1hr(1:365), dlm_10hr(1:365)  
      real  dlm_100hr(1:365), dlm_1000hr (1:365)
      logical present(1:npft),tree(1:npft)

c     local variables
       real moist_lg_d1hr,moist_1hr,moist_10_100hr
      real net_moist_lg, net_moist_1hr, net_moist_10hr
      real net_moist_100hr, net_moist_1000hr
      real pot_fc_1hr_total,mw_weight
      real cf_lg, cf_1hr,cf_10hr,cf_100hr,cf_1000hr

c     initialise
      fuel_consum=0.0
      moist_lg_d1hr=0.0
      moist_1hr=0.0
      moist_10_100hr=0.0
      mw_weight=0.0
      cf_1hr=0.0
      cf_10hr=0.0
      cf_100hr=0.0

       cf_1000hr=0
       cf_lg=0
c     dlm_lg as a function of upper soil moisture
       net_moist_lg = dlm_lg(d) / moistfactor_livegrass 
       net_moist_1hr = dlm_1hr(d) / moistfactor_1hr
       net_moist_10hr = dlm_10hr(d) / moistfactor_10hr
       net_moist_100hr = dlm_100hr(d) / moistfactor_100hr
c       if (d.gt.200) then
c          print*,'day= ',d,dlm_1000hr(d), moistfactor_1000hr
c       endif
       net_moist_1000hr = dlm_1000hr(d) / moistfactor_1000hr
c       print*,'moistfactor',moistfactor_livegrass,moistfactor_1hr,
c     *   moistfactor_10hr,moistfactor_100hr,moistfactor_1000hr
c       pause
c     fuel consumption depending on fuel moisture
c     ACHTUNG: This is using upper soil moisture until fuel moisture can be calculated
c     per dead fuel class
       

c        if (fuel_1hr(pft).lt.0) print*,'error line 7245 1hr ',
c     *   'fuel_1hr(pft),fc_1hr(pft)',fuel_1hr(pft),fc_1hr(pft)
c        if (fuel_10hr(pft).lt.0) print*,'error line 7245 10hr ',
c     *   'fuel_10hr(pft),fc_10hr(pft)',fuel_10hr(pft),fc_10hr(pft)
c        if (fuel_100hr(pft).lt.0) print*,'error line 7245 100hr ',
c     *   'fuel_100hr(pft),fc_100hr(pft)',fuel_100hr(pft),fc_100hr(pft)
c        if (fuel_1000hr(pft).lt.0) print*,'error line 7245 1000hr ',
c     *   'fuel_1000hr(pft),fc_1000hr(pft)',
c     *   fuel_1000hr(pft),fc_1000hr(pft)

      do pft=1,npft
        fc_1hr(pft)=0.0
        fc_lg(pft)=0.0
        fc_10hr(pft)=0.0
        fc_100hr(pft)=0.0
        fc_1000hr(pft)=0.0
        tau_l(pft)=0.0 

       if (present(pft)) then

c     1hr fuel consumption
      if(net_moist_1hr.le.0.18) then
         cf_1hr=1.0 
         fc_1hr(pft)= cf_1hr *(1.0-MINER_TOT)*(fuel_1hr(pft)/0.45)
     *     *fire_frac
      else
         if(net_moist_1hr.gt.0.18.and.net_moist_1hr.le.0.73) then
             cf_1hr=1.10-0.62*net_moist_1hr
             fc_1hr(pft)=cf_1hr *(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft)/0.45)*fire_frac
         else
           if(net_moist_1hr.gt.0.73.and.net_moist_1hr.le.1.0) then
            cf_1hr=2.45-2.45*net_moist_1hr
            fc_1hr(pft)=cf_1hr*(1.0-MINER_TOT)*
     *                   (fuel_1hr(pft)/0.45)*fire_frac
           else 
            fc_1hr(pft)=0.0
            cf_1hr=0.0
           endif
         endif
      endif

c    time required for cambial kill, 1hr fuel [gBiomass/cm��]
c	/1e4 to get from 1m2 to cm2, the units used by P&R(1986) 
         tau_l(pft)=tau_l(pft)+((fuel_1hr(pft)/0.45/1e4)*
     *                     (1.0-((1.0-cf_1hr)**0.5)))

*      print*,'tau_l, 1hr',d,dlm(d),pft,(fuel_1hr(pft)/0.45),cf_1hr,
*     * tau_l(pft)

c     livegrass consumption. Same slopes as 1hr dead. A.
      if (.not.tree(pft)) then
      if(net_moist_lg.le.0.18) then
         fc_lg(pft)=1.0*(1.0-MINER_TOT)*
     *              pot_fc_lg(pft)*fire_frac
      else
         if(net_moist_lg.gt.0.18.and.net_moist_lg.le.0.73) then
           fc_lg(pft)=(1.10-0.62*net_moist_lg)*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft)*fire_frac
         else
           if(net_moist_lg.gt.0.73.and.net_moist_lg.le.1.0) then
            fc_lg(pft)=(2.45-2.45*net_moist_lg)*(1.0-MINER_TOT)*
     *                 pot_fc_lg(pft)*fire_frac
           else  
            fc_lg(pft)=0.0
           endif
         endif
      endif
       fc_lg(pft)= 0.05 * pot_fc_lg(pft) 
      fc_lg(pft)=0.0
      endif

c     10hr fuel consumption
      if(net_moist_10hr.le.0.12) then
        cf_10hr=1.0 
       fc_10hr(pft)=cf_10hr*(1.0-MINER_TOT)*
     *              (fuel_10hr(pft)/0.45)*fire_frac
      else
       if(net_moist_10hr.gt.0.12.and.net_moist_10hr.le.0.51) then
          cf_10hr=1.09-0.72*net_moist_10hr
          fc_10hr(pft)=cf_10hr*(1.0-MINER_TOT)*
     *           (fuel_10hr(pft)/0.45)*fire_frac
         else
           if(net_moist_10hr.gt.0.51.and.net_moist_10hr.le.1.0) then
             cf_10hr=1.47-1.47*net_moist_10hr
             fc_10hr(pft)=(1.47-1.47*net_moist_10hr)*
     *             (1.0-MINER_TOT)*(fuel_10hr(pft)/0.45)*fire_frac
           else 
             fc_10hr(pft)=0.0
             cf_10hr=0.0
           endif
         endif
      endif

c    time required for cambial kill, 10hr fuel
         tau_l(pft)=tau_l(pft)+((fuel_10hr(pft)/0.45/1e4)*
     *                     (1.0-((1.0-cf_10hr)**0.5)))

c     100hr fuel consumption
      if(net_moist_100hr.le.0.38) then
        cf_100hr=0.98-0.85*net_moist_100hr
        fc_100hr(pft)=cf_100hr*(1.0-MINER_TOT)*
     *               (fuel_100hr(pft)/0.45)*fire_frac
      else
        if(net_moist_100hr.gt.0.38.and.net_moist_100hr.le.1.0) then
           cf_100hr=1.06-1.06*net_moist_100hr
           fc_100hr(pft)=cf_100hr*(1.0-MINER_TOT)
     *               *(fuel_100hr(pft)/0.45)*fire_frac
        else 
           fc_100hr(pft)=0.0
           cf_100hr=0.0
        endif
      endif


c    time required for cambial kill, 100hr fuel
         tau_l(pft)=tau_l(pft)+((fuel_100hr(pft)/0.45/1e4)*
     *                     (1.0-((1.0-cf_100hr)**0.5)))
c    Peterson& Ryan: 39.4*tau_l, assuming particle density of 510, we have 513, thus 39.09*tau_l
         tau_l(pft)=39.09*tau_l(pft)

c    1000hr fuel consumption, not influencing rate of spread or I_surface (Rothermel 1972)

        if (net_moist_1000hr.le.1) then
         fc_1000hr(pft)=(-0.8*net_moist_1000hr+0.8)*(1.0-MINER_TOT)*
     *                    fuel_1000hr(pft)*fire_frac 
        else
	fc_1000hr(pft) = 0
        endif
c	Allan: Approximate form. No data.

c    total fuel consumption (without 1000hr fuel) in g biomass per m��
c    Used to calculate fire intensity in the FLAMING FRONT.
      fuel_consum=fuel_consum+fc_1hr(pft)+fc_10hr(pft)+fc_100hr(pft)
                 

*       print*,'fuel classes',pft,pot_fc_1hr(pft),pot_fc_10hr(pft),
*     *           pot_fc_100hr(pft)
c      print*,'fuel_consum',d,pft,fuel_consum,
c     *      fc_1hr(pft),fc_10hr(pft),fc_100hr(pft)
*      pause
       endif  

      enddo 

c        if (fuel_1hr(pft).lt.0) print*,'error line 7472 1hr ',
c     *   'fuel_1hr(pft),fc_1hr(pft)',fuel_1hr(pft),fc_1hr(pft)
c        if (fuel_10hr(pft).lt.0) print*,'error line 7472 10hr ',
c     *   'fuel_10hr(pft),fc_10hr(pft)',fuel_10hr(pft),fc_10hr(pft)
c        if (fuel_100hr(pft).lt.0) print*,'error line 7472 100hr ',
c     *   'fuel_100hr(pft),fc_100hr(pft)',fuel_100hr(pft),fc_100hr(pft)
c        if (fuel_1000hr(pft).lt.0) print*,'error line 7472 1000hr ',
c     *   'fuel_1000hr(pft),fc_1000hr(pft)',
c     *   fuel_1000hr(pft),fc_1000hr(pft)

      return
      end

c -----------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE ESTABLISHMENT
c     Establishment of new individuals (saplings) of woody PFTs,
c     grass establishment, removal of PFTs not adapted to current climate,
c     update of individual structure and FPC.

      subroutine establishment(pftpar,present,survive,estab,nind,lm_ind,
     *  sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,hm_sapl,
     *  crownarea,fpc_grid,lai_ind,height,dbh,tau_c,cl_t,sla,wooddens,
     *  latosa,mprec,reinickerp,
     *  litter_ag,litter_bg,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,
     *  tree,allom1,allom2,allom3,acflux_estab,
     *  leafondays,leafoffdays,leafon,mnpp,anpp,mnpp_add,anpp_add,year)

      implicit none

c     PARAMETERS
      integer npft,npftpar,nsoilpar,year
        parameter (npft=9,npftpar=50,nsoilpar=7)
      real pi
        parameter (pi=3.14159265)
      real aprec_min_estab                
        parameter (aprec_min_estab=100.0) 
      real estab_max                      
*        parameter (estab_max=0.12)        
        parameter (estab_max=0.24)
      real nind_min                       
        parameter (nind_min=1.0E-10)      

c     ARGUMENTS
      logical present(1:npft),survive(1:npft),estab(1:npft)
      logical tree(1:npft)
      real pftpar(1:npft,1:npftpar)
      real nind(1:npft)
      real lm_ind(1:npft),sm_ind(1:npft)
      real rm_ind(1:npft),hm_ind(1:npft)
      real lm_sapl(1:npft),sm_sapl(1:npft)
      real rm_sapl(1:npft),hm_sapl(1:npft)
      real crownarea(1:npft)
      real fpc_grid(1:npft)
      real lai_ind(1:npft)
      real height(1:npft)
      real dbh(1:npft),tau_c(0:4,1:npft),cl_t(0:4,1:npft)
      real sla(1:npft)
      real wooddens,latosa,reinickerp
      real mprec(1:12)
      real litter_ag(1:npft),litter_bg(1:npft)
      real fuel_1hr(1:npft),fuel_10hr(1:npft),fuel_100hr(1:npft)
      real fuel_1000hr(1:npft)
      real allom1,allom2,allom3
      real acflux_estab
      integer leafondays(1:npft),leafoffdays(1:npft)
      logical leafon(1:npft)
      real anpp(1:npft),mnpp(1:12,1:npft)
      real mnpp_add(1:12),anpp_add        


c     LOCAL VARIABLES
      integer pft,m
      integer npft_estab  
      real aprec          
      real fpc_tree_total,fpc_tree_new 
      real estab_rate     
                          
      real estab_grid     
      real nind_old       
      real stemdiam       
      real sm_ind_temp    
      real fpc_ind        
      real crownarea_max  
      real fpc_total,fpc_new      
      real bare           
      integer ngrass
      real fpc_grass_total,fpc_grass_new
      real bare_max
      integer class
      real param1(1:npft),param2(1:npft)  
      real bt(0:4,1:npft),crown(1:npft)      
      real dbh_class(0:4,1:npft),height_class(0:4,1:npft)

c     Kill PFTs not adapted to current climate, introduce newly "adapted" PFTs

c     Initialize buffers for npps of pfts that go not present
      anpp_add=0.
      do m=1,12
         mnpp_add(m)=0.
      enddo


c     Calculate annual precipitation
      bare_max=0.0
      aprec=0.0
      do m=1,12
        aprec=aprec+mprec(m)
      enddo

      do pft=1,npft

c         Kirsten: parameter for bark thickness
          crown(pft)=pftpar(pft,44)
          param1(pft)=pftpar(pft,47)
          param2(pft)=pftpar(pft,48)

        if (present(pft).and.
     *    (.not.survive(pft).or.nind(pft).lt.nind_min)) then     

          present(pft)=.false.

c         Add up NPP of PFTs that are killed in extra balance
          do m=1,12
             mnpp_add(m)=mnpp_add(m)+mnpp(m,pft)
             anpp_add=anpp_add+mnpp(m,pft)
          enddo

c         Add killed biomass to litter

          if (tree(pft)) then
            litter_ag(pft)=litter_ag(pft)+(lm_ind(pft)+sm_ind(pft)+
     *        hm_ind(pft))*nind(pft)
c         KIRSTEN: per fuel class
            fuel_1hr(pft)=fuel_1hr(pft)+(lm_ind(pft)+0.045*sm_ind(pft)+
     *        0.045*hm_ind(pft))*nind(pft)
            fuel_10hr(pft)=fuel_10hr(pft)+(0.075*hm_ind(pft)+
     *                   0.075*sm_ind(pft))*nind(pft)
            fuel_100hr(pft)=fuel_100hr(pft)+(0.21*hm_ind(pft)+
     *                   0.21*sm_ind(pft))*nind(pft)
            fuel_1000hr(pft)=fuel_1000hr(pft)+(0.67*hm_ind(pft)+
     *                   0.67*sm_ind(pft))*nind(pft)

          else  
            litter_ag(pft)=litter_ag(pft)+lm_ind(pft)*nind(pft)
c         KIRSTEN: 1hr fuel class
            fuel_1hr(pft)=fuel_1hr(pft)+lm_ind(pft)*nind(pft)
          endif

          litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*nind(pft)
c          print*,'litter_bg pres',pft,litter_bg(pft)

        elseif (.not.present(pft).and.survive(pft).and.estab(pft)
     *    .and.aprec.ge.aprec_min_estab) then

c         Introduce PFT if conditions suitable for establishment

          present(pft)=.true.
          if (tree(pft)) then
            nind(pft)=0.0
          else
            nind(pft)=1.0   
          endif

          lm_ind(pft)=0.0
          sm_ind(pft)=0.0
          rm_ind(pft)=0.0
          hm_ind(pft)=0.0
          fpc_grid(pft)=0.0
          anpp(pft)=0.0
          do m=1,12
             mnpp(m,pft)=0.0
          enddo


          if(.not.tree(pft)) crownarea(pft)=1.0
          leafon(pft)=.true.
          leafondays(pft)=0.0
          leafoffdays(pft)=0.0

        endif

      enddo


c     SAPLING AND GRASS ESTABLISHMENT

c     Calculate total woody FPC and number of woody PFTs present and
c     able to establish

      fpc_tree_total=0.0
      fpc_tree_new=0.0
      fpc_grass_total=0.0
      fpc_grass_new=0.0
      npft_estab=0
      fpc_total=0.0
      fpc_new=0.0
      ngrass=0

      do pft=1,npft
         if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
              if (estab(pft)) npft_estab=npft_estab+1
            else
              ngrass=ngrass+1
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)
            endif
            fpc_total=fpc_total+fpc_grid(pft)
         endif
      enddo 


      acflux_estab=0.0
c     Prohibit establishment under extreme temperature or water stress.

      do pft=1,npft

      if (aprec.ge.aprec_min_estab.and.
     *  npft_estab.gt.0) then


c       Calculate establishment rate over available space, per tree PFT
c       Maximum establishment rate reduced by shading as tree FPC approaches 1
c       Total establishment rate partitioned equally among regenerating woody
c       PFTs

        estab_rate=estab_max*(1.0-exp(5.0*(fpc_tree_total-1.0)))/
     *    real(npft_estab)

c       Calculate grid-level establishment rate per woody PFT
c       Space available for woody PFT establishment is proportion of grid cell
c       not currently occupied by woody PFTs.

        estab_grid=estab_rate*(1.0-fpc_tree_total)

      else  

        estab_grid=0.0

      endif

        if (present(pft).and.tree(pft).and.estab(pft)) then
          crownarea_max=pftpar(pft,20)


c         Add new saplings to current population

          nind_old=nind(pft)
          nind(pft)=nind_old+estab_grid

c          if (year.gt.1050) 
c     *    print*,'nind(pft) establ',pft,nind(pft),estab_grid


*          if (nind(pft).gt.0.0) then
            lm_ind(pft)=(lm_ind(pft)*nind_old+lm_sapl(pft)*estab_grid)/
     *        nind(pft)
            sm_ind_temp=(sm_ind(pft)*nind_old+sm_sapl(pft)*estab_grid)/
     *        nind(pft)
            hm_ind(pft)=(hm_ind(pft)*nind_old+hm_sapl(pft)*estab_grid)/
     *        nind(pft)
            rm_ind(pft)=(rm_ind(pft)*nind_old+rm_sapl(pft)*estab_grid)/
     *        nind(pft)

*          else
*           lm_ind(pft)=0.0
*           sm_ind(pft)=0.0
*           hm_ind(pft)=0.0
*           rm_ind(pft)=0.0
*          endif


c         Accumulate biomass increment due to sapling establishment

          acflux_estab=acflux_estab+(lm_sapl(pft)+sm_sapl(pft)+
     *      hm_sapl(pft)+rm_sapl(pft))*estab_grid

c         Calculate height, diameter and crown area for new average
c         individual such that the basic allometric relationships (A-C below)
c         are satisfied.

c         (A) (leaf area) = latosa * (sapwood xs area)
c                (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c         (B) (leaf mass) = lmtorm * (root mass)
c         (C) height = allom2 * (stem diameter)**allom3
c                (source?)
c         (D) (crown area) = min (allom1 * (stem diameter)**reinickerp,
c                                 crownarea_max)

c         From (A),
c          (1) sap_xsa = lm_ind * sla / latosa
c          (2) wooddens = (sm_ind + hm_ind) / stemvolume
c          (3) stemvolume = stem_xsa * height
c         From (1), (2) & (3),
c          (4) stem_xsa = (sm_ind + hm_ind) / wooddens / height
c          (5) stem_xsa = pi * (stemdiam**2) / 4
c         From (5),
c          (6) stemdiam = ( 4 * stem_xsa / pi )**0.5
c         From (4) & (6),
c          (7) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens / height /
c                pi )**0.5
c         From (C) & (7),
c          (8) stemdiam = ( 4 * (sm_ind + hm_ind) / wooddens /
c                ( allom2 * stemdiam**allom3 ) / pi )**0.5
c         From (8),
c          (9) stemdiam = ( 4 * (sm_ind + hm_ind ) / wooddens / pi /
c                allom2 )**( 1 / (2 + allom3) )

          stemdiam=(4.0*(sm_ind_temp+hm_ind(pft))/wooddens/pi/allom2)**
     *      (1.0/(2.0+allom3))                 
          height(pft)=allom2*stemdiam**allom3  
          crownarea(pft)=min(crownarea_max,
     *      allom1*stemdiam**reinickerp)       

c         Kirsten recalculate dbh, bt and tau_c for average individual
c          ACHTUNG: Is this needed here? prob. no influence on fire for the next year...
c    KIRSTEN: put in uniform distribution of dbh = 5 classes
            dbh(pft)=stemdiam*100.0 

           do class=0,4
             dbh_class(class,pft)=(2.0*dbh(pft)+
     *                  (class*0.25*2.0*dbh(pft))- 0.125*2.0*dbh(pft))

            bt(class,pft)=(param1(pft)*dbh_class(class,pft)+param2(pft))
              tau_c(class,pft)=2.9*(bt(class,pft)**2.0)

              height_class(class,pft)=2*height(pft)-
     *              (class*0.25*2.0*height(pft))+ 0.125*2.0*height(pft)   
              cl_t(class,pft)=height_class(class,pft)*crown(pft)  
           enddo 

c            bt(pft)=(param1(pft)*dbh(pft)+param2(pft))
c            tau_c(pft)=2.9*(bt(pft)**2.0)


c         Recalculate sapwood mass, transferring excess sapwood to heartwood
c         compartment, if necessary to satisfy Eqn A

          sm_ind(pft)=lm_ind(pft)*height(pft)*wooddens*sla(pft)/latosa
          hm_ind(pft)=hm_ind(pft)+(sm_ind_temp-sm_ind(pft))

c          Update LAI and FPC

           if (crownarea(pft).gt.0.0) then
             lai_ind(pft)=(lm_ind(pft)*sla(pft))/crownarea(pft)
           else
             lai_ind(pft)=0.0
           endif

           fpc_ind=(1.0-exp(-0.5*lai_ind(pft)))
           fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
       endif
      enddo 

      fpc_tree_total=0.0
      fpc_grass_total=0.0
      fpc_total=0.0

      do pft=1,npft
         if (present(pft)) then
            if (tree(pft)) then
              fpc_tree_total=fpc_tree_total+fpc_grid(pft)
            else
              fpc_grass_total=fpc_grass_total+fpc_grid(pft)
            endif
            fpc_total=fpc_total+fpc_grid(pft)
         endif
      enddo 
      if (fpc_tree_total.gt.0.95) then
        do pft=1,npft
          if(tree(pft).and.present(pft)) then
              nind_old=nind(pft)
              nind(pft)=nind(pft)/(fpc_tree_total/0.95)
              fpc_grid(pft)=fpc_grid(pft)/(fpc_tree_total/0.95)
              litter_ag(pft)=litter_ag(pft)+(lm_ind(pft)+sm_ind(pft)+
     *              hm_ind(pft))*(nind_old-nind(pft))

c      Kirsten: fuel class
             fuel_1hr(pft)=fuel_1hr(pft)+(lm_ind(pft)+0.045*sm_ind(pft)+
     *              0.045*hm_ind(pft))*(nind_old-nind(pft))
              fuel_10hr(pft)=fuel_10hr(pft)+((0.075*hm_ind(pft)+
     *                   0.075*sm_ind(pft))*(nind_old-nind(pft)))
              fuel_100hr(pft)=fuel_100hr(pft)+((0.21*hm_ind(pft)+
     *                   0.21*sm_ind(pft))*(nind_old-nind(pft)))
              fuel_1000hr(pft)=fuel_1000hr(pft)+((0.67*hm_ind(pft)+
     *                  0.67*sm_ind(pft))*(nind_old-nind(pft)))

              litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*(nind_old-
     *              nind(pft))
c              print*,'litter_bg, fpc0.95',pft,litter_bg(pft)
          endif

        enddo

        fpc_total=fpc_total-(fpc_tree_total-0.95)
        if (fpc_total.gt.1.0) fpc_total=1.0
        fpc_tree_total=0.95

      endif 	


c     SECTION FOR GRASSES

      do pft=1,npft
	
        if (present(pft).and..not.tree(pft)) then
           if (estab(pft)) then
c          Grasses can establish in non-vegetated areas
              if (ngrass.gt.0) then
                 bare=(1.0-fpc_total)/real(ngrass)
              else
                 bare=0.0
              endif

              bare_max=( (-2.0*crownarea(pft)*alog(1.0-bare-
     *             fpc_grid(pft))/sla(pft))-
     *             lm_ind(pft) )/lm_sapl(pft)

c         Accumulate biomass increment due to grass establishment

             if (bare.le.0.0) then
              litter_bg(pft)=litter_bg(pft)+(bare*rm_sapl(pft))*(-1)
c              print*,'litter_bg, grass',pft,litter_bg(pft)
              litter_ag(pft)=litter_ag(pft)+(bare*lm_sapl(pft))*(-1)

c            KIRSTEN: per fuel class
              fuel_1hr(pft)=fuel_1hr(pft)+(bare*lm_sapl(pft))*(-1)
             else
               if(bare.gt.bare_max) bare=bare_max
               lm_ind(pft)=lm_ind(pft)+bare*lm_sapl(pft)
               rm_ind(pft)=rm_ind(pft)+bare*rm_sapl(pft)
	     endif
              acflux_estab=acflux_estab+bare*(lm_sapl(pft)+
     *             rm_sapl(pft))*crownarea(pft)
           endif


           if (lm_ind(pft).le.0.0) then
              present(pft)=.false.
              do m=1,12
                 mnpp_add(m)=mnpp_add(m)+mnpp(m,pft)
                 anpp_add=anpp_add+mnpp(m,pft)
              enddo
              litter_bg(pft)=litter_bg(pft)+rm_ind(pft)*nind(pft)
c              print*,'not present',pft,litter_bg(pft)
           endif
        endif
      enddo

c     recalculate fpc's

      do pft=1,npft
         if (present(pft)) then
           if (crownarea(pft).gt.0.0) then
             lai_ind(pft)=(lm_ind(pft)*sla(pft))/crownarea(pft)
           else
             lai_ind(pft)=0.0
           endif

           fpc_ind=(1.0-exp(-0.5*lai_ind(pft)))
           fpc_grid(pft)=crownarea(pft)*nind(pft)*fpc_ind
         endif
      enddo
c           if(year.ge.1000) print*,year,(litter_bg(pft),pft=1,npft)
c (fpc_grid(pft),pft=1,npft)

      return
      end

***************************************************************************      c
c         subroutine calculates the area of the pixel(0.5x0.5 or other)degress
c         alberte's method

      subroutine pixelarea(lat,area)
cf2py inten(in) lat
cf2py intent(out) area
      implicit none
      real lat,ilat
      double precision pie,re,dip,equ,
     >     sum,earth
       real area

*       re = 6.378E6
       pie = 4. * ATAN(1.)
       dip = pie/180.

*       area=((111.0*10**3*0.167)**2)*cos(lat*dip)
        area=((111.0*10**3*0.5)**2)*cos(lat*dip)



       return
       end
c -----------------------------------------------------------------------------
c                                 REFERENCES
c -----------------------------------------------------------------------------

c Carslaw, HS & Jaeger JC 1959 Conduction of Heat in Solids, Oxford University
c   Press, London
c Collatz, GJ, Ball, JT, Grivet C & Berry, JA 1991 Physiological and
c   environmental regulation of stomatal conductance, photosynthesis and
c   transpiration: a model that includes a laminar boundary layer. Agricultural
c   and Forest Meteorology 54: 107-136
c Collatz, GJ, Ribas-Carbo, M & Berry, JA 1992 Coupled photosynthesis-stomatal
c   conductance models for leaves of C4 plants. Australian Journal of Plant
c   Physiology 19: 519-538
c Farquhar GD & von Caemmerer 1982 Modelling of photosynthetic response to
c   environmental conditions. In: Lange, OL, Nobel PS, Osmond CB, Ziegler H
c   (eds) Physiological Plant Ecology II: Water Relations and Carbon
c   Assimilation, Vol 12B. Springer, Berlin, pp 549-587.
c Foley J A 1995 An equilibrium model of the terrestrial carbon budget
c   Tellus (1995), 47B, 310-319
c Harper JL 1977 Population Biology of Plants, Academic Press, London
c Haxeltine A & Prentice IC 1996 BIOME3: an equilibrium terrestrial biosphere
c   model based on ecophysiological constraints, resource availability, and
c   competition among plant functional types. Global Biogeochemical Cycles 10:
c    693-709
c Haxeltine A & Prentice IC 1996a A general model for the light-use efficiency
c   of primary production. Functional Ecology 10: 551-561
c Henderson-Sellers, A & Robinson, PJ 1986 Contemporary Climatology. Longman,
c   Essex.
c Jarvis, PG & McNaughton KG 1986 Stomatal control of transpiration: scaling up
c   from leaf to region. Advances in Ecological Research 15: 1-49
c Jury WA, Gardner WR & Gardner WH 1991 Soil Physics 5th ed, John Wiley, NY
c Larcher W 1983 Physiological Plant Ecology, 2nd ed, Springer-Verlag, Berlin
c Lloyd, J & Taylor JA 1994 On the temperature dependence of soil respiration
c   Functional Ecology 8: 315-323
c Monsi, M & Saeki, T 1953 Ueber den Lichtfaktor in den Pflanzengesellschaften
c   und seine Bedeutung fuer die Stoffproduktion. Japanese Journal of Botany
c   14: 22-52
c Monteith, JL & Unsworth, MH 1990 Principles of Environmental Physics, 2nd ed,
c   Arnold, London
c Prentice, IC, Sykes, MT & Cramer W 1993 A simulation model for the transient
c   effects of climate change on forest landscapes. Ecological Modelling 65:
c   51-70.
c Press, WH, Teukolsky, SA, Vetterling, WT & Flannery, BT. 1986. Numerical
c   Recipes in FORTRAN, 2nd ed. Cambridge University Press, Cambridge
c Reich, PB, Walters, MB & Ellsworth, DS, 1997. From tropics to tundra: global
c   convergence in plant functioning. Proceedings of the National Academy of
c   Sciences USA 94: 13730-13734.
c Ryan, GR 1991 Effects of climate change on plant respiration. Ecological
c   applications 1: 157-167
c Shinozaki, K, Yoda, K, Hozumi, K & Kira, T 1964 A quantitative analysis of
c   plant form - the pipe model theory. I. basic analyses. Japanese Journal of
c   Ecology 14: 97-105
c Shinozaki, K, Yoda, K, Hozumi, K & Kira, T 1964 A quantitative analysis of
c   plant form - the pipe model theory. II. further evidence of the theory and
c   its application in forest ecology. Japanese Journal of Ecology 14: 133-139
c Sprugel, DG, Ryan, MG, Brooks, JR, Vogt, KA, Martin, TA, Respiration from the
c   organ level to the stand (in press '93 ---> look up)
c van Duin, RHA 1963 The influence of soil management on the temperature
c   wave near the surface. Tech Bull 29 Inst for Land and Water Management
c   Research, Wageningen, Netherlands
c Waring, RH Schroeder, PE & Oren, R 1982 Application of the pipe model theory
c   to predict canopy leaf area. Canadian Journal of Forest Research 12:
c   556-560

c -----------------------------------------------------------------------------
c//////////////////////////////////////////////////////////////////////////////
c******************************************************************************
c     SUBROUTINE PFTPARAMETERS
c     Assignment of PFT-specific parameters and bioclimatic limits
c     Definition of initial sapling and grass mass structure
c     Calculation of maximum crown area for woody PFTs

      subroutine pftparameters(pftpar,sla,tree,evergreen,
     *  summergreen,raingreen,needle,boreal,lm_sapl,sm_sapl,hm_sapl,
     *  rm_sapl,latosa,allom1,allom2,allom3,
     *  allom4,wooddens,reinickerp)
cf2py intent(out) pftpar, sla, tree, evergreen
cf2py intent(out) summergreen, raingreen, needle, boreal
cf2py intent(out) lm_sapl, sm_sapl, hm_sapl
cf2py intent(out) rm_sapl, latosa, allom1, allom2, allom3, allom4
cf2py intent(out) wooddens, reinickerp

      implicit none

c     PARAMETERS:
      integer npft,npftpar,nsoilpar
        parameter (npft=9,npftpar=50,nsoilpar=7)
      real pi
        parameter (pi=3.14159265)

c     ARGUMENTS:
      real pftpar(1:npft,1:npftpar),sla(1:npft)
      logical tree(1:npft),evergreen(1:npft)
      logical summergreen(1:npft),raingreen(1:npft),needle(1:npft)
      logical boreal(1:npft)
      real lm_sapl(1:npft),sm_sapl(1:npft),hm_sapl(1:npft)
      real rm_sapl(1:npft)
      real latosa,wooddens,reinickerp
      real allom1,allom2,allom3,allom4

c     LOCAL VARIABLES:
      integer n,pft
      real table(1:npft,1:npftpar)
      real lai_sapl       !sapling or initial grass LAI
      real x
      real lmtorm         !non-waterstressed leafmass to rootmass ratio
      real stemdiam       !sapling stem diameter
      real height_sapl    !sapling height
      !From main LPJ code, where the vars are defined
      allom1 = 100.0
      allom2 = 40.0
      allom3 = 0.67
      allom4 = 0.3
      reinickerp = 1.6
      latosa = 6.0e3
      wooddens = 2.0e5


c-----------------------------------------------------------------------------

c     PFT PARAMETERS

c      1  fraction of roots in upper soil layer
c      2  plants with C4 (1) or C3 (0) photosynthetic pathway
c      3  water scalar value at which leaves shed by drought deciduous PFT
c      4  canopy conductance component (gmin, mm/s) not associated with
c         photosynthesis (Haxeltine & Prentice 1996, Table 4)
c      5  maintenance respiration coefficient
c      6  flammability threshold
c      7  maximum foliar N content (mg/g)
c         (Haxeltine & Prentice 1996a, Fig 4)
c      8  fire resistance nindex
c      9  leaf turnover period (years)
c     10  leaf longevity (years)
c     11  sapwood turnover period (sapwood converted to heartwood) (years)
c     12  root turnover period (years)
c     13  leaf C:N mass ratio
c     14  sapwood C:N mass ratio
c     15  root C:N mass ratio
c     16  leaf type: broadleaved (1), needleleaved (2) or grass (3)
c     17  phenology type: evergreen (1), summergreen (2), raingreen (3),
c         any type (4)
c     18  leaf to root ratio under non-water stressed conditions
c     19  summergreen phenology ramp, GDD requirement to grow full leaf canopy
c     20  tree maximum crown area (m2)
c     21  sapling (or grass on initialisation) LAI
c     22  sapling [(heartwood mass) + (sapwood mass)] / (sapwood mass)
c     23  boreal pft (1), non-boreal pft (0)
c     24  low temperature limit for CO2 uptake
c     25  lower range of temperature optimum for photosynthesis
c     26  upper range of temperature optimum for photosynthesis
c     27  high temperature limit for CO2 unptake

c     BIOCLIMATIC LIMITS

c     28 minimum coldest monthly mean temperature
c     29 maximum coldest monthly mean temperature
c     30 minimum growing degree days (at or above 5 deg C)
c     31 upper limit of temperature of the warmest month
c     32 lower limit of growth efficiency (g/m2)
c
c     PARAMETERS ADDED LATER
c
c     33 GDD base
c     34 20-year average min warmest - coldest month temperature range
c     35 emax (mm/day)
c     36 intc (Interception storage parameter, unitless)
c     37 fuel bulk density  (kg/m��)
c     38 emission factor CO2
c	39 emission factor CO
c	40 emission factor CH4
c     41 emission factor VOC
c     42 emission factor TPM
c     43 emission factor NOx
c     44 proportion of tree height as crown
c     45 f parameter for flame length
c     46 g parameter for flame length  Kirsten, 4.11.04: not used at the moment
c     47 param1 for bark thickness
c     48 param2 for bark thickness
c     49 r(ck) crown damage: postfire mortality
c     50 p   crown damage: postfire mortality

c NOTE ON THE INTRODUCTION OF THE LARCH PFT:
c assume that Larch canopy conductance is higher than for other boreals
c assume that phenology ramp is quicker than for broadleaved trees
c also that ramp is applied to a GDD base of 2, not 5
      data ((table(pft,n),n=1,8),pft=1,npft) /

c     ---------------------------------------------------------------------
c          1      2      3      4      5      6      7      8          PFT
c     ---------------------------------------------------------------------

     *  0.85,   0.0,  0.00,   0.5,  0.20,  0.20, 100.0,  0.60,        !  1
     *  0.50,   0.0,  0.10,   0.5,  0.20,  0.30, 100.0,  0.70,        !  2
     *  0.60,   0.0,  0.00,   0.3,  1.20,  0.30, 100.0,  0.12,        !  3
     *  0.50,   0.0,  0.10,   0.3,  1.20,  0.30, 100.0,  0.50,        !  4
     *  0.70,   0.0,  0.00,   0.5,  1.20,  0.30, 120.0,  0.12,        !  5
     *  0.90,   0.0,  0.00,   0.3,  1.20,  0.35, 100.0,  0.12,        !  6
c     *  0.90,   0.0,  0.00,   0.3,  1.30,  0.35, 100.0,  0.12,        !  7
     *  0.90,   0.0,  0.00,   0.5,  1.20,  0.35, 100.0,  0.12,        !  8
     *  0.90,   0.0,  0.20,   0.5,  1.30,  0.20, 100.0,  0.01,        !  9
     *  0.70,   1.0,  0.20,   0.5,  0.70,  0.20, 100.0,  0.01/        ! 10

      data ((table(pft,n),n=9,17),pft=1,npft) /

c     ---------------------------------------------------------------------
c          9     10     11     12     13     14     15     16     17   PFT
c     ---------------------------------------------------------------------

     *   2.0,  2.00,  20.0,   2.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  1
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   3.0, !  2
     *   4.0,  4.00,  20.0,   4.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  3
     *   1.0,  1.00,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   1.0, !  4
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, !  5
     *   4.0,  4.00,  20.0,   4.0,  29.0, 330.0,  29.0,   2.0,   1.0, !  6
c     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   2.0,   2.0, !  7
     *   1.0,  0.50,  20.0,   1.0,  29.0, 330.0,  29.0,   1.0,   2.0, !  8
     *   1.0,  0.50,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0, !  9
     *   1.0,  0.50,   1.0,   2.0,  29.0,   0.0,  29.0,   3.0,   4.0/ ! 10

      data ((table(pft,n),n=18,23),pft=1,npft) /

c     ---------------------------------------------------
c           18      19     20      21     22     23     pft
c     ---------------------------------------------------

     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  1
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  2
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  3
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    0.0,    !  4
     *     1.0,   200.0,  15.0,  1.500,  1.2,    0.0,    !  5
     *     1.0,  1000.0,  15.0,  1.500,  1.2,    1.0,    !  6
c     *     1.0,   100.0,  15.0,  1.500,  1.2,    1.0,    !  7
     *     1.0,   200.0,  15.0,  1.500,  1.2,    1.0,    !  8
     *    0.65,   100.0,   0.0,  0.001,  1.2,    1.0,    !  9
     *    0.75,   100.0,   0.0,  0.001,  1.2,    0.0/    ! 10

      data ((table(pft,n),n=24,27),pft=1,npft) /
c     -------------------------------------
c          24     25     26      27    PFT
c     -------------------------------------
     *    2.0,   25.0,  30.0,   55.0, !  1
     *    2.0,   25.0,  30.0,   55.0, !  2
     *   -4.0,   20.0,  30.0,   42.0, !  3
     *   -4.0,   20.0,  30.0,   42.0, !  4
     *   -4.0,   20.0,  25.0,   38.0, !  5
     *   -4.0,   15.0,  25.0,   38.0, !  6
c     *   -4.0,   15.0,  25.0,   38.0, !  7
     *   -4.0,   15.0,  25.0,   38.0, !  8
     *   -4.0,   10.0,  30.0,   45.0, !  9
     *    6.0,   20.0,  45.0,   55.0/ ! 10


      data ((table(pft,n),n=28,32),pft=1,npft) /

c     ------------------------------------------------------------
c          28       29       30       31       32       PFT
c     ------------------------------------------------------------
     *    15.5,  1000.0,    0.0,    1000.0,    0.0,     !  1
     *    13.5,  1000.0,    0.0,    1000.0,    0.0,     !  2
     *    -2.0,    15.0,  900.0,    1000.0,    0.0,     !  3
     *     1.0,    15.5, 1200.0,    1000.0,    0.0,     !  4
     *   -17.0,    15.5, 1200.0,    1000.0,    0.0,     !  5
     *   -32.5,    -4.0,  600.0,    1000.0,    0.0,     !  6
c     * -1000.0,    -2.0,  350.0,      23.0,    0.0,     !  7
     * -1000.0,    -4.0,  350.0,    1000.0,    0.0,     !  8
     *   -32.0,  1000.0,    0.0,    1000.0,    0.0,     !  9
     *    10.0,  1000.0,    0.0,    1000.0,    0.0/     ! 10


      data ((table(pft,n),n=33,37),pft=1,npft) /

c     ---------------------------------------------
c          33       34        35    36     37      PFT
c     ---------------------------------------------
     *     5.0,  -1000.0,     7.0,  0.02, 14.0,   !  1
     *     5.0,  -1000.0,     7.0,  0.02, 12.0,   !  2
     *     5.0,  -1000.0,     5.0,  0.06, 23.0,   !  3
     *     5.0,  -1000.0,     5.0,  0.02, 28.0,   !  4
     *     5.0,  -1000.0,     5.0,  0.02, 22.0,   !  5
     *     5.0,  -1000.0,     5.0,  0.06, 18.0,   !  6
c     *     2.0,   43.0,       5.0,  0.06, 16.0,    !  7
     *     5.0,  -1000.0,     6.0,  0.06, 16.0,   !  8
     *     5.0,  -1000.0,     5.0,  0.01,  4.0,   !  9
     *     5.0,  -1000.0,     7.0,  0.01,  2.0/   ! 10


      data ((table(pft,n),n=38,45),pft=1,npft) /

c     --------------------------------------------------------------------
c          38       39      40      41      42   43       44      45    PFT
c     --------------------------------------------------------------------
     *    1580.0, 103.0,   6.80,   8.10,   8.50, 1.999, 0.3334, 0.160, !  1
     *    1664.0,  63.0,   2.20,   3.40,   8.50, 2.540, 0.10,   0.351, !  2
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, ! 3
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.070, !  4
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, ! 5
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.6667, 0.094, !  6
c     *     2.0,   43.0,       5.0,  0.06,   !  7
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.3334, 0.094, !  8
     *    1568.0, 106.0,   4.80,   5.70,  17.60, 3.240, 0.0,    0.0, !  9
     *    1664.0,  63.0,   2.20,   3.40,   8.50, 2.540, 0.0,    0.0/! 10
c----------------------------------------------------------------------------

      data ((table(pft,n),n=46,npftpar),pft=1,npft) /

c     -------------------------------------------------
c         46       47      48     49   50     PFT
c     -------------------------------------------------
     *    0.630, 0.0301, 0.0281, 0.95, 3.00,  !  1
     *    0.460, 0.1085, 0.2120, 0.05, 3.00,  !  2
     *    0.667, 0.0670, 0.5590, 0.95, 3.75,  !  3
     *    0.560, 0.0451, 0.1412, 0.95, 3.00,  !  4
     *    0.667, 0.0347, 0.1086, 0.95, 3.00,  !  5
     *    0.667, 0.0292, 0.2632, 0.95, 3.00,  !  6
c     *     2.0,   43.0,       5.0,  0.06,   !  7
     *    0.667, 0.0347, 0.1086, 0.95, 3.00,  !  8
     *    0.0,      0.0,    0.0, 0.00, 0.00,  !  9
     *    0.0,      0.0,    0.0, 0.00, 0.00 / ! 10
c----------------------------------------------------------------------------

      do pft=1,npft

c       Transfer parameter values to array pftpar

        do n=1,npftpar
          pftpar(pft,n)=table(pft,n)
        enddo

c       Assign leaf and phenology logicals

        if (pftpar(pft,16).le.2.0) then
          tree(pft)=.true.
          if (pftpar(pft,16).eq.2.0) then
            needle(pft)=.true.
          else
            needle(pft)=.false.
          endif
        else
          tree(pft)=.false.
          needle(pft)=.false.
        endif

        if (pftpar(pft,17).eq.1.0) then
          evergreen(pft)=.true.
          summergreen(pft)=.false.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.2.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.true.
          raingreen(pft)=.false.
        elseif (pftpar(pft,17).eq.3.0) then
          evergreen(pft)=.false.
          summergreen(pft)=.false.
          raingreen(pft)=.true.
        else
          evergreen(pft)=.true.
          summergreen(pft)=.true.
          raingreen(pft)=.true.
        endif

        if (pftpar(pft,23).eq.1.0) then
          boreal(pft)=.true.
        else
          boreal(pft)=.false.
        endif

c       Calculate specific leaf area (SLA) for each PFT from leaf longevity
c       Include conversion (multiplier of 2.0) from m2/g(dry wt) to m2/gC
c       Equation based on Reich et al 1997, Fig 1f:

c       SLA = 2e-4 * exp(6.15 - 0.46 ln (leaf_longevity * 12))

c       SLA in m2/gC, leaf_longevity in years

        sla(pft)=2.0e-4*exp(6.15-0.46*log(pftpar(pft,10)*12.0))

c       Define initial mass structure

        lai_sapl=pftpar(pft,21)

        if (tree(pft)) then  !woody PFTs

c         Calculate leafmass for a sapling individual
c          (1) lai = leafmass * sla / (crown area)
c          (2) (leaf area) = latosa * (sapwood xs area)
c                 (Pipe Model, Shinozaki et al. 1964a,b; Waring et al 1982)
c          (3) (crown area) = allom1 * (stem diameter) ** reinickerp
c                 (Reinickes theory)
c         From (1),
c          (4) leafmass = lai * (crown area) / sla
c         From (1) & (3),
c          (5) leafmass = lai * allom1 * (stem diameter)**reinickerp / sla
c         From (2),
c          (6) leafmass = latosa * (sapwood xs area) / sla
c          (7) (sapwood xs area) = pi * (sapwood diameter)**2 / 4
c         From (6) and (7),
c          (8) leafmass = latosa * pi * (sapwood diameter)**2 / 4 / sla
c         From (8),
c          (9) (sapwood diameter) = [ 4 * leafmass * sla / pi / latosa ]**0.5
c         (10) (stem diameter) = (sapwood diameter) + (heartwood diameter)
c         Define x,
c         (11) x = [ (sapwood diameter)+(heartwood diameter) ] /
c                  (sapwood diameter)
c         From (10) & (11),
c         (12) (stem diameter) = x * (sapwood diameter)
c         From (5), (9) & (12),
c         (13) leafmass = lai * allom1 * x**reinickerp *
c                       (4*leafmass*sla/pi/latosa)**(reinickerp*0.5) / sla
c         From (13),
c         (14) leafmass = [ lai * allom1 * x**reinickerp *
c              (4*sla/pi/latosa)**(reinickerp*0.5) / sla ]**(2/(2-reinickerp))

          x=pftpar(pft,22)

          lm_sapl(pft)=(lai_sapl*allom1*x**reinickerp*
     *      (4.0*sla(pft)/pi/latosa)**(reinickerp*0.5) / sla(pft))**
     *      (2.0/(2.0-reinickerp))  !eqn 14

c         Calculate sapling stem diameter
c         From (9) & (12),
c         (15) (stem diameter) = x * [ 4 * leafmass * sla / pi / latosa ]**0.5

          stemdiam=x*(4.0*lm_sapl(pft)*sla(pft)/pi/latosa)**0.5  !Eqn 15

c         Calculate sapling height
c         (16) height = allom2 * (stem diameter)**allom3 (source?)

          height_sapl=allom2*stemdiam**allom3   !Eqn 16

c         Calculate sapling sapwood mass
c         (17) (sapwood volume) = height * (sapwood xs area)
c         (18) (sapwood xs area) = leafmass * sla / latosa
c         From (17) & (18),
c         (19) (sapwood volume) = height * leafmass * sla / latosa
c         (20) (sapwood mass) = (wood density) * (sapwood volume)
c         From (19) & (20),
c         (21) (sapwood mass) = (wood density) * height * leafmass * sla /
c                latosa

          sm_sapl(pft)=wooddens*height_sapl*lm_sapl(pft)*sla(pft)/
     *      latosa   !Eqn 21

c         Calculate sapling heartwood mass
c         From (11),
c         (22) (heartwood mass) = (x-1) * (sapwood mass)

          hm_sapl(pft)=(x-1.0)*sm_sapl(pft)  !Eqn 22

        else ! grass PFT

          lm_sapl(pft)=lai_sapl/sla(pft)

        endif

c       Calculate sapling or initial grass rootmass
c       (23) lmtorm = (leafmass) / (rootmass)

        lmtorm=pftpar(pft,18)
        rm_sapl(pft)=(1.0/lmtorm)*lm_sapl(pft)  !From Eqn 23

      enddo ! pft loop

      return
      end