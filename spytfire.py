!    -*- f90 -*-
! Note: the context of this file is case sensitive.

python module spytfire ! in 
    interface  ! in :spytfire
        subroutine fire(year,pftpar,dtemp,dtemp_min,dtemp_max,dprec,dwindsp,lightn,dphen,litter_ag,litter_bg,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,acflux_fire,mcflux_fire,afire_frac,lm_ind,rm_ind,sm_ind,hm_ind,nind,dw1,present,tree,lat,mw1,fpc_grid,popden,a_nd_array,height,height_class,dbh,tau_c,cl_t,num_fire,annum_fire,area_burnt,an_areafires,mfdi,an_fdi,an_fseason,mcflux_trace,acflux_trace,m_fc_crown,an_fc_crown,m_i_surface,an_i_surface) ! in :spytfire:spitfire.f
            integer intent(in) :: year
            real dimension(9,50),intent(in) :: pftpar
            real dimension(365),intent(in) :: dtemp
            real dimension(365),intent(in) :: dtemp_min
            real dimension(365),intent(in) :: dtemp_max
            real dimension(365),intent(in) :: dprec
            real dimension(365),intent(in) :: dwindsp
            real dimension(12),intent(in) :: lightn
            real dimension(365,9),intent(in) :: dphen
            real dimension(9),intent(inout) :: litter_ag
            real dimension(9),intent(inout) :: litter_bg
            real dimension(9),intent(inout) :: fuel_1hr
            real dimension(9),intent(inout) :: fuel_10hr
            real dimension(9),intent(inout) :: fuel_100hr
            real dimension(9),intent(inout) :: fuel_1000hr
            real intent(out) :: acflux_fire
            real dimension(12) :: mcflux_fire
            real intent(out) :: afire_frac
            real dimension(9),intent(inout) :: lm_ind
            real dimension(9),intent(inout) :: rm_ind
            real dimension(9),intent(inout) :: sm_ind
            real dimension(9),intent(inout) :: hm_ind
            real dimension(9),intent(inout) :: nind
            real dimension(365),intent(in) :: dw1
            logical dimension(9),intent(in) :: present
            logical dimension(9),intent(in) :: tree
            real intent(in) :: lat
            real dimension(12),intent(in) :: mw1
            real dimension(9),intent(inout) :: fpc_grid
            real intent(in) :: popden
            real dimension(12) :: a_nd_array
            real dimension(9),intent(in) :: height
            real dimension(5,9),intent(in) :: height_class
            real dimension(9),intent(in) :: dbh
            real dimension(5,9),intent(inout) :: tau_c
            real dimension(5,9),intent(in) :: cl_t
            real dimension(12),intent(out) :: num_fire
            real intent(out) :: annum_fire
            real dimension(12),intent(out) :: area_burnt
            real intent(out) :: an_areafires
            real dimension(12),intent(out) :: mfdi
            real intent(out) :: an_fdi
            real intent(out) :: an_fseason
            real dimension(12,6),intent(out) :: mcflux_trace
            real dimension(6),intent(out) :: acflux_trace
            real dimension(12),intent(out) :: m_fc_crown
            real intent(out) :: an_fc_crown
            real dimension(12),intent(out) :: m_i_surface
            real intent(out) :: an_i_surface
        end subroutine fire
        subroutine fire_danger_index(d_fdi,dlm,dlm_lg,dtemp_min,dtemp_max,dprec,d,moistfactor,fuel_1hr_total,fuel_10hr_total,fuel_100hr_total,dead_fuel,char_moistfactor,ratio_dead_fuel,ratio_live_fuel,dlm_1hr,dlm_10hr,dlm_100hr,dlm_1000hr,year) ! in :spytfire:spitfire.f
            real dimension(365) :: d_fdi
            real dimension(365) :: dlm
            real dimension(365) :: dlm_lg
            real dimension(365) :: dtemp_min
            real dimension(365) :: dtemp_max
            real dimension(365) :: dprec
            integer :: d
            real :: moistfactor
            real :: fuel_1hr_total
            real :: fuel_10hr_total
            real :: fuel_100hr_total
            real :: dead_fuel
            real :: char_moistfactor
            real :: ratio_dead_fuel
            real :: ratio_live_fuel
            real dimension(365) :: dlm_1hr
            real dimension(365) :: dlm_10hr
            real dimension(365) :: dlm_100hr
            real dimension(365) :: dlm_1000hr
            integer :: year
        end subroutine fire_danger_index
        function human_ign(popden,a_nd,year) ! in :spytfire:spitfire.f
            real :: popden
            real :: a_nd
            integer :: year
            real :: human_ign
        end function human_ign
        subroutine rate_of_spread(u_front,base_wind,dens_fuel_ave,sigma,dlm,d,net_fuel,moistfactor,h,char_dens_fuel_ave,char_sigma,char_net_fuel,char_moistfactor,gamma) ! in :spytfire:spitfire.f
            real :: u_front
            real :: base_wind
            real :: dens_fuel_ave
            real :: sigma
            real dimension(365) :: dlm
            integer :: d
            real :: net_fuel
            real :: moistfactor
            real :: h
            real :: char_dens_fuel_ave
            real :: char_sigma
            real :: char_net_fuel
            real :: char_moistfactor
            real :: gamma
        end subroutine rate_of_spread
        subroutine fuel_consumption(npft,present,fuel_consum,fire_frac,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,pot_fc_lg,tree,fuel_1hr_total,moistfactor,dlm,d,miner_tot,fc_1hr,fc_lg,fc_10hr,fc_100hr,fc_1000hr,cf,char_moistfactor,dlm_lg,dlm_1hr,dlm_10hr,dlm_100hr,dlm_1000hr,moistfactor_livegrass,moistfactor_1hr,moistfactor_10hr,moistfactor_100hr,moistfactor_1000hr) ! in :spytfire:spitfire.f
            integer optional,check(len(present)>=npft),depend(present) :: npft=len(present)
            logical dimension(npft) :: present
            real :: fuel_consum
            real :: fire_frac
            real dimension(npft),depend(npft) :: fuel_1hr
            real dimension(npft),depend(npft) :: fuel_10hr
            real dimension(npft),depend(npft) :: fuel_100hr
            real dimension(npft),depend(npft) :: fuel_1000hr
            real :: livegrass
            real dimension(npft),depend(npft) :: pot_fc_lg
            logical dimension(npft),depend(npft) :: tree
            real :: fuel_1hr_total
            real :: moistfactor
            real dimension(365) :: dlm
            integer :: d
            real :: miner_tot
            real dimension(npft),depend(npft) :: fc_1hr
            real dimension(npft),depend(npft) :: fc_lg
            real dimension(npft),depend(npft) :: fc_10hr
            real dimension(npft),depend(npft) :: fc_100hr
            real dimension(npft),depend(npft) :: fc_1000hr
            real dimension(npft),depend(npft) :: cf
            real :: char_moistfactor
            real dimension(365) :: dlm_lg
            real dimension(365) :: dlm_1hr
            real dimension(365) :: dlm_10hr
            real dimension(365) :: dlm_100hr
            real dimension(365) :: dlm_1000hr
            real :: moistfactor_livegrass
            real :: moistfactor_1hr
            real :: moistfactor_10hr
            real :: moistfactor_100hr
            real :: moistfactor_1000hr
        end subroutine fuel_consumption
        subroutine fuel_consumption_netmoist(npft,present,fuel_consum,fire_frac,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,livegrass,pot_fc_lg,tree,fuel_1hr_total,moistfactor,dlm,d,miner_tot,fc_1hr,fc_lg,fc_10hr,fc_100hr,fc_1000hr,tau_l,char_moistfactor,dlm_lg,dlm_1hr,dlm_10hr,dlm_100hr,dlm_1000hr,moistfactor_livegrass,moistfactor_1hr,moistfactor_10hr,moistfactor_100hr,moistfactor_1000hr) ! in :spytfire:spitfire.f
            integer optional,check(len(present)>=npft),depend(present) :: npft=len(present)
            logical dimension(npft) :: present
            real :: fuel_consum
            real :: fire_frac
            real dimension(npft),depend(npft) :: fuel_1hr
            real dimension(npft),depend(npft) :: fuel_10hr
            real dimension(npft),depend(npft) :: fuel_100hr
            real dimension(npft),depend(npft) :: fuel_1000hr
            real :: livegrass
            real dimension(npft),depend(npft) :: pot_fc_lg
            logical dimension(npft),depend(npft) :: tree
            real :: fuel_1hr_total
            real :: moistfactor
            real dimension(365) :: dlm
            integer :: d
            real :: miner_tot
            real dimension(npft),depend(npft) :: fc_1hr
            real dimension(npft),depend(npft) :: fc_lg
            real dimension(npft),depend(npft) :: fc_10hr
            real dimension(npft),depend(npft) :: fc_100hr
            real dimension(npft),depend(npft) :: fc_1000hr
            real dimension(npft),depend(npft) :: tau_l
            real :: char_moistfactor
            real dimension(365) :: dlm_lg
            real dimension(365) :: dlm_1hr
            real dimension(365) :: dlm_10hr
            real dimension(365) :: dlm_100hr
            real dimension(365) :: dlm_1000hr
            real :: moistfactor_livegrass
            real :: moistfactor_1hr
            real :: moistfactor_10hr
            real :: moistfactor_100hr
            real :: moistfactor_1000hr
        end subroutine fuel_consumption_netmoist
        subroutine establishment(pftpar,present,survive,estab,nind,lm_ind,sm_ind,rm_ind,hm_ind,lm_sapl,sm_sapl,rm_sapl,hm_sapl,crownarea,fpc_grid,lai_ind,height,dbh,tau_c,cl_t,sla,wooddens,latosa,mprec,reinickerp,litter_ag,litter_bg,fuel_1hr,fuel_10hr,fuel_100hr,fuel_1000hr,tree,allom1,allom2,allom3,acflux_estab,leafondays,leafoffdays,leafon,mnpp,anpp,mnpp_add,anpp_add,year) ! in :spytfire:spitfire.f
            real dimension(9,50) :: pftpar
            logical dimension(9) :: present
            logical dimension(9) :: survive
            logical dimension(9) :: estab
            real dimension(9) :: nind
            real dimension(9) :: lm_ind
            real dimension(9) :: sm_ind
            real dimension(9) :: rm_ind
            real dimension(9) :: hm_ind
            real dimension(9) :: lm_sapl
            real dimension(9) :: sm_sapl
            real dimension(9) :: rm_sapl
            real dimension(9) :: hm_sapl
            real dimension(9) :: crownarea
            real dimension(9) :: fpc_grid
            real dimension(9) :: lai_ind
            real dimension(9) :: height
            real dimension(9) :: dbh
            real dimension(5,9) :: tau_c
            real dimension(5,9) :: cl_t
            real dimension(9) :: sla
            real :: wooddens
            real :: latosa
            real dimension(12) :: mprec
            real :: reinickerp
            real dimension(9) :: litter_ag
            real dimension(9) :: litter_bg
            real dimension(9) :: fuel_1hr
            real dimension(9) :: fuel_10hr
            real dimension(9) :: fuel_100hr
            real dimension(9) :: fuel_1000hr
            logical dimension(9) :: tree
            real :: allom1
            real :: allom2
            real :: allom3
            real :: acflux_estab
            integer dimension(9) :: leafondays
            integer dimension(9) :: leafoffdays
            logical dimension(9) :: leafon
            real dimension(12,9) :: mnpp
            real dimension(9) :: anpp
            real dimension(12) :: mnpp_add
            real :: anpp_add
            integer :: year
        end subroutine establishment
        subroutine pixelarea(lat,area) ! in :spytfire:spitfire.f
            real :: lat
            real :: area
        end subroutine pixelarea
    end interface 
end python module spytfire

! This file was auto-generated with f2py (version:2_5222).
! See http://cens.ioc.ee/projects/f2py2e/
