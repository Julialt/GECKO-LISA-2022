PROGRAM postproc
  USE IO
  USE COMPUTE
  USE UPDATE_SPEC
  USE NCUTILS
  USE MACHINE_LEARNING, ONLY: OUTPUT_ENVIRON, OUTPUT_BINS
  USE CONC, ONLY: numsp
  implicit none

  integer    :: ncid
  integer    :: ibox,tskip

  ! get_precursor code, seed mass, seed molw
  call get_userinfo()

!---------------------------------------
! set skip_time  = (precursor lifetime)/e
! NB: if a SMALLER skip_time is requested by user,
! that will take precedence

  ! save input skip_time, override it below for calculation
  tskip = skip_time 
  !skip_time = 1
  boxname = "bottom"  
  !pressure = 1013.25

  if (input_type == "binary") then

    filename = "dictionary"
    call read_dict()
    filename = "outdat.ppf"
    call read_ppf()
    
    filename = "pvap.dat"
    if (flag_pvap .or. flag_cstar .or. flag_bubble) call read_pvap()
    filename = "henry.dat"
    if (flag_henry) call read_henry()
    
  else if (input_type == "netcdf") then
    filename = "outdat.nc"
    CALL open_ncfile_readonly(filename,ncid)
    call read_ncdf_dict(ncid)
    call read_ncdf_ppf(ncid, 1)
    if (flag_pvap .or. flag_cstar .or. flag_bubble .or. flag_ml) call read_ncdf_pvap(ncid)
    if (flag_henry .or. flag_ml) call read_ncdf_henry(ncid)
    if (flag_ml .OR. flag_OHR .OR. flag_NO3R .OR. flag_O3R) call read_ncdf_kreac(ncid)
    if (flag_rates) call read_ncdf_reacrate(ncid)
  else
    write(6,*) "input_type not recognised: ", input_type
    stop
  endif
 
!   sort species array to accelerate things later
  write(6,*)  'sorting species array...'
  call sort_species()  
  write(6,*)  'sorted'
  
  ! update info on all species in ppf filename
  call update_all_species()

  ! calculate times when precursor has been consumed e-times
  call calc_precursor_lifetime()

  skip_time = MIN(tskip,MAX( 1, NINT(REAL(efold_indices(2)-efold_indices(1))/(EXP(1.)))))

! endif

 PRINT*,"skip time =", skip_time
!---------------------------------------
! re-read, retaining only the lifetime values

  DO ibox = 1,nbox  

  SELECT CASE(ibox)
    CASE(1)
! bottombox
      boxname = "bottom"
      ! already read input, above

    CASE(2)
      boxname = "top"

      if (input_type == "binary") then
        filename = "outdat.pff"
        call read_ppf()
      else if (input_type == "netcdf") then
        filename = "outdat.nc"
        call read_ncdf_ppf(ncid, ibox)
      else
        write(6,*) "input_type not recognised: ", input_type
        stop
      endif

    ! sort & update info on all species in ppf filename
      call sort_species()  
      call update_all_species()

    END SELECT

  ! (re)calculate times when precursor has been consumed e-times
  call calc_precursor_lifetime()
  
  ! print environmental parameters (rh, Temp, pbl_height)
  call calc_environmental_param()

  ! ! O/C, H/C, N/C ratios
  if (flag_atomratios)  then
    call calc_atom_ratios()
    if (flag_potential_atomratios) then
      call calc_potential_frag_dimer_atom_ratio()
    endif
  endif

  ! ! AMS factors
  if (flag_amsfactors) call calc_ams_factors()

  ! ! count functions  and ROF/C and output their time evolution
  if (flag_functions) call count_functions()

  ! ! calculate carbon chain distribution
  if (flag_carbonchain) call carbon_chain()

  ! print selected species (default units = molec/cc)
  if (flag_selected) call calc_selected_species()

  ! ! calculate distributions depending on phases
  !! carbon_distribution, mass distribution
  if (flag_phasedist) call phase_distribution()

  ! ! number of species explaining X% of the mass
  if (flag_contributingspecs) call calc_contributing_species()

  ! ! calculate soa_yield
  if (flag_soayield) call calc_soa_yield()

  ! ! find most important reaction rates for input species
  if (flag_rates) call calc_species_rates(ncid)

  ! ! find most important species (default = molec/cc, optional = ppb 
  ! ! Alternative units = ppbC (only)
  if (flag_topspec) then
    if (flag_ppbC) then
      call calc_top_species_ppbC()
    else
      call calc_top_species()
! JTML test !
      call calc_top_sp_fgrp()
    endif
  endif

  ! ! most abundant N-containing species (default = molec/cc, optional = ppb 
  if (flag_topNspc) then
      call calc_top_Nspc()
  endif

  ! find most important elemental formulas
  if (flag_chon) call calc_chon()

    ! ! find most important CHON
    ! already covered in calc_chon
    ! call calc_top_chon()

  ! ! pvap distribution
  !! if (flag_pvap) call calc_pvap_distribution()
  if (flag_pvap) call calc_pvap_distribution_lifetime()
  if (flag_cstar) call calc_cstar_distribution_lifetime()

  ! ! Henry coefficients distribution
  if (flag_henry) call calc_henry_distribution_lifetime()
!    call calc_henry_distribution()

  ! ! distribution of concentration in (pvap/cstar)/OSC space   
  if (flag_bubble)then
    if(flag_cstar) then
      call calc_bubble_cstar_lifetime()
    else
      call calc_bubble_pvap_lifetime()
    endif
  endif
!    call calc_bubble()

  ! ! simulated mass spectrum at different stages
  if (flag_massspectrum) call calc_mass_spectrum_lifetime()
!    call calc_mass_spectrum()

  if (flag_elementscontrib) call calc_elements_mass_contribution()

! calculates glass transition temperature
  if (flag_phasestate) call calc_phasestate()

  if (flag_entropy) call calc_entropy()
  
  if (flag_nitrates) call calc_nitrates()

  if (flag_chochonfreq) call calc_chochonfreq()
  
  if (flag_dyn_filter) call calc_dyn_filters()

  if (flag_dbeai) then
    call calc_dbe()
    call calc_ai()
  endif
  
  if (flag_ohexposure) call calc_oh_exposure()
  if (flag_no3exposure) call calc_no3_exposure()

  if (flag_OHR) then
     call calc_tot_reactivity("OH ")
     call calc_prod_reactivity("OH ")
     call calc_spec_reactivity("OH ")
  endif
  if (flag_NO3R) then
     call calc_tot_reactivity("NO3")
     call calc_prod_reactivity("NO3")
     call calc_spec_reactivity("NO3")
  endif
  if (flag_O3R) then
     call calc_tot_reactivity("O3 ")
     call calc_prod_reactivity("O3 ")
     call calc_spec_reactivity("O3 ")
  endif

  if (flag_ml) then
     call output_environ()
     call output_bins()
  endif 

  ! ! OSc vs nC
! SUBROUTINE COMMENTED OUT
  !if (flag_oscnc) call calc_osc_nc()

  ! ! adjacency of functional groups
! NO SUBROUTINE EXISTS
  ! call functions_adjacency()

  ! ! find most important CHON
! NO SUBROUTINE EXISTS
  ! call calc_top_chon()

  ! ! find most important cofunctionalities (NN, KD, etc..)
! NO SUBROUTINE EXISTS
  ! call cofunctionalities()

  ! wsy
  if (flag_vbs_param) call calc_vbs_parameterization()

  ENDDO ! ibox

END PROGRAM
!======================================================================
