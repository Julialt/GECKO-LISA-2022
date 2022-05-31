PROGRAM postproc
  USE IO
  USE COMPUTE
  USE UPDATE_SPEC
  USE netcdf

  implicit none

  integer ncid

  ! get_precursor code, seed mass, seed molw
  call get_userinfo()

  filename = "outdat.nc"

  CALL open_ncfile_readonly(filename,ncid)

  call read_ncdf_dict(ncid)

  call read_ncdf_ppf(ncid)

  if (flag_pvap) call read_ncdf_pvap(ncid)

  if (flag_henry) call read_ncdf_henry(ncid)


! bottombox
  boxname = "bottom"

  ! update info on all species in ppf filename
  call update_all_species()

  ! calculate times when precursor has been consumed e-times
  call calc_precursor_lifetime()

  ! print selected species
  if (flag_selected) call calc_selected_species()

  ! find most important elemental formulas
  if (flag_chon) call calc_chon()


  ! ! calculate distributions depending on phases
  !! carbon_distribution, mass distribution
  if (flag_phasedist) call phase_distribution()

  ! ! number of species explaining X% of the mass
  if (flag_contributingspecs) call calc_contributing_species()


  ! ! calculate soa_yield
  if (flag_soayield) call calc_soa_yield()

  ! ! count functions  and ROF/C and output their time evolution
  if (flag_functions) call count_functions()

  ! ! calculate carbon chain distribution
  if (flag_carbonchain) call carbon_chain()

  ! ! find most important species
  if (flag_topspec) call calc_top_species()

  ! ! find most important CHON
  ! call calc_top_chon()


  ! ! find most important cofunctionalities (NN, KD, etc..)
  ! call cofunctionalities()

  ! ! O/C, H/C, N/C ratios
  if (flag_atomratios) call calc_atom_ratios()

  ! ! pvap distribution
  if (flag_pvap) call calc_pvap_distribution()
  !! if (flag_pvap) call calc_vbs_parameterization()

  if (flag_henry) call calc_henry_distribution()

  if (flag_elementscontrib) call calc_elements_mass_contribution()

  if (flag_phasestate) call calc_phasestate()

  if (flag_entropy) call calc_entropy()

  if (flag_chochonfreq) call calc_chochonfreq()

  if (flag__dbeai) then
    call calc_dbe()
    call calc_ai()
  endif
  ! !! simulated mass spectru at different stages
  if (flag_massspectrum) call calc_mass_spectrum()
  ! ! OSc vs nC
  !if (flag_oscnc) call calc_osc_nc()

  ! ! adjacency of functional groups
  ! call functions_adjacency()


  ! wsy
  if (flag_vbs_param) call calc_vbs_parameterization()



  if (nbox == 2) then

! topbox
  boxname = "top"

    call read_ncdf_ppf(ncid)

  ! update info on all species in pff filename
    call update_all_species()


    ! print selected species
    if (flag_selected) call calc_selected_species()

    ! find most important elemental formulas
    if (flag_chon) call calc_chon()

    ! ! calculate distributions depending on phases
    !! carbon_distribution, mass distribution
    if (flag_phasedist) call phase_distribution()

    ! ! number of species explaining X% of the mass
    if (flag_contributingspecs) call calc_contributing_species()

    ! ! count functions  and ROF/C and output their time evolution
    if (flag_functions) call count_functions()

    ! ! calculate carbon chain distribution
    if (flag_carbonchain) call carbon_chain()

    ! ! find most important species
    if (flag_topspec) call calc_top_species()

    ! ! find most important CHON
    ! call calc_top_chon()


    ! ! find most important cofunctionalities (NN, KD, etc..)
    ! call cofunctionalities()

    ! ! O/C, H/C, N/C ratios
    if (flag_atomratios) call calc_atom_ratios()

    ! ! pvap distribution at different stages
    if (flag_pvap) call calc_pvap_distribution()

    if (flag_henry) call calc_henry_distribution()

    if (flag_elementscontrib) call calc_elements_mass_contribution()

    if (flag_phasestate) call calc_phasestate()

    if (flag_entropy) call calc_entropy()

    if (flag_chochonfreq) call calc_chochonfreq()

    if (flag__dbeai) then
      call calc_dbe()
      call calc_ai()
    endif
    ! !! simulated mass spectru at different stages
    if (flag_massspectrum) call calc_mass_spectrum()

    !if (flag_oscnc) call calc_osc_nc()
  endif

  CALL close_ncfile(ncid)

END PROGRAM
