MODULE USER_INPUT
  USE PARAMETERS
  IMPLICIT NONE

  character*(maxlsp)         :: precursor_codes(mxprecu), selected_species(mxselspe), rate_species
  character*(6)              :: input_type, output_type
  integer                    :: nbox, skip_time, n_topspecies, n_toprates
  real                       :: seed_mass, seed_molw
  character*(mfilt), dimension(100) :: gas_dyn_filter, aer_dyn_filter
  LOGICAL    :: flag_chon, flag_selected, flag_phasedist, flag_contributingspecs, flag_soayield, &
      flag_functions, flag_carbonchain, flag_topspec, flag_atomratios, flag_pvap, flag_cstar, flag_henry, &
      flag_elementscontrib, flag_phasestate, flag_entropy,flag_chochonfreq,flag_dbeai,flag_massspectrum, &
      flag_amsfactors, flag_bubble, flag_nitrates, flag_smiles, flag_vbs_param, flag_dyn_filter, flag_ohexposure, &
      flag_no3exposure, flag_potential_atomratios, flag_ml, flag_ppb, flag_ppbC, flag_OHR, flag_NO3R, flag_O3R, &
      flag_rates, flag_topNspc


  NAMELIST  /userinput/  input_type, output_type, nbox, precursor_codes, seed_mass, seed_molw, selected_species, rate_species, &
    flag_chon, flag_selected, flag_phasedist, flag_contributingspecs, flag_soayield, &
    flag_functions, flag_carbonchain, flag_topspec, flag_atomratios, flag_pvap, flag_cstar, flag_henry, &
    flag_elementscontrib, flag_phasestate, flag_entropy,flag_chochonfreq,flag_dbeai,flag_massspectrum, &
    skip_time, flag_amsfactors, n_topspecies, flag_bubble, flag_nitrates, flag_smiles, flag_vbs_param, &
    gas_dyn_filter, aer_dyn_filter, flag_ohexposure, flag_no3exposure, flag_potential_atomratios, flag_ml, &
    flag_ppb, flag_ppbC, flag_OHR, flag_NO3R, flag_O3R, flag_rates, n_toprates, flag_topNspc
END MODULE
