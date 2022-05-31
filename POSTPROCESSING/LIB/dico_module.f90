MODULE DICO
  USE PARAMETERS
  USE USER_INPUT
  USE CONC
  IMPLICIT NONE

  integer, parameter                 :: GAS = 1, AEROSOL = 2, PRECURSOR = 3, CO_CO2 = 4,&
      AQUEOUS = 5, WALL = 6, DIMER = 7, SPINUP_SPEC = 8, INORG = 9, VBS = 10, LAST_PHASE = 11
                                        ! LAST_PHASE is only used to count how many phases we have
  character(len=*), parameter        :: phase_letter = "GAV"

  integer, parameter                 :: function_types = 22
  
  character(len=1), dimension(function_types) :: func_codes=(/'A', 'B', 'D', 'E', 'F', 'G', 'H', 'K', 'L', &
                                                'N', 'O', 'P', 'R', 'T', 'U', 'V', 'X', '1', '2', '3', '4', 'S' /)

  type dico_info
    character*(lco)                  ::  code
    character*(lfo)                  ::  formula
    character*(lfl)                  ::  functions
    integer                          ::  igen
  end type dico_info

  type henry_data
    real                             ::  cf, Keff, rf
  end type henry_data
  
  type koh_data
    real                             :: A, n, Ea, k_298
  end type koh_data
  type kno3_data
    real                             :: A, n, Ea, k_298
  end type kno3_data
  type ko3_data
    real                             :: A, n, Ea, k_298
  end type ko3_data

  type pvap_data
    real                             ::  Tb, dB
    real                             ::  pvap_atm_298, pvap_Cstar_298, delta_H_vap_J_per_mol
  end type pvap_data


  type species_detailed_info
    character*(maxlsp)               ::  code
    character*(lfo)                  ::  formula
    character*(lfl)                  ::  functions
    integer                          ::  phase
    integer                          ::  nfunc

    integer                          ::  dict_indx
    integer                          ::  nc, no, nn, nh, nr, ns, nfl, nbr, ncl
    real                             ::  osc
    real                             ::  nox
    real                             ::  molw
    real                             ::  kendrick_mass
    real                             ::  kendrick_mass_defect
    real                             ::  dbe, ai
    type(henry_data)                 ::  henry
    type(pvap_data)                  ::  pvap
    type(koh_data)                   ::  koh
    type(kno3_data)                  ::  kno3
    type(ko3_data)                   ::  ko3
    character*(16)                   ::  chon
    integer                          ::  igen
    
    logical                          ::  is_nitrate
    integer                          ::  n_nitrates
    integer                          ::  nitrate_substitution(max_nitrates)
  end type species_detailed_info

! number of species in the dictionary
  integer                                                ::  ndic
  type(dico_info),  dimension(:), allocatable            ::  dictionary
  type(species_detailed_info), dimension(:), allocatable :: species

  CONTAINS

  subroutine get_number(chem, molmass, ic,ih,in,io,ir,is,ifl,ibr,icl)
    character*(lfo), intent(in)        :: chem
    integer, intent(out)               :: ic,ih,in,io,ir,is,ifl,ibr,icl
    integer                            :: nc
    character*(lfo)                    :: tchem
    logical                            :: dummy_fg
    real                               :: molmass
    ic= 0
    ih= 0
    in= 0
    io= 0
    ir= 0
    is= 0
    ifl= 0
    ibr= 0
    icl = 0
    tchem = ''
    dummy_fg = .true.

    tchem = chem
    if (chem(1:3) .eq. "#mm" .or. chem(1:3) .eq. "#bb") then
      tchem = ''
      tchem = chem(4:)
    elseif (chem(1:1) .eq. "#") then
      tchem = ''
      tchem = chem(2:)
    endif

    if (tchem(1:1) /= 'C' .and. tchem(1:1) /= 'c' .and. tchem(1:2) /= '-O') return

    nc = index(tchem, ' ') -1
    CALL molarw_number_withexceptions(tchem,molmass, dummy_fg,     &
                    nc,ic,ih,in,io,ir,is,ifl,ibr,icl)

  end subroutine

  elemental function update_species_info(species) result(new_species)
    type(species_detailed_info), intent(in)    :: species
    type(species_detailed_info)                :: new_species

    new_species = species



  END FUNCTION update_species_info

  elemental function write_chon(species) result(new_species)
    type(species_detailed_info), intent(in)    :: species
    type(species_detailed_info)                :: new_species

    character*(4)                              :: cn, hn, on, nn, sn

    cn = ''
    hn = ''
    on = ''
    nn = ''
    sn = ''
    new_species = species
    if(new_species%nc > 0) then
      if(new_species%nc > 9) then
        write(cn, '(A1,I2)') 'C',new_species%nc
      else
        write(cn, '(A1,I1)') 'C',new_species%nc
      endif
    endif
    if(new_species%nh > 0) then
      if(new_species%nh > 9) then
        write(hn, '(A1,I2)') 'H', new_species%nh
      else
        write(hn, '(A1,I1)') 'H',new_species%nh
      endif
    endif
    if(new_species%no > 0) then
      if(new_species%no > 9) then
       write(on, '(A1,I2)') 'O',new_species%no
      else
        write(on, '(A1,I1)') 'O',new_species%no
      endif
    endif
    if(new_species%nn > 0) then
      if(new_species%nn > 9) then
       write(nn, '(A1,I2)') 'N',new_species%nn
      else
        write(nn, '(A1,I1)') 'N',new_species%nn
      endif
    endif
    if(new_species%ns > 0) then
      if(new_species%ns > 9) then
       write(nn, '(A1,I2)') 'S',new_species%ns
      else
        write(nn, '(A1,I1)') 'S',new_species%ns
      endif
    endif


    write(new_species%chon, '(A)') trim(cn)//trim(hn)//trim(on)//trim(nn)


  end function write_chon


END MODULE DICO
