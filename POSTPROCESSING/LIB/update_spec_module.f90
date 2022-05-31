MODULE UPDATE_SPEC
  USE DICO
  USE SORTING, only : find_species_index
  USE PARAMETERS, only : real_kind
  USE USER_INPUT, only: flag_nitrates
  CONTAINS

  SUBROUTINE update_all_species()
    integer           :: i,n, j, ispec

    write(6,*) '-- getting info about all species --'
    ! here we do everything that can be parallelized
!    species = update_species_info(species)

    ! look up dictionary for formulas and functions
    do i = 1, ndic
      do j=1, len_trim(phase_letter)
        ispec = find_species_index(phase_letter(j:j)//dictionary(i)%code)
        !if (ispec .ne. -1) then
        if (ispec .gt. 0) then
          species(ispec)%formula = dictionary(i)%formula
          species(ispec)%functions = dictionary(i)%functions
          species(ispec)%igen = dictionary(i)%igen
          species(ispec)%dict_indx = i
          select case (phase_letter(j:j))
            case ("G")
              species(ispec)%phase = GAS
            case ("A")
              species(ispec)%phase = AEROSOL
            case ("W")
              species(ispec)%phase = WALL
            case ("D")
              species(ispec)%phase = DIMER
            case ("V")
              species(ispec)%phase = VBS
          end select
          if (species(ispec)%code(1:5) == "GCO2 " .or. &
              species(ispec)%code(1:5) == "GCO  " ) then
            species(ispec)%phase = CO_CO2
          endif
          if (species(ispec)%code(1:5) == "GCH4 ") then
            species(ispec)%phase = INORG
          endif
! SPECIAL FOR OSULFA TYPO: Sept 2020
          IF (species(ispec)%formula(1:32) == "CH3C(OH)CH2(OH))CH(OH)CH2(OSO3) ") then
              species(ispec)%formula(1:32)  = "CH3C(OH)(CH2(OH))CH(OH)CH2(OSO3)"
          ENDIF

          if ( any(precursor_codes == trim(species(ispec)%code))) then
            species(ispec)%phase = PRECURSOR
          endif
          if (species(ispec)%code(1:2) .eq. "GX" .or. species(ispec)%code(1:2) .eq. "GY") then
            species(ispec)%phase = SPINUP_SPEC
          endif
!        else 
!          write(6,*) " could not find: ", phase_letter(j:j)//dictionary(i)%code
!          write(6,*) " in subroutine update_all_species"

        endif
      enddo
    enddo

    ! and the rest
    do i = 1, numsp
      call get_number(species(i)%formula, &
                      species(i)%molw,&
                      species(i)%nc,&
                      species(i)%nh,&
                      species(i)%nn,&
                      species(i)%no,&
                      species(i)%nr,&
                      species(i)%ns,&
                      species(i)%nfl,&
                      species(i)%nbr,&
                      species(i)%ncl)
 ! now we can estimate carbon oxidation state
      species(i)%osc = 0.
      if (species(i)%nc > 0) then
        species(i)%osc = (2*real(species(i)%no, kind = real_kind) - &
                        real(species(i)%nh, kind = real_kind)) / &
                        real(species(i)%nc, kind = real_kind)
      endif
    enddo

    !calculate kendrick mass and kendrick mass defect,
    ! defining m(CH2)  = 14.000
    species%kendrick_mass = species%molw * 14./14.027
    ! assuming nominal mass must the rounded value of the exact mass
    species%kendrick_mass_defect = anint(species%molw) - species%kendrick_mass

    ! calculate double bond equivalent
    species%dbe = species%nc - species%nh/2 + species%nn/2 + 1
    ! calculate aromaticity index
    where (species%dbe <= 0 .or. (species%nc - species%no - species%ns - species%nn) <= 0)
      species%ai = 0
    elsewhere
      species%ai = (1 + species%nc - species%no - species%ns - 0.5*species%nh)/ &
                   (species%nc - species%no - species%ns - species%nn)
    end where

! DO NOT CHANGE THE PHASE OF LUMPED SPECIES
    where (species%formula(1:4) == "Lump")
      species%phase = species%phase
    else where (species%nc == 0)
      species%phase = INORG
    end where

    species = write_chon(species)

    precursor_indx = 0
    n=1
    do i=lbound(species, 1), ubound(species, 1)
      if(species(i)%phase == PRECURSOR) then
        precursor_indx(n) = i
        n=n+1
      endif
    enddo

    if (all(precursor_indx == 0)) then
      write(6,*) ' --error--  precursors ',precursor_codes, ' not found in species list'
      stop
    endif
    
    ! do nitrates
    if (flag_nitrates) then
      write(6,*) ' ... updating nitrates'
      call init_nitrates()
      write(6,*) ' ... done updating nitrates'
    endif

  END SUBROUTINE
  

  SUBROUTINE init_nitrates()
    integer :: ispec, j
  
  ! grbond
    integer   :: nc, bond(mca,mca), dbflg, nring
    character*(lgr)  :: group(mca)
    character*(lfo) :: chem, tchem
  
    
    do ispec = 1, numsp
      if(index(species(ispec)%functions, 'N') .eq. 0) cycle
      species(ispec)%is_nitrate = .true.
      chem = ""
      tchem = ""
      tchem = species(ispec)%formula
      if (tchem(1:3) .eq. "#mm") then 
        chem = tchem(4:lfo)
      else if (tchem(1:1) .eq. "#" ) then 
        chem = tchem(2:lfo)
      else
        chem = tchem
      endif
      if (chem(1:1) .ne. "C" .and. chem(1:1) .ne. "c") cycle
      nc = len(trim(chem))
      group = ""
      bond = 0
      dbflg = 0
      nring = 0      
      species(ispec)%n_nitrates = 0
      call grbond(chem,nc,group,bond,dbflg,nring)
      
      do j = 1, mca
        if (group(j)(1:1) .eq. ' ') exit
        if (index(group(j), '(ONO2)') .eq. 0) cycle
        species(ispec)%n_nitrates = species(ispec)%n_nitrates + 1   
        if(species(ispec)%n_nitrates > max_nitrates) then
          write(6, *) "Error, too many nitrates on molecule"
          write(6,*) chem
          write(6,*) "increase max_nitrates"
          exit
        endif
        if (group(j)(1:3) == "CH2") then
          species(ispec)%nitrate_substitution(species(ispec)%n_nitrates) = 1
        else if (group(j)(1:2) == "CH") then
          species(ispec)%nitrate_substitution(species(ispec)%n_nitrates) = 2
        else if (group(j)(1:1) == "C") then
          species(ispec)%nitrate_substitution(species(ispec)%n_nitrates) = 3
        else
          species(ispec)%nitrate_substitution(species(ispec)%n_nitrates) = 0          
        endif
      enddo
    enddo
  END SUBROUTINE
  


  END MODULE

