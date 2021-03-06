MODULE compute
  USE CONC
  USE USER_INPUT
  USE IO
  USE BINS
  USE dynamic_filter
  IMPLICIT NONE

  CONTAINS

!-----------------------------------------------------------------------
  SUBROUTINE calc_environmental_param()
!-----------------------------------------------------------------------
! always called
! output filename = "environmental_parameters"
!-----------------------------------------------------------------------
    real, dimension(:,:), allocatable      :: environ_params
    integer, parameter                     :: ind_time = 1, ind_rh = 2, ind_temp = 3, ind_pres = 4
    character*(header_length)              :: header(ind_pres)

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(environ_params(ndat, ind_pres))

    header(ind_time) = "Time [s]"
    header(ind_rh)   = "RH [%]"
    header(ind_temp) = "Temperature [K]"
    header(ind_pres) = "Pressure [hPa]"

    environ_params(:, ind_time) = time
    environ_params(:, ind_rh)   = rh
    environ_params(:, ind_temp) = temperature
    environ_params(:, ind_pres) = pressure

    filename = "environmental_parameters"
    call write_2D_array(environ_params, header , filename)

    deallocate(environ_params)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  END SUBROUTINE calc_environmental_param
!-----------------------------------------------------------------------
  SUBROUTINE calc_precursor_lifetime()
!-----------------------------------------------------------------------
! always called
! output filename = "precu_lifetime"
!-----------------------------------------------------------------------
    integer                               :: i
    character*(header_length)             :: header(2)

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    if(nprecu .eq. 1) then

      write(6,*) '-- calculating precursors lifetime --'

      allocate(precursor_lifetime(ndat, 2))

      precursor_lifetime(:,1) = time
      precursor_lifetime(:,2) =  concentrations(1, precursor_indx(1))/concentrations(:, precursor_indx(1))
      where (precursor_lifetime(:,2) > 0)
        precursor_lifetime(:,2) = log(precursor_lifetime(:,2))
      elsewhere
        precursor_lifetime(:,2) = -1
      end where

      header(1) = "Time [s]"
      header(2) = "Precursor Lifetimes"
      filename = "precu_lifetime"
      call write_2D_array(precursor_lifetime, header , filename)

      ! fill array with closest times indices corresponding to tau = 1, 2, 3, 4, 5
      do i = 1, 5
        efold_indices(i) = minloc(abs(precursor_lifetime(:,2) - i), dim = 1)
      enddo

      deallocate(precursor_lifetime)
    else ! if more than one precursor, just split in 5 equal parts
      do i = 1, 5
        efold_indices(i) = i * ndat / 5
      enddo
    endif

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  END SUBROUTINE calc_precursor_lifetime
!-----------------------------------------------------------------------
  SUBROUTINE phase_distribution()
!-----------------------------------------------------------------------
! called with flag_phasedist
! output filename(1) = "carb_dist_time"
! output filename(2) = "mass_dist_time"
!-----------------------------------------------------------------------

    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 9, ncols = nphase+1
    integer                               :: itime, ispec, iphase
    real, dimension(:,:), allocatable     :: carb_dist, mass_dist
    character*(header_length)             :: header(ncols)
    real, dimension(:), allocatable       :: ppbfac

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating carbon phase distribution --'

    allocate(carb_dist(ndat, ncols), mass_dist(ndat, ncols))
    header(1) = "Time [s]"
    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase+1) = "Gas Phase"
        case (AEROSOL)
          header(iphase+1) = "Aerosol Phase"
        case (PRECURSOR)
          header(iphase+1) = "Precursor"
        case (CO_CO2)
          header(iphase+1) = "CO+CO2"
        case (AQUEOUS)
          header(iphase+1) = "Aqueous Phase"
        case (WALL)
          header(iphase+1) = "Wall"
        case (DIMER)
          header(iphase+1) = "Dimers"
        case (SPINUP_SPEC)
          header(iphase+1) = "Spinup Species"
        case (INORG)
          header(iphase+1) = "Inorganics+CH4"
        case default
          header(iphase+1) = ""
      end select
    enddo

! factor to convert to ppbC: time dependent because depends on pressure and temp
    allocate(ppbfac(ndat))
    ppbfac = 7.245461056e+09*pressure/temperature

    do itime = 1, ndat
      carb_dist(itime, :) = 0.
      mass_dist(itime, :) = 0.
      do ispec = 1, numsp
        do iphase = 1, nphase
          if (species(ispec)%phase == iphase) then
            carb_dist(itime, iphase+1) = carb_dist(itime, iphase+1) + &
              concentrations(itime, ispec) * species(ispec)%nC/ ppbfac(itime)
            mass_dist(itime, iphase+1) = mass_dist(itime, iphase+1) + &
              concentrations(itime, ispec) * species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
            exit
          endif
        enddo
      enddo
    enddo

    carb_dist(:, 1) = time
    mass_dist(:, 1) = time
    filename = "carb_dist_time"
    call write_2D_array(carb_dist, header , filename)
    filename = "mass_dist_time"
    call write_2D_array(mass_dist, header , filename)

    deallocate(ppbfac, carb_dist, mass_dist)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  END SUBROUTINE phase_distribution
!-----------------------------------------------------------------------
  subroutine calc_soa_yield()
!-----------------------------------------------------------------------
! called with flag_soayield
! output filename = "soa_yield_time"
!-----------------------------------------------------------------------

    real, dimension(:), allocatable    :: mass_precu, mass_soa
    real, dimension(:,:), allocatable  :: soa_yield
    character*(header_length)          :: header(2)
    integer                            :: itime, ispec

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating soa yield --'
    allocate(mass_precu(ndat), mass_soa(ndat), soa_yield(ndat, 2))

    do itime = 1, ndat
      mass_precu(itime) = 0
      mass_soa(itime) = 0
      do ispec = 1, numsp ! precaution in case there are several precursors
        if (species(ispec)%phase == PRECURSOR) then
          mass_precu(itime) = mass_precu(itime) + &
            concentrations(itime, ispec) * species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
        elseif (species(ispec)%phase == AEROSOL) then
          mass_soa(itime) = mass_soa(itime) + &
            concentrations(itime, ispec) * species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
        endif
      enddo
    enddo

    header(1) = 'Time [s]'
    header(2) = 'SOA yield [%]'
    soa_yield(:, 1) = time
    soa_yield(:, 2) = 100*mass_soa(:) / (mass_precu(1) - mass_precu(:))
    where (mass_precu == mass_precu(1))
      soa_yield(:,2) = 0.0
    end where

    filename = "soa_yield_time"
    call write_2D_array(soa_yield, header , filename)

    deallocate(mass_precu, mass_soa, soa_yield)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_soa_yield
!-----------------------------------------------------------------------
  subroutine count_functions()
!-----------------------------------------------------------------------
! called with flag_functions
! output filename(1) = "rofc_gas_time"
! output filename(2) = "rofc_aer_time"
! output filename(3) = "functions_mass_gas_time"
! output filename(4) = "functions_mass_aer_time"
! output filename(5) = "rofctot_gas_time"
! output filename(5) = "rofctot_aer_time"
!-----------------------------------------------------------------------
    character*(15)                    :: func_descs(function_types)
    real, dimension(:,:), allocatable :: rof_c_aer, rof_c_gas
    real, dimension(:,:), allocatable :: rof_cp_aer, rof_cp_gas
    real, dimension(:,:), allocatable :: massfunc_aer, massfunc_gas
    real, dimension(:), allocatable   :: ncarbon_aer, ncarbon_gas
    real, dimension(:), allocatable   :: ncarbon_tot_aer, ncarbon_tot_gas

    character*(header_length)          :: header(function_types+1)

    integer                 :: itime, ispec, ifunc, ifl

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating functions distribution --'
    write(6,*) '------- for MULTI-CARBON SPECIES -------'


    func_descs(1)  =  "-CO(OH)"
    func_descs(2)  =  "-Br"
    func_descs(3)  =  "-CHO"
    func_descs(4)  =  "R-O-R"
    func_descs(5)  =  "-F"
    func_descs(6)  =  "-CO(OOH)"
    func_descs(7)  =  "-OOH"
    func_descs(8)  =  '>CO'
    func_descs(9)  =  '-Cl'
    func_descs(10) =  '-ONO2'
    func_descs(11) =  '-OH'
    func_descs(12) =  '-CO(OONO2)'
    func_descs(13) =  'aromatic rings'
    func_descs(14) =  "aliphatic rings"
    func_descs(15) =  '>C=C<'
    func_descs(16) =  "-NO2"
    func_descs(17) =  "-C=C=O"
    func_descs(18) =  "->C(O.)"
    func_descs(19) =  "->C(OO.)"
    func_descs(20) =  "-CO(OO.)"
    func_descs(21) =  ">C.(OO.)"
    func_descs(22) =  "-OSO3"
    allocate(rof_c_aer(ndat, function_types+1), &
             rof_c_gas(ndat, function_types+1), &
             rof_cp_gas(ndat, function_types+1), &
             rof_cp_aer(ndat, function_types+1), &
             massfunc_aer(ndat, function_types+1), &
             massfunc_gas(ndat, function_types+1), &
             ncarbon_aer(ndat), &
             ncarbon_gas(ndat), &
             ncarbon_tot_gas(ndat), &
             ncarbon_tot_aer(ndat))
    ncarbon_gas = 0.
    ncarbon_tot_gas = 0.
    ncarbon_tot_aer = 0.
    ncarbon_aer = 0.
    do itime = 1, ndat
      rof_c_aer(itime, :) = 0.
      rof_c_gas(itime, :) = 0.
      rof_cp_gas(itime, :) = 0.
      rof_cp_aer(itime, :) = 0.
      massfunc_gas(itime, :) = 0.
      massfunc_aer(itime, :) = 0.
      do ispec = 1, numsp

        if (species(ispec)%phase == GAS) then
          ncarbon_gas(itime) = ncarbon_gas(itime)+ &
            species(ispec)%nc*concentrations(itime, ispec)
          do ifl = 1, lfl
            do ifunc = 1, function_types
              if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc)) then
              !if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc) .AND. &
              !    species(ispec)%nc.GE.1 ) then
                rof_c_gas(itime, ifunc) = rof_c_gas(itime, ifunc) +  &
                  concentrations(itime, ispec)
                massfunc_gas(itime, ifunc) = massfunc_gas(itime, ifunc) + &
                  concentrations(itime, ispec)* species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
              endif
            enddo
          enddo
        elseif (species(ispec)%phase == AEROSOL) then
          ncarbon_aer(itime) = ncarbon_aer(itime)+ &
            species(ispec)%nc*concentrations(itime, ispec)
          do ifl = 1, lfl
            do ifunc = 1, function_types
              if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc)) then
              !if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc) .AND. &
              !    species(ispec)%nc.GE.1 ) then
                rof_c_aer(itime, ifunc) = rof_c_aer(itime, ifunc) +  &
                  concentrations(itime, ispec)
                massfunc_aer(itime, ifunc) = massfunc_aer(itime, ifunc) + &
                  concentrations(itime, ispec)* species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
              endif
            enddo
          enddo
        endif
        if (species(ispec)%phase == GAS .or. &
            species(ispec)%phase == PRECURSOR) then
          ncarbon_tot_gas(itime) = ncarbon_tot_gas(itime)+ &
            species(ispec)%nc*concentrations(itime, ispec)
          do ifl = 1, lfl
            do ifunc = 1, function_types
              if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc) .AND. &
                  species(ispec)%nc.GE.1 ) then
                rof_cp_gas(itime, ifunc) = rof_cp_gas(itime, ifunc) +  &
                  concentrations(itime, ispec)
              endif
            enddo
          enddo
        endif
        if (species(ispec)%phase == AEROSOL) then 
          ncarbon_tot_aer(itime) = ncarbon_tot_aer(itime)+ &
            species(ispec)%nc*concentrations(itime, ispec)
          do ifl = 1, lfl
            do ifunc = 1, function_types
              if (species(ispec)%functions(ifl:ifl) == func_codes(ifunc) .AND. &
                  species(ispec)%nc.GE.1 ) then
                rof_cp_aer(itime, ifunc) = rof_cp_aer(itime, ifunc) +  &
                  concentrations(itime, ispec)
              endif
            enddo
          enddo
        endif
      enddo
    enddo

    do ifunc = 1, function_types
      where (ncarbon_gas > 0)
        rof_c_gas(:, ifunc) = rof_c_gas(:, ifunc)/ncarbon_gas(:)
      elsewhere
        rof_c_gas(:, ifunc) = 0.
      end where
      where (ncarbon_aer > 0)
        rof_c_aer(:, ifunc) = rof_c_aer(:, ifunc)/ncarbon_aer(:)
      elsewhere
        rof_c_aer(:, ifunc) = 0.
      end where
      where (ncarbon_tot_gas > 0)
        rof_cp_gas(:, ifunc) = rof_cp_gas(:, ifunc)/ncarbon_tot_gas(:)
      elsewhere
        rof_cp_gas(:, ifunc) = 0.
      end where
      where (ncarbon_tot_aer > 0)
        rof_cp_aer(:, ifunc) = rof_cp_aer(:, ifunc)/ncarbon_tot_aer(:)
      elsewhere
        rof_cp_aer(:, ifunc) = 0.
      end where
    enddo

    header(function_types+1) = 'Time [s]'
    header(1:function_types) = func_descs(:)

    filename = "rofc_gas_time"
    rof_c_gas(:, function_types+1) = time
    call write_2D_array(rof_c_gas, header , filename)
    filename = "rofc_aer_time"
    rof_c_aer(:, function_types+1) = time
    call write_2D_array(rof_c_aer, header , filename)

    filename = "rofctot_gas_time"
    rof_cp_gas(:, function_types+1) = time
    call write_2D_array(rof_cp_gas, header , filename)
    filename = "rofctot_aer_time"
    rof_cp_aer(:, function_types+1) = time
    call write_2D_array(rof_cp_aer, header , filename)


    filename = "functions_mass_gas_time"
    massfunc_gas(: , function_types+1) = time
    call write_2D_array(massfunc_gas, header , filename)
    filename = "functions_mass_aer_time"
    massfunc_aer(: , function_types+1) = time
    call write_2D_array(massfunc_aer, header , filename)

    deallocate(rof_c_aer,    rof_c_gas, &
               rof_cp_gas, ncarbon_tot_gas,&
               rof_cp_aer, ncarbon_tot_aer,&
               ncarbon_aer,  ncarbon_gas, &
               massfunc_aer, massfunc_gas)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine count_functions
!-----------------------------------------------------------------------
  subroutine carbon_chain()
!-----------------------------------------------------------------------
! called with flag_carbonchain
! output filename(1) = "carbon_chain_gas_ppbC_time"
! output filename(2) = "carbon_chain_aer_ppbC_time"
!-----------------------------------------------------------------------

    real, dimension(:,:), allocatable :: carbon_chain_dist_gas, carbon_chain_dist_aer
    integer                           :: ic, itime, maxc
    character*(header_length), dimension(:), allocatable   :: header
    real                              :: ppbfac

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating chain size distribution --'
    ! find max number of carbon atoms in a species
    maxc = maxval(species(:)%nc)

    allocate(carbon_chain_dist_gas(ndat, maxc + 1), &
             carbon_chain_dist_aer(ndat, maxc + 1), &
             header(maxc + 1))

    do itime=1, ndat
      carbon_chain_dist_gas(itime, :) = 0
      carbon_chain_dist_aer(itime, :) = 0
      ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
      do ic = 1, maxc
        write(header(ic),*) ic
        carbon_chain_dist_gas(itime, ic) = sum( &
          concentrations(itime, :)*species%nc/ppbfac, &
           mask = (species%nc == ic) .and. (species%phase == GAS))
        carbon_chain_dist_aer(itime, ic) = sum( &
          concentrations(itime, :)*species%nc/ppbfac, &
          mask = (species%nc == ic) .and. (species%phase == AEROSOL))

      enddo
    enddo
    header(maxc+1) = "Time [s]"
    carbon_chain_dist_gas(:, maxc+1) = time
    carbon_chain_dist_aer(:, maxc+1) = time
    filename = "carbon_chain_gas_ppbC_time"
    call write_2D_array(carbon_chain_dist_gas, header , filename)
    filename = "carbon_chain_aer_ppbC_time"
    call write_2D_array(carbon_chain_dist_aer, header , filename)

    deallocate(carbon_chain_dist_gas, carbon_chain_dist_aer, header)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine carbon_chain
!-----------------------------------------------------------------------
  subroutine calc_chon()
!-----------------------------------------------------------------------
! called with flag_chon
! output filename(1) = "top_chon_ppbc_gas_time"
! output filename(2) = "top_chon_ppbc_aer_time"
!-----------------------------------------------------------------------

    integer                                         :: topN, n_chon
    character*(16), dimension(:), allocatable       :: topN_chon
    character*(16), dimension(:), allocatable       :: chon_name_table
    logical, dimension(:), allocatable              :: unique_chon_mask, phase_mask
    integer, dimension(:), allocatable              :: index_vector, phase_index_vector, top_chon_indx

    integer                                         :: itime, ispec, iphase, n_phase, ichon, itop
    real, dimension(:,:), allocatable               :: top_chon_concentrations
    real, dimension(:), allocatable                 :: time_integrated_sum
    character*(header_length), dimension(:), allocatable   :: header
    real                              :: ppbfac
    logical                                         :: found_fg

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    topN = 10
    allocate(top_chon_concentrations(ndat, topN+2 ), &
             header(topN + 2), &
             unique_chon_mask(numsp), phase_mask(numsp), &
             top_chon_indx(topN), topN_chon(topN))
  ! find all unique chon
    write(6, *) ' -- calculating top contributing elemental formulas --'
    do iphase = GAS, AEROSOL, AEROSOL-GAS
      if (iphase == GAS) then
        write(6,*) '      in GAS phase'
      else
        write(6,*) '      in AEROSOL phase'
      endif

      phase_mask = species%phase == iphase
      n_phase = count(phase_mask)
      allocate(phase_index_vector(n_phase))
      phase_index_vector = pack((/ (ispec, ispec = 1, numsp) /), &
                                 phase_mask)

      unique_chon_mask = .FALSE.
      unique_chon_mask(phase_index_vector) = .TRUE.
      do ispec = n_phase,2,-1
         unique_chon_mask(phase_index_vector(ispec)) =  &
           .not.(any(species(phase_index_vector(1:ispec-1))%chon== &
                     species(phase_index_vector(ispec))%chon))
      end do
      n_chon = count(unique_chon_mask)

      write(6, *) '       found ', n_chon, ' unique elemental formulas '
      ! Make an index vector
      allocate(index_vector(n_chon), chon_name_table(n_chon), &
               time_integrated_sum(n_chon))
      index_vector = pack( (/ (ispec, ispec=1,numsp) /) ,unique_chon_mask)
      chon_name_table = species(index_vector)%chon



      top_chon_concentrations=0
      time_integrated_sum = 0

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ichon = 1, n_chon

          time_integrated_sum(ichon) = time_integrated_sum(ichon) + &
            sum(concentrations(itime,:)*species%nc/ppbfac, &
                mask = species%phase == iphase .and. &
                species%chon == chon_name_table(ichon))

        enddo
      enddo

      top_chon_indx(1) = maxloc(time_integrated_sum, dim = 1)
      topN_chon(1) = chon_name_table(top_chon_indx(1))
      time_integrated_sum(top_chon_indx(1)) = 0
      do itop =2, topN
        top_chon_indx(itop) = maxloc(time_integrated_sum, dim = 1)
        topN_chon(itop) = chon_name_table(top_chon_indx(itop))
        time_integrated_sum(top_chon_indx(itop)) = 0
      enddo

      header(1:topN) = topN_chon
      header(topN+1) = "Others [ppbC]"
      header(topN+2) = "Time [s]"

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec = 1, numsp
          if (species(ispec)%phase /= iphase) CYCLE
          found_fg = .false.

          do ichon = 1, topN
            if (species(ispec)%chon == topN_chon(ichon)) then
              found_fg = .true.

              top_chon_concentrations(itime, ichon) = &
                top_chon_concentrations(itime, ichon) + concentrations(itime, ispec) * &
                  species(ispec)%nc/ppbfac
              exit
            endif
          enddo
          if (.not. found_fg) then
            top_chon_concentrations(itime, topN+1) = &
                top_chon_concentrations(itime, topN+1) + concentrations(itime, ispec) * &
                  species(ispec)%nc/ppbfac
          endif
        enddo
      enddo

      top_chon_concentrations(:, topN+2) = time
      if (iphase == GAS) then
        filename = "top_chon_ppbc_gas_time"
      else if (iphase == AEROSOL) then
        filename = "top_chon_ppbc_aer_time"
      endif

      call write_2D_array(top_chon_concentrations, header , filename)



      ! find the most concentrated chon formulas
      deallocate(index_vector, chon_name_table, phase_index_vector, &
                 time_integrated_sum)
    enddo
    deallocate(top_chon_concentrations, header, unique_chon_mask, phase_mask, &
               top_chon_indx, topN_chon)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_chon
!-----------------------------------------------------------------------
  subroutine calc_atom_ratios()
!-----------------------------------------------------------------------
! called with flag_atom_ratios
! output filename = "atom_ratios_gas"
! output filename = "atom_ratios_aer"
!-----------------------------------------------------------------------

    real, dimension(:,:), allocatable    :: atom_ratios
    character*(header_length), dimension(:), allocatable   :: header
    real, dimension(:), allocatable   :: sumnc
    integer                              :: itime, i

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating O/C, H/C and N/C ratios --'
    allocate(atom_ratios(ndat, 4), header(4), sumnc(ndat))

    header(1) = "O/C"
    header(2) = "H/C"
    header(3) = "N/C"
    header(4) = "Time [s]"

! first in the gas phase
    do itime = 1, ndat
      atom_ratios(itime, 1) = sum( &
            concentrations(itime,:)*species%no, &
            mask = species%phase == GAS)
      atom_ratios(itime, 2) = sum( &
            concentrations(itime,:)*species%nh, &
            mask = species%phase == GAS)
      atom_ratios(itime, 3) = sum( &
            concentrations(itime,:)*species%nn, &
            mask = species%phase == GAS)
      sumnc(itime) = sum( &
            concentrations(itime, :)*species%nc, &
            mask = species%phase == GAS)
    enddo
    do i = 1,3
      atom_ratios(:,i) = atom_ratios(:,i)/sumnc(:)
    enddo

    atom_ratios(:,4) = time
    filename = "atom_ratios_gas"
    call write_2D_array(atom_ratios, header , filename)

! then in the aerosol phase
    do itime = 1, ndat
      atom_ratios(itime, 1) = sum( &
            concentrations(itime,:)*species%no, &
            mask = species%phase == AEROSOL)
      atom_ratios(itime, 2) = sum( &
            concentrations(itime,:)*species%nh, &
            mask = species%phase == AEROSOL)
      atom_ratios(itime, 3) = sum( &
            concentrations(itime,:)*species%nn, &
            mask = species%phase == AEROSOL)
      sumnc(itime) = sum( &
            concentrations(itime, :)*species%nc, &
            mask = species%phase == AEROSOL)
    enddo
    do i = 1,3
      atom_ratios(:,i) = atom_ratios(:,i)/sumnc(:)
    enddo

    filename = "atom_ratios_aer"
    call write_2D_array(atom_ratios, header , filename)


    deallocate(atom_ratios, header, sumnc)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_atom_ratios
!-----------------------------------------------------------------------
  subroutine calc_ams_factors()
!-----------------------------------------------------------------------
! called with flag calc_amsfactors
! output filename = "ams_factors_aer_ug"
!-----------------------------------------------------------------------

    character*(header_length), dimension(:), allocatable   :: header
    real, dimension(:,:), allocatable :: ams_factors
    integer, parameter :: nfactors = 6
    integer, parameter :: ind_mo_ooa = 1, ind_lo_ooa = 2, ind_iepox_soa = 3
    integer, parameter :: ind_adoa = 4, ind_bboa = 5, ind_hoa = 6, ind_time = 7
    integer :: itime, ifac, ispec

    real, dimension(nfactors, 4) :: limits
    integer, parameter :: ind_looc = 1, ind_hioc = 2, ind_lohc = 3, ind_hihc = 4

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    limits(ind_mo_ooa, ind_looc) = 1.09 - 0.17
    limits(ind_mo_ooa, ind_hioc) = 1.09 + 0.17
    limits(ind_mo_ooa, ind_lohc) = 1.27 - 0.12
    limits(ind_mo_ooa, ind_hihc) = 1.27 + 0.12

    limits(ind_lo_ooa, ind_looc) = 0.72 - 0.10
    limits(ind_lo_ooa, ind_hioc) = 0.72 + 0.10
    limits(ind_lo_ooa, ind_lohc) = 1.49 - 0.07
    limits(ind_lo_ooa, ind_hihc) = 1.49 + 0.07

    limits(ind_iepox_soa, ind_looc) = 0.93 - 0.10
    limits(ind_iepox_soa, ind_hioc) = 0.93 + 0.10
    limits(ind_iepox_soa, ind_lohc) = 1.39 - 0.07
    limits(ind_iepox_soa, ind_hihc) = 1.39 + 0.07

    limits(ind_adoa, ind_looc) = 0.40 - 0.05
    limits(ind_adoa, ind_hioc) = 0.40 + 0.05
    limits(ind_adoa, ind_lohc) = 1.63 - 0.02
    limits(ind_adoa, ind_hihc) = 1.63 + 0.02

    limits(ind_bboa, ind_looc) = 0.61 - 0.08
    limits(ind_bboa, ind_hioc) = 0.61 + 0.08
    limits(ind_bboa, ind_lohc) = 1.57 - 0.04
    limits(ind_bboa, ind_hihc) = 1.57 + 0.04

    limits(ind_hoa, ind_looc) = 0.18 - 0.02
    limits(ind_hoa, ind_hioc) = 0.18 + 0.02
    limits(ind_hoa, ind_lohc) = 1.94 - 0.02
    limits(ind_hoa, ind_hihc) = 1.94 + 0.02

    allocate(ams_factors(ndat, nfactors + 1), &
             header(nfactors + 1))

    header(ind_mo_ooa) = "MO-OOA [??g m-3]"
    header(ind_lo_ooa) = "LO-OOA [??g m-3]"
    header(ind_iepox_soa) = "IEPOX-SOA [??g m-3]"
    header(ind_adoa) = "ADOA [??g m-3]"
    header(ind_bboa) = "BBOA [??g m-3]"
    header(ind_hoa) = "HOA [??g m-3]"
    header(ind_time) = "Time [s]"

    do itime = 1, ndat
      ams_factors(itime, ind_time) = time(itime)
      do ifac = 1, nfactors
        ams_factors(itime, ifac) = 0.
        do ispec = 1, numsp
          if(species(ispec)%phase == AEROSOL .and. &
            species(ispec)%nc > 0 .and. &
           (species(ispec)%no / species(ispec)%nc) > limits(ifac, ind_looc) .and. &
           (species(ispec)%no / species(ispec)%nc) < limits(ifac, ind_hioc) .and. &
           (species(ispec)%nh / species(ispec)%nc) > limits(ifac, ind_lohc) .and. &
           (species(ispec)%nh / species(ispec)%nc) < limits(ifac, ind_hihc)) then

            ams_factors(itime, ifac) = ams_factors(itime, ifac) + &
              concentrations(itime, ispec)*species(ispec)%molw*1.660578e-12
          endif
        enddo
      enddo
    enddo

    filename = "ams_factors_aer_ug"
    call write_2D_array(ams_factors, header , filename)

    deallocate(ams_factors, header)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_ams_factors
!-----------------------------------------------------------------------
  subroutine calc_pvap_distribution(opt_resolution)
!-----------------------------------------------------------------------
! called with flag_pvap
! output filename(1) ="pvap_atm_distribution_ppbC_time"
! output filename(2) ="pvap_atm_distribution_ug_time"
!-----------------------------------------------------------------------

  ! optional resolution of the pvap spectrum, in log10(atm) units
  ! default is 2
    integer, optional :: opt_resolution
    integer   :: resolution
    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 3
    INTEGER,PARAMETER :: ncols = nphase+3
    INTEGER,PARAMETER :: ind_lowerlim = ncols-2
    INTEGER,PARAMETER :: ind_upperlim = ncols-1
    INTEGER,PARAMETER :: ind_time = ncols
    real, dimension(:,:), allocatable   :: pvap_distribution_ppbC, pvap_distribution_ug
    character*(header_length), dimension(:), allocatable   :: header_ppbC, header_ug
    integer                             :: nbins, ibin, i, iphase, minpvap, maxpvap
    integer, dimension(:), allocatable  :: upperlim, lowerlim
    real                                :: ppbfac, ugfac
    logical, dimension(:), allocatable  :: spec_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(spec_mask(numsp))

    if (present(opt_resolution)) then
      resolution = opt_resolution
    else
      resolution = 2
    endif

    write(6,*) '-- calculating pvap distribution --'
    ! compute number of bins
    minpvap = floor(log10(minval(species%pvap%pvap_atm_298,  &
                     mask = species%pvap%pvap_atm_298 > 0)))
    maxpvap = ceiling(log10(maxval(species%pvap%pvap_atm_298)))
    nbins = floor(real((maxpvap - minpvap)/resolution, kind = 8))
! constrain bin center to be EVEN if resolution = 2
    if(resolution.EQ.2.AND.MOD(minpvap,2).EQ.0)minpvap = minpvap -1
    if(resolution.EQ.2.AND.MOD(maxpvap,2).EQ.0)maxpvap = maxpvap +1

! NB: formula truncates top bin if (maxpvap-minpvap)/resolution is non-integral
    if(MOD((maxpvap - minpvap),resolution).gt.0) nbins = nbins+1

    allocate(upperlim(nbins), lowerlim(nbins), &
             pvap_distribution_ppbC(ndat*nbins, ncols), &
             pvap_distribution_ug(ndat*nbins, ncols), &
             header_ppbC(ncols), header_ug(ncols))


    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header_ppbC(iphase) = "CO+CO2 [ppbC]"
          header_ug(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          header_ug(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header_ppbC(iphase) = "Wall [ppbC]"
          header_ug(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header_ppbC(iphase) = "Dimers [ppbC]"
          header_ug(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header_ppbC(iphase) = "Spinup Species [ppbC]"
          header_ug(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header_ppbC(iphase) = "Inorganics [ppbC]"
          header_ug(iphase) = "Inorganics [ug m-3]"
      end select
    enddo

    header_ppbC(ind_lowerlim) = "Pvap bin lower limit [log10(atm)]"
    header_ppbC(ind_upperlim) = "Pvap bin upper limit [log10(atm)]"
    header_ppbC(ind_time) = "Time [s]"
    header_ug(ind_lowerlim) = "Pvap bin lower limit [log10(atm)]"
    header_ug(ind_upperlim) = "Pvap bin upper limit [log10(atm)]"
    header_ug(ind_time) = "Time [s]"


    do ibin = 1, nbins
      lowerlim(ibin) = minpvap + resolution*(ibin-1)
      upperlim(ibin) = minpvap + resolution*ibin
    enddo
    ugfac = 1.660578e-12

    pvap_distribution_ppbC = 0.
    pvap_distribution_ug = 0.
    do i=1,ndat
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_lowerlim) = lowerlim
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_upperlim) = upperlim
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_time) = time(i)
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,ind_lowerlim) = lowerlim
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,ind_upperlim) = upperlim
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,ind_time) = time(i)

       ppbfac = 7.245461056e+09*pressure(i)/temperature(i)
       do iphase = 1, nphase
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%pvap%pvap_atm_298) >= lowerlim(ibin) .and. &
                     log10(species%pvap%pvap_atm_298) < upperlim(ibin)
           pvap_distribution_ppbC(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%nc/ppbfac, &
              mask = spec_mask)
           pvap_distribution_ug(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%molw * ugfac, &
              mask = spec_mask)
         enddo
      enddo
    enddo
    filename="pvap_atm_distribution_ppbC_time"
    call write_2D_array(pvap_distribution_ppbC, header_ppbC , filename)
    filename="pvap_atm_distribution_ug_time"
    call write_2D_array(pvap_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               pvap_distribution_ppbC, pvap_distribution_ug, &
               header_ppbC, header_ug)
  
    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_pvap_distribution
!-------------------------------------------------------------------
  subroutine calc_pvap_distribution_lifetime(opt_resolution)
!-------------------------------------------------------------------
! called with flag_pvap
! output filename(1) ="pvap_atm_distribution_ppbC_lifetime"
! output filename(2) ="pvap_atm_distribution_ug_lifetime"
!-------------------------------------------------------------------

  ! optional resolution of the pvap spectrum, in log10(atm) units
  ! default is 2
    integer, optional :: opt_resolution
    integer   :: resolution
    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 3
    INTEGER,PARAMETER :: ncols = nphase+3
    INTEGER,PARAMETER :: ind_lowerlim = ncols-2
    INTEGER,PARAMETER :: ind_upperlim = ncols-1
    INTEGER,PARAMETER :: ind_time = ncols
    real, dimension(:,:), allocatable   :: pvap_distribution_ppbC, pvap_distribution_ug
    character*(header_length), dimension(:), allocatable   :: header_ppbC, header_ug
    integer                             :: nbins, ibin, ilt, iphase, minpvap, maxpvap, nlt
    integer                             :: itime
    integer, dimension(:), allocatable  :: upperlim, lowerlim
    real                                :: ppbfac, ugfac
    logical, dimension(:), allocatable  :: spec_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(spec_mask(numsp))

    if (present(opt_resolution)) then
      resolution = opt_resolution
    else
      resolution = 2
    endif

    nlt = size(efold_indices,1)

    write(6,*) '-- calculating pvap distribution by lifetime --'
    ! compute number of bins
    minpvap = floor(log10(minval(species%pvap%pvap_atm_298,  &
                     mask = species%pvap%pvap_atm_298 > 0)))
    maxpvap = ceiling(log10(maxval(species%pvap%pvap_atm_298)))
   ! constrain bin center to be EVEN if resolution = 2
    if(resolution.EQ.2.AND.MOD(minpvap,2).EQ.0)minpvap = minpvap -1
    if(resolution.EQ.2.AND.MOD(maxpvap,2).EQ.0)maxpvap = maxpvap +1

    nbins = floor(real((maxpvap - minpvap)/resolution, kind = 8))
! NB: formula truncates top bin if (maxpvap-minpvap)/resolution is non-integral
    if(MOD((maxpvap - minpvap),resolution).gt.0) nbins = nbins+1

    allocate(upperlim(nbins), lowerlim(nbins), &
             pvap_distribution_ppbC(nlt*nbins, ncols), &
             pvap_distribution_ug(nlt*nbins, ncols), &
             header_ppbC(ncols), header_ug(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header_ppbC(iphase) = "CO+CO2 [ppbC]"
          header_ug(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          header_ug(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header_ppbC(iphase) = "Wall [ppbC]"
          header_ug(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header_ppbC(iphase) = "Dimers [ppbC]"
          header_ug(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header_ppbC(iphase) = "Spinup Species [ppbC]"
          header_ug(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header_ppbC(iphase) = "Inorganics [ppbC]"
          header_ug(iphase) = "Inorganics [ug m-3]"
      end select
    enddo

    header_ppbC(ind_lowerlim) = "Pvap bin lower limit [atm]"
    header_ppbC(ind_upperlim) = "Pvap bin upper limit [atm]"
    header_ppbC(ind_time  ) = "Time [precursor lifetimes]"
    header_ug(ind_lowerlim) = "Pvap bin lower limit [atm]"
    header_ug(ind_upperlim) = "Pvap bin upper limit [atm]"
    header_ug(ind_time) = "Time [precursor lifetimes]"


    do ibin = 1, nbins
      lowerlim(ibin) = minpvap + resolution*(ibin-1)
      upperlim(ibin) = minpvap + resolution*ibin
    enddo
    ugfac = 1.660578e-12

    pvap_distribution_ppbC = 0.
    pvap_distribution_ug = 0.
    do ilt=1,nlt
       itime = efold_indices(ilt)
       pvap_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_time) = ilt
       pvap_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       pvap_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim
       
       pvap_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_time) = ilt
       pvap_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       pvap_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim

       ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
       do iphase = 1, nphase
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%pvap%pvap_atm_298) >= lowerlim(ibin) .and. &
                     log10(species%pvap%pvap_atm_298) < upperlim(ibin)
           pvap_distribution_ppbC(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%nc/ppbfac, &
                   mask = spec_mask)
           pvap_distribution_ug(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%molw * ugfac, &
                   mask = spec_mask)
         enddo
      enddo
    enddo
    filename="pvap_atm_distribution_ppbC_lifetime"
    call write_2D_array(pvap_distribution_ppbC, header_ppbC , filename)
    filename="pvap_atm_distribution_ug_lifetime"
    call write_2D_array(pvap_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               pvap_distribution_ppbC, pvap_distribution_ug, &
               header_ppbC, header_ug)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_pvap_distribution_lifetime
!-------------------------------------------------------------------
  subroutine calc_cstar_distribution_lifetime(opt_resolution)
!-------------------------------------------------------------------
! called with flag_cstar
! output filename="cstar_atm_distribution_ppbC_lifetime"
! output filename="cstar_atm_distribution_ug_lifetime"
!-------------------------------------------------------------------

  ! optional resolution of the cstar spectrum, in log10(C*) units
  ! default is 2
    integer, optional :: opt_resolution
    integer   :: resolution
    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 3
    INTEGER,PARAMETER :: ncols = nphase+3
    INTEGER,PARAMETER :: ind_lowerlim = ncols-2
    INTEGER,PARAMETER :: ind_upperlim = ncols-1
    INTEGER,PARAMETER :: ind_time = ncols
    real, dimension(:,:), allocatable   :: cstar_distribution_ppbC, cstar_distribution_ug
    character*(header_length), dimension(:), allocatable   :: header_ppbC, header_ug
    integer                             :: nbins, ibin, ilt, iphase, mincstar, maxcstar, nlt
    integer                             :: itime, ispec
    integer, dimension(:), allocatable  :: upperlim, lowerlim
    real                                :: ppbfac, ugfac
    logical, dimension(:), allocatable  :: spec_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(spec_mask(numsp))

    if (present(opt_resolution)) then
      resolution = opt_resolution
    else
      resolution = 2
    endif

    nlt = size(efold_indices,1)

    write(6,*) '-- calculating cstar  --'
    do ispec = 1, numsp
      if(species(ispec)%pvap%pvap_atm_298.gt.0)then
           species(ispec)%pvap%pvap_cstar_298 = &
            (1e6*species(ispec)%molw*species(ispec)%pvap%pvap_atm_298)/ &
            (8.2e-5*298)
      endif
    enddo

    write(6,*) '-- calculating cstar distribution by lifetime --'
    ! compute number of bins
    mincstar = floor(log10(minval(species%pvap%pvap_cstar_298,  &
                     mask = species%pvap%pvap_cstar_298 > 0)))
    maxcstar = ceiling(log10(maxval(species%pvap%pvap_cstar_298)))
! constrain bin center to be EVEN if resolution = 2
    if(resolution.EQ.2.AND.MOD(mincstar,2).EQ.0)mincstar = mincstar -1
    if(resolution.EQ.2.AND.MOD(maxcstar,2).EQ.0)maxcstar = maxcstar +1

    nbins = floor(real((maxcstar - mincstar)/resolution, kind = 8))
! NB: formula truncates top bin if (maxCstar-minCstar)/resolution is non-integral
    if(MOD((maxcstar - mincstar),resolution).gt.0) nbins = nbins+1


    allocate(upperlim(nbins), lowerlim(nbins), &
             cstar_distribution_ppbC(nlt*nbins, ncols), &
             cstar_distribution_ug(nlt*nbins, ncols), &
             header_ppbC(ncols), header_ug(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header_ppbC(iphase) = "CO+CO2 [ppbC]"
          header_ug(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          header_ug(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header_ppbC(iphase) = "Wall [ppbC]"
          header_ug(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header_ppbC(iphase) = "Dimers [ppbC]"
          header_ug(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header_ppbC(iphase) = "Spinup Species [ppbC]"
          header_ug(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header_ppbC(iphase) = "Inorganics [ppbC]"
          header_ug(iphase) = "Inorganics [ug m-3]"
        case default
          header_ppbC(iphase) = ""
          header_ug(iphase) = ""
      end select
    enddo

    header_ppbC(ind_lowerlim) = "Cstar bin lower limit [log_10(ug.m-3)]"
    header_ppbC(ind_upperlim) = "Cstar bin upper limit [log_10(ug.m-3)]"
    header_ppbC(ind_time )= "Time [precursor lifetimes]"
    header_ug(ind_lowerlim) = "Cstar bin lower limit [log_10(ug.m-3)]"
    header_ug(ind_upperlim) = "Cstar bin upper limit [log_10(ug.m-3)]"
    header_ug(ind_time )= "Time [precursor lifetimes]"


    do ibin = 1, nbins
      lowerlim(ibin) = mincstar + resolution*(ibin-1)
      upperlim(ibin) = mincstar + resolution*ibin
    enddo
    ugfac = 1.660578e-12

    cstar_distribution_ppbC = 0.
    cstar_distribution_ug = 0.
    do ilt=1,nlt
       itime = efold_indices(ilt)
       cstar_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_time )= ilt
       cstar_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       cstar_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim
       
       cstar_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_time )= ilt
       cstar_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       cstar_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim

       ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
       do iphase = 1, nphase
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%pvap%pvap_cstar_298) >= lowerlim(ibin) .and. &
                     log10(species%pvap%pvap_cstar_298) < upperlim(ibin)
           cstar_distribution_ppbC(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%nc/ppbfac, &
                   mask = spec_mask)
           cstar_distribution_ug(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%molw * ugfac, &
                   mask = spec_mask)
         enddo
      enddo
    enddo
    filename="cstar_atm_distribution_ppbC_lifetime"
    call write_2D_array(cstar_distribution_ppbC, header_ppbC , filename)
    filename="cstar_atm_distribution_ug_lifetime"
    call write_2D_array(cstar_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               cstar_distribution_ppbC, cstar_distribution_ug, &
               header_ppbC, header_ug)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_cstar_distribution_lifetime
!-----------------------------------------------------------------------
  subroutine calc_henry_distribution(opt_resolution)
!-----------------------------------------------------------------------
! called with flag_henry (unless commented out in main.f90)
! output filename(1) = "kh_Matm_distribution_ppbC_time"
! output filename(2) = "kh_Matm_distribution_ug_time"
!-----------------------------------------------------------------------

  ! optional resolution of the Kh spectrum, in log10(M atm-1) units
  ! default is 2
    integer, optional :: opt_resolution
    integer   :: resolution
    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 3
    INTEGER,PARAMETER :: ncols = nphase+3
    INTEGER,PARAMETER :: ind_lowerlim=ncols-2
    INTEGER,PARAMETER :: ind_upperlim=ncols-1
    INTEGER,PARAMETER :: ind_time=ncols
    real, dimension(:,:), allocatable   :: kh_distribution_ppbC, kh_distribution_ug
    character*(header_length), dimension(:), allocatable   :: header_ppbC, header_ug
    integer                             :: nbins, ibin, i, iphase, minkh, maxkh
    integer, dimension(:), allocatable  :: upperlim, lowerlim
    real                                :: ppbfac, ugfac
    logical, dimension(:), allocatable  :: spec_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(spec_mask(numsp))

    if (present(opt_resolution)) then
      resolution = opt_resolution
    else
      resolution = 2
    endif

    write(6,*) '-- calculating henry distribution --'
    ! compute number of bins
    minkh = floor(log10(minval(species%henry%Keff,  &
                     mask = species%henry%Keff > 0)))
    maxkh = ceiling(log10(maxval(species%henry%Keff)))
! constrain bin center to be EVEN if resolution = 2
    if(resolution.EQ.2.AND.MOD(minkh,2).EQ.0)minkh = minkh -1
    if(resolution.EQ.2.AND.MOD(maxkh,2).EQ.0)maxkh = maxkh +1

    nbins = floor(real((maxkh - minkh)/resolution, kind = 8))
! NB: formula truncates top bin if (maxkh-minkh)/resolution is non-integral
    if(MOD((maxkh - minkh),resolution).gt.0) nbins = nbins+1


    allocate(upperlim(nbins), lowerlim(nbins), &
             kh_distribution_ppbC(ndat*nbins, ncols), &
             kh_distribution_ug(ndat*nbins, ncols), &
             header_ppbC(ncols), header_ug(ncols))


    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header_ppbC(iphase) = "CO+CO2 [ppbC]"
          header_ug(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          header_ug(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header_ppbC(iphase) = "Wall [ppbC]"
          header_ug(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header_ppbC(iphase) = "Dimers [ppbC]"
          header_ug(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header_ppbC(iphase) = "Spinup Species [ppbC]"
          header_ug(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header_ppbC(iphase) = "Inorganics [ppbC]"
          header_ug(iphase) = "Inorganics [ug m-3]"
        case default
          header_ppbC(iphase) = ""
          header_ug(iphase) = ""
      end select
    enddo

    header_ppbC(ind_lowerlim) = "Kh bin lower limit [log10(M/atm)]"
    header_ppbC(ind_upperlim) = "Kh bin upper limit [log10(M/atm)]"
    header_ppbC(ind_time) = "Time [s]"
    header_ug(ind_lowerlim) = "Kh bin lower limit [log10(M/atm)]"
    header_ug(ind_upperlim) = "Kh bin upper limit [log10(M/atm)]"
    header_ug(ind_time) = "Time [s]"


    do ibin = 1, nbins
      lowerlim(ibin) = minkh +resolution*(ibin-1)
      upperlim(ibin) = minkh + resolution*ibin
    enddo
    ugfac = 1.660578e-12

    kh_distribution_ppbC = 0.
    kh_distribution_ug = 0.
    do i=1,ndat
       kh_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_time) = time(i)
       kh_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_upperlim) = upperlim
       kh_distribution_ppbC((i-1)*nbins+1:i*nbins,ind_lowerlim) = lowerlim
       kh_distribution_ug((i-1)*nbins+1:i*nbins,ind_time) = time(i)
       kh_distribution_ug((i-1)*nbins+1:i*nbins,ind_upperlim) = upperlim
       kh_distribution_ug((i-1)*nbins+1:i*nbins,ind_lowerlim) = lowerlim

       ppbfac = 7.245461056e+09*pressure(i)/temperature(i)
       do iphase = 1, nphase
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%henry%Keff) >= lowerlim(ibin) .and. &
                     log10(species%henry%Keff) < upperlim(ibin)
           kh_distribution_ppbC(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%nc/ppbfac, &
              mask = spec_mask)
           kh_distribution_ug(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%molw * ugfac, &
              mask = spec_mask)
         enddo
      enddo
    enddo
    filename="kh_Matm_distribution_ppbC_time"
    call write_2D_array(kh_distribution_ppbC, header_ppbC , filename)
    filename="kh_Matm_distribution_ug_time"
    call write_2D_array(kh_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               kh_distribution_ppbC, kh_distribution_ug, &
               header_ppbC, header_ug)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_henry_distribution
!------------------------------------------------------------------------
  subroutine calc_henry_distribution_lifetime(opt_resolution)
!------------------------------------------------------------------------
! called with flag_henry (unless commented out in main.f90)
! output filename(1) = "Henry_Matm_distribution_ppbC_lifetime"
! output filename(2) = "Henry_Matm_distribution_ug_lifetime"
!------------------------------------------------------------------------

  ! optional resolution of the Kh spectrum, in log10(M atm-1) units
  ! default is 2
    integer, optional :: opt_resolution
    integer   :: resolution
    !INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: nphase = 3
    INTEGER,PARAMETER :: ncols = nphase+3
    INTEGER,PARAMETER :: ind_lowerlim = ncols-2
    INTEGER,PARAMETER :: ind_upperlim = ncols-1
    INTEGER,PARAMETER :: ind_time = ncols
    real, dimension(:,:), allocatable   :: kh_distribution_ppbC, kh_distribution_ug
    character*(header_length), dimension(:), allocatable   :: header_ppbC, header_ug
    integer                             :: nbins, ibin, ilt, iphase, minkh, maxkh
    integer                             :: nlt, itime
    integer, dimension(:), allocatable  :: upperlim, lowerlim
    real                                :: ppbfac, ugfac
    logical, dimension(:), allocatable  :: spec_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate(spec_mask(numsp))

    if (present(opt_resolution)) then
      resolution = opt_resolution
    else
      resolution = 2
    endif

    write(6,*) '-- calculating henry distribution by lifetime --'
    ! compute number of bins
    minkh = floor(log10(minval(species%henry%Keff,  &
                     mask = species%henry%Keff > 0)))
    maxkh = ceiling(log10(maxval(species%henry%Keff)))
   ! constrain bin center to be EVEN if resolution = 2
    if(resolution.EQ.2.AND.MOD(minkH,2).EQ.0)minkH = minkH -1
    if(resolution.EQ.2.AND.MOD(maxkH,2).EQ.0)maxkH = maxkH +1

    nbins = floor(real((maxkh - minkh)/resolution, kind = 8))
! NB: formula truncates top bin if (maxkh-minkh)/resolution is non-integral
    if(MOD((maxkh - minkh),resolution).gt.0) nbins = nbins+1

    nlt = size(efold_indices,1)

    allocate(upperlim(nbins), lowerlim(nbins), &
             kh_distribution_ppbC(nlt*nbins, ncols), &
             kh_distribution_ug(nlt*nbins, ncols), &
             header_ppbC(ncols), header_ug(ncols))


    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header_ppbC(iphase) = "CO+CO2 [ppbC]"
          header_ug(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          header_ug(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header_ppbC(iphase) = "Wall [ppbC]"
          header_ug(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header_ppbC(iphase) = "Dimers [ppbC]"
          header_ug(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header_ppbC(iphase) = "Spinup Species [ppbC]"
          header_ug(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header_ppbC(iphase) = "Inorganics [ppbC]"
          header_ug(iphase) = "Inorganics [ug m-3]"
      end select
    enddo

    header_ppbC(ind_lowerlim)   = "Kh bin lower limit [atm]"
    header_ppbC(ind_upperlim) = "Kh bin upper limit [atm]"
    header_ppbC(ind_time) = "Time [precursor lifetime]"
    header_ug(ind_lowerlim)   = "Kh bin lower limit [atm]"
    header_ug(ind_upperlim) = "Kh bin upper limit [atm]"
    header_ug(ind_time) = "Time [precursor lifetime]"


    do ibin = 1, nbins
      lowerlim(ibin) = minkh +resolution*(ibin-1)
      upperlim(ibin) = minkh + resolution*ibin
    enddo
    ugfac = 1.660578e-12

    kh_distribution_ppbC = 0.
    kh_distribution_ug = 0.
    do ilt=1,nlt
       itime = efold_indices(ilt)
       kh_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_time) = ilt
       kh_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       kh_distribution_ppbC((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim
       
       kh_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_time) = ilt
       kh_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_upperlim) = upperlim
       kh_distribution_ug((ilt-1)*nbins+1:ilt*nbins,ind_lowerlim) = lowerlim

       ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
       do iphase = 1, nphase
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%henry%Keff) >= lowerlim(ibin) .and. &
                     log10(species%henry%Keff) < upperlim(ibin)
           kh_distribution_ppbC(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%nc/ppbfac, &
                   mask = spec_mask )
           kh_distribution_ug(nbins*(ilt-1)+ibin, iphase) =  &
              sum( concentrations(itime, :)*species%molw * ugfac, &
                   mask = spec_mask )
         enddo
      enddo
    enddo
    filename="Henry_Matm_distribution_ppbC_lifetime"
    call write_2D_array(kh_distribution_ppbC, header_ppbC , filename)
    filename="Henry_Matm_distribution_ug_lifetime"
    call write_2D_array(kh_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               kh_distribution_ppbC, kh_distribution_ug, &
               header_ppbC, header_ug)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_henry_distribution_lifetime
!-------------------------------------------------------------------
  subroutine calc_top_sp_fgrp()
!-------------------------------------------------------------------
! EXPERIMENTAL !
! calculates top species with a PARTICULAR FUNCTIONAL GROUP
! called with flag_topspec (WITHOUT flag_ppbc)
! output filename(1a) = "top_species_molec_gas_time"
! output filename(1b) = "top_species_molec_aer_time"
! if also called with flag_ppb
! output filename(3a) = "top_species_ppb_gas_time"
! output filename(3b) = "top_species_ppb_aer_time"
!-------------------------------------------------------------------

    integer, dimension(:), allocatable  :: top_species_indx
    character*(maxlsp), dimension(:), allocatable :: top_species_names
    real, dimension(:,:), allocatable   :: top_species_molec
    real, dimension(:,:), allocatable   :: top_species_ppb
    character*(header_length) :: smile
    character*(header_length), dimension(:), allocatable :: header1,header2, &
                                                            header3,header4
    real, dimension(:), allocatable     :: time_integrated_sum
    real                                :: ppbfac
    integer                             :: N, itime, itop, iphase, ispec, ifl
    logical                             :: found_fg

    character(len=1),parameter :: target_func_code = "3"

    CHARACTER*(10)    :: date,time1,time2
    CHARACTER*(filenames_length) :: filename1,filename2,filename3,filename4
    CALL date_and_time(date,time1)

    N = n_topspecies
    allocate(top_species_names(N), top_species_indx(N), &
             top_species_molec(ndat, N +2), &
             top_species_ppb(ndat, N +2), &
             time_integrated_sum(numsp), &
             header1(N+2),header2(N+2),header3(N+2),header4(N+2))


    write(6,*) '-- finding top ', N, ' species (molec/cc) with chosen fgrp --'
    if (flag_ppb) &
    write(6,*) '-- finding top ', N, ' species (ppb) with chosen fgrp --'
    ! start in the gas phase
    ! first we need to find the topN-1 species (Nth species is "Others")
    ! with respect to time integrated mixing ratio in molec/cc (ppb)

    do iphase = GAS, AEROSOL, AEROSOL-GAS

      top_species_molec = 0
      top_species_ppb = 0
      time_integrated_sum = 0
      
      do itime = 1, ndat
        do ispec = 1, numsp
          do ifl = 1, lfl
            if (species(ispec)%functions(ifl:ifl) == target_func_code) then
              if(species(ispec)%phase == iphase) then
                time_integrated_sum(ispec) = time_integrated_sum(ispec)+ &
                concentrations(itime, ispec)
              endif
            endif
          enddo
        enddo
      enddo

      top_species_indx(1) = maxloc(time_integrated_sum, dim = 1)
      top_species_names(1) = species(top_species_indx(1))%code
      time_integrated_sum(top_species_indx(1)) = 0
      do itop = 2, N
        top_species_indx(itop) = maxloc(time_integrated_sum, dim = 1)
        top_species_names(itop) = species(top_species_indx(itop))%code
        time_integrated_sum(top_species_indx(itop)) = 0
      enddo

      do itop = 1, N
        header1(itop) = trim(species(top_species_indx(itop))%formula)// &
                   ":"//trim(species(top_species_indx(itop))%code)
        header3(itop) = header1(itop)
        if (flag_smiles) then
          CALL smiles(species(top_species_indx(itop))%formula,smile)
          ispec=INDEX(smile," ")-1
          header2(itop) = trim(smile(1:ispec))// &
                     ":"//trim(species(top_species_indx(itop))%code)
          header4(itop) = header2(itop)
        endif
      enddo
      header1(N+1) = "Others [molec/cc]"
      header1(N+2) = "Time [s]"
      header3(N+1) = "Others [ppb]"
      header3(N+2) = "Time [s]"
      if (flag_smiles) then
        header2(N+1) = "Others [molec/cc]"
        header2(N+2) = "Time [s]"
        header4(N+1) = "Others [ppb]"
        header4(N+2) = "Time [s]"
      endif

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec=  1, numsp
          if (species(ispec)%phase /= iphase) CYCLE
          do ifl = 1, lfl
            if (species(ispec)%functions(ifl:ifl) == target_func_code) then
              found_fg = .false.
              do itop = 1, N
                if (ispec == top_species_indx(itop)) then
                  found_fg = .true.
                  top_species_molec(itime,itop) = concentrations(itime, ispec)
                  top_species_ppb(itime,itop) = concentrations(itime, ispec) &
                                              / ppbfac
                  exit
                endif
              enddo
              if (.not. found_fg) then
                top_species_molec(itime, N+1) = top_species_molec(itime, N+1) + &
                                      concentrations(itime, ispec)
               top_species_ppb(itime, N+1) = top_species_ppb(itime, N+1) + &
                                      concentrations(itime, ispec) / ppbfac
              endif
            endif
         enddo
        enddo
      enddo

      top_species_molec(:,N+2) = time
      top_species_ppb(:,N+2) = time

! chemical formulae headers
      if (iphase == GAS) then
        filename1 = "top_sp_fgrp_molec_gas_time"
        filename3 = "top_sp_fgrp_ppb_gas_time"
      else if (iphase == AEROSOL) then
        filename1 = "top_sp_fgrp_molec_aer_time"
        filename3 = "top_sp_fgrp_ppb_aer_time"
      endif
      call write_2D_array(top_species_molec, header1 , filename1)
      if(flag_ppb) &
      call write_2D_array(top_species_ppb, header3 , filename3)
! SMILES headers
      if (flag_smiles) then
        if (iphase == GAS) then
          filename2 = "top_smiles_molec_gas_time"
          filename4 = "top_smiles_ppb_gas_time"
        else if (iphase == AEROSOL) then
          filename2 = "top_smiles_molec_aer_time"
          filename4 = "top_smiles_ppb_aer_time"
        endif
        call write_2D_array(top_species_molec, header2 , filename2)
        if(flag_ppb) &
        call write_2D_array(top_species_ppb, header4 , filename4)
      endif
    enddo


    deallocate(top_species_names, top_species_indx, &
               top_species_molec, top_species_ppb, &
               time_integrated_sum, header1,header2,header3,header4)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_top_sp_fgrp
!-------------------------------------------------------------------
  subroutine calc_top_species()
!-------------------------------------------------------------------
! called with flag_topspec (unless flag_ppbc is active)
! output filename(1a) = "top_species_molec_gas_time"
! output filename(1b) = "top_species_molec_aer_time"
! if also called with flag_smiles
! output filename(2a) = "top_smiles_molec_gas_time"
! output filename(2b) = "top_smiles_molec_aer_time"
! if also called with flag_ppb
! output filename(3a) = "top_species_ppb_gas_time"
! output filename(3b) = "top_species_ppb_aer_time"
! if also called with flag_ppb AND flag_smiles
! output filename(4a) = "top_smiles_ppb_gas_time"
! output filename(4b) = "top_smiles_ppb_aer_time"
!-------------------------------------------------------------------

    integer, dimension(:), allocatable  :: top_species_indx
    character*(maxlsp), dimension(:), allocatable :: top_species_names
    real, dimension(:,:), allocatable   :: top_species_molec
    real, dimension(:,:), allocatable   :: top_species_ppb
    character*(header_length) :: smile
    character*(header_length), dimension(:), allocatable :: header1,header2, &
                                                            header3,header4
    real, dimension(:), allocatable     :: time_integrated_sum
    real                                :: ppbfac
    integer                             :: N, itime, itop, iphase, ispec
    logical                             :: found_fg

    CHARACTER*(10)    :: date,time1,time2
    CHARACTER*(filenames_length) :: filename1,filename2,filename3,filename4
    CALL date_and_time(date,time1)

    N = n_topspecies
    allocate(top_species_names(N), top_species_indx(N), &
             top_species_molec(ndat, N +2), &
             top_species_ppb(ndat, N +2), &
             time_integrated_sum(numsp), &
             header1(N+2),header2(N+2),header3(N+2),header4(N+2))


    write(6,*) '-- finding top ', N, ' species (molec/cc) --'
    if (flag_ppb) &
    write(6,*) '-- finding top ', N, ' species (ppb) --'
    ! start in the gas phase
    ! first we need to find the topN-1 species (Nth species is "Others")
    ! with respect to time integrated mixing ratio in molec/cc (ppb)

    do iphase = GAS, AEROSOL, AEROSOL-GAS

      top_species_molec = 0
      top_species_ppb = 0
      time_integrated_sum = 0
      
      do itime = 1, ndat
        do ispec = 1, numsp
          if(species(ispec)%phase == iphase) then
            time_integrated_sum(ispec) = time_integrated_sum(ispec)+ &
              concentrations(itime, ispec)
          endif
        enddo
      enddo

      top_species_indx(1) = maxloc(time_integrated_sum, dim = 1)
      top_species_names(1) = species(top_species_indx(1))%code
      time_integrated_sum(top_species_indx(1)) = 0
      do itop = 2, N
        top_species_indx(itop) = maxloc(time_integrated_sum, dim = 1)
        top_species_names(itop) = species(top_species_indx(itop))%code
        time_integrated_sum(top_species_indx(itop)) = 0
      enddo

      do itop = 1, N
        header1(itop) = trim(species(top_species_indx(itop))%formula)// &
                   ":"//trim(species(top_species_indx(itop))%code)
        header3(itop) = header1(itop)
        if (flag_smiles) then
          CALL smiles(species(top_species_indx(itop))%formula,smile)
          ispec=INDEX(smile," ")-1
          header2(itop) = trim(smile(1:ispec))// &
                     ":"//trim(species(top_species_indx(itop))%code)
          header4(itop) = header2(itop)
        endif
      enddo
      header1(N+1) = "Others [molec/cc]"
      header1(N+2) = "Time [s]"
      header3(N+1) = "Others [ppb]"
      header3(N+2) = "Time [s]"
      if (flag_smiles) then
        header2(N+1) = "Others [molec/cc]"
        header2(N+2) = "Time [s]"
        header4(N+1) = "Others [ppb]"
        header4(N+2) = "Time [s]"
      endif

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec=  1, numsp
          if (species(ispec)%phase /= iphase) CYCLE
          found_fg = .false.
          do itop = 1, N
            if (ispec == top_species_indx(itop)) then
              found_fg = .true.
              top_species_molec(itime,itop) = concentrations(itime, ispec)
              top_species_ppb(itime,itop) = concentrations(itime, ispec) &
                                          / ppbfac
              exit
            endif
          enddo
          if (.not. found_fg) then
            top_species_molec(itime, N+1) = top_species_molec(itime, N+1) + &
                                  concentrations(itime, ispec)
            top_species_ppb(itime, N+1) = top_species_ppb(itime, N+1) + &
                                  concentrations(itime, ispec) / ppbfac
          endif
        enddo
      enddo

      top_species_molec(:,N+2) = time
      top_species_ppb(:,N+2) = time

! chemical formulae headers
      if (iphase == GAS) then
        filename1 = "top_species_molec_gas_time"
        filename3 = "top_species_ppb_gas_time"
      else if (iphase == AEROSOL) then
        filename1 = "top_species_molec_aer_time"
        filename3 = "top_species_ppb_aer_time"
      endif
      call write_2D_array(top_species_molec, header1 , filename1)
      if(flag_ppb) &
      call write_2D_array(top_species_ppb, header3 , filename3)
! SMILES headers
      if (flag_smiles) then
        if (iphase == GAS) then
          filename2 = "top_smiles_molec_gas_time"
          filename4 = "top_smiles_ppb_gas_time"
        else if (iphase == AEROSOL) then
          filename2 = "top_smiles_molec_aer_time"
          filename4 = "top_smiles_ppb_aer_time"
        endif
        call write_2D_array(top_species_molec, header2 , filename2)
        if(flag_ppb) &
        call write_2D_array(top_species_ppb, header4 , filename4)
      endif
    enddo


    deallocate(top_species_names, top_species_indx, &
               top_species_molec, top_species_ppb, &
               time_integrated_sum, header1,header2,header3,header4)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_top_species
!-------------------------------------------------------------------
  subroutine calc_top_species_ppbc()
!-------------------------------------------------------------------
! called with flag_topspec (unless commented out in main.f90)
! output filename(1a) = "top_species_ppbC_gas_time"
! output filename(1b) = "top_species_ppbC_aer_time"
! if also called with flag_smiles
! output filename(2a) = "top_smiles_ppbC_gas_time"
! output filename(2b) = "top_smiles_ppbC_aer_time"
!-------------------------------------------------------------------

    integer, dimension(:), allocatable  :: top_species_indx
    character*(maxlsp), dimension(:), allocatable :: top_species_names
    real, dimension(:,:), allocatable   :: top_species_concentrations
    character*(header_length) :: smile
    character*(header_length), dimension(:), allocatable   :: header,header2
    real, dimension(:), allocatable     :: time_integrated_sum
    real                                :: ppbfac

    integer                             :: N, itime, itop, iphase, ispec
    logical                             :: found_fg

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    N = n_topspecies
    allocate(top_species_names(N), top_species_indx(N), &
             top_species_concentrations(ndat, N +2), &
             time_integrated_sum(numsp), header(N+2),header2(N+2))


    write(6,*) '-- finding top ', N, ' species (ppbC) --'
    ! start in the gas phase
    ! first we need to find the topN-1 species (Nth species is "Others")
    ! with respect to time integrated mixing ratio in ppbC

    do iphase = GAS, AEROSOL, AEROSOL-GAS

      top_species_concentrations=0
      time_integrated_sum = 0
      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec = 1, numsp
          if(species(ispec)%phase == iphase) then
            time_integrated_sum(ispec) = time_integrated_sum(ispec)+ &
              concentrations(itime, ispec)*species(ispec)%nc/ppbfac
          endif
        enddo
      enddo

      top_species_indx(1) = maxloc(time_integrated_sum, dim = 1)
      top_species_names(1) = species(top_species_indx(1))%code
      time_integrated_sum(top_species_indx(1)) = 0
      do itop = 2, N
        top_species_indx(itop) = maxloc(time_integrated_sum, dim = 1)
        top_species_names(itop) = species(top_species_indx(itop))%code
        time_integrated_sum(top_species_indx(itop)) = 0
      enddo

      do itop = 1, N
        header(itop) = trim(species(top_species_indx(itop))%formula)// &
          ":"//trim(species(top_species_indx(itop))%code)
        if (flag_smiles) then
          CALL smiles(species(top_species_indx(itop))%formula,smile)
          ispec=INDEX(smile," ")-1
          header2(itop) = trim(smile(1:ispec))// &
            ":"//trim(species(top_species_indx(itop))%code)
        endif
      enddo
      header(N+1) = "Others [ppbC]"
      header(N+2) = "Time [s]"
      if (flag_smiles) then
        header2(N+1) = "Others [ppbC]"
        header2(N+2) = "Time [s]"
      endif

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec=  1, numsp
          if (species(ispec)%phase /= iphase) CYCLE
          found_fg = .false.
          do itop = 1, N
            if (ispec == top_species_indx(itop)) then
              found_fg = .true.
              top_species_concentrations(itime,itop) = concentrations(itime, ispec)* &
                                  species(ispec)%nc/ppbfac
              exit
            endif
          enddo
          if (.not. found_fg) then
            top_species_concentrations(itime, N+1) = top_species_concentrations(itime, N+1) + &
                                  concentrations(itime, ispec)*species(ispec)%nc/ppbfac
          endif
        enddo
      enddo

      top_species_concentrations(:,N+2) = time
! chemical formulae headers
      if (iphase == GAS) then
        filename = "top_species_ppbc_gas_time"
      else if (iphase == AEROSOL) then
        filename = "top_species_ppbc_aer_time"
      endif
      call write_2D_array(top_species_concentrations, header , filename)
! SMILES headers
      if (flag_smiles) then
        if (iphase == GAS) then
          filename = "top_smiles_ppbc_gas_time"
        else if (iphase == AEROSOL) then
          filename = "top_smiles_ppbc_aer_time"
        endif
        call write_2D_array(top_species_concentrations, header2 , filename)
      endif
    enddo


    deallocate(top_species_names, top_species_indx, &
               top_species_concentrations, &
               time_integrated_sum, header,header2)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_top_species_ppbc
!-------------------------------------------------------------------
  subroutine calc_top_Nspc()
!-------------------------------------------------------------------
! called with flag_topNspc
! output filename(1a) = "top_Nspc_molec_gas_time"
! output filename(1b) = "top_Nspc_molec_aer_time"
! if also called with flag_smiles
! output filename(2a) = "top_Nsmiles_molec_gas_time"
! output filename(2b) = "top_Nsmiles_molec_aer_time"
! if also called with flag_ppb
! output filename(3a) = "top_Nspc_ppb_gas_time"
! output filename(3b) = "top_Nspc_ppb_aer_time"
! if also called with flag_ppb AND flag_smiles
! output filename(4a) = "top_Nsmiles_ppb_gas_time"
! output filename(4b) = "top_Nsmiles_ppb_aer_time"
!-------------------------------------------------------------------

    integer, dimension(:), allocatable  :: top_species_indx
    character*(maxlsp), dimension(:), allocatable :: top_species_names
    real, dimension(:,:), allocatable   :: top_species_molec
    real, dimension(:,:), allocatable   :: top_species_ppb
    character*(header_length) :: smile
    character*(header_length), dimension(:), allocatable :: header1,header2, &
                                                            header3,header4
    real, dimension(:), allocatable     :: time_integrated_sum
    real                                :: ppbfac
    integer                             :: N, itime, itop, iphase, ispec
    logical                             :: found_fg

    CHARACTER*(10)    :: date,time1,time2
    CHARACTER*(filenames_length) :: filename1,filename2,filename3,filename4
    CALL date_and_time(date,time1)

    N = n_topspecies
    allocate(top_species_names(N), top_species_indx(N), &
             top_species_molec(ndat, N +2), &
             top_species_ppb(ndat, N +2), &
             time_integrated_sum(numsp), &
             header1(N+2),header2(N+2),header3(N+2),header4(N+2))


    write(6,*) '-- finding top ', N, ' N-containing species (molec/cc) --'
    if (flag_ppb) &
    write(6,*) '-- finding top ', N, ' N-containing species (ppb) --'
    ! start in the gas phase
    ! first we need to find the topN-1 species (Nth species is "Others")
    ! with respect to time integrated mixing ratio in molec/cc (ppb)

    do iphase = GAS, AEROSOL, AEROSOL-GAS

      top_species_molec = 0
      top_species_ppb = 0
      time_integrated_sum = 0
      
      do itime = 1, ndat
        do ispec = 1, numsp
          if(species(ispec)%phase == iphase .AND. &
             species(ispec)%nn .GT. 0) then
            time_integrated_sum(ispec) = time_integrated_sum(ispec)+ &
              concentrations(itime, ispec)
          endif
        enddo
      enddo

      top_species_indx(1) = maxloc(time_integrated_sum, dim = 1)
      top_species_names(1) = species(top_species_indx(1))%code
      time_integrated_sum(top_species_indx(1)) = 0
      do itop = 2, N
        top_species_indx(itop) = maxloc(time_integrated_sum, dim = 1)
        top_species_names(itop) = species(top_species_indx(itop))%code
        time_integrated_sum(top_species_indx(itop)) = 0
      enddo

      do itop = 1, N
        header1(itop) = trim(species(top_species_indx(itop))%formula)// &
                   ":"//trim(species(top_species_indx(itop))%code)
        header3(itop) = header1(itop)
        if (flag_smiles) then
          CALL smiles(species(top_species_indx(itop))%formula,smile)
          ispec=INDEX(smile," ")-1
          header2(itop) = trim(smile(1:ispec))// &
                     ":"//trim(species(top_species_indx(itop))%code)
          header4(itop) = header2(itop)
        endif
      enddo
      header1(N+1) = "Others [molec/cc]"
      header1(N+2) = "Time [s]"
      header3(N+1) = "Others [ppb]"
      header3(N+2) = "Time [s]"
      if (flag_smiles) then
        header2(N+1) = "Others [molec/cc]"
        header2(N+2) = "Time [s]"
        header4(N+1) = "Others [ppb]"
        header4(N+2) = "Time [s]"
      endif

      do itime = 1, ndat
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        do ispec=  1, numsp
          if (species(ispec)%phase /= iphase .OR. &
              species(ispec)%nn .eq. 0) CYCLE
          found_fg = .false.
          do itop = 1, N
            if (ispec == top_species_indx(itop)) then
              found_fg = .true.
              top_species_molec(itime,itop) = concentrations(itime, ispec)
              top_species_ppb(itime,itop) = concentrations(itime, ispec) &
                                          / ppbfac
              exit
            endif
          enddo
 ! "others" should ONLY include N-containing species
          if (.not. found_fg .and. species(ispec)%nn .gt. 0) then
            top_species_molec(itime, N+1) = top_species_molec(itime, N+1) + &
                                  concentrations(itime, ispec)
            top_species_ppb(itime, N+1) = top_species_ppb(itime, N+1) + &
                                  concentrations(itime, ispec) / ppbfac
          endif
        enddo
      enddo

      top_species_molec(:,N+2) = time
      top_species_ppb(:,N+2) = time

! chemical formulae headers
      if (iphase == GAS) then
        filename1 = "top_Nspc_molec_gas_time"
        filename3 = "top_Nspc_ppb_gas_time"
      else if (iphase == AEROSOL) then
        filename1 = "top_Nspc_molec_aer_time"
        filename3 = "top_Nspc_ppb_aer_time"
      endif
      call write_2D_array(top_species_molec, header1 , filename1)
      if(flag_ppb) &
      call write_2D_array(top_species_ppb, header3 , filename3)
! SMILES headers
      if (flag_smiles) then
        if (iphase == GAS) then
          filename2 = "top_Nsmiles_molec_gas_time"
          filename4 = "top_Nsmiles_ppb_gas_time"
        else if (iphase == AEROSOL) then
          filename2 = "top_Nsmiles_molec_aer_time"
          filename4 = "top_Nsmiles_ppb_aer_time"
        endif
        call write_2D_array(top_species_molec, header2 , filename2)
        if(flag_ppb) &
        call write_2D_array(top_species_ppb, header4 , filename4)
      endif
    enddo

    deallocate(top_species_names, top_species_indx, &
               top_species_molec, top_species_ppb, &
               time_integrated_sum, header1,header2,header3,header4)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_top_Nspc
!-------------------------------------------------------------------
  subroutine calc_contributing_species()
!-------------------------------------------------------------------
! called with flag_contributingspecs
! output filename(1) = "contributing_species_molec_gas_time"
! output filename(2) = "contributing_species_molec_aer_time"
!-------------------------------------------------------------------

    real, dimension(:,:), allocatable  :: contributing_species
    character*(header_length), dimension(:), allocatable :: header
    integer                            :: itime, ispec, iphase, itop
    real, dimension(:), allocatable    :: temp_conc
    real                               :: cumul, total
    logical                            :: fg_50

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)


    write(6,*) '-- calculating number of species contributing to the majority of concentration --'

    allocate(contributing_species(ndat, 3), header(3), temp_conc(numsp))
    header(1) = "Time [s]"
    header(2) = "50% species"
    header(3) = "90% species"
    do iphase = GAS, AEROSOL, AEROSOL-GAS
      contributing_species = 0
      do itime = 2, ndat
        temp_conc = concentrations(itime, :)
        where (species%phase /= iphase)
          temp_conc = 0.
        end where
        total = sum(temp_conc)
        fg_50 = .true.
        cumul = 0
        do ispec = 1, numsp
          itop = maxloc(temp_conc, dim = 1)
          cumul = cumul + temp_conc(itop)
          if(fg_50 .and. cumul > 0.5*total) then
            contributing_species(itime, 2) = ispec
            fg_50 = .false.
          endif
          if(cumul > 0.9*total) then
            contributing_species(itime, 3) = ispec
            exit
          endif
          temp_conc(itop) = 0
        enddo
      enddo

      contributing_species(:,1) = time
      if (iphase == GAS) then
        filename = "contributing_species_molec_gas_time"
      else if (iphase == AEROSOL) then
        filename = "contributing_species_molec_aer_time"
      endif
      call write_2D_array(contributing_species, header , filename)
    enddo

    deallocate(contributing_species, header,temp_conc)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_contributing_species
!-------------------------------------------------------------------
  subroutine calc_selected_species()
!-------------------------------------------------------------------
! called with flag_selected
! output filename = "selected_species_molec_time"
! output filename = "selected_species_ppb_time"
!-------------------------------------------------------------------

    real, allocatable, dimension(:,:)   ::  selected_species_conc
    real, allocatable, dimension(:,:)   ::  selected_species_ppb
    character*(header_length), dimension(:), allocatable :: header
    integer                             :: nsel, isel, ispec, itime
    integer, allocatable, dimension(:)  :: sel_indices
    real                                :: ppbfac

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    nsel = 0
    do isel = 1, mxselspe
      if(selected_species(isel) .ne. '') then
        nsel = nsel+1
      else
        exit
      endif
    enddo

    allocate(header(nsel+1), selected_species_conc(ndat, nsel +1 ),selected_species_ppb(ndat, nsel +1 ),sel_indices(nsel))

    selected_species_conc = 0.
    sel_indices = 0
    header = ''
    do isel =1, nsel
      ispec = find_species_index(selected_species(isel))
      !if (ispec .ne. -1) then
      if (ispec .gt. 0) then
        sel_indices(isel) = ispec
        selected_species_conc(:, isel) = concentrations(:, ispec)
        if(flag_ppb) then
          do itime = 1, ndat
            ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
            selected_species_ppb(itime, isel) = concentrations(itime, ispec)/ppbfac
          enddo
        endif
        if (species(ispec)%formula == "") then
          header(isel) = selected_species(isel)
        else
          header(isel) = trim(species(ispec)%formula)// &
                    ":"//trim(species(ispec)%code)
        endif
      else
        WRITE(6,*) "Selected species could not be found"
        WRITE(6,*) selected_species(isel)
      endif
    enddo

    selected_species_conc(:, nsel+1) = time
    selected_species_ppb(:, nsel+1) = time
    header(nsel+1) = "Time [s]"

    filename = "selected_species_molec_time"
    call write_2D_array(selected_species_conc, header, filename)

    if(flag_ppb)then
      filename = "selected_species_ppb_time"
      call write_2D_array(selected_species_ppb, header, filename)
    endif

    deallocate(header, selected_species_conc, selected_species_ppb, sel_indices)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_selected_species
!-------------------------------------------------------------------
  subroutine calc_species_rates(ncid)
!-------------------------------------------------------------------
! called with flag_rates
! also requires a species name assigned to "rate_species"
! output filename = "top_rates_{rate_species}_time"
!-------------------------------------------------------------------
    USE NCUTILS

    integer, dimension(:), allocatable  :: rxn_indx
    character*(header_length), dimension(:), allocatable :: rxn_name
    character*(header_length), dimension(:), allocatable :: header
    real, dimension(:,:), allocatable   :: top_rates_flux
    real, dimension(:), allocatable     :: time_integrated_prod
    real, dimension(:), allocatable     :: time_integrated_loss
    integer                             :: id_rate_sp, ind_time
    integer                             :: L,N,M, itime, itop, ispec, ireac, ipos
    integer                             :: ncid, mbox
    logical                             :: found_fg

    CHARACTER*(10)               :: date,time1,time2
    CHARACTER*(filenames_length) :: filename
    CALL date_and_time(date,time1)

    N = n_toprates
    L = N+1+1
    M = N+1+N
    ind_time = M+2

    allocate(rxn_indx(M),&
             rxn_name(M), &
             top_rates_flux(ndat,ind_time), &
             time_integrated_prod(numre), &
             time_integrated_loss(numre), &
             header(ind_time))

    write(6,*) '-- finding top ', N, ' reaction rates (molec/cc/s) --'
    write(6,*) '-- for species ',rate_species,' --'

    ! first we need to find the top (N-1)/2 reactions in each direction
    ! with respect to time integrated flux in molec/cc
    ! (Nth rate is "Others")

    id_rate_sp = 0
    rxn_indx = 0
    rxn_name = "none"
    top_rates_flux = 0.
    time_integrated_prod = 0.
    time_integrated_loss = 0.
    header=""

! allocate this here instead of in io_module.f90, to save memory

    ALLOCATE(reacrate(numre,1,1))
    CALL eznc_get_dimension(ncid,"mbox",mbox)

! find target species index in alphabetically sorted species list
    do ispec = 1, numsp
      if (species(ispec)%code == rate_species) then
        id_rate_sp = ispec
        print*,species(ispec)%code ,id_rate_sp
        exit
      endif
    enddo

    do itime = 1,ndat
      !print*,"itime = ",itime
      !CALL eznc_get_3Dreal(ncid,"reacrate", maxre,  mbox,  ndat, &
      !                           reacrate,1,numre,1,1,itime,itime)
      CALL eznc_get_3Dreal(ncid,"reacrate",reacrate, &
                                           1,numre,1, &
                                           1,1,1, &
                                           itime,itime,1)
     
      do ireac = 1, numre
! production
        if(species(id_rate_sp)%phase .ne. PRECURSOR) then
          do ipos = 1, mxright
            if(idpdstoi(ireac,ipos) == id_rate_sp) then
       !print*,"found product: ",ireac,species(idrestoi(ireac,1))%code&
       !                       //"+"//species(idrestoi(ireac,2))%code
             time_integrated_prod(ireac) = time_integrated_prod(ireac)+ &
                 reacrate(ireac,1,1) * pdstoicf(ireac,ipos)
            endif
          enddo
        endif
! loss
        do ipos = 1, mxleft
          if(idrestoi(ireac,ipos) == id_rate_sp) then
       !print*,"found reactant: ",ireac,species(idrestoi(ireac,1))%code&
       !                        //"+"//species(idrestoi(ireac,2))%code
            time_integrated_loss(ireac) = time_integrated_loss(ireac)- &
                reacrate(ireac,1,1) * restoicf(ireac,ipos)
          endif
        enddo
      if(time_integrated_prod(ireac).GT.1e-20.OR.&
         time_integrated_loss(ireac).LT.-1e-20)  &
          print*,time_integrated_prod(ireac) ,time_integrated_loss(ireac)
      enddo
    enddo

! top prod reaction
    if(species(id_rate_sp)%phase .ne. PRECURSOR) then
    if(maxval(time_integrated_prod).gt.0)then
      rxn_indx(1) = maxloc(time_integrated_prod, dim = 1)
      time_integrated_prod(rxn_indx(1)) = 0
      if(numstoi(rxn_indx(1),1).EQ.1)then
        rxn_name(1) = species(idrestoi(rxn_indx(1),1))%code&
                    (1:LEN_TRIM(species(idrestoi(rxn_indx(1),1))%code))//&
                            '_UNIMOLEC                         '
      ELSE
        rxn_name(1) = species(idrestoi(rxn_indx(1),1))%code&
                    (1:LEN_TRIM(species(idrestoi(rxn_indx(1),1))%code))//&
                           '+'//species(idrestoi(rxn_indx(1),2))%code&
                    (1:LEN_TRIM(species(idrestoi(rxn_indx(1),2))%code))
      endif
      !print*,"top prod",rxn_indx(1),rxn_name(1)
    endif
    endif

! top loss reaction (index L=N+2)
    rxn_indx(L) = minloc(time_integrated_loss, dim = 1)
    time_integrated_loss(rxn_indx(L)) = 0
    if(numstoi(rxn_indx(L),1).EQ.1)then
      rxn_name(L) = species(idrestoi(rxn_indx(L),1))%code&
                  (1:LEN_TRIM(species(idrestoi(rxn_indx(L),1))%code))//&
                           '_UNIMOLEC'
    ELSE
      rxn_name(L) = species(idrestoi(rxn_indx(L),1))%code&
                  (1:LEN_TRIM(species(idrestoi(rxn_indx(L),1))%code))//&
                         '+'//species(idrestoi(rxn_indx(L),2))%code&
                  (1:LEN_TRIM(species(idrestoi(rxn_indx(L),2))%code))
    endif
    !print*,"top loss",rxn_indx(L),rxn_name(L)

! sorted production reactions
    if(species(id_rate_sp)%phase .ne. PRECURSOR) then
    do itop = 2, N
      if(maxval(time_integrated_prod).gt.0)then
        rxn_indx(itop) = maxloc(time_integrated_prod, dim = 1)
        time_integrated_prod(rxn_indx(itop)) = 0
        if(numstoi(rxn_indx(itop),1).EQ.1)then
        rxn_name(itop) = &
                        species(idrestoi(rxn_indx(itop),1))%code&
            (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),1))%code))//&
                   '_UNIMOLEC'
        ELSE
        rxn_name(itop) = &
                          species(idrestoi(rxn_indx(itop),1))%code&
              (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),1))%code))//&
                     '+'//species(idrestoi(rxn_indx(itop),2))%code&
              (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),2))%code))
        endif
        print*,"prod",rxn_indx(itop),rxn_name(itop)
      endif
    enddo
    endif


! sorted loss reactions
    do itop = L+1,M
      rxn_indx(itop) = max(minloc(time_integrated_loss, dim = 1),0)
      time_integrated_loss(rxn_indx(itop)) = 0
      if(minval(time_integrated_loss).LT.0)then
        if(numstoi(rxn_indx(itop),1).EQ.1)then
        rxn_name(itop) = &
                        species(idrestoi(rxn_indx(itop),1))%code&
            (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),1))%code))//&
                        '_UNIMOLEC'
        ELSE
        rxn_name(itop) = &
                        species(idrestoi(rxn_indx(itop),1))%code&
            (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),1))%code))//&
                   '+'//species(idrestoi(rxn_indx(itop),2))%code&
            (1:LEN_TRIM(species(idrestoi(rxn_indx(itop),2))%code))
        endif
        print*,"loss",rxn_indx(itop),rxn_name(itop)
      endif
    enddo

    do itop = 1, N
      header(itop) = rxn_name(itop)
    enddo
    header(N+1) = "Others [molec/cc/s]"
    do itop = L,M
      header(itop) = rxn_name(itop)
    enddo
    header(M+1) = "Others [molec/cc/s]"
    header(M+2) = "Time [s]"

!    PRINT*,"================="
    do itime = 1, ndat
      top_rates_flux(itime, ind_time) = time(itime)
      !CALL eznc_get_3Dreal(ncid,"reacrate", maxre,  mbox,  ndat, &
      !                           reacrate,1,numre,1,1,itime,itime)
      CALL eznc_get_3Dreal(ncid,"reacrate",reacrate, &
                                           1,numre,1, &
                                           1,1,1, &
                                           itime,itime,1)

      do ireac=  1, numre
! top prod rates (not precursor)
        if(species(id_rate_sp)%phase .ne. PRECURSOR) then
        found_fg = .false.
        do itop = 1, N
          if (ireac == rxn_indx(itop)) then
            found_fg = .true.
            do ipos=1,numstoi(ireac,2)
              if(idpdstoi(ireac,ipos) == id_rate_sp) then
                top_rates_flux(itime,itop) = reacrate(ireac,1,1) &
                                           * pdstoicf(ireac,ipos)
                !print*,"P",itime,rxn_indx(itop),top_rates_flux(itime,itop)
              endif
            enddo
            exit
          endif
        enddo
        if (.not. found_fg) then
          do ipos=1,numstoi(ireac,2)
            if(idpdstoi(ireac,ipos) == id_rate_sp) then
              top_rates_flux(itime, N+1) = top_rates_flux(itime, N+1) + &
                 reacrate(ireac,1,1) * pdstoicf(ireac,ipos)
            endif
          enddo
          !print*,"O",itime,ireac,top_rates_flux(itime, N+1)
        endif
        endif
! top loss rates
        found_fg = .false.
        do itop = L,M
          if (ireac == rxn_indx(itop)) then
            found_fg = .true.
            do ipos=1,numstoi(ireac,1)
              if(idrestoi(ireac,ipos) == id_rate_sp) then
                !print*,itime,ireac,ipos,restoicf(ireac,ipos)
                top_rates_flux(itime,itop) = - reacrate(ireac,1,1) &
                                             * restoicf(ireac,ipos)
                !print*,"L",itime,ireac,top_rates_flux(itime,itop)
              endif
            enddo
            exit
          endif
        enddo
        if (.not. found_fg) then
          do ipos=1,numstoi(ireac,1)
            if(idrestoi(ireac,ipos) == id_rate_sp) then
              !print*,itime,ireac,ipos,restoicf(ireac,ipos)
              top_rates_flux(itime, M+1) = top_rates_flux(itime, M+1) - &
                 reacrate(ireac,1,1) * restoicf(ireac,ipos)
            endif
          enddo
          !print*,"T",itime,ireac,top_rates_flux(itime, M+1)
        endif
      enddo
    enddo

    filename = "top_rates_"//rate_species(1:LEN_TRIM(rate_species))//"_time"
    call write_2D_array(top_rates_flux, header , filename)

    deallocate(header)
    deallocate(idrestoi, idpdstoi, numstoi, &
               restoicf, pdstoicf)
    deallocate(reacrate)
    deallocate(top_rates_flux)
    deallocate(time_integrated_prod, time_integrated_loss)
    deallocate(rxn_indx, rxn_name)

!    CALL date_and_time(date,time2)
!    print*,time1," ",time2

  end subroutine calc_species_rates

!-------------------------------------------------------------------
  subroutine calc_elements_mass_contribution()
!-------------------------------------------------------------------
! called with flag_elementscontrib
! output filename = "elements_contribution_gas_time"
! output filename = "elements_contribution_aer_time"
!-------------------------------------------------------------------

    real, allocatable, dimension(:, :)   :: elements_contrib
    character*(header_length), dimension(:), allocatable :: header
    integer, parameter  :: ind_c = 1, ind_n = 2, ind_o = 3, ind_h = 4, ind_time = 5
    integer  :: iphase, itime

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    ! we want to count C, H, O, N atoms number contributions to both AEROSOL and GAS phases
    write(6,*) '-- calculating elements contributions to organic matter --'

    allocate(elements_contrib(ndat, 5), header(5))

    header(ind_c) = "C [atoms cm-3]"
    header(ind_n) = "N [atoms cm-3]"
    header(ind_o) = "O [atoms cm-3]"
    header(ind_h) = "H [atoms cm-3]"
    header(ind_time) = "Time [s]"

    do iphase = GAS, AEROSOL, AEROSOL-GAS
      elements_contrib = 0.
      do itime = 1, ndat
        elements_contrib(itime, ind_time) = time(itime)
        elements_contrib(itime, ind_c) = sum(concentrations(itime, :)*species%nc, &
            mask = (species%phase == iphase))
        elements_contrib(itime, ind_n) = sum(concentrations(itime, :)*species%nn, &
            mask = (species%phase == iphase))
        elements_contrib(itime, ind_o) = sum(concentrations(itime, :)*species%no, &
            mask = (species%phase == iphase))
        elements_contrib(itime, ind_h) = sum(concentrations(itime, :)*species%nh, &
            mask = (species%phase == iphase))
      enddo

      if (iphase == GAS) then
        filename = "elements_contribution_gas_time"
      else if (iphase == AEROSOL) then
        filename = "elements_contribution_aer_time"
      endif
      call write_2D_array(elements_contrib, header, filename)
    enddo

    deallocate(header, elements_contrib)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_elements_mass_contribution
!-------------------------------------------------------------------
  subroutine calc_phasestate()
!-------------------------------------------------------------------
! called with flag_phasestate
! output filename = "aerosol_glass_transition_temp"
!-------------------------------------------------------------------

    ! this subroutine estimates dry glass transition temperature in the aeroosol phase
    ! Tg (K) following Shiraiwa et al (2017).
    ! updated following DeRieux et al., (2018)
    ! This is a crude estimate but it gives an idea of the phase state of the aerosol particles
    ! (Tg/T)> 1 => solid
    ! (Tg/T) < 1 => semi-solid or liquid

    real, allocatable, dimension(:,:)  :: tg
    character*(header_length), dimension(:), allocatable :: header

    real :: mass_frac, total_mass, mw, oc, temp_tg
    real :: no, nc, nh
    real, parameter :: A=-21.57,B=1.51,C=-1.7E-3,D=131.4,E=-0.25
    real, parameter :: n0C_CH = 1.93, bC_CH = 61.99, bH_CH = -113.33, bCH_CH = 28.74
    real, parameter :: n0C_CHO = 12.13, bC_CHO = 10.95, bH_CHO = -41.82, bCH_CHO = 21.61, bO_CHO = 118.96, bCO_CHO = -24.38

! hypotheses for humid aerosol
    real, parameter :: kappa = 0.1, low_kappa = 0.05, high_kappa = 0.15 ! hygroscopicity parameter e.g. Gunthe et al. 2009
    real, parameter :: rho_water = 1. !  water density [g cm-3]
    real, parameter :: rho_soa   = 1.4 ! soa organic density [g cm-3] (shiraiwa et al 2017)
    real, parameter :: tgw = 136 ! glass transition temperature for water [K]
    real, parameter :: kgt = 2.5, low_kgt = 1.5, high_kgt = 3.5 ! gordon taylor constant (shiraiwa et al 2017)

! estimate viscosity
    real, parameter :: n_inf = 1e-5 ! viscosity at infinite temperature [Pa s]
    real            :: t0 ! Vogel Temperature [K]
    real, parameter :: frag = 10 ! fragility parameter
    real    :: org_frac, mh2o
    real    :: mh2o_low, mh2o_high, org_frac_low, org_frac_high

    integer :: itime, ispec
    integer, parameter :: ind_time = 1, ind_temp = 2, ind_tg_dry = 3, ind_tg_humid = 4
    integer, parameter :: ind_low_kgt = 5, ind_high_kgt = 6
    integer, parameter :: ind_low_kappa = 7, ind_high_kappa = 8, ind_viscosity = 9

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating dry/humid glass transition temperature --'

    allocate(tg(ndat,9), header(9))
    header(ind_time) = "Time [s]"
    header(ind_temp) = "Temperature [K]"
    header(ind_tg_dry)   = "Dry Glass Transition Temperature [K]"
    header(ind_tg_humid)   = "Humid Glass Transition Temperature [K]"
    header(ind_low_kgt)   = "Low Gordon Taylor Glass Transition Temperature [K]"
    header(ind_high_kgt)   = "High Gordon Taylor Glass Transition Temperature [K]"
    header(ind_low_kappa)   = "Low Kappa Glass Transition Temperature [K]"
    header(ind_high_kappa)   = "High Kappa Glass Transition Temperature [K]"
    header(ind_viscosity) = "Viscosity [Pa s]"

    tg(:, ind_time) = time
    tg(:, ind_temp) = temperature

    do itime=1, ndat
      tg(itime, ind_tg_dry) = 0.
      tg(itime, ind_tg_humid) = 0.
      total_mass = sum(concentrations(itime,:) * species%molw, mask = species%phase == AEROSOL)
      mh2o = kappa*rho_water*total_mass/(rho_soa*(100/rh(itime) - 1))
      mh2o_low = low_kappa*rho_water*total_mass/(rho_soa*(100/rh(itime) - 1))
      mh2o_high = high_kappa*rho_water*total_mass/(rho_soa*(100/rh(itime) - 1))
      org_frac   = total_mass / (mh2o + total_mass)
      org_frac_low = total_mass / (mh2o_low + total_mass)
      org_frac_high = total_mass / (mh2o_high + total_mass)
      do ispec = 1, numsp
        if(species(ispec)%phase .ne. AEROSOL .or. species(ispec)%nc == 0) cycle
        mass_frac = concentrations(itime,ispec) * species(ispec)%molw / total_mass
        if(mass_frac <= 1e-4) cycle
        mw = species(ispec)%molw
        oc = real(species(ispec)%no, kind = 8) / real(species(ispec)%nc, kind = 8)
        no = real(species(ispec)%no, kind = 8)
        nc = real(species(ispec)%nc, kind = 8)
        nh = real(species(ispec)%nh, kind = 8)
        if (species(ispec)%no == 0) then ! CH case
          temp_tg = (n0C_CH + log(nc))*bC_CH + log(nh)*bH_CH + log(nc)*log(nh)*bCH_CH
        else
          temp_tg = (n0C_CHO + log(nc))*bC_CHO + log(nh)*bH_CHO + log(nc)*log(nh)*bCH_CHO + &
                    log(no)*bO_CHO + log(nc)*log(no)*bCO_CHO
        endif
        tg(itime, ind_tg_dry) = tg(itime, ind_tg_dry)+ &
           mass_frac * temp_tg
      enddo
      tg(itime, ind_tg_humid) = ((1-org_frac)*tgw+org_frac*tg(itime, ind_tg_dry)/kgt)/&
            ((1-org_frac) +org_frac/kgt)
      tg(itime, ind_low_kappa) = ((1-org_frac_low)*tgw+org_frac_low*tg(itime, ind_tg_dry)/kgt)/&
            ((1-org_frac_low) +org_frac_low/kgt)
      tg(itime, ind_high_kappa) = ((1-org_frac_high)*tgw+org_frac_high*tg(itime, ind_tg_dry)/kgt)/&
            ((1-org_frac_high) +org_frac_high/kgt)
      tg(itime, ind_low_kgt) = ((1-org_frac)*tgw+org_frac*tg(itime, ind_tg_dry)/low_kgt)/&
            ((1-org_frac) +org_frac/low_kgt)
      tg(itime, ind_high_kgt) = ((1-org_frac)*tgw+org_frac*tg(itime, ind_tg_dry)/high_kgt)/&
            ((1-org_frac) +org_frac/high_kgt)
      t0 = 39.17*tg(itime, ind_tg_humid)/(frag + 39.17)
      tg(itime, ind_viscosity) = n_inf*exp(T0*frag/(tg(itime, ind_temp) - t0))
    enddo

    filename = "aerosol_glass_transition_temp"
    call write_2D_array(tg, header, filename)

    deallocate(header, tg)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_phasestate
!-------------------------------------------------------------------
  subroutine calc_entropy()
!-------------------------------------------------------------------
! called with flag_entropy
! output filename = "aerosol_entropy"
!-------------------------------------------------------------------

!  this subroutine calculates the entropy of aerosol particles Hi
! according to Riemer & West (2013)
    real, allocatable, dimension(:,:)  :: h
    character*(header_length), dimension(:), allocatable :: header
    real :: mass_frac, total_mass
    real, dimension(:), allocatable :: cumulative_total_mass
    integer :: itime, ispec, itime2
    integer, parameter :: ind_time = 1, ind_h = 2, ind_h2 = 3, ind_avrg = 4, ind_bulk = 5

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating aerosols entropy --'

    allocate(h(ndat, 5), header(5), cumulative_total_mass(ndat))
    header(ind_time) = "Time [s]"
    header(ind_h) = "1st Order Mixing Entropy"
    header(ind_h2) = "2nd Order Mixing Entropy"
    header(ind_avrg) = "Time Average mixing entropy"
    header(ind_bulk) = "Cumulative bulk mixing entropy"


    h(:,ind_time) = time
    cumulative_total_mass = 0.

    do itime = 1, ndat
      h(itime, ind_h) = 0.
      h(itime, ind_h2) = 0.
      total_mass = sum(concentrations(itime,:) * species%molw, mask = species%phase == AEROSOL)
      if (itime > 1) then
        cumulative_total_mass(itime) = cumulative_total_mass(itime-1) + total_mass
      else
        cumulative_total_mass(itime) = total_mass
      endif

      do ispec = 1, numsp
        if(species(ispec)%phase .ne. AEROSOL) cycle
        mass_frac   = concentrations(itime,ispec) * species(ispec)%molw / total_mass
        if(mass_frac <= 0.) cycle
        h(itime, ind_h) = h(itime, ind_h) - &
          mass_frac*log(mass_frac)
        h(itime, ind_h2) = h(itime, ind_h2) + mass_frac**2
      enddo
      h(itime, ind_h2) = 1- h(itime, ind_h2)
    enddo

    do itime = 1, ndat
      h(itime,ind_avrg) = 0.
      do itime2 = 1, itime
        total_mass = sum(concentrations(itime,:) * species%molw, mask = species%phase == AEROSOL)
        mass_frac = total_mass / cumulative_total_mass(itime)
        if(mass_frac <= 0.) cycle
        h(itime,ind_avrg) = h(itime,ind_avrg) - &
          mass_frac*log(mass_frac)
      enddo
    enddo


    do itime = 1, ndat
      h(itime, ind_bulk) = 0.
      do ispec = 1, numsp
        if(species(ispec)%phase .ne. AEROSOL) cycle
        total_mass = 0.
        do itime2 = 1, itime
          total_mass = total_mass + concentrations(itime2, ispec) * species(ispec)%molw
        enddo
        mass_frac = total_mass / cumulative_total_mass(itime)
        if(mass_frac <= 0.) cycle
        h(itime, ind_bulk) = h(itime, ind_bulk) - &
          mass_frac*log(mass_frac)
      enddo
    enddo

    filename = "aerosol_entropy"
    call write_2D_array(h, header, filename)

    deallocate(header, h, cumulative_total_mass)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_entropy
!-------------------------------------------------------------------
  subroutine calc_chochonfreq()
!-------------------------------------------------------------------
! called with flag_chochonfreq
! output filename = "chochonfreq_"//trim(formulas(iformula))
!-------------------------------------------------------------------

    ! calculate for each carbon number, each phase, the number of CHO, CHNO, CHOS, CHNOS containing species
    real, dimension(:,:), allocatable :: chochonfreq_molec
    character*(header_length), dimension(:), allocatable :: header_molec
    integer, dimension(:), allocatable :: nc_bins
    logical, dimension(:), allocatable  :: species_mask
    integer :: min_nc, max_nc, nbins
    INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: ncols = nphase+2
    integer, parameter :: icho =1, ichno = 2, ichos = 3, ichnos = 4, nformulas = 4
    integer, parameter :: ind_nc = nphase+1, ind_time = nphase+2
    character*(5), dimension(nformulas) :: formulas
    integer  :: itime, iphase, ibin, iformula
    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating CHO, CHNO, CHOS, CHNOS frequency --'

    min_nc = minval(species%nc)
    max_nc = maxval(species%nc)
    nbins = max_nc - min_nc + 1

    formulas(icho) = "CHO"
    formulas(ichno) = "CHNO"
    formulas(ichos) = "CHOS"
    formulas(ichnos) = "CHNOS"

    allocate(chochonfreq_molec(ndat*nbins, ncols), &
             header_molec(ncols), nc_bins(nbins), species_mask(numsp))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header_molec(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header_molec(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header_molec(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header_molec(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header_molec(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header_molec(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header_molec(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header_molec(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header_molec(iphase) = "Inorganics [molec cm-3]"
        case default
          header_molec(iphase) = ""
      end select
    enddo

    do ibin = 1, nbins
      nc_bins(ibin) = min_nc + ibin -1
    enddo

    header_molec(ind_nc)   = "Carbon Number"
    header_molec(ind_time) = "Time [s]"

    do iformula=1,nformulas
      write(6,*) "begin ",formulas(iformula)
      chochonfreq_molec = 0.
      do itime=1,ndat
        chochonfreq_molec((itime-1)*nbins+1:itime*nbins,ind_time) = time(itime)
        chochonfreq_molec((itime-1)*nbins+1:itime*nbins,ind_nc) = real(nc_bins,kind = 8)
        do iphase = 1, nphase
          do ibin = 1, nbins
            species_mask = species%phase == iphase .and. &
                      species%nc == nc_bins(ibin)
            chochonfreq_molec(nbins*(itime-1)+ibin, iphase) = sum( &
               concentrations(itime, :), mask = species_mask)
          enddo
        enddo
      enddo

      write(6,*) "end ", formulas(iformula)
      filename = "chochonfreq_"//trim(formulas(iformula))
      call write_2D_array(chochonfreq_molec, header_molec, filename)
    enddo

    deallocate(chochonfreq_molec, header_molec, nc_bins, species_mask)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_chochonfreq
!-------------------------------------------------------------------
  subroutine calc_dbe()
!-------------------------------------------------------------------
! called with flag_dbeai
! output filename = "dbe"
!-------------------------------------------------------------------

    ! calculate mass_weighted average DBE and AI, for each phase
    INTEGER,PARAMETER:: nphase = LAST_PHASE-2
    INTEGER,PARAMETER:: ncols = nphase+1
    
    real, dimension(:,:), allocatable :: dbe
    character*(header_length), dimension(:), allocatable :: header
    logical, dimension(:), allocatable  :: species_mask

    integer, parameter :: ind_time = ncols
    real, parameter :: ugfac = 1.660578e-12
    real :: total_mass
    integer :: itime, iphase

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)


    write(6,*) " -- calculating DBE -- "

    allocate(dbe(ndat, ncols), &
             header(ncols), &
             species_mask(numsp))
    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header(iphase) = "Inorganics [molec cm-3]"
        case default
          header(iphase) = ""
      end select
    enddo
    header(ind_time) = "Time [s]"

    dbe = 0.
    species_mask = .FALSE.
    do itime = 1, ndat
      dbe(itime,ind_time) = time(itime)
      do iphase = 1, nphase
        species_mask = species%phase == iphase
        total_mass = sum( concentrations(itime, :)*species%molw * ugfac, &
              mask = species_mask)
        if (total_mass <= 0.) then
          dbe(itime, iphase) = 0.
        else
          dbe(itime, iphase) = sum( species%dbe*concentrations(itime, :)*species%molw * ugfac, &
              mask = species_mask)/total_mass
        endif

      enddo
    enddo

    filename = "dbe"
    call write_2D_array(dbe, header, filename)

    deallocate(dbe, header, species_mask)
  
    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_dbe
!-------------------------------------------------------------------
  subroutine calc_ai()
!-------------------------------------------------------------------
! called with flag_dbeai
! output filename = "ai"
! calculate mass_weighted average DBE and AI, for each phase
!-------------------------------------------------------------------

    real, dimension(:,:), allocatable :: ai
    character*(header_length), dimension(:), allocatable :: header
    logical, dimension(:), allocatable  :: species_mask

    INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: ncols = nphase+1
    integer, parameter :: ind_time = ncols
    real, parameter :: ugfac = 1.660578e-12
    real :: total_mass
    integer :: itime, iphase


    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)


    write(6,*) " -- calculating AI -- "

    allocate(ai(ndat, ncols), &
             header(ncols), &
             species_mask(numsp))
    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header(iphase) = "Inorganics [molec cm-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_time) = "Time [s]"

    ai = 0.
    species_mask = .FALSE.
    do itime = 1, ndat
      ai(itime,ind_time) = time(itime)
      do iphase = 1, nphase
        species_mask = species%phase == iphase
        total_mass = sum( concentrations(itime, :)*species%molw * ugfac, &
              mask = species_mask)
        if(total_mass <= 0.) then
          ai(itime, iphase) = 0.
        else
          ai(itime, iphase) = sum( species%ai*concentrations(itime, :)*species%molw * ugfac, &
              mask = species_mask)/total_mass
        endif
      enddo
    enddo

    filename = "ai"
    call write_2D_array(ai, header, filename)

    deallocate(ai,header,species_mask)

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_ai
!---------------------------------------------------------------------
  subroutine calc_mass_spectrum()
!---------------------------------------------------------------------
! called with flag_massspectrum (unless commented out in main.f90)
!---------------------------------------------------------------------
! output filename = "mass_spectrum_time"
! calculate a mass spectrum (no fragmentation)
!---------------------------------------------------------------------

    real, dimension(:,:), allocatable :: mass_spec
    character*(header_length), dimension(:), allocatable :: header
    INTEGER,PARAMETER :: nphase = LAST_PHASE-2
    INTEGER,PARAMETER :: ncols = nphase+3
    integer, parameter :: ind_lowerlim=ncols-2
    integer, parameter :: ind_upperlim=ncols-1
    integer, parameter :: ind_time=ncols
    integer :: iphase, itime, ibin, ispec

    type(bins_class) :: mass_bins

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) " -- calculating mass spectrum -- "
    call mass_bins%initialize(real(floor(minval(species%molw)), kind = 8), & ! min mass,
                              real(ceiling(maxval(species%molw)), kind = 8), & ! max_mass
                              1.0)

    allocate(mass_spec(ndat*mass_bins%nbins, ncols), header(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header(iphase) = "Inorganics [molec cm-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_time) = "Time [s]"
    header(ind_lowerlim) = "Mass bin lower limit"
    header(ind_upperlim) = "Mass bin upper limit"

    mass_spec = 0.
    do itime=1,ndat
      mass_spec((itime-1)*mass_bins%nbins + 1:itime*mass_bins%nbins, ind_lowerlim) = mass_bins%lowerlim
      mass_spec((itime-1)*mass_bins%nbins + 1:itime*mass_bins%nbins, ind_upperlim) = mass_bins%upperlim
      mass_spec((itime-1)*mass_bins%nbins + 1:itime*mass_bins%nbins, ind_time) = time(itime)

      do ispec =1,numsp
        ibin = mass_bins%find_index(species(ispec)%molw)
        if (ibin == -1 .or. ibin < 1 .or. ibin > mass_bins%nbins) then
          write(6,*) "error in finding bin index for value:", species(ispec)%molw
          write(6,*) "ibin=",ibin
          STOP
        endif
        mass_spec(mass_bins%nbins*(itime-1)+ibin,species(ispec)%phase) = &
          mass_spec(mass_bins%nbins*(itime-1)+ibin,species(ispec)%phase)+concentrations(itime,ispec)
      enddo

    enddo

    filename = "mass_spectrum_time"
    call write_2D_array(mass_spec, header, filename)

    deallocate(mass_spec,header)
    call mass_bins%destroy()

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_mass_spectrum
!---------------------------------------------------------------------
  subroutine calc_mass_spectrum_lifetime()
!---------------------------------------------------------------------
! called with flag_massspectrum 
! output filename = "mass_spectrum_lifetime"
! calculate a "mass spectrum" at each precursor lifetime (no fragmentation)
!---------------------------------------------------------------------

    real, dimension(:,:), allocatable :: mass_spec
    character*(header_length), dimension(:), allocatable :: header
    integer, parameter :: nphase = LAST_PHASE-2
    integer, parameter :: ncols = nphase+3
    integer, parameter :: ind_lowerlim=nphase+1
    integer, parameter :: ind_upperlim=nphase+2
    integer, parameter :: ind_lifetime=nphase+3
    integer :: ilt, iphase, itime, ibin, ispec, nlt

    type(bins_class) :: mass_bins

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) " -- calculating mass spectrum by lifetime -- "
    call mass_bins%initialize(real(floor(minval(species%molw)), kind = 8), & ! min mass,
                              real(ceiling(maxval(species%molw)), kind = 8), & ! max_mass
                              1.0)
    nlt = size(efold_indices,1)

    allocate(mass_spec(nlt*mass_bins%nbins, ncols), header(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header(iphase) = "Inorganics [molec cm-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_lowerlim) = "Mass bin lower limit"
    header(ind_upperlim) = "Mass bin upper limit"
    header(ind_lifetime) = "Time [precursor lifetimes]"

    mass_spec = 0.
    do ilt=1,nlt
      itime = efold_indices(ilt)
      mass_spec((ilt-1)*mass_bins%nbins + 1:ilt*mass_bins%nbins, ind_lowerlim) = mass_bins%lowerlim
      mass_spec((ilt-1)*mass_bins%nbins + 1:ilt*mass_bins%nbins, ind_upperlim) = mass_bins%upperlim
      mass_spec((ilt-1)*mass_bins%nbins + 1:ilt*mass_bins%nbins, ind_lifetime) = ilt

      do ispec =1,numsp
        !print*,"ispec = ",ispec
        ibin = mass_bins%find_index(species(ispec)%molw)
        if (ibin == -1 .or. ibin < 1 .or. ibin > mass_bins%nbins) then
          write(6,*) "error in finding bin index for value:", species(ispec)%molw
          write(6,*) "ibin=",ibin
          STOP
        endif
        mass_spec(mass_bins%nbins*(ilt-1)+ibin,species(ispec)%phase) = &
        mass_spec(mass_bins%nbins*(ilt-1)+ibin,species(ispec)%phase)+concentrations(itime,ispec)
      enddo

    enddo

    filename = "mass_spectrum_lifetime"
    call write_2D_array(mass_spec, header, filename)

    deallocate(mass_spec,header)
    call mass_bins%destroy()

    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_mass_spectrum_lifetime
!-----------------------------------------------------------------------
  subroutine calc_bubble()
!-----------------------------------------------------------------------
! called with flag_bubble (unless commented out in main.f90)
! output filename = "bubble_molec_time"
!-----------------------------------------------------------------------

    real, dimension(:,:), allocatable :: bubble
    character*(header_length), dimension(:), allocatable :: header
    integer, parameter :: nphase = LAST_PHASE-2
    integer, parameter :: ncols = nphase+5
    integer, parameter :: ind_pvaplowerlim=nphase+1
    integer, parameter :: ind_pvapupperlim=nphase+2
    integer, parameter :: ind_osclowerlim=nphase+3
    integer, parameter :: ind_oscupperlim=nphase+4
    integer, parameter :: ind_time=ncols
    integer :: iphase, itime, ipvap, iosc, ispec, iline
    real, parameter :: ugfac = 1.660578e-12
    real :: minpvap, maxpvap
    type(bins_class) :: pvap_bins, osc_bins
    
    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) " -- calculating bubble plot -- "
    minpvap = real(floor(minval(log10(species%pvap%pvap_atm_298), mask = species%pvap%pvap_atm_298 > 0)), kind = 8)
    maxpvap = real(ceiling(maxval(log10(species%pvap%pvap_atm_298), mask = species%pvap%pvap_atm_298 > 0)), kind = 8)
    call pvap_bins%initialize(minpvap, maxpvap, 1.0)
    print *, pvap_bins%nbins                
    call osc_bins%initialize(real(floor(minval(species%osc)), kind = 8), & 
                              real(ceiling(maxval(species%osc)), kind = 8), & 
                              0.2)
    print *, osc_bins%nbins           
                              
    allocate(bubble(ndat*pvap_bins%nbins*osc_bins%nbins, ncols), &
             header(ncols))
    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [molec cm-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [molec cm-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [molec cm-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [molec cm-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [molec cm-3]"
        case (WALL)
          header(iphase) = "Wall [molec cm-3]"
        case (DIMER)
          header(iphase) = "Dimers [molec cm-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [molec cm-3]"
        case (INORG)
          header(iphase) = "Inorganics [molec cm-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_time) = "Time [s]"
    header(ind_pvaplowerlim) = "Pvap bin lower limit [log10(atm)]"
    header(ind_pvapupperlim) = "Pvap bin upper limit [log10(atm)]"
    header(ind_osclowerlim) = "OSc bin lower limit"
    header(ind_oscupperlim) = "OSc bin upper limit"
    
    bubble = 0.
    iline = 0
    do itime = 1, ndat
      do ipvap = 1, pvap_bins%nbins
        do iosc = 1, osc_bins%nbins    
          iline = iline + 1    
          bubble(iline, ind_time) = time(itime)
          bubble(iline, ind_pvaplowerlim) = pvap_bins%lowerlim(ipvap)
          bubble(iline, ind_pvapupperlim) = pvap_bins%upperlim(ipvap)
          bubble(iline, ind_osclowerlim) = osc_bins%lowerlim(iosc)
          bubble(iline, ind_oscupperlim) = osc_bins%upperlim(iosc)          
        enddo
      enddo
    enddo
    
    do itime = 1, ndat
      do ispec = 1, numsp
        if (species(ispec)%pvap%pvap_atm_298 <= 0) cycle
      
        ipvap = pvap_bins%find_index(log10(species(ispec)%pvap%pvap_atm_298))      
        iosc = osc_bins%find_index(species(ispec)%osc)
        
        iline = osc_bins%nbins*pvap_bins%nbins*(itime-1) + &
                osc_bins%nbins*(ipvap-1) + &
                iosc
      !print*,itime,ispec,iline,species(ispec)%phase, &
      !       species(ispec)%molw,concentrations(itime,ispec)
      !print*,'-bubble-' ,bubble(iline, species(ispec)%phase)
        bubble(iline, species(ispec)%phase) = bubble(iline, species(ispec)%phase) + &
            concentrations(itime,ispec)*species(ispec)%molw*ugfac
      enddo
    enddo

    filename = "bubble_molec_time"
    call write_2D_array(bubble, header, filename)


    deallocate(bubble,header)
    call pvap_bins%destroy()
    call osc_bins%destroy()
    
    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_bubble
!-----------------------------------------------------------------------
  subroutine calc_bubble_pvap_lifetime()
!-----------------------------------------------------------------------
! called with BOTH flag_bubble and flag_pvap
! output filename = "bubble_mass_pvap_lifetime"
!-----------------------------------------------------------------------

    real, dimension(:,:), allocatable :: bubble
    character*(header_length), dimension(:), allocatable :: header
    integer, parameter :: nphase = 3
    integer, parameter :: ncols = nphase+5
    integer, parameter :: ind_pvaplowerlim=nphase+1
    integer, parameter :: ind_pvapupperlim=nphase+2
    integer, parameter :: ind_osclowerlim=nphase+3
    integer, parameter :: ind_oscupperlim=nphase+4
    integer, parameter :: ind_lifetime=ncols
    integer :: iphase, itime, ilt, ipvap, iosc, ispec, iline, nlt
    real, parameter :: ugfac = 1.660578e-12
    real, parameter :: limit = 0.01
    real :: minpvap, maxpvap
    real :: minosc, maxosc
    type(bins_class) :: pvap_bins, osc_bins
    
    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) " -- calculating bubble plot (OSc v pvap) by lifetime -- "
    minpvap = real(floor(minval(log10(species%pvap%pvap_atm_298), mask = species%pvap%pvap_atm_298 > 0)), kind = 8)
    maxpvap = real(ceiling(maxval(log10(species%pvap%pvap_atm_298), mask = species%pvap%pvap_atm_298 > 0)), kind = 8)

    minosc = real(floor(minval(species%osc, mask = species%osc > -10)), kind = 8)
    maxosc = real(ceiling(maxval(species%osc, mask = species%osc > -10)), kind = 8)

    call pvap_bins%initialize(minpvap, maxpvap, 1.0)
    print *, pvap_bins%nbins                

    call osc_bins%initialize(minosc, maxosc, 0.2)
    !call osc_bins%initialize(real(floor(minval(species%osc)), kind = 8), &
    !                         real(ceiling(maxval(species%osc)), kind = 8), &
    !                         0.2)
    print *, osc_bins%nbins           
                    
    nlt = size(efold_indices,1)
          
    allocate(bubble(nlt*pvap_bins%nbins*osc_bins%nbins, ncols), &
             header(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header(iphase) = "Inorganics [ug m-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_pvaplowerlim) = "Pvap bin lower limit [log10(atm)]"
    header(ind_pvapupperlim) = "Pvap bin upper limit [log10(atm)]"
    header(ind_osclowerlim) = "OSc bin lower limit"
    header(ind_oscupperlim) = "OSc bin upper limit"
    header(ind_lifetime) = "Time [precursor lifetimes]"
    
    bubble = 0.
    iline = 0
    do ilt = 1, nlt
      itime = efold_indices(ilt)
      do ipvap = 1, pvap_bins%nbins
        do iosc = 1, osc_bins%nbins    
          iline = iline + 1    
          bubble(iline, ind_pvaplowerlim) = pvap_bins%lowerlim(ipvap)
          bubble(iline, ind_pvapupperlim) = pvap_bins%upperlim(ipvap)
          bubble(iline, ind_osclowerlim) = osc_bins%lowerlim(iosc)
          bubble(iline, ind_oscupperlim) = osc_bins%upperlim(iosc)          
          bubble(iline, ind_lifetime) = ilt
        enddo
      enddo
    enddo
   
    do ilt = 1, nlt
      itime = efold_indices(ilt)
      do ispec = 1, numsp
        !if (species(ispec)%pvap%pvap_atm_298 <= 0) cycle
        if (species(ispec)%pvap%pvap_atm_298.LE.0.or. &
            species(ispec)%nc.LE.0) cycle
        if (species(ispec)%phase.gt.3) cycle

        ipvap = pvap_bins%find_index(log10(species(ispec)%pvap%pvap_atm_298))      
        iosc = osc_bins%find_index(species(ispec)%osc)
        
        iline = osc_bins%nbins*pvap_bins%nbins*(ilt-1) + &
                osc_bins%nbins*(ipvap-1) + &
                iosc
        bubble(iline, species(ispec)%phase) = bubble(iline, species(ispec)%phase) + &
            concentrations(itime,ispec)*species(ispec)%molw*ugfac

! output top few species (diagnostic)
        if(concentrations(itime,ispec)*species(ispec)%molw*ugfac.gt.limit)then
          WRITE(98,*)pvap_bins%lowerlim(ipvap),pvap_bins%upperlim(ipvap), &
                     osc_bins%lowerlim(iosc),osc_bins%upperlim(iosc), &
                     concentrations(itime,ispec)*species(ispec)%molw*ugfac, &
                     species(ispec)%code, &
                     species(ispec)%formula
        endif
      !print*,'-bubble-',ilt,itime,ispec,iline,species(ispec)%phase,&
      !         bubble(iline,species(ispec)%phase)
      enddo
    enddo

    filename = "bubble_mass_pvap_lifetime"
    call write_2D_array(bubble, header, filename)

    deallocate(bubble,header)
    call pvap_bins%destroy()
    call osc_bins%destroy()
    
    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_bubble_pvap_lifetime
!-----------------------------------------------------------------------
  subroutine calc_bubble_cstar_lifetime()
!-----------------------------------------------------------------------
! called with BOTH flag_bubble and flag_cstar
! output filename = "bubble_mass_cstar_lifetime"
!-----------------------------------------------------------------------

    real, dimension(:,:), allocatable :: bubble
    character*(header_length), dimension(:), allocatable :: header
    !integer, parameter :: nphase = LAST_PHASE-2
    integer, parameter :: nphase = 3
    integer, parameter :: ncols = nphase+5
    integer, parameter :: ind_cstarlowerlim=nphase+1
    integer, parameter :: ind_cstarupperlim=nphase+2
    integer, parameter :: ind_osclowerlim=nphase+3
    integer, parameter :: ind_oscupperlim=nphase+4
    integer, parameter :: ind_lifetime=ncols
    integer :: iphase, itime, ilt, icstar, iosc, ispec, iline, nlt
    real, parameter :: ugfac = 1.660578e-12
    real, parameter :: limit = 0.01
    real :: mincstar, maxcstar
    real :: minosc, maxosc
    type(bins_class) :: cstar_bins, osc_bins

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) " -- calculating bubble plot (OSc v Cstar) by lifetime -- "

    do ispec = 1, numsp
      if(species(ispec)%pvap%pvap_atm_298.gt.0)then
! DEBUG
!        PRINT*,ispec,species(ispec)%code,species(ispec)%pvap%pvap_atm_298 ! DBG
           species(ispec)%pvap%pvap_cstar_298 = &
            (1e6*species(ispec)%molw*species(ispec)%pvap%pvap_atm_298)/ &
            (8.2e-5*298)
!      ELSE ! DBG
!        PRINT*,ispec,species(ispec)%code ! DBG
      endif
    enddo

    mincstar = real(floor(minval(log10(species%pvap%pvap_cstar_298), mask = species%pvap%pvap_cstar_298 > 0)), kind = 8)
    maxcstar = real(ceiling(maxval(log10(species%pvap%pvap_cstar_298), mask = species%pvap%pvap_cstar_298 > 0)), kind = 8)

    print *, "cstar_bins%nbins"
    call cstar_bins%initialize(mincstar, maxcstar, 1.0)
    print *, cstar_bins%nbins                

    minosc = real(floor(minval(species%osc, mask = species%osc > -10)), kind = 8)
    maxosc = real(ceiling(maxval(species%osc, mask = species%osc > -10)), kind = 8)

    print *, "osc_bins%nbins"
    call osc_bins%initialize(minosc, maxosc, 0.2)
!    call osc_bins%initialize(real(floor(minval(species%osc)), kind = 8), &
!                             real(ceiling(maxval(species%osc)), kind = 8), &
!                             0.2)
    print *, osc_bins%nbins           
                    
    nlt = size(efold_indices,1)
          
    allocate(bubble(nlt*cstar_bins%nbins*osc_bins%nbins, ncols), &
             header(ncols))

    do iphase = 1, nphase
      select case (iphase)
        case (GAS)
          header(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header(iphase) = "Precursor [ug m-3]"
        case (CO_CO2)
          header(iphase) = "CO+CO2 [ug m-3]"
        case (AQUEOUS)
          header(iphase) = "Aqueous Phase [ug m-3]"
        case (WALL)
          header(iphase) = "Wall [ug m-3]"
        case (DIMER)
          header(iphase) = "Dimers [ug m-3]"
        case (SPINUP_SPEC)
          header(iphase) = "Spinup Species [ug m-3]"
        case (INORG)
          header(iphase) = "Inorganics [ug m-3]"
        case default
          header(iphase) = ""
      end select
    enddo

    header(ind_cstarlowerlim) = "Cstar bin lower limit [log10(ug.m-3)]"
    header(ind_cstarupperlim) = "Cstar bin upper limit [log10(ug.m-3)]"
    header(ind_osclowerlim) = "OSc bin lower limit"
    header(ind_oscupperlim) = "OSc bin upper limit"
    header(ind_lifetime) = "Time [precursor lifetimes]"
    
    bubble = 0.
    iline = 0
    do ilt = 1, nlt
      itime = efold_indices(ilt)
      do icstar = 1, cstar_bins%nbins
        do iosc = 1, osc_bins%nbins    
          iline = iline + 1    
          bubble(iline, ind_cstarlowerlim) = cstar_bins%lowerlim(icstar)
          bubble(iline, ind_cstarupperlim) = cstar_bins%upperlim(icstar)
          bubble(iline, ind_osclowerlim) = osc_bins%lowerlim(iosc)
          bubble(iline, ind_oscupperlim) = osc_bins%upperlim(iosc)          
          bubble(iline, ind_lifetime) = ilt
        enddo
      enddo
    enddo
 
    do ilt = 1, nlt
      itime = efold_indices(ilt)
      do ispec = 1, numsp
      if(species(ispec)%pvap%pvap_atm_298.gt.0.and.species(ispec)%nc.gt.0)then
          print*,ilt,ispec
! here's the problem. where ispc > nsat

        icstar = cstar_bins%find_index(log10(species(ispec)%pvap%pvap_cstar_298))      
        iosc = osc_bins%find_index(species(ispec)%osc)
        
        iline = osc_bins%nbins*cstar_bins%nbins*(ilt-1) + &
                osc_bins%nbins*(icstar-1) + &
                iosc

        do iphase = 1, nphase
          if (species(ispec)%phase == iphase) then

            bubble(iline, iphase) = bubble(iline, iphase) + &
                concentrations(itime,ispec)*species(ispec)%molw*ugfac

! output top few species (diagnostic)
            if(concentrations(itime,ispec)*species(ispec)%molw*ugfac &
                                                        .gt.limit)then
              WRITE(98,*)cstar_bins%lowerlim(icstar), &
                         cstar_bins%upperlim(icstar), &
                         osc_bins%lowerlim(iosc), &
                         osc_bins%upperlim(iosc), &
               concentrations(itime,ispec)*species(ispec)%molw*ugfac, &
                         species(ispec)%code, &
                         species(ispec)%formula
            endif
          endif
        enddo
        endif
      enddo
    enddo

    filename = "bubble_mass_Cstar_lifetime"
    call write_2D_array(bubble, header, filename)

    deallocate(bubble,header)
    call cstar_bins%destroy()
    call osc_bins%destroy()
    
    CALL date_and_time(date,time2)
    !print*,time1," ",time2

  end subroutine calc_bubble_cstar_lifetime
!-----------------------------------------------------------------------
  subroutine calc_nitrates()
!-----------------------------------------------------------------------
! called with flag_nitrates
! output filename = "nitrate_func_aer"
! output filename = "nitrate_subs_aer"
!-----------------------------------------------------------------------
  
    character*(header_length), dimension(:), allocatable :: header_func
    real, dimension(:,:), allocatable :: functionalization
    integer, parameter :: IND_TIME = 1, IND_MFUNC = 2, IND_DFUNC = 3, IND_TFUNC = 4, IND_PFUNC = 5, IND_FTOTAL = 6
    character*(header_length), dimension(:), allocatable :: header_sub
    real, dimension(:,:), allocatable :: substitution
    integer, parameter :: IND_PRIM = 2, IND_SEC = 3, IND_TERT = 4, IND_STOTAL = 5 
    
    real, parameter :: ugfac = 1.660578e-12
    
    integer  :: itime, ispec, nfunc, initrate, sub
    real     :: mass, molec
    
    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    allocate( &
      functionalization(ndat, IND_FTOTAL), &
      substitution(ndat, IND_STOTAL), &
      header_func(IND_FTOTAL), &
      header_sub(IND_STOTAL))
      
    header_func(IND_TIME) = "Time [s]"
    header_func(IND_MFUNC) = "Monofunctional Nitrates [ug m-3]"
    header_func(IND_DFUNC) = "Difunctional Nitrates [ug m-3]"
    header_func(IND_TFUNC) = "Trifunctional Nitrates [ug m-3]"
    header_func(IND_PFUNC) = "Morefunctional Nitrates [ug m-3]"
    header_func(IND_FTOTAL) = "Total Nitrates [ug m-3]"
    
    header_sub(IND_TIME)  = "Time [s]"
    header_sub(IND_PRIM)  = "Primary Nitrates [nitrate m-3]"
    header_sub(IND_SEC)   = "Secondary Nitrates [nitrate m-3]"
    header_sub(IND_TERT)  = "Tertiary Nitrates [nitrate m-3]"
    header_sub(IND_STOTAL)= "Total Nitrates [nitrate m-3]"
    
    functionalization = 0.
    substitution = 0.
    
    do itime = 1, ndat
      do ispec = 1, numsp
        if (species(ispec)%phase .ne. AEROSOL) cycle
        if (.not. species(ispec)%is_nitrate) cycle
        molec = concentrations(itime, ispec)
        mass  = molec*species(ispec)%molw*ugfac
        nfunc = species(ispec)%nfunc
        functionalization(itime, IND_TIME) = time(itime)
        functionalization(itime, IND_FTOTAL) = functionalization(itime, IND_FTOTAL) + &
                                               mass                
        if (nfunc == 1) then
          functionalization(itime, IND_MFUNC) = functionalization(itime, IND_MFUNC) + &
                                               mass     
        else if (nfunc == 2) then   
          functionalization(itime, IND_DFUNC) = functionalization(itime, IND_DFUNC) + &
                                               mass
        else if (nfunc == 3) then   
          functionalization(itime, IND_TFUNC) = functionalization(itime, IND_TFUNC) + &
                                               mass
        else if (nfunc > 4) then
          functionalization(itime, IND_PFUNC) = functionalization(itime, IND_PFUNC) + &
                                               mass
        endif
        
        substitution(itime, IND_TIME) = time(itime)
        do initrate = 1, species(ispec)%n_nitrates
          sub = species(ispec)%nitrate_substitution(initrate)
          if (sub == 1) then
            substitution(itime, IND_PRIM) = substitution(itime, IND_PRIM) + molec
          else if (sub == 2) then
            substitution(itime, IND_SEC) = substitution(itime, IND_SEC) + molec
          else if (sub == 3) then
            substitution(itime, IND_TERT) = substitution(itime, IND_TERT) + molec
          endif
        enddo      
        substitution(itime, IND_STOTAL) = substitution(itime, IND_PRIM) + &
                                          substitution(itime, IND_SEC)  + &
                                          substitution(itime, IND_TERT)
      enddo    
    enddo
    
    filename = "nitrate_func_aer"
    call write_2D_array(functionalization, header_func, filename)
    
    filename = "nitrate_subs_aer"
    call write_2D_array(substitution, header_sub, filename)
    
    deallocate(functionalization, substitution, header_func, header_sub)
    
    CALL date_and_time(date,time2)
    PRINT*,time1," ",time2
  
  end subroutine calc_nitrates

!-----------------------------------------------------------------------
  subroutine calc_osc_nc()
!-----------------------------------------------------------------------
! called with flag_OSc
! output filename = "osc_nc"
!-----------------------------------------------------------------------

!    real, dimension(:,:), allocatable :: osc_nc
!    integer, dimension(:), allocatable :: osc_bins, nc_bins
!    character*(header_length), dimension(:), allocatable :: header
!    integer :: min_osc, max_osc, nosc_bin
!    integer :: min_nc, max_nc, nnc_bin

!    integer :: itime, iphase, ispec, iosc, inc
!    integer, parameter :: ind_time = 1, ind_osc = 2, ind_nc = 3, ind_phase = 4, ind_mass = 5
    !integer, dimension(1), parameter :: phases = (/ GAS, AEROSOL, PRECURSOR, CO_CO2, AQUEOUS /)
    !integer, parameter :: nphase = size(phases)
!    write(6,*) '-- calculating oxidation state vs nc --'

!    min_osc = floor(min(species%osc))
!    max_osc = floor(max(species%osc))
 !   nosc_bin = max_osc - min_osc + 1
!
!    min_nc = floor(min(species%nc))
!    max_nc = floor(max(species%nc))
!    nnc_bin = max_nc - min_nc + 1

!    allocate(osc_bins(nosc_bin), nc_bins(nnc_bin))
!    do iosc = 1, nosc_bin
!      osc_bins(iosc) = min_osc + iosc -1
!    enddo
!
!    do inc = 1, nnc_bin
!      nc_bins(inc) = min_nc + inc - 1
!    enddo
!
!    allocate(osc_nc(ndat*nphase*nosc_bin*nnc_bin, 5))
!    osc_nc = 0.
!
!    do itime = 1, ndat
!      do ispec = 1, numsp
!        iphase = species(ispec)%phase
!        iosc = ceiling(species(ispec)%osc-min_osc)
!        inc = ceiling(species(ispec)%nc-min_nc)
!        imass = (itime-1)*(iphase-1)*(iosc-1)*(inc-1)
!        osc_nc(imass, ind_mass) = osc_nc(imass, ind_mass) + &
!          concentrations(itime, ispec) * species(ispec)%molw * 1.660578e-12
!      enddo
!    enddo



    !allocate(osc_nc(ndat, 5))

    !filename = "osc_nc"
    !call write_2D_array(osc_nc, header, filename)

    !deallocate(header, osc_nc)
  end subroutine calc_osc_nc
!-----------------------------------------------------------------------
  subroutine calc_vbs_parameterization()
!-----------------------------------------------------------------------
! called with flag_vbs_param
! output filename="vbs_parameterization_ppbC_time"
! output filename="vbs_parameterization_ug_time"
!-----------------------------------------------------------------------
  
  integer           :: resolution
  integer           :: n_VBS_param
  integer           :: nbins, ibin, i, j, ispec, iphase, minpvap, maxpvap, LAST_PHASE_shrink
  
  character*(header_length), dimension(:), allocatable :: header_ppbC, header_ug
  integer, dimension(:), allocatable :: upperlim, lowerlim
  real, dimension(:,:), allocatable  :: pvap_distribution_ppbC, pvap_distribution_ug
  real                               :: ppbfac, ugfac, molec_cm3_to_ppb

  real, dimension(:,:), allocatable  :: atom_ratios
  real, dimension(:), allocatable    :: mass_precu, conc_precu
  real, dimension(:,:), allocatable  :: mass_soa, mass_cg, soa_yield, cg_yield
  
  logical, dimension(:), allocatable :: spec_mask, spec_oa_mask, spec_cg_mask

    CHARACTER*(10)    :: date,time1,time2
    CALL date_and_time(date,time1)

    write(6,*) '-- calculating VBS parameterizations --'
  
  ! -----------------------  
  ! adding 6 VBS parameters
  ! -----------------------
  n_VBS_param = 8
  LAST_PHASE_shrink = LAST_PHASE - 6
  
  resolution = 1
    minpvap = -6
    maxpvap = 6
  ! two more bins added: smaller than minpvap, and larger than maxpvap
    nbins = floor(real((maxpvap - minpvap)/resolution, kind = 8)) + 2

    allocate(spec_mask(numsp), spec_oa_mask(numsp), spec_cg_mask(numsp))
    allocate(upperlim(nbins), lowerlim(nbins), &
             pvap_distribution_ppbC(ndat*nbins, LAST_PHASE_shrink+2+n_VBS_param), &
             pvap_distribution_ug(ndat*nbins, LAST_PHASE_shrink+2+n_VBS_param), &
             header_ppbC(LAST_PHASE_shrink+2+n_VBS_param), header_ug(LAST_PHASE_shrink+2+n_VBS_param))    
  allocate(atom_ratios(ndat, 4))
    allocate(mass_precu(ndat), conc_precu(ndat), mass_soa(ndat, nbins), mass_cg(ndat, nbins))
  allocate(soa_yield(ndat, nbins), cg_yield(ndat, nbins))

    do iphase = 1, LAST_PHASE_shrink -1
      select case (iphase)
        case (GAS)
          header_ppbC(iphase) = "Gas Phase [ppbC]"
          header_ug(iphase) = "Gas Phase [ug m-3]"
        case (AEROSOL)
          header_ppbC(iphase) = "Aerosol Phase [ppbC]"
          header_ug(iphase) = "Aerosol Phase [ug m-3]"
        case (PRECURSOR)
          header_ppbC(iphase) = "Precursor [ppbC]"
          header_ug(iphase) = "Precursor [ug m-3]"
        case default
          header_ppbC(iphase) = ""
          header_ug(iphase) = ""
        ! case (CO_CO2)
          ! header_ppbC(iphase) = "CO+CO2 [ppbC]"
          ! header_ug(iphase) = "CO+CO2 [ug m-3]"
        ! case (AQUEOUS)
          ! header_ppbC(iphase) = "Aqueous Phase [ppbC]"
          ! header_ug(iphase) = "Aqueous Phase [ug m-3]"
        ! case (WALL)
          ! header_ppbC(iphase) = "Wall [ppbC]"
          ! header_ug(iphase) = "Wall [ug m-3]"
        ! case (DIMER)
          ! header_ppbC(iphase) = "Dimers [ppbC]"
          ! header_ug(iphase) = "Dimers [ug m-3]"
        ! case (SPINUP_SPEC)
          ! header_ppbC(iphase) = "Spinup Species [ppbC]"
          ! header_ug(iphase) = "Spinup Species [ug m-3]"
        ! case (INORG)
          ! header_ppbC(iphase) = "Inorganics [ppbC]"
          ! header_ug(iphase) = "Inorganics [ug m-3]"
      end select
    enddo
  
    header_ppbC(LAST_PHASE_shrink)   = "Pvap bin lower limit [log10(Cstar)]"
    header_ppbC(LAST_PHASE_shrink+1) = "Pvap bin upper limit [log10(Cstar)]"
    header_ppbC(LAST_PHASE_shrink+2) = "Time [s]"
    header_ug(LAST_PHASE_shrink)   = "Pvap bin lower limit [log10(Cstar)]"
    header_ug(LAST_PHASE_shrink+1) = "Pvap bin upper limit [log10(Cstar)]"
    header_ug(LAST_PHASE_shrink+2) = "Time [s]"

  ! -----------  
  ! VBS headers
  ! -----------
  header_ppbC(LAST_PHASE_shrink+3) = "Kh [M/atm]"
  header_ppbC(LAST_PHASE_shrink+4) = "MW [g/mol]"
  header_ppbC(LAST_PHASE_shrink+5) = "Yield (SOA)"
  header_ppbC(LAST_PHASE_shrink+6) = "Yield (total)"
  header_ppbC(LAST_PHASE_shrink+7) = "Atomic O/C ratio"
  header_ppbC(LAST_PHASE_shrink+8) = "Atomic H/C ratio"
  header_ppbC(LAST_PHASE_shrink+9) = "Atomic N/C ratio"
  header_ppbC(LAST_PHASE_shrink+10) = "Vaporization Enthalpy [J/mol]"
  
  header_ug(LAST_PHASE_shrink+3) = "Kh [M/atm]"
  header_ug(LAST_PHASE_shrink+4) = "MW [g/mol]"
  header_ug(LAST_PHASE_shrink+5) = "Yield (SOA)"
  header_ug(LAST_PHASE_shrink+6) = "Yield (total)"
  header_ug(LAST_PHASE_shrink+7) = "Atomic O/C ratio"
  header_ug(LAST_PHASE_shrink+8) = "Atomic H/C ratio"
  header_ug(LAST_PHASE_shrink+9) = "Atomic N/C ratio"
  header_ug(LAST_PHASE_shrink+10) = "Vaporization Enthalpy [J/mol]"

    do ibin = 1, nbins
     if (ibin==1) then
          lowerlim(ibin) = minpvap + resolution*(ibin-1) - 10
          upperlim(ibin) = minpvap + resolution*(ibin-1)     
     elseif (ibin==nbins) then
          lowerlim(ibin) = minpvap + resolution*(ibin-2)
          upperlim(ibin) = minpvap + resolution*(ibin-2) + 10
     else
       lowerlim(ibin) = minpvap + resolution*(ibin-2)
       upperlim(ibin) = minpvap + resolution*(ibin-1)
     endif
    enddo
  
    molec_cm3_to_ppb = 8.314 * 298.0/101325.0/6.0232e+8
  ugfac = 1.660578e-12

    pvap_distribution_ppbC = 0.
    pvap_distribution_ug = 0.
  conc_precu = 0.0
  mass_precu = 0.0
    mass_soa = 0.0
    mass_cg = 0.0
    soa_yield = 0.0
    cg_yield = 0.0
    atom_ratios = 0.0
    
    
    do i=1,ndat
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink+2) = time(i)
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink+1) = upperlim
       pvap_distribution_ppbC((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink) = lowerlim
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink+2) = time(i)
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink+1) = upperlim
       pvap_distribution_ug((i-1)*nbins+1:i*nbins,LAST_PHASE_shrink) = lowerlim

       ppbfac = 7.245461056e+09*pressure(i)/temperature(i)
       do iphase = 1, LAST_PHASE_shrink-1
         do ibin = 1, nbins
           spec_mask = species%phase == iphase .and. &
                     log10(species%pvap%pvap_Cstar_298) >= lowerlim(ibin) .and. &
                     log10(species%pvap%pvap_Cstar_298) < upperlim(ibin)
           pvap_distribution_ppbC(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%nc/ppbfac, &
              mask = spec_mask)
           pvap_distribution_ug(nbins*(i-1)+ibin, iphase) = sum( &
              concentrations(i, :)*species%molw * ugfac, &
              mask = spec_mask)
         enddo
       enddo
    
     do ispec = 1, numsp
        if (species(ispec)%phase == PRECURSOR) then
         ! write (6,*) species(ispec)%formula, species(ispec)%henry%Keff, species(ispec)%pvap%delta_H_vap_J_per_mol
         ! write (6,*) species(ispec)%formula, species(ispec)%molw
             mass_precu(i) = mass_precu(i) + &
                             concentrations(i, ispec) * &
               species(ispec)%molw * 1.660578e-12 ! conversion to ug m-3
             conc_precu(i) = conc_precu(i) + concentrations(i, ispec) * molec_cm3_to_ppb
      endif
     end do
       ! write (6,*) conc_precu(i), mass_precu(i)


      do ibin = 1, nbins
         spec_mask = (species%phase == GAS .or. &
                  species%phase == AEROSOL) .and. &
                   species%henry%Keff > 0.0 .and. &
                   log10(species%pvap%pvap_Cstar_298) >= lowerlim(ibin) .and. &
               log10(species%pvap%pvap_Cstar_298) < upperlim(ibin)
       ! ---------------------------------------------------------------------------------
         ! calculate mass-weighted mean Henry's law constant & molec weight for each VBS bin
       ! ---------------------------------------------------------------------------------
     if (sum(concentrations(i, :),  mask = spec_mask) > 0.0) then
           pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+3) = sum( &
                  concentrations(i, :)*species%henry%Keff, &
                  mask = spec_mask) / sum(concentrations(i, :),  mask = spec_mask)    
         pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+4) = sum( &
                  concentrations(i, :)*species%molw, &
                  mask = spec_mask) / sum(concentrations(i, :),  mask = spec_mask)  
     endif

       ! -------------------------------------------------------------------
         ! calculate mass-weighted mean Vaporization Enthalpy for each VBS bin
       ! -------------------------------------------------------------------
         ! spec_mask = (species%phase == GAS .or. &
                  ! species%phase == AEROSOL) .and. &
                    ! log10(species%pvap%pvap_atm_298) >= lowerlim(ibin) .and. &
                ! log10(species%pvap%pvap_atm_298) < upperlim(ibin)  
     if (sum(concentrations(i, :),  mask = spec_mask) > 0.0) then
           pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+10) = sum( &
                  concentrations(i, :)*species%pvap%delta_H_vap_J_per_mol, &
                  mask = spec_mask) / sum(concentrations(i, :),  mask = spec_mask)    
     endif
     
       ! -----------------------------------------------------------
         ! calculate mass-weighted mean atomic ratios for each VBS bin
       ! -----------------------------------------------------------
         atom_ratios(i, 1) = sum(concentrations(i,:)*species%no, mask = spec_mask)
         atom_ratios(i, 2) = sum(concentrations(i,:)*species%nh, mask = spec_mask)
         atom_ratios(i, 3) = sum(concentrations(i,:)*species%nn, mask = spec_mask)     
         atom_ratios(i, 4) = sum(concentrations(i,:)*species%nc, mask = spec_mask)     
       if (atom_ratios(i, 4) > 0.0) then
          pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+7) = atom_ratios(i, 1)/atom_ratios(i, 4)
          pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+8) = atom_ratios(i, 2)/atom_ratios(i, 4)
          pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+9) = atom_ratios(i, 3)/atom_ratios(i, 4)
     endif


       ! --------------------------------------------
         ! calculate mass yield of SOA for each VBS bin
       ! --------------------------------------------
     spec_oa_mask = species%phase == AEROSOL .and. &
                   log10(species%pvap%pvap_Cstar_298) >= lowerlim(ibin) .and. &
               log10(species%pvap%pvap_Cstar_298) < upperlim(ibin)
     mass_soa(i, ibin) = sum(concentrations(i, :) * species%molw * 1.660578e-12, &
                             mask = spec_oa_mask)  ! conversion to ug m-3
     ! write (6,*) mass_soa(i, ibin)
     if (i > 1 .and. abs(mass_soa(i, ibin)) > 0.000001) then  ! threshhold: 0.001 ng/m3
        soa_yield(i, ibin) = mass_soa(i, ibin) / (mass_precu(1) - mass_precu(i))
        ! soa_yield(i, ibin) = (mass_soa(i, ibin) - mass_soa(i-1, ibin)) / (mass_precu(i-1) - mass_precu(i))
     endif
         pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+5) = 100.0 * soa_yield(i, ibin)
         ! write (6,*) mass_soa(i, ibin)

       ! ----------------------------------------------------------------
         ! calculate mass yield of total condensable gases for each VBS bin
       ! ----------------------------------------------------------------
     spec_cg_mask = species%phase /= PRECURSOR .and. &
                    log10(species%pvap%pvap_Cstar_298) >= lowerlim(ibin) .and. &
                  log10(species%pvap%pvap_Cstar_298) < upperlim(ibin)
     mass_cg(i, ibin) = sum(concentrations(i, :) * species%molw * 1.660578e-12, &
                             mask = spec_cg_mask)  ! conversion to ug m-3
     if (i > 1 .and. abs(mass_cg(i, ibin)) > 0.000001) then  ! threshhold: 0.001 ng/m3
        cg_yield(i, ibin) = mass_cg(i, ibin) / (mass_precu(1) - mass_precu(i))
        ! cg_yield(i, ibin) = (mass_cg(i, ibin) - mass_cg(i-1, ibin)) / (mass_precu(i-1) - mass_precu(i))
     endif
         pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+6) = 100.0 * cg_yield(i, ibin)
     ! write (6,*) mass_cg(i, ibin)
     
         ! -------------------------------    
     ! VBS outputs will be the same...
         ! -------------------------------     
     do j = 3, 10
        pvap_distribution_ug(nbins*(i-1)+ibin, LAST_PHASE_shrink+j) = pvap_distribution_ppbC(nbins*(i-1)+ibin, LAST_PHASE_shrink+j)
     end do
        
      enddo
    
    enddo
    filename="vbs_parameterization_ppbC_time"
    call write_2D_array(pvap_distribution_ppbC, header_ppbC , filename)
    filename="vbs_parameterization_ug_time"
    call write_2D_array(pvap_distribution_ug, header_ug , filename)

    deallocate(upperlim, lowerlim, spec_mask, &
               pvap_distribution_ppbC, pvap_distribution_ug, &
               header_ppbC, header_ug, &
         atom_ratios, mass_precu, conc_precu, mass_soa, mass_cg, &
         soa_yield, cg_yield, spec_oa_mask, spec_cg_mask)
  
     CALL date_and_time(date,time2)
    ! print*,time1," ",time2

  end subroutine calc_vbs_parameterization

!-------------------------------------------------------------------------- 
  subroutine calc_dyn_filters()
!-------------------------------------------------------------------------- 
! called with flag_dyn_filter
! output filename = "aer_filtered_"//trim(aer_dyn_filter(i))  
! output filename = "gas_filtered_"//trim(gas_dyn_filter(i))  
!-------------------------------------------------------------------------- 
    integer     :: i
    real, allocatable, dimension(:,:) :: output_array
    character*(header_length), dimension(4) :: header
    real :: ppbfac
    real, parameter :: ugfac = 1.660578e-12
    integer, parameter :: ind_time = 1, ind_mass = 2, ind_ppbc = 3, ind_num = 4
    integer  :: itime
    logical, dimension(:), allocatable  :: species_mask
    
    header(ind_time) = "Time [s]"
    header(ind_mass) = "Mass Concentration [ug m-3]"
    header(ind_ppbc) = "Mixing Ratio [ppbC]"
    header(ind_num) = "Number Concentration [molec cm-3]"
    
    allocate(output_array(ndat, 4), &
             species_mask(numsp))
    
    
    ppbfac = 7.245461056e+09*pressure(i)/temperature(i)
    do i = lbound(aer_dyn_filter, 1), ubound(aer_dyn_filter, 1)
      if (aer_dyn_filter(i)(1:1) == ' ') cycle
      write(6,*) " -- treating aerosol dynamic filter #", i, ": ", trim(aer_dyn_filter(i))
      call calc_anyphase_filters(AEROSOL,aer_dyn_filter(i), species_mask)
    
      do itime = 1, ndat
        output_array(itime, ind_time) = time(itime)
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        output_array(itime, ind_num) = sum(concentrations(itime,:), mask = species_mask)
        output_array(itime, ind_ppbc) = sum(concentrations(itime,:) * species%nC/ ppbfac, &
                                          mask = species_mask)
        output_array(itime, ind_mass) = sum(concentrations(itime,:) * species%molw * ugfac, &
                                          mask = species_mask)      
      enddo
    
      filename = "aer_filtered_"//trim(aer_dyn_filter(i))  
      call write_2D_array(output_array, header , filename)
    enddo
    
    do i = lbound(gas_dyn_filter, 1), ubound(gas_dyn_filter, 1)
      if (gas_dyn_filter(i)(1:1) == ' ') cycle
      write(6,*) " -- treating gas dynamic filter #", i, ": ", trim(gas_dyn_filter(i))
      call calc_anyphase_filters(GAS,gas_dyn_filter(i), species_mask)
    
      do itime = 1, ndat
        output_array(itime, ind_time) = time(itime)
        ppbfac = 7.245461056e+09*pressure(itime)/temperature(itime)
        output_array(itime, ind_num) = sum(concentrations(itime,:), mask = species_mask)
        output_array(itime, ind_ppbc) = sum(concentrations(itime,:) * species%nC/ ppbfac, &
                                          mask = species_mask)
        output_array(itime, ind_mass) = sum(concentrations(itime,:) * species%molw * ugfac, &
                                          mask = species_mask)      
      enddo
    
      filename = "gas_filtered_"//trim(gas_dyn_filter(i))  
      call write_2D_array(output_array, header , filename)
    enddo
    
    end subroutine calc_dyn_filters
!------------------------------------------------------------------------    
    subroutine calc_oh_exposure()
!------------------------------------------------------------------------    
! called with flag_ohexposure
! output filename = "oh_exposure"
! output total oh exposure until time t
!------------------------------------------------------------------------    
     
    real, allocatable, dimension(:,:) :: oh_exp
    character*(header_length), dimension(2) :: header
    
    integer, parameter :: ind_time = 1, ind_ohexp = 2
    integer            :: itime, ioh
    real               :: dt, old_oh, new_oh
    
    write(6,*) " -- calculating cumulative OH exposure -- "
    
    ioh = find_species_index("GHO")
    if (ioh == 0) then
      write(6,*) " ! cannot find OH, OH exposure calculation cancelled "
      return
    endif
    
    allocate(oh_exp(ndat,2))
    header(ind_time) = "Time [s]"
    header(ind_ohexp) = "cumulative OH exposure [molec cm-3 s]"
    oh_exp(:,ind_time) = time
    oh_exp(1,ind_ohexp) = 0.
    old_oh = concentrations(1, ioh)
    
    do itime = 2, ndat
      ! integrate oh concentration between itime -1 and itime
      dt = time(itime) - time(itime - 1)
      new_oh = concentrations(itime, ioh)
      oh_exp(itime, ind_ohexp) = oh_exp(itime - 1, ind_ohexp) + dt * (old_oh + new_oh) / 2
      old_oh = new_oh
    enddo
    
    filename = "oh_exposure"
    call write_2D_array(oh_exp, header, filename)
    deallocate(oh_exp)   
   
    
    end subroutine calc_oh_exposure 
!--------------------------------------------------------------------------
    subroutine calc_no3_exposure()
!------------------------------------------------------------------------    
! called with flag_no3exposure
! output filename = "no3_exposure"
! output total no3 exposure until time t
!------------------------------------------------------------------------    
     
    real, allocatable, dimension(:,:) :: no3_exp
    character*(header_length), dimension(2) :: header
    
    integer, parameter :: ind_time = 1, ind_no3exp = 2
    integer            :: itime, ino3
    real               :: dt, old_no3, new_no3
    
    write(6,*) " -- calculating cumulative NO3 exposure -- "
    
    ino3 = find_species_index("GHO")
    if (ino3 == 0) then
      write(6,*) " ! cannot find NO3, NO3 exposure calculation cancelled "
      return
    endif
    
    
    allocate(no3_exp(ndat,2))
    header(ind_time) = "Time [s]"
    header(ind_no3exp) = "cumulative NO3 exposure [molec cm-3 s]"
    no3_exp(:,ind_time) = time
    no3_exp(1,ind_no3exp) = 0.
    old_no3 = concentrations(1, ino3)
    
    do itime = 2, ndat
      ! integrate no3 concentration between itime -1 and itime
      dt = time(itime) - time(itime - 1)
      new_no3 = concentrations(itime, ino3)
      no3_exp(itime, ind_no3exp) = no3_exp(itime - 1, ind_no3exp) + dt * (old_no3 + new_no3) / 2
      old_no3 = new_no3
    enddo
    
    filename = "no3_exposure"
    call write_2D_array(no3_exp, header, filename)
    deallocate(no3_exp)   
    
    end subroutine calc_no3_exposure 

!--------------------------------------------------------------------------
    subroutine calc_tot_reactivity(ksp)
!------------------------------------------------------------------------    
! called with flag_OHR/NO3R/O3R
! output filename = "OH/NO3/O3_reactivity"
! output total OH/NO3/O3 reactivity time series
!------------------------------------------------------------------------    
     
    real, allocatable, dimension(:,:) :: totr
    character*(header_length), dimension(2) :: header
    
    integer, parameter :: ind_time = 1, ind_totr = 2
    integer            :: itime, ispec
    real               :: k_298, totr_i 
    character(LEN=3)   :: ksp
    
    write(6,*) " -- calculating instantaneous net reactivity -- "
    
    allocate(totr(ndat,2))
    header(ind_time) = "Time [s]"
    totr(:,ind_time) = time

    PRINT*,"ksp = ",ksp

    SELECT CASE(ksp)
      CASE("OH ")
        header(ind_totr) = "OHR [s-1]"
      CASE("NO3")
        header(ind_totr) = "NO3R [s-1]"
      CASE("O3 ")
        header(ind_totr) = "O3R [s-1]"
     end SELECT
  
    do itime = 1, ndat
      totr(itime, ind_totr) = 0.
      ! integrate reaction rate over all active species
      do ispec = 1, numsp
        totr_i = 0

        if (species(ispec)%phase .NE. GAS .and. &
            species(ispec)%phase .NE. PRECURSOR ) cycle
        SELECT CASE(ksp)
          CASE("OH ")
            k_298 = species(ispec)%koh%k_298
          CASE("NO3")
            k_298 = species(ispec)%kno3%k_298
          CASE("O3 ")
            k_298 = species(ispec)%ko3%k_298
          CASE DEFAULT
            print*,"calc_tot_reactivity: cannot assign k_298"
            STOP
        end SELECT
        if (k_298 .LT. 0.) cycle

        totr_i = concentrations(itime,ispec) * k_298
        totr(itime, ind_totr) = totr(itime, ind_totr) + totr_i
!        print*,ispec,species(ispec)%code, &
!                     concentrations(itime,ispec), k_298, totr_i
      enddo
    enddo
    
    filename = ksp(1:LEN_TRIM(ksp))//"_reactivity"
    call write_2D_array(totr, header, filename)
    deallocate(totr)   
 
    end subroutine calc_tot_reactivity 

!--------------------------------------------------------------------------
    subroutine calc_prod_reactivity(ksp)
!------------------------------------------------------------------------    
! called with flag_OHR/NO3R/O3R
! output filename = "OH/NO3/O3_reac_prods"
! output product OH/NO3/O3 reactivity time series
!------------------------------------------------------------------------    
     
    real, allocatable, dimension(:,:) :: prodr
    character*(header_length), dimension(2) :: header
    
    integer, parameter :: ind_time = 1, ind_prodr = 2
    integer            :: itime, ispec
    real               :: k_298, prodr_i 
    character(LEN=3)   :: ksp
    
    write(6,*) " -- calculating instantaneous net reactivity -- "
    
    allocate(prodr(ndat,2))
    header(ind_time) = "Time [s]"
    prodr(:,ind_time) = time

    PRINT*,"ksp = ",ksp

    SELECT CASE(ksp)
      CASE("OH ")
        header(ind_prodr) = "OHR [s-1]"
      CASE("NO3")
        header(ind_prodr) = "NO3R [s-1]"
      CASE("O3 ")
        header(ind_prodr) = "O3R [s-1]"
     end SELECT
  
    do itime = 1, ndat
      prodr(itime, ind_prodr) = 0.
      ! integrate reaction rate over all active species
      do ispec = 1, numsp
        prodr_i = 0

        if (species(ispec)%phase .NE. GAS ) cycle
        SELECT CASE(ksp)
          CASE("OH ")
            k_298 = species(ispec)%koh%k_298
          CASE("NO3")
            k_298 = species(ispec)%kno3%k_298
          CASE("O3 ")
            k_298 = species(ispec)%ko3%k_298
          CASE DEFAULT
            print*,"calc_prod_reactivity: cannot assign k_298"
            STOP
        end SELECT
        if (k_298 .LT. 0.) cycle

        prodr_i = concentrations(itime,ispec) * k_298
        prodr(itime, ind_prodr) = prodr(itime, ind_prodr) + prodr_i
!        print*,ispec,species(ispec)%code, &
!                     concentrations(itime,ispec), k_298, prodr_i
      enddo
    enddo
    
    filename = ksp(1:LEN_TRIM(ksp))//"_reac_prods"
    call write_2D_array(prodr, header, filename)
    deallocate(prodr)   
 
    end subroutine calc_prod_reactivity 

!--------------------------------------------------------------------------
!------------------------------------------------------------------------    
    subroutine calc_spec_reactivity(ksp)
!------------------------------------------------------------------------    
! called with flag_OHR/NO3R/O3R
! output filename = "OH/NO3/O3_top_spc_reac_time"
! output OH/NO3/O3 reactivity time series, speciated by top species
!------------------------------------------------------------------------    
     
    integer, dimension(:), allocatable  :: top_species_indx
    character*(maxlsp), dimension(:), allocatable :: top_species_names
    real, dimension(:,:), allocatable   :: top_species_reac

    real               :: k_298
    character(LEN=3)   :: ksp

    character*(header_length) :: smile
    character*(header_length), dimension(:), allocatable :: header1,header2
    real, dimension(:), allocatable     :: time_integrated_reac
    integer                             :: N, itime, itop, ispec
    logical                             :: found_fg

    CHARACTER*(10)    :: date,time1,time2
    CHARACTER*(filenames_length) :: filename1,filename2
    CALL date_and_time(date,time1)

    N = n_topspecies
    allocate(top_species_names(N), top_species_indx(N), &
             top_species_reac(ndat, N +2), &
             time_integrated_reac(numsp), &
             header1(N+2),header2(N+2))


    write(6,*) "-- calculating reactivity of top ", N, " species (s-1) --"

    ! gas phase + precursor only
    ! first we need to find the top N-1 species (Nth species is "Others")
    ! with respect to time integrated mixing ratio in molec/cc (ppb)
    
    PRINT*,"ksp = ",ksp
    top_species_reac = 0
   time_integrated_reac = 0

      ! integrate reaction rate over all active species
    do ispec = 1, numsp

      if (species(ispec)%phase .NE. GAS .and. &
          species(ispec)%phase .NE. INORG .and. &
          species(ispec)%phase .NE. PRECURSOR ) cycle

      SELECT CASE(ksp)
        CASE("OH ")
          k_298 = species(ispec)%koh%k_298
        CASE("NO3")
          k_298 = species(ispec)%kno3%k_298
        CASE("O3 ")
          k_298 = species(ispec)%ko3%k_298
        CASE DEFAULT
          print*,"calc_spec_reactivity: cannot assign k_298"
          STOP
      end SELECT
      if (k_298 .LT. 0.) cycle

      do itime = 1, ndat
!        print*,ispec,species(ispec)%code, &
!                     concentrations(itime,ispec), k_298, totr_i

        time_integrated_reac(ispec) = time_integrated_reac(ispec)+ &
                                      concentrations(itime,ispec)*k_298

      enddo
    enddo

    top_species_indx(1) = maxloc(time_integrated_reac, dim = 1)
    top_species_names(1) = species(top_species_indx(1))%code
    time_integrated_reac(top_species_indx(1)) = 0
    do itop = 2, N
      top_species_indx(itop) = maxloc(time_integrated_reac, dim = 1)
      top_species_names(itop) = species(top_species_indx(itop))%code
      time_integrated_reac(top_species_indx(itop)) = 0
    enddo

    do itop = 1, N
      header1(itop) = trim(species(top_species_indx(itop))%formula)// &
                 ":"//trim(species(top_species_indx(itop))%code)
      if (flag_smiles) then
        CALL smiles(species(top_species_indx(itop))%formula,smile)
        ispec=INDEX(smile," ")-1
        header2(itop) = trim(smile(1:ispec))// &
                   ":"//trim(species(top_species_indx(itop))%code)
      endif
    enddo
    header1(N+1) = "Others [s-1]"
    header1(N+2) = "Time [s]"
    header2(N+1) = "Others [s-1]"
    header2(N+2) = "Time [s]"

    do itime = 1, ndat
      do ispec=  1, numsp
        found_fg = .false.
        if (species(ispec)%phase .NE. GAS .and. &
            species(ispec)%phase .NE. PRECURSOR ) cycle

        SELECT CASE(ksp)
          CASE("OH ")
            k_298 = species(ispec)%koh%k_298
          CASE("NO3")
            k_298 = species(ispec)%kno3%k_298
          CASE("O3 ")
            k_298 = species(ispec)%ko3%k_298
          CASE DEFAULT
            print*,"calc_tot_reactivity: cannot assign k_298"
            STOP
        end SELECT
        if (k_298 .LT. 0.) cycle

        do itop = 1, N
          if (ispec == top_species_indx(itop)) then
            found_fg = .true.
            top_species_reac(itime,itop) = concentrations(itime, ispec)*k_298
            exit
          endif
        enddo
        if (.not. found_fg) then
          top_species_reac(itime, N+1) = top_species_reac(itime, N+1) + &
                                concentrations(itime, ispec)*k_298
        endif
      enddo
    enddo

    top_species_reac(:,N+2) = time

! filename chars
    SELECT CASE(ksp)
      CASE("OH ")
        k_298 = species(ispec)%koh%k_298
      CASE("NO3")
        k_298 = species(ispec)%kno3%k_298
      CASE("O3 ")
        k_298 = species(ispec)%ko3%k_298
      CASE DEFAULT
        print*,"calc_tot_reactivity: cannot assign k_298"
        STOP
    end SELECT

! chemical formulae headers
    filename1 = ksp(1:LEN_TRIM(ksp))//"_top_spc_reac_time"
    call write_2D_array(top_species_reac, header1 , filename1)
! SMILES headers
    if (flag_smiles) then
        filename2 = ksp(1:LEN_TRIM(ksp))//"_top_sml_reac_time"
      call write_2D_array(top_species_reac, header2 , filename2)
    endif

    deallocate(top_species_names, top_species_indx, &
               top_species_reac, time_integrated_reac, &
               header1,header2)
 
    end subroutine calc_spec_reactivity 

!--------------------------------------------------------------------------
    subroutine calc_potential_frag_dimer_atom_ratio()
!--------------------------------------------------------------------------
! called with flag_potential_atomratios
! output filename = "atom_ratios_aer_frag"
! output filename = "atom_ratios_aer_dimer"
! output what would O/C, H/C be if all C10 species fragmented or dimerized in aerosol phase
!--------------------------------------------------------------------------
      character*(header_length), dimension(4) :: header

      real, dimension(:,:), allocatable    :: atom_ratios_frag, atom_ratios_dimer
      real, dimension(:), allocatable      :: sumnc
      integer                              :: itime, ispec
      integer, parameter                   :: ind_oc = 1, ind_hc = 2, ind_nc = 3, ind_time = 4

      CHARACTER*(10)    :: date,time1,time2
      CALL date_and_time(date,time1)

      write(6,*) '-- calculating O/C, H/C and N/C ratios for frag and dimers --'
      allocate(atom_ratios_frag(ndat, 4), atom_ratios_dimer(ndat,4), sumnc(ndat))

      header(ind_oc)   = "O/C"
      header(ind_hc)   = "H/C"
      header(ind_nc)   = "N/C"
      header(ind_time) = "Time [s]"

      sumnc = 0.

     atom_ratios_dimer = 0.
     atom_ratios_frag = 0.
  ! in the aerosol phase
      do itime = 1, ndat
       atom_ratios_dimer(itime, ind_time) = time(itime)
       atom_ratios_frag(itime, ind_time) = time(itime)       
        do ispec = 1, numsp
          if (species(ispec)%phase /= AEROSOL) cycle
          sumnc(itime) = sumnc(itime) + concentrations(itime, ispec)*species(ispec)%nc
          atom_ratios_dimer(itime, ind_nc) = atom_ratios_dimer(itime, ind_nc) + &
                                             concentrations(itime, ispec)*species(ispec)%nn
          atom_ratios_frag(itime, ind_nc) = atom_ratios_frag(itime, ind_nc) + &
                                            concentrations(itime, ispec)*species(ispec)%nn
          if (species(ispec)%nc == 10 .and. species(ispec)%no > 4 .and. species(ispec)%nh > 1) then
            atom_ratios_dimer(itime, ind_oc) = atom_ratios_dimer(itime, ind_oc) + &
                                               concentrations(itime, ispec)*(2*species(ispec)%no - 4.)/2.
            atom_ratios_dimer(itime, ind_hc) = atom_ratios_dimer(itime, ind_hc) + &
                                               concentrations(itime, ispec)*(2*species(ispec)%nh - 1.)/2.
            atom_ratios_frag(itime, ind_oc) = atom_ratios_frag(itime, ind_oc) + &
                                              concentrations(itime, ispec)*(species(ispec)%no + 2.)
            atom_ratios_frag(itime, ind_hc) = atom_ratios_frag(itime, ind_hc) + &
                                              concentrations(itime, ispec)*(species(ispec)%nh - 2.)
          else
            atom_ratios_dimer(itime, ind_oc) = atom_ratios_dimer(itime, ind_oc) + &
                                               concentrations(itime, ispec)*species(ispec)%no
            atom_ratios_dimer(itime, ind_hc) = atom_ratios_dimer(itime, ind_hc) + &
                                               concentrations(itime, ispec)*species(ispec)%nh
            atom_ratios_frag(itime, ind_oc) = atom_ratios_frag(itime, ind_oc) + &
                                              concentrations(itime, ispec)*species(ispec)%no
            atom_ratios_frag(itime, ind_hc) = atom_ratios_frag(itime, ind_hc) + &
                                              concentrations(itime, ispec)*species(ispec)%nh
          endif
        enddo
        atom_ratios_dimer(itime, ind_oc) = atom_ratios_dimer(itime, ind_oc) / sumnc(itime)
        atom_ratios_dimer(itime, ind_hc) = atom_ratios_dimer(itime, ind_hc) / sumnc(itime)
        atom_ratios_dimer(itime, ind_nc) = atom_ratios_dimer(itime, ind_nc) / sumnc(itime)
        atom_ratios_frag(itime, ind_oc) = atom_ratios_frag(itime, ind_oc) / sumnc(itime)
        atom_ratios_frag(itime, ind_hc) = atom_ratios_frag(itime, ind_hc) / sumnc(itime)
        atom_ratios_frag(itime, ind_nc) = atom_ratios_frag(itime, ind_nc) / sumnc(itime)
      enddo

      filename = "atom_ratios_aer_frag"
      call write_2D_array(atom_ratios_frag, header , filename)
      filename = "atom_ratios_aer_dimer"
      call write_2D_array(atom_ratios_dimer, header , filename)


      deallocate(atom_ratios_dimer, atom_ratios_frag, sumnc)
      CALL date_and_time(date,time2)

    end subroutine calc_potential_frag_dimer_atom_ratio

END MODULE compute
