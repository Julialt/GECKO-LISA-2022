MODULE IO
  USE PARAMETERS
  USE DICO
  USE CONC
  USE RATE
  USE USER_INPUT
  USE SORTING
  USE NCUTILS
  IMPLICIT NONE
  integer, parameter                           :: read_unit = 19, write_unit = 20
  integer, parameter                           :: filenames_length = 400
  integer, parameter                           :: header_length = lfo+lco
  character*(filenames_length)                 :: filename
  character*(10)                               :: boxname

  CONTAINS
  SUBROUTINE get_userinfo()
    integer                                    :: i

      ! default values if not provided
    nbox = 1
    nprecu = 0
    n_topspecies = 10
    n_toprates = 5
    input_type = "netcdf"
    !input_type = "binary"
    output_type = "ascii"
    !output_type = "netcdf"
    precursor_codes = " "
    selected_species = ''
    rate_species = ''

    gas_dyn_filter = " "
    aer_dyn_filter = " "
    seed_mass = 10.
    seed_molw = 250.
    skip_time = 1
    flag_amsfactors = .false.
    flag_atomratios = .false.
    flag_bubble = .false.
    flag_carbonchain = .false.
    flag_chochonfreq = .false.
    flag_chon = .false.
    flag_contributingspecs = .false.
    flag_cstar = .false.
    flag_dbeai = .false.
    flag_dyn_filter = .false.
    flag_elementscontrib = .false.
    flag_entropy = .false.
    flag_functions = .false.
    flag_henry = .false.
    flag_massspectrum = .false.
    flag_ml = .false.
    flag_nitrates = .false.
    flag_NO3R = .false.
    flag_O3R = .false.
    flag_OHR = .false.
    flag_ohexposure = .false.
    flag_no3exposure = .false.
    flag_phasedist = .false.
    flag_phasestate = .false.
    flag_potential_atomratios = .false.
    flag_ppb = .false.
    flag_ppbC = .false.
    flag_pvap = .false.
    flag_rates = .false.
    flag_selected = .false.
    flag_smiles = .false.
    flag_soayield = .false.
    flag_topspec = .false.
    flag_topNspc = .false.

    OPEN(read_unit, file = "userinput.nml")
    READ(read_unit, nml = userinput)
    CLOSE(read_unit)

    nprecu = count(precursor_codes .ne. " ")
    if (nprecu == 0) then
      write(6,*) ' --error-- no precursor code was initialized in '
      write(6,*) ' userinput.nml'
      stop
    endif
    do i = 1, nprecu
      precursor_codes(i) = trim(precursor_codes(i))
      PRINT*,"PRECU= ",precursor_codes(i)
    enddo
    
    if (count(gas_dyn_filter .ne. " ") + count(aer_dyn_filter .ne. " ") > 0) then
      flag_dyn_filter = .true.
    endif

! make sure at least one primary selector is present
! NB: flag_smiles, flag_ppb, flag_ppbC are all secondary selectors
    if (.not.(flag_chon .or. flag_selected .or. flag_phasedist .or. flag_contributingspecs .or. flag_soayield .or. &
      flag_functions .or. flag_carbonchain .or. flag_topspec .or. flag_atomratios .or. flag_pvap.or.flag_cstar .or. &
      flag_henry .or. flag_elementscontrib .or. flag_phasestate .or. flag_entropy .or. flag_chochonfreq .or. flag_dbeai .or. &
      flag_massspectrum .or. flag_amsfactors .or. flag_bubble .or. flag_nitrates .or. flag_dyn_filter .or. flag_ohexposure .or. &
      flag_no3exposure .or. flag_potential_atomratios .or. flag_ml .or. &
      flag_OHR .or. flag_NO3R .or. flag_O3R .or. flag_rates .or. flag_topNspc)) then
      write(6,*) ' --error-- no posttreatment flag was activated in'
      write(6,*) ' userinput.nml'
      write(6,*) 'we should stop here'
      stop
    endif

  END SUBROUTINE

  SUBROUTINE read_dict()

    integer                                                   :: i
    character*(ldi)                                           :: line
  !find number of species

    write(6,*) ' opening: ', trim(filename)
    open(read_unit, file = filename, status = 'OLD')
    ndic = 0
    do i = 1, maxsp
      line = ""
      read(read_unit, '(a)') line
      if (line(1:4) == "****") exit
      if (line(1:5) == '     ' .or. line(1:5) == 'HV   ' .or. &
          line(1:5) == 'M    ' .or. line(1:5) == '(M)  ') cycle
      ndic = ndic +1
    enddo

    write(6,*) 'found ', ndic, ' dictionary species in ', trim(filename)
    if (ndic >= maxsp) then
      write(6,*) ' -- error -- too many species in dictionary '
      STOP
    endif

    write(6,*) 'reading ', trim(filename)
    if( allocated(dictionary)) deallocate(dictionary)
    allocate(dictionary(ndic))
    rewind(read_unit)
    do i = 1, maxsp
      line = ""
      read(read_unit,'(a)') line
      if (line(1:4) == "****") exit
      if (line(1:5) == '     ' .or. line(1:5) == 'HV   ' .or. &
          line(1:5) == 'M    ' .or. line(1:5) == '(M)  ') cycle
      read(line, '(A6, 3X, A120, 2X, A15)') dictionary(i)%code, dictionary(i)%formula, dictionary(i)%functions
    enddo

    close(read_unit)

  END SUBROUTINE read_dict


  SUBROUTINE read_ppf()

    integer                                      :: numextra, lensp, dummy1, dummy2, i, j
    integer                                      :: counter
    character*(maxlsp)                           :: namextra1, namextra2, namextra3
    real                                         :: t0, dummy

    write(6,*) ' opening: ', trim(filename)
    open(read_unit, file = filename, status = 'OLD', form = 'UNFORMATTED')

    ! read parameters of the file and check the size
    read(read_unit) numsp, numextra, lensp, dummy1, dummy2
    write(6,*) 'found ', numsp, ' mechanism species in ', trim(filename)

    if (numsp >= maxsp) then
      write(6,*) ' -- error -- too many species in ppf file '
      STOP
    endif

    if (numextra /= 3) then
      WRITE(6,*) '--error--, length of the header not correct'
      STOP
    endif

    if(.not. ALLOCATED(species)) allocate(species(numsp))

    ! read header
    read(read_unit) namextra1, namextra2, namextra3, (species(i)%code, i=1, numsp)


    IF (namextra1(1:5).ne.'TIME ') THEN
      WRITE(6,*) '--error--, first header expected to be TIME'
      STOP
    ENDIF

    IF (namextra2(1:6).ne.'TEMPER') THEN
      WRITE(6,*) '--error--, second header expected to be TEMPERATURE'
      STOP
    ENDIF

    IF (namextra3(1:6).ne.'HUMIDI') THEN
      WRITE(6,*) '--error--, second header expected to be HUMIDITY'
      STOP
    ENDIF

    ndat = 0
    counter = 0
    t0 = 0.
    dummy = 0.
    ! we read the full first hour and every skip_time-th timestep after that
    do
      read(read_unit, end = 50) dummy
      counter = counter + 1
      if ( counter == 1) then 
        t0 = dummy
      endif      
      if ( dummy <= (t0 + 3600) .or. modulo(counter,skip_time) == 0) then
        ndat = ndat+1
      endif
    enddo

50 write(6,*) 'found ', ndat, ' data points in ', trim(filename)

   if (ndat >= mdat) then
      write(6,*)  ' -- error -- too many data points in ppf file '
      STOP
    endif

    if(allocated(concentrations)) then
      deallocate(concentrations, time, temperature, pressure, rh)
    endif
    allocate(concentrations(ndat, numsp), time(ndat), temperature(ndat), pressure(ndat), rh(ndat))

    rewind(read_unit)
    write(6,*) 'reading ', trim(filename)
    read(read_unit) ! skip header
    read(read_unit) ! skip header
    counter = 0
    i = 1
    do
      counter = counter + 1
      read(read_unit) time(i), temperature(i), rh(i), (concentrations(i, j), j =1, numsp)
      if(time(i) <= (t0 + 3600) .or. modulo(counter,skip_time) == 0) then
        i = i + 1
        if (i == (ndat + 1)) exit
      endif
    enddo
    

    pressure = 1013.25

    close(read_unit)
  END SUBROUTINE read_ppf


  SUBROUTINE read_pvap()
    real, dimension(:), allocatable    :: Tb, dB
    character*(maxlsp), dimension(:), allocatable :: name
    integer                            :: npvap, i,j, ispec
  real                               :: pvap_atm_280_temp

! initialize
    species%pvap%pvap_atm_298 = 0.
    species%pvap%pvap_Cstar_298 = 0.
    species%pvap%delta_H_vap_J_per_mol = 0.

    open(read_unit, file = filename, status = 'OLD')
    npvap = 0
    do
      read(read_unit, fmt=*, end = 60)
      npvap = npvap+1
    enddo
60  npvap = npvap - 2 ! account for header and END
    write(6,*) 'found ', npvap, ' pvap values in ', trim(filename)

    if (npvap == 0) return
    allocate(name(npvap), Tb(npvap), dB(npvap))

    rewind(read_unit)
    read(read_unit, fmt=*)
    do i =1, npvap
    !three columns to read, name, Tb, dB
      read(read_unit,'(A7, 2x, f6.1, 2x, f8.4)') name(i), Tb(i), dB(i)
! DEBUG
      !print*, name(i), Tb(i), dB(i)
! END DEBUG
    enddo

    do i = 1, npvap
      do j = 1, len_trim(phase_letter)
        ispec = find_species_index(phase_letter(j:j)//name(i)(2:7))
        !if (ispec .ne. -1) then
        if (ispec .gt. 0) then
          species(ispec)%pvap%Tb = Tb(i)
          species(ispec)%pvap%dB = dB(i)
          species(ispec)%pvap%pvap_atm_298 = &
            exp((4.1012 + dB(i)) * ((298/Tb(i)-1)/(298/Tb(i) - 0.125))*log(10.))
          !species(ispec)%pvap%pvap_cstar_298 = &
          !  (1e6*seed_molw*seed_mass*species(ispec)%pvap%pvap_atm_298)/(8.2e-5*298)
! Cstar corrected per CU work
          call get_number(species(ispec)%formula, species(ispec)%molw,&
                      species(ispec)%nc, species(ispec)%nh, species(ispec)%nn, &
                      species(ispec)%no, species(ispec)%nr, species(ispec)%ns, &
                      species(ispec)%nfl, species(ispec)%nbr, species(ispec)%ncl)

          species(ispec)%pvap%pvap_cstar_298 = &
            (1e6*species(ispec)%molw*species(ispec)%pvap%pvap_atm_298)/ &
            (8.2e-5*298)

    ! I'm too lazy to solve the clausius-clapeyron equation so here's a quick and dirty way
    ! Siyuan Wang (siyuan@ucar.edu)
    pvap_atm_280_temp = exp((4.1012 + dB(i)) * ((280./Tb(i)-1.)/(280./Tb(i) - 0.125))*log(10.))
    species(ispec)%pvap%delta_H_vap_J_per_mol = & 
                        log(pvap_atm_280_temp/species(ispec)%pvap%pvap_atm_298) * &
                        8.314 / (1./298. - 1./280.)
        else
          write(6,*) " could not find: ", phase_letter(j:j)//name(i)(2:7)
          write(6,*) " in subroutine read_pvap"
        endif
      enddo

    enddo

  END SUBROUTINE read_pvap
!-------------------------------------------------------------------
  SUBROUTINE read_henry
    real, dimension(:), allocatable    :: cf, Keff, rf
    character*(maxlsp), dimension(:), allocatable :: name
    integer                            :: nhenry, i,j, ispec

    open(read_unit, file = filename, status = 'OLD')
    nhenry = 0
    do
      read(read_unit, fmt=*, end = 70)
      nhenry = nhenry+1
    enddo
70  nhenry = nhenry -11 ! account for END and 10 lines preamble
   write(6,*) 'found ', nhenry, ' henry values in ', trim(filename)

    if (nhenry == 0) return
    allocate(name(nhenry), cf(nhenry), Keff(nhenry), rf(nhenry) )

    rewind(read_unit)
    ! skip the first 10 lines
    do i =1, 10
      read(read_unit,*)
    enddo
    do i =1, nhenry
    !four columns to cf, Keff, rf, name
      read(read_unit,'(E7.1,2x,E7.1,2x,E7.1,2x,a7,1x)') cf(i), Keff(i), rf(i), name(i)
    enddo


    do i = 1, nhenry
      do j=1, len_trim(phase_letter)
        ispec = find_species_index(phase_letter(j:j)// &
                name(i)(2:len_trim(name(i))))
        if (ispec .gt. 0) then
          species(ispec)%henry%cf = cf(i)
          species(ispec)%henry%Keff = Keff(i)
          species(ispec)%henry%rf = rf(i)
        else if (ispec .lt. 0) then
          write(6,*) " could not find: ", phase_letter(j:j)//name(i)(2:7)
          write(6,*) " in subroutine read_henry"
        endif
      enddo
    enddo

    deallocate(name, cf, Keff, rf )

  END SUBROUTINE read_henry

!-------------------------------------------------------------------
  SUBROUTINE read_ncdf_dict(ncid)
!! subroutine needs to know values of: {maxlsp,lfl,lfo} = {6,120,15} ! OK
!! subroutine needs to access routines in ncutil.f !OK
!! subroutine DOES NOT need to "use netcdf" !OK

    integer                                  :: ncid,varid,status

  !read number of species
    CALL eznc_get_dimension(ncid,"mxdic",ndic)

    write(6,*) 'found ', ndic, ' dict species in ', trim(filename)
    if (ndic >= maxsp) then
      write(6,*) ' -- error -- too many species in dictionary '
      STOP
    endif

    write(6,*) 'reading ', trim(filename)
    if( allocated(dictionary)) deallocate(dictionary)

    allocate(dictionary(ndic))

   !read dictionary names, formulae, functional group codes
   
!    CALL eznc_get_1Dchar(ncid,"dicnam",lco,ndic,dictionary(:)%code,1,ndic)
!    CALL eznc_get_1Dchar(ncid,"chem",lfo,ndic,dictionary(:)%formula,1,ndic)
!    CALL eznc_get_1Dchar(ncid,"code",lfl,ndic,dictionary(:)%functions,1,ndic)
    CALL eznc_get_1Dchar(ncid,"dicnam",lco,dictionary(:)%code,1,ndic,1)
    CALL eznc_get_1Dchar(ncid,"chem",lfo,dictionary(:)%formula,1,ndic,1)
    CALL eznc_get_1Dchar(ncid,"code",lfl,dictionary(:)%functions,1,ndic,1)

    status = NF90_INQ_VARID(ncid,"igen",varid)
    IF (status==NF90_NOERR) &
    !CALL eznc_get_1Dint(ncid,"igen",ndic,dictionary(:)%igen,1,ndic,1)
    CALL eznc_get_1Dint(ncid,"igen",dictionary(:)%igen,1,ndic,1)

  END SUBROUTINE read_ncdf_dict

!-------------------------------------------------------------------
  SUBROUTINE read_ncdf_ppf(ncid, ibox)
   USE netcdf
!! subroutine needs to know values of: {nbox} = {?} !OK
!! subroutine needs to access routines in ncutil.f ! OK
!! subroutine DOES NOT need to "use netcdf" ! OK

    integer, intent(in)                          :: ibox
    integer                                      :: ncid,varid,status
    integer                                      :: mbox
    integer                                      :: sumcfix
    integer                                      :: isp,idat
    !integer                                      :: ndat
    integer                                      :: numconc
    real    :: seedval
!    real(kind=8),dimension(:,:),allocatable      :: temp_conc
    real(kind=8),dimension(:,:,:),allocatable    :: temp_conc
    real(kind=8),dimension(:,:),allocatable      :: t_temperature
    real(kind=8),dimension(:),allocatable        :: sumc
!    integer                                      :: numextra
!    character*(maxlsp)                           :: namextra1, namextra2, namextra3

    ! read parameters of the file and check the size
    write(6,*) 'reading ', trim(filename)

    CALL eznc_get_0Dint(ncid,"numsp",numsp)
!    CALL eznc_get_0Dint(ncid,"numextra",numextra)
    !?? we are not reading lensp (it's not currently used in NetCDF file - should it be?) or dummy1/2.

    write(6,*) 'found ', numsp, ' mech species in ', trim(filename)

    if (numsp >= maxsp) then
      write(6,*) ' -- error -- too many species in output file '
      STOP
    endif

    IF(.not. ALLOCATED(species)) allocate(species(numsp))

    ! read header
!    read(read_unit) namextra1, namextra2, namextra3, (species(i)%code, i=1, numsp)
!    IF (namextra1(1:5).ne.'TIME ') (etc) => "time"
!    IF (namextra2(1:6).ne.'TEMPER') (etc) => "tempbot" or "temptop"
    !! humidity is not currently recorded in NetCDF file - TO DO !!
!    IF (namextra3(1:6).ne.'HUMIDI') (etc) => ""

    !CALL eznc_get_1Dchar(ncid,"chrsp",maxlsp,numsp,species(:)%code,1,numsp)
    CALL eznc_get_1Dchar(ncid,"chrsp",maxlsp,species(:)%code,1,numsp,1)

    ! read number of model output times ("ntout") => rddat
    CALL eznc_get_dimension(ncid,"ntout",rddat)

    ! read number of boxes
    call eznc_get_dimension(ncid,"mbox",mbox)
    call eznc_get_0Dint(ncid, "nbox", nbox)

   write(6,*) 'found ', rddat, ' data points in ', trim(filename)

   if (rddat >= mdat) then
      write(6,*)  ' -- error -- too many data points in output file '
      STOP
    endif

   ! define ndat (number of times to use in postproc) according to 
   ! user input variable skip_time
   ! where skip_time = 1 means "use all output times"
   IF(skip_time.gt.1)THEN
     ndat = INT((rddat-1)/skip_time)+1
   ELSE
     ndat = rddat
   ENDIF
   numconc = ndat*numsp
   write(6,*) 'calculating',ndat,"*",numsp,"=",numconc,'conc points in ',trim(filename)

!   if (numconc >= mconc) then
!      write(6,*)  ' -- error -- too many conc points in output file '
!      STOP
!    endif

    if (ibox > mbox) then
      write(6, *) ' -- error -- trying to read from a box that does not exist'
      write(6, *) ibox, '>', mbox
      stop
    endif

    if(allocated(rh)) deallocate(rh)
    if(allocated(sza)) deallocate(sza)
    if(allocated(sumc)) deallocate(sumc)
    if(allocated(time)) deallocate(time)
    if(allocated(seed)) deallocate(seed)
    if(allocated(pressure)) deallocate(pressure)
    if(allocated(temp_conc)) deallocate(temp_conc)
    if(allocated(temperature)) deallocate(temperature)
    if(allocated(t_temperature)) deallocate(t_temperature)
    if(allocated(concentrations)) deallocate(concentrations)

    allocate(concentrations(ndat,numsp), time(ndat), temperature(ndat), pressure(ndat),&
             rh(ndat),  sza(ndat), seed(ndat))
             !rh(ndat), temp_conc(numsp,1,ndat), sza(ndat), seed(ndat))

    PRINT*,"time"
    !CALL eznc_get_1Dreal(ncid,"time",ndat,time(1:ndat),1,ndat)
      CALL eznc_get_1Dreal(ncid,"time",time(1:ndat),1,rddat,skip_time)
    PRINT*,time
!
    PRINT*,"temperature"
    !call eznc_get_2Dreal(ncid,"temp",mbox,ndat,temperature(1:ndat),ibox,ibox,1,ndat)
      allocate(t_temperature(1,ndat))
      call eznc_get_2Dreal(ncid,"temp",t_temperature(1,1:ndat),&
                                     ibox,ibox,1,1,rddat,skip_time)
      temperature(:) = t_temperature(1,:)
      deallocate(t_temperature)
    !PRINT*,temperature(1:ndat)

!    call eznc_get_3Dreal(ncid, "conc",numsp, ndat, nbox, &
!                              temp_conc(1:numsp,          1:ndat), &
!                                        1,numsp,ibox,ibox,1,ndat)
    PRINT*,"conc"
      allocate(temp_conc(numsp,1,ndat))
      call eznc_get_3Dreal(ncid,"conc", temp_conc, & !(1:numsp,1,1:ndat), &
                                        1,numsp,1, &
                                        ibox,ibox,1, &
                                        1,rddat,skip_time)
      !DO isp=1,numsp
      DO idat=1,ndat
        concentrations(idat,1:numsp)=temp_conc(1:numsp,1,idat)
      ENDDO
      !ENDDO
      deallocate(temp_conc)
    PRINT*,concentrations(1:ndat,1)

    call eznc_get_0Dreal(ncid, "cnv", seedval)
    seed = seedval

    status = NF90_INQ_VARID(ncid,"rh",varid)
    IF (status==NF90_NOERR) &
    !call eznc_get_2Dreal(ncid,"rh",mbox,ndat,rh(1:ndat),ibox,ibox,1,ndat)
    call eznc_get_2Dreal(ncid,"rh",rh(1:ndat),ibox,ibox,1,1,rddat,skip_time)

    status = NF90_INQ_VARID(ncid,"sza",varid)
    IF (status==NF90_NOERR) &
    !call eznc_get_2Dreal(ncid, "sza",  mbox, ndat, sza(1:ndat),      ibox, ibox, 1, ndat)
    call eznc_get_2Dreal(ncid, "sza",sza(1:ndat),ibox,ibox,1,1,rddat,skip_time)

    status = NF90_INQ_VARID(ncid,"sumcfix",varid)
    IF (status==NF90_NOERR) &
    CALL eznc_get_0Dint(ncid,"sumcfix",sumcfix)


! pressure

! If supplies as pres(bar), convert atm to hPa
    status = NF90_INQ_VARID(ncid,"pres",varid)
    IF (status==NF90_NOERR) THEN
      !call eznc_get_2Dreal(ncid, "pres", mbox, ndat, pressure(1:ndat), ibox, ibox, 1, ndat)
      call eznc_get_2Dreal(ncid,"pres",pressure(1:ndat),ibox,ibox,1,1,rddat,skip_time)
      pressure = pressure * 1000
    ELSE

! Else supplied as sumc(molec/cc), 
! calculate pressure in bar & convert to mbar(hPa)
      IF(sumcfix.EQ.1)THEN
        allocate (sumc(1))
        !call eznc_get_1Dreal(ncid, "sumc", mbox, sumc(ibox:ibox), ibox, ibox )
        call eznc_get_1Dreal(ncid,"sumc",sumc(ibox:ibox),ibox,ibox,skip_time)
        !pressure(1:ndat) = sumc(1)*8.31446*temperature(1:ndat)/6.022E+22
        pressure(1:ndat) = sumc(1)*8.31446*temperature(1:ndat)/6.022E+22
        pressure = pressure*1000
      ELSE
        !allocate (sumc(ndat))
        allocate (sumc(ndat))
        !call eznc_get_2Dreal(ncid, "sumc", mbox, ndat, sumc(1:ndat), ibox, ibox, 1, ndat)
        !pressure(1:ndat) = sumc(1:ndat)*8.31446*temperature(1:ndat)/6.022E+22
        call eznc_get_2Dreal(ncid,"sumc",sumc,ibox,ibox,1,1,rddat,skip_time)
        pressure(1:ndat) = sumc(1:ndat)*8.31446*temperature(1:ndat)/6.022E+22
        pressure = pressure*1000
      ENDIF
    ENDIF

    !concentrations = transpose(temp_conc)

  END SUBROUTINE read_ncdf_ppf

!-------------------------------------------------------------------

  SUBROUTINE read_ncdf_pvap(ncid)
!! subroutine needs to know values of: {maxlsp, mxsat} ! OK
!! subroutine needs to access routines in ncutil.f ! OK
!! subroutine DOES NOT need to "use netcdf" !OK
    integer, dimension(:), allocatable :: satid
    real, dimension(:), allocatable    :: Tb, dB
    real, dimension(:,:), allocatable  :: satdat
!    character*(maxlsp), dimension(:), allocatable :: satname
    integer                            :: npvap, i, ispec !, j
    integer                            :: ncid,varid,status
    integer                            :: mxsat
    real                               :: pvap_atm_280_temp

    species%pvap%Tb = -9999.
    species%pvap%dB = -9999.
    species%pvap%pvap_atm_298 = -9999.

! read max # of pvap ("mxsat")
    CALL eznc_get_dimension(ncid,"mxsat",mxsat)
! read npvap ("nsat")
    CALL eznc_get_0Dint(ncid,"nsat",npvap)
    write(6,*) 'found ', npvap, ' pvap values in ', trim(filename)

    if (npvap == 0) return
    allocate(satid(numsp),satdat(npvap,2))
    !allocate(satname(npvap))
    allocate(Tb(npvap), dB(npvap))

! read name !! CAUSES PROG TO CRASH !!
!    CALL eznc_get_1Dchar(ncid,"namsat",maxlsp,mxsat,satname,1,npvap)
!    PRINT*,"CAUTION: READING NAMNAN NOT NAMSAT: TEMPORARY HARDWIRE"
!    CALL eznc_get_1Dchar(ncid,"namnan",maxlsp,mxsat,satname,1,npvap)

! read name, Tb, dB
    !CALL eznc_get_2Dreal(ncid,"nandat",mxsat,2,satdat(1:npvap,1:2),1,npvap,1,2)
    CALL eznc_get_2Dreal(ncid,"nandat",satdat(1:npvap,1:2),1,npvap,1,1,2,1)
      Tb=satdat(:,1)
      dB=satdat(:,2)

! find index for species...
!... or you could read it from the NetCDF file ...
    status = NF90_INQ_VARID(ncid,"satid",varid)
    IF (status==NF90_NOERR) THEN
      !CALL eznc_get_1Dint(ncid,"satid",maxsp,satid(1:numsp),1,numsp)
      CALL eznc_get_1Dint(ncid,"satid",satid(1:numsp),1,numsp,1)
      do ispec =  numsp,1,-1 ! fill in the "AER" values too for these old files
        i = satid(ispec)
         if(i .ne. 0 .and. satid(ispec+1).EQ.0) THEN
           satid(ispec+1)=satid(ispec)
        ENDIF
      enddo
    ELSE
      !CALL eznc_get_1Dint(ncid,"nanid",maxsp,satid(1:numsp),1,numsp)
      CALL eznc_get_1Dint(ncid,"nanid",satid(1:numsp),1,numsp,1)
      do ispec =  numsp,1,-1 ! fill in the "AER" values too for these old files
        i = satid(ispec)
        if (i .ne. 0 ) satid(ispec+1)=satid(ispec)
      enddo
    ENDIF
    
    do ispec = 1, numsp
      i = satid(ispec)
      if (i .ne. 0 ) then
        species(ispec)%pvap%Tb = Tb(i)
        species(ispec)%pvap%dB = dB(i)
        species(ispec)%pvap%pvap_atm_298 = &
          exp((4.1012 + dB(i)) * ((298/Tb(i)-1)/(298/Tb(i) - 0.125))*log(10.))
          !species(ispec)%pvap%pvap_cstar_298 = &
          !  (1e6*seed_molw*seed_mass*species(ispec)%pvap%pvap_atm_298)/(8.2e-5*298)
! Cstar corrected per CU work
      call get_number(species(ispec)%formula, species(ispec)%molw,&
                      species(ispec)%nc, species(ispec)%nh, species(ispec)%nn, &
                      species(ispec)%no, species(ispec)%nr, species(ispec)%ns, &
                      species(ispec)%nfl, species(ispec)%nbr, species(ispec)%ncl)
          species(ispec)%pvap%pvap_cstar_298 = &
            (1e6*species(ispec)%molw*species(ispec)%pvap%pvap_atm_298)/ &
            (8.2e-5*298)

    ! I'm too lazy to solve the clausius-clapeyron equation so here's a quick and dirty way
    ! Siyuan Wang (siyuan@ucar.edu)
    pvap_atm_280_temp = exp((4.1012 + dB(i)) * ((280./Tb(i)-1.)/(280./Tb(i) - 0.125))*log(10.))
    species(ispec)%pvap%delta_H_vap_J_per_mol = & 
                        log(pvap_atm_280_temp/species(ispec)%pvap%pvap_atm_298) * &
                        8.314 / (1./298. - 1./280.)
      endif
    enddo

    DEALLOCATE (satdat,satid,Tb,dB)
   
  END SUBROUTINE read_ncdf_pvap
!-----------------------------------------------------------
  SUBROUTINE read_ncdf_henry(ncid)
!! subroutine needs to know values of: {maxlsp, mxdep}
!! subroutine needs to access routines in ncutil.f
!! subroutine DOES NOT need to "use netcdf"
    integer, dimension(:), allocatable :: depid,iddep
    real, dimension(:), allocatable    :: cf, Keff, rf
    real, dimension(:,:), allocatable  :: depdat
    integer                            :: nhenry, i, ispec
    integer                            :: ncid, status, varid
    integer                            :: mxdep
    
    species%henry%cf = -9999.
    species%henry%Keff = -9999.
    species%henry%rf   = -9999.
    
! read max # of henry ("mxdep")
    CALL eznc_get_dimension(ncid,"mxdep",mxdep)
! read nhenry ("ndep")
    CALL eznc_get_0Dint(ncid,"ndep",nhenry)

    if (nhenry == 0) return
    allocate(cf(nhenry), Keff(nhenry), rf(nhenry) )
    allocate(depid(numsp),iddep(nhenry),depdat(nhenry,3))

!read names
!    CALL eznc_get_1Dchar(ncid,"depnam",maxlsp,mxdep,name,1,nhenry)
!read depdat(cf, Keff, rf) and distribute into constituent arrays
    !CALL eznc_get_2Dreal(ncid,"depdat",mxdep,3,depdat(1:nhenry,1:3),1,nhenry,1,3)
    CALL eznc_get_2Dreal(ncid,"depdat",depdat(1:nhenry,1:3),1,nhenry,1,1,3,1)
    cf=depdat(:,1)
    Keff=depdat(:,2)
    rf=depdat(:,3)

! find index for species
!... or you could read it from the NetCDF file ...
! NB: older files have iddep, newer ones have depid
    status = NF90_INQ_VARID(ncid,"depid",varid)
    IF (status==NF90_NOERR) THEN
      !CALL eznc_get_1Dint(ncid,"depid",maxsp,depid(1:numsp),1,numsp)
      CALL eznc_get_1Dint(ncid,"depid",depid(1:numsp),1,numsp,1)
      do ispec = 1, numsp
        i = depid(ispec)
        if (i == 0) cycle
        species(ispec)%henry%cf = cf(i)
        species(ispec)%henry%Keff = Keff(i)
        species(ispec)%henry%rf = rf(i)
      enddo
    ELSE ! older files only contain iddep
      !CALL eznc_get_1Dint(ncid,"iddep",mxdep,iddep(1:nhenry),1,nhenry)
      CALL eznc_get_1Dint(ncid,"iddep",iddep(1:nhenry),1,nhenry,1)
      do i = 1, nhenry
        ispec = iddep(i)
        species(ispec)%henry%cf = cf(i)
        species(ispec)%henry%Keff = Keff(i)
        species(ispec)%henry%rf = rf(i)
        species(ispec+1)%henry%cf = cf(i)
        species(ispec+1)%henry%Keff = Keff(i)
        species(ispec+1)%henry%rf = rf(i)
      enddo
    ENDIF

    DEALLOCATE (depdat,depid,iddep, cf, Keff, rf)

  END SUBROUTINE read_ncdf_henry
  
!-------------------------------------------------------------------
  SUBROUTINE read_ncdf_kreac(ncid)
!! subroutine needs to know values of: {maxlsp, mxkOH/NO3/O3 etc}
!! subroutine needs to access routines in ncutil.f
!! subroutine DOES NOT need to "use netcdf"
    integer, dimension(:), allocatable :: kreacid
    real, dimension(:), allocatable    :: A, n, Ea
    real, dimension(:,:), allocatable  :: kreacdat
    integer                            :: nkreac, i,k, ispec
    integer                            :: ncid
    integer                            :: mxkreac
    character(LEN=7) :: name_numk, name_mxk, name_kdat, name_kid
    character(LEN=3) :: ksp
   
    species%koh%A = -9999.
    species%koh%n = -9999.
    species%koh%Ea   = -9999.
    species%koh%k_298 = -9999.
    species%kno3%A = -9999.
    species%kno3%n = -9999.
    species%kno3%Ea   = -9999.
    species%kno3%k_298 = -9999.
    species%ko3%A = -9999.
    species%ko3%n = -9999.
    species%ko3%Ea   = -9999.
    species%ko3%k_298 = -9999.
   
    allocate(kreacid(maxsp))

    do k = 1,3
      SELECT CASE (k)
      CASE (1) 
        ksp = "OH "
      CASE (2) 
        ksp = "NO3"
      CASE (3) 
        ksp = "O3 "
      END SELECT

      name_numk = "nk"//ksp(1:LEN_TRIM(ksp))
      name_mxk  = "mxk"//ksp(1:LEN_TRIM(ksp))
      name_kdat = "k"//ksp(1:LEN_TRIM(ksp))//"dat"
      name_kid  = "k"//ksp(1:LEN_TRIM(ksp))//"id"

      !print*,k,name_numk, name_mxk, name_kdat, name_kid

! read max & actual # of ksp species ("mxkreac" "nkreac")
      CALL eznc_get_dimension(ncid,name_mxk(1:LEN_TRIM(name_mxk)),mxkreac)
      CALL eznc_get_0Dint(ncid,name_numk(1:LEN_TRIM(name_numk)),nkreac)

      print*,mxkreac,nkreac

      if (nkreac == 0) cycle
      allocate(A(nkreac), n(nkreac), Ea(nkreac) )
      allocate(kreacdat(nkreac,3))

!read kreacdat(cf, Keff, rf) and distribute into constituent arrays
      !CALL eznc_get_2Dreal(ncid,name_kdat(1:LEN_TRIM(name_kdat)), &
      !                     mxkreac,3,kreacdat(1:nkreac,1:3),1,nkreac,1,3)
      CALL eznc_get_2Dreal(ncid,name_kdat(1:LEN_TRIM(name_kdat)), &
                           kreacdat(1:nkreac,1:3),1,nkreac,1,1,3,1)
      A=kreacdat(:,1)
      n=kreacdat(:,2)
      Ea=kreacdat(:,3)

! find index for species
!... or you could read it from the NetCDF file ...
      !CALL eznc_get_1Dint(ncid,name_kid(1:LEN_TRIM(name_kid)), &
      !                    maxsp,kreacid(1:numsp),1,numsp)
      CALL eznc_get_1Dint(ncid,name_kid(1:LEN_TRIM(name_kid)), &
                               kreacid(1:numsp),1,numsp,1)

      do ispec = 1, numsp
        i = kreacid(ispec)
        if (i == 0) cycle


        SELECT CASE (k)
        CASE (1)
          species(ispec)%koh%A = A(i)
          species(ispec)%koh%n = n(i)
          species(ispec)%koh%Ea = Ea(i)
          species(ispec)%koh%k_298 = A(i)*(298**n(i))*exp(-Ea(i)/298.)
          !PRINT*,ispec,i,species(ispec)%code,species(ispec)%koh%k_298
        CASE (2)
          species(ispec)%kno3%A = A(i)
          species(ispec)%kno3%n = n(i)
          species(ispec)%kno3%Ea = Ea(i)
          species(ispec)%kno3%k_298 = A(i)*(298**n(i))*exp(-Ea(i)/298.)
          !PRINT*,ispec,i,species(ispec)%code,species(ispec)%kno3%k_298
        CASE (3)
          species(ispec)%ko3%A = A(i)
          species(ispec)%ko3%n = n(i)
          species(ispec)%ko3%Ea = Ea(i)
          species(ispec)%ko3%k_298 = A(i)*(298**n(i))*exp(-Ea(i)/298.)
          !PRINT*,ispec,i,species(ispec)%code,species(ispec)%ko3%k_298
        CASE DEFAULT
          PRINT*,"read_ncdf_kreac: cannot assign species%...%k_298"
          
        END SELECT

      enddo ! ispec = 1,numsp

      DEALLOCATE (kreacdat, A, n, Ea)
    enddo ! k = 1,3
    
    DEALLOCATE (kreacid)

  END SUBROUTINE read_ncdf_kreac
!-------------------------------------------------------------------
  SUBROUTINE read_ncdf_reacrate(ncid)
!! subroutine needs to know values of: {maxlsp, mxkOH/NO3/O3 etc}
!! subroutine needs to access routines in ncutil.f
!! subroutine needs to know working dimensions ndat, numre
!! subroutine DOES NOT need to "use netcdf"
    integer :: ncid,mbox
    integer,parameter :: ibox = 1

    CALL eznc_get_dimension(ncid,"mbox",mbox)
    CALL eznc_get_dimension(ncid,"maxre",maxre)
    CALL eznc_get_dimension(ncid,"mxleft",mxleft)
    CALL eznc_get_dimension(ncid,"mxright",mxright)
    CALL eznc_get_0Dint(ncid,"numre",numre)

    ALLOCATE(numstoi(numre,2),        &
             idrestoi(numre,mxleft),  &
             idpdstoi(numre,mxright), &
             restoicf(numre,mxleft),  &
             pdstoicf(numre,mxright))

    CALL eznc_get_2Dint(ncid,"numstoi", &! maxre,  mxleft, &
                              numstoi,1,numre,1,1,mxleft,1)
    CALL eznc_get_2Dint(ncid,"idpdstoi", &!maxre,  mxright, &
                              idpdstoi,1,numre,1,1,mxright,1)
    CALL eznc_get_2Dint(ncid,"idrestoi", &! maxre,  mxleft, &
                              idrestoi,1,numre,1,1,mxleft,1)
    CALL eznc_get_2Dreal(ncid,"pdstoicf", &! maxre,  mxright, &
                               pdstoicf,1,numre,1,1,mxright,1)
    CALL eznc_get_2Dreal(ncid,"restoicf", &! maxre,  mxleft, &
                               restoicf,1,numre,1,1,mxleft,1)

  END SUBROUTINE read_ncdf_reacrate

!-------------------------------------------------------------------

  SUBROUTINE write_2D_array(array, header, filename)
    character*(filenames_length), intent(in)                    :: filename
    character*(filenames_length)                                :: full_filename
    real, dimension(:,:), intent(in)                            :: array
    character*(header_length), dimension(:),   intent(in)       :: header
    integer                                                     :: i, j

    if(size(header) /= size(array,2)) then
      write(6,*) " --error-- header size doesn't fit the data in write_2D_array "
      stop
    endif

! write Netcdf output
    IF(output_type.eq."netcdf") then
      CALL write_2D_nc_op(array, header, filename)
    ELSE

! write csv output
    IF(nbox.EQ.1)THEN
    full_filename = trim(filename)//".csv"
    ELSE
    full_filename = trim(boxname)//"_"//trim(filename)//".csv"
    ENDIF
    open(write_unit, file =full_filename)

    write(write_unit,*) (trim(header(i)), ", ", i =lbound(header,1), ubound(header, 1)-1), header(ubound(header, 1))

    do i=lbound(array,1), ubound(array,1)
      write(write_unit,*) (array(i,j), ", ", j = lbound(array,2), ubound(array,2)-1), array(i, ubound(array, 2))
    enddo

    close(write_unit)

    ENDIF
  END SUBROUTINE
!-------------------------------------------------------------------
  SUBROUTINE write_2D_nc_op(array, header, filename)
    character*(filenames_length), intent(in)                    :: filename
    character*(filenames_length)                                :: full_filename
    real, dimension(:,:), intent(in)                            :: array
    character*(header_length), dimension(:),   intent(in)       :: header
    integer                                                     :: i
    integer                                       :: ncid,dim1,dim2,nchars

! find array sizes
    if(size(header) /= size(array,2)) then
      write(6,*) " --error-- header size doesn't fit the data in write_2D_array "
      stop
    endif

    dim1 = size(array,1)
    dim2 = size(array,2)
    nchars = 0
    DO i = 1,dim2
      nchars = MAX(nchars,LEN_TRIM(header(i)))
    ENDDO

! define output file
    IF(nbox.EQ.1)THEN
    full_filename = trim(filename)//".nc"
    ELSE
    full_filename = trim(boxname)//"_"//trim(filename)//".nc"
    ENDIF

! open file
    CALL open_ncfile_new(full_filename,ncid)

! define dimensions
    CALL eznc_def_dim(ncid,"dim1",dim1)
    CALL eznc_def_dim(ncid,"dim2",dim2)
    CALL eznc_def_dim(ncid,"nchars",nchars)

! define variables
    CALL eznc_def_1Dchar(ncid,"header","nchars","dim2")
    CALL eznc_def_2Dreal(ncid,"array","dim1","dim2")

! switch to data mode for data writing
    CALL switch_ncfile_to_data_mode(ncid)

! write ouput variables
    DO i = 1,dim2
      CALL eznc_put_1Dchar(ncid,"header",header(i),nchars,i,i)
    ENDDO
    CALL eznc_put_2Dreal(ncid,"array",array,1,dim1,1,dim2)

! close file
    CALL close_ncfile(ncid)

  END SUBROUTINE
END MODULE
