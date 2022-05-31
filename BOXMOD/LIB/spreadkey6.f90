      SUBROUTINE readkey6

      USE flags_module,ONLY: isopsoa_fg,emis_fg,noem_fg,mixing_fg, &
                             isopsoa_fac,iscape, lightsonoff_fg, &
                             inject_fg
      USE akparameter_module
      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: nskip,ntprint,ntstep,tstart,tstop
      USE constraints_module,ONLY : constrained_thing
      USE solver_params_module,ONLY: atol,rtol,dtmin
      USE forcing_params_module
      USE module_data_gecko_main
      USE module_chamber_tools,ONLY: set_light, set_injection, &
                            init_lights, init_injections, A_V_ratio 

! local
      LOGICAL lokerr
      CHARACTER(4)  keyword
      CHARACTER(76) line
      INTEGER  i_val, ierr
      INTEGER  nrtol, natol, ndtmin
      INTEGER  ncons,ntcon
      INTEGER  nemi,ntem
      INTEGER,PARAMETER :: nchar=76
      INTEGER,PARAMETER :: lin=12
      REAL     tlen
      REAL     r_val

! only for chamber setup
      integer :: ninject, nlights
      real    :: inject_tstart,inject_tstop,inject_conc
      real    :: light_ton, light_toff
      character(maxlsp) :: inject_code

      PRINT*,"in spreadkey6"
!***************************
! initialise               *
!***************************
! LOGICAL
      lokerr=.false.
      
      natol = 0
      nrtol = 0
      ndtmin = 0

!***************************
!  reading input           *
!***************************

      PRINT*,"open indat.key"
      OPEN(lin, file='indat.key', form='formatted',  status='old')

      WRITE(lout,*)
      WRITE(lout,*)'    Keyword input:'
      WRITE(lout,*)

! read next input line
90    CONTINUE
      READ(lin,'(a4,(a))',end=99)keyword,line

! DEBUG
      !PRINT*,line

! check for end of file - skip to end and return
      IF(keyword(1:3).eq.'END')THEN
        WRITE(lout,*)"END OF INDAT.KEY ENCOUNTERED"
        GOTO 100
      ENDIF

! check if the line is a comment - skip to next read
      IF(keyword(1:1).eq.'.'.or.keyword(1:1).eq.'/'.or. &
         keyword(1:1).eq.'!') GOTO 90

!-----------------------------------------------------------
! FLAGS FOR SIMULATION OPTIONS ARE READ BY readflags.f !
! list of flags: if any found, skip to next read
!-----------------------------------------------------------
      IF(keyword(1:4).EQ.'CHAM'.OR. &
         keyword(1:4).EQ.'CROM'.OR. &
         keyword(1:4).EQ.'DEPO'.OR. &
         keyword(1:4).EQ.'DIMR'.OR. &
         keyword(1:4).EQ.'DYNF'.OR. &
         keyword(1:4).EQ.'IFMT'.OR. &
         keyword(1:4).EQ.'LAGR'.OR. &
         keyword(1:3).EQ.'OFR'.OR. &
         keyword(1:4).EQ.'NVID'.OR. &
         keyword(1:4).EQ.'PREV'.OR. &
         keyword(1:4).EQ.'PSSF'.OR. &
         keyword(1:4).EQ.'PVAP'.OR. &
         keyword(1:4).EQ.'RRAT'.OR. &
         keyword(1:4).EQ.'RO2F'.OR. &
         keyword(1:4).EQ.'SOAF'.OR. &
         keyword(1:4).EQ.'ISOA'.OR. &
         keyword(1:4).EQ.'VBSF'.OR. &
         keyword(1:4).EQ.'VBAG'.OR. &
         keyword(1:4).EQ.'WALL') THEN
!        PRINT*,"FLAG: ",keyword
        GOTO 90
      ENDIF

! check for end of flags section
      IF(keyword(1:4).EQ.'DATA')THEN
        WRITE(lout,*)"END OF FLAGS SECTION ENCOUNTERED"
        PRINT*,"END OF FLAGS SECTION ENCOUNTERED"
        GOTO 90
      ENDIF

!==================================================================
! data inputs:
!---------------

! fix number of boxes (now also in readkeyflags.f)
      IF(keyword.eq.'NBOX')THEN
        nbox=i_val(line,1,nchar,1,ierr)
        PRINT*,"nbox ",nbox
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        IF(nbox.GT.2)THEN
          WRITE(lout,*)' --error-- NBOX > 2'
          WRITE(lout,*)' change NBOX or reconfigure model'
          STOP
        ENDIF
        GOTO 90

! temperature (sine function)
      ELSE IF(keyword.eq.'TEMP')THEN
        IF (nbox.eq.0) THEN
          WRITE(lout,*) '--error--, TEMP was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        DO ibox=1,nbox
           READ (lin,'(a76)') line
           WRITE (lout,'(a76)') line
           tempm(ibox)=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading first parameter'
             WRITE(lout,*)'            of TEMP for box :',ibox
           ENDIF
           tempa(ibox)=r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading second parameter'
             WRITE(lout,*)'            of TEMP for box :',ibox
           ENDIF
           temptm(ibox)=r_val(line,1,nchar,3,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading third parameter'
             WRITE(lout,*)'            of TEMP for box :',ibox
           ENDIF
        ENDDO
        GOTO 90

! temperature (tabular input)
      ELSE IF(keyword.eq.'TKTB')THEN
        IF (nbox.eq.0) THEN
          WRITE(lout,*) '--error--, TKTB was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        ntk=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword TKTB'
        ENDIF
        IF(ntk.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword TKTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,ntk
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         tktim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of TKTB for line :',i
         ENDIF
         DO ibox=1,nbox
           tkval(ibox,i)=r_val(line,1,nchar,ibox+1,ierr)
         ENDDO
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of TKTB for line :',i
         ENDIF
       ENDDO
       GOTO 90

! pressure (tabular input)
      ELSE IF(keyword.eq.'PRTB')THEN
        IF (nbox.eq.0) THEN
          WRITE(lout,*) '--error--, PRTB was READ before'
          WRITE(lout,*) '           setting the number of BOX'
        ENDIF
        npr=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword PRTB'
        ENDIF
        IF(npr.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword PRTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,npr
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         prtim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of PRTB for line :',i
         ENDIF
         DO ibox=1,nbox
           prval(ibox,i)=r_val(line,1,nchar,ibox+1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of PRTB for line :',i
           ENDIF
         ENDDO
       ENDDO
       GOTO 90

! "vertical" dilution rate (constant)
      ELSE IF(keyword.eq.'DILF') THEN
        dilfix = 1
        dilconst = r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90
! dilution rate (tabular input)
      ELSE IF(keyword.eq.'DLTB')THEN
        IF (nbox.eq.0) THEN
          WRITE(lout,*) '--error--, DLTB was READ before'
          WRITE(lout,*) '           setting the number of BOX'
        ENDIF
        ndil=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword DLTB'
        ENDIF
        IF(ndil.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword DLTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,ndil
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         diltim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of DLTB for line :',i
         ENDIF
         dilval(i)=r_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of DLTB for line :',i
         ENDIF
       ENDDO
       GOTO 90

! calcul des equilibre thermo par scape
      ELSE IF(keyword.eq.'SCAP')THEN
        iscape=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        print*,"spreadkey: iscape =",iscape
        GOTO 90

! set atmospheric number density
      ELSE IF(keyword.eq.'NDEN')THEN
        IF (nbox.eq.0) THEN
          lokerr=.true.
          WRITE(lout,*) '--error--,',keyword,' was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        sumcfix = 1
        DO ibox=1,nbox
          sumc(ibox)=r_val(line,1,nchar,ibox,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            for box :',ibox
          ENDIF
        ENDDO
        GOTO 90

! set volume mixing ratio of water
      ELSE IF(keyword.eq.'WATR')THEN
        IF (nbox.eq.0) THEN
          lokerr=.true.
          WRITE(lout,*) '--error--,',keyword,' was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        waterfix = 1
        DO ibox=1,nbox
          water(ibox)=r_val(line,1,nchar,ibox,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            for box :',ibox
          ENDIF
        ENDDO
        GOTO 90
! set value of pressure if fixed
      ELSE IF(keyword.eq.'PRES') THEN
        presfix = 1
        DO ibox=1,nbox
          prconst(ibox) = r_val(line,1,nchar,ibox,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            for box :',ibox
          ENDIF
        ENDDO
        GOTO 90

! set value of jO2
      ELSE IF(keyword.eq.'JO2 ')THEN
        jo2=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! set value of jH2O
      ELSE IF(keyword.eq.'JH2O')THEN
        jh2o=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! set value of 185nm photon flux
      ELSE IF(keyword.eq.'F185')THEN
        f185=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! set value of 254nm photon flux
      ELSE IF(keyword.eq.'F254')THEN
        f254=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! regarde quel type de saison (two seasons allowed : 1=winter
! 2=summer)
      ELSE IF(keyword.eq.'SEAS')THEN
        iseas=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! humidite relative (sine function)
      ELSE IF(keyword.eq.'RHUM')THEN
        IF (nbox.eq.0) THEN
          WRITE(lout,*) '--error--, RHUM was read before'
          WRITE(lout,*) '           setting NBOX'
          stop
        ENDIF
        DO ibox=1,nbox
           READ (lin,'(a76)') line
           WRITE (lout,'(a76)') line
           rhm(ibox)=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading first parameter'
             WRITE(lout,*)'            of RHUM for box :',ibox
           ENDIF
           rha(ibox)=r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading second parameter'
             WRITE(lout,*)'            of RHUM for box :',ibox
           ENDIF
           rhtm(ibox)=r_val(line,1,nchar,3,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading third parameter'
             WRITE(lout,*)'            of RHUM for box :',ibox
           ENDIF
        ENDDO
        GOTO 90

! relative humidity (tabular input)
      ELSE IF(keyword.eq.'RHTB')THEN
        IF (nbox.eq.0) THEN
          lokerr=.true.
          WRITE(lout,*) '--error--, RHTB was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        nrh=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword RHTB'
        ENDIF
        IF(nrh.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword RHTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,nrh
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         rhtim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of RHTB for line :',i
         ENDIF
         DO ibox=1,nbox
           rhval(ibox,i)=r_val(line,1,nchar,ibox+1,ierr)
         ENDDO
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of RHTB for line :',i
         ENDIF
       ENDDO
       GOTO 90

! vitesse du vent
      ELSE IF(keyword.eq.'WIND')THEN
           windm=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading first parameter'
             WRITE(lout,*)'            of WIND for box :',ibox
           ENDIF
           winda=r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading second parameter'
             WRITE(lout,*)'            of WIND for box :',ibox
           ENDIF
           windtm=r_val(line,1,nchar,3,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading third parameter'
             WRITE(lout,*)'            of WIND for box :',ibox
           ENDIF
        GOTO 90

! wind speed (tabular input)
      ELSE IF(keyword.eq.'WSTB')THEN
        nws=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword WSTB'
        ENDIF
        IF(nws.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword WSTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,nws
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         wstim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of WSTB for line :',i
         ENDIF
         wsval(i)=r_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of WSTB for line :',i
         ENDIF
       ENDDO
       GOTO 90
       
       ELSE IF (keyword .eq. 'IFAC') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining isopsoa_fac when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif     
         isopsoa_fac=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading ',keyword
         ENDIF           
       GOTO 90
       
       ELSE IF (keyword .eq. 'NAER') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining naer when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif
         if(index(line, 'TOP ') .ne. 0) then
           ibox = 2
         else if (index(line, 'BOT ') .ne. 0) then
           ibox = 1
         else
           write(lout,*) ' --error--'
           write(lout,*) ' must specify which box the NAER'
           write(lout,*) ' keyword is applied to'
           write(lout,*) ' e.g. NAER BOT 8 or NAER TOP 6'
           stop
         endif       
       
         naer_const(ibox)%npoints = i_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword NAER'
         ENDIF
         IF(naer_const(ibox)%npoints.GT.maxinput)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword NAER'
           WRITE(lout,*)'            maximum size is', maxinput
         ENDIF
         DO i = 1, naer_const(ibox)%npoints
           READ (lin,'(a76)') line
           naer_const(ibox)%time(i) = r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of NAER for line :',i
           ENDIF
           naer_const(ibox)%values(i) = r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of NAER for line :',i
           ENDIF
         ENDDO
         goto 90

       
       ELSE IF (keyword .eq. 'PHAR') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining ph when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif
         if(index(line, 'TOP ') .ne. 0) then
           ibox = 2
         else if (index(line, 'BOT ') .ne. 0) then
           ibox = 1
         else
           write(lout,*) ' --error--'
           write(lout,*) ' must specify which box the PHAR'
           write(lout,*) ' keyword is applied to'
           write(lout,*) ' e.g. PHAR BOT 8 or PHAR TOP 6'
           stop
         endif       
       
         ph_const(ibox)%npoints = i_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword PHAR'
         ENDIF
         IF(ph_const(ibox)%npoints.GT.maxinput)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword PHAR'
           WRITE(lout,*)'            maximum size is', maxinput
         ENDIF
         DO i = 1, ph_const(ibox)%npoints
           READ (lin,'(a76)') line
           ph_const(ibox)%time(i) = r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of PHAR for line :',i
           ENDIF
           ph_const(ibox)%values(i) = r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of PHAR for line :',i
           ENDIF
         ENDDO
         goto 90
         
       ELSE IF (keyword .eq. 'SULF') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining sulfate when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif
         if(index(line, 'TOP ') .ne. 0) then
           ibox = 2
         else if (index(line, 'BOT ') .ne. 0) then
           ibox = 1
         else
           write(lout,*) ' --error--'
           write(lout,*) ' must specify which box the SULF'
           write(lout,*) ' keyword is applied to'
           write(lout,*) ' e.g. SULF BOT 8 or SULF TOP 6'
           stop
         endif       
       
         sulfate_const(ibox)%npoints = i_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword SULF'
         ENDIF
         IF(sulfate_const(ibox)%npoints.GT.maxinput)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword SULF'
           WRITE(lout,*)'            maximum size is', maxinput
         ENDIF
         DO i = 1, sulfate_const(ibox)%npoints
           READ (lin,'(a76)') line
           sulfate_const(ibox)%time(i) = r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of SULF for line :',i
           ENDIF
           sulfate_const(ibox)%values(i) = r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of SULF for line :',i
           ENDIF
         ENDDO
         goto 90
         
       ELSE IF (keyword .eq. 'NITR') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining nitrates when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif
         if(index(line, 'TOP ') .ne. 0) then
           ibox = 2
         else if (index(line, 'BOT ') .ne. 0) then
           ibox = 1
         else
           write(lout,*) ' --error--'
           write(lout,*) ' must specify which box the NITR'
           write(lout,*) ' keyword is applied to'
           write(lout,*) ' e.g. NITR BOT 8 or NITR TOP 6'
           stop
         endif       
       
         nitrate_const(ibox)%npoints = i_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword NITR'
         ENDIF
         IF(nitrate_const(ibox)%npoints.GT.maxinput)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword NITR'
           WRITE(lout,*)'            maximum size is', maxinput
         ENDIF
         DO i = 1, nitrate_const(ibox)%npoints
           READ (lin,'(a76)') line
           nitrate_const(ibox)%time(i) = r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of NITR for line :',i
           ENDIF
           nitrate_const(ibox)%values(i) = r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of NITR for line :',i
           ENDIF
         ENDDO
         goto 90
       
       ELSE IF (keyword .eq. 'KAPA') then
         if (isopsoa_fg .eq. 0) then
           write(lout,*) ' --warning--'
           write(lout,*) 'constraining kappa when not '
           write(lout,*) 'making isoprene soa '
           write(lout,*) 'isopsoa_fg = 0'
         endif
         if(index(line, 'TOP ') .ne. 0) then
           ibox = 2
         else if (index(line, 'BOT ') .ne. 0) then
           ibox = 1
         else
           write(lout,*) ' --error--'
           write(lout,*) ' must specify which box the KAPA'
           write(lout,*) ' keyword is applied to'
           write(lout,*) ' e.g. KAPA BOT 8 or KAPA TOP 6'
           stop
         endif       
       
         kappa_const(ibox)%npoints = i_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword KAPA'
         ENDIF
         IF(kappa_const(ibox)%npoints.GT.maxinput)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword KAPA'
           WRITE(lout,*)'            maximum size is', maxinput
         ENDIF
         DO i = 1, kappa_const(ibox)%npoints
           READ (lin,'(a76)') line
           kappa_const(ibox)%time(i) = r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of KAPA for line :',i
           ENDIF
           kappa_const(ibox)%values(i) = r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading value parameter'
             WRITE(lout,*)'            of KAPA for line :',i
           ENDIF
         ENDDO
         goto 90

! calcul des equilibre thermo par scape
! definie la proportion des differentes surfaces
!      ELSE IF(keyword.eq.'URB1')THEN
!        psurf(1)=r_val(line,1,nchar,1,ierr)
!        IF(ierr.NE.0)THEN
!          lokerr=.true.
!          WRITE(lout,*)' --error--  while reading ',keyword
!        ENDIF
!        GOTO 90
!      ELSE IF(keyword.eq.'CULT')THEN
!        psurf(2)=r_val(line,1,nchar,1,ierr)
!        IF(ierr.NE.0)THEN
!          lokerr=.true.
!          WRITE(lout,*)' --error--  while reading ',keyword
!        ENDIF
!        GOTO 90
!      ELSE IF(keyword.eq.'FLEA')THEN
!        psurf(3)=r_val(line,1,nchar,1,ierr)
!        IF(ierr.NE.0)THEN
!          lokerr=.true.
!          WRITE(lout,*)' --error--  while reading ',keyword
!        ENDIF
!        GOTO 90
!      ELSE IF(keyword.eq.'FCON')THEN
!        psurf(4)=r_val(line,1,nchar,1,ierr)
!        IF(ierr.NE.0)THEN
!          lokerr=.true.
!          WRITE(lout,*)' --error--  while reading ',keyword
!        ENDIF
!       GOTO 90
       ELSE IF(keyword.eq.'SURF')THEN
         nsd=i_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword SURF'
         ENDIF
         IF(nsd.GT.mhd)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword SURF'
           WRITE(lout,*)'            maximum size is', mhd
         ENDIF
         DO i=1,nsd
           READ (lin,'(a76)') line
           WRITE (lout,'(a76)') line
           surft(i)=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of SURF for line :',i
             WRITE(lout,*)'            parameter 1'
           ENDIF
           psurf(i,1)=r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading height parameter'
             WRITE(lout,*)'            of SURF for line :',i
             WRITE(lout,*)'            parameter 2'
           ENDIF
           psurf(i,2)=r_val(line,1,nchar,3,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading height parameter'
             WRITE(lout,*)'            of SURF for line :',i
             WRITE(lout,*)'            parameter 3'
           ENDIF
           psurf(i,3)=r_val(line,1,nchar,4,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading height parameter'
             WRITE(lout,*)'            of SURF for line :',i
             WRITE(lout,*)'            parameter 4'
           ENDIF
           psurf(i,4)=r_val(line,1,nchar,5,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading height parameter'
             WRITE(lout,*)'            of SURF for line :',i
             WRITE(lout,*)'            parameter 5'
           ENDIF
         ENDDO
         GOTO 90

! emission depending on surface type
      ELSEIF (keyword .eq. 'SEMI') THEN
        write(6,*) line
        isurf = i_val(line,1,nchar,1,ierr)
        write(6,*) isurf
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword,'(1)'
        ENDIF

        if (isurf .gt. msur) then
          WRITE(lout,*)' --error-- isurf > msurf: update akparameter'
        ENDIF

        nemi = i_val(line,1,nchar,2,ierr)
        IF(nemi.GT.0) emis_fg=1
        WRITE(6,*)nemi," EMITTED SPECIES"
        WRITE(lout,*)nemi," EMITTED SPECIES"
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword,'(1)'
        ENDIF
        IF(nemi.GT.maxem)THEN
          lokerr=.true.
          WRITE(lout,*)' --error-- nemi > maxem: update akparameter'
        ENDIF
        surface_emi(isurf)%nemis = nemi

! read individual species headers
        DO i = 1,nemi
93        READ (lin,'(a76)') line

! check if the line is a comment
          IF(line(1:1).eq.'.'.or.line(1:1).eq.'/'.or. &
             line(1:1).eq.'!') GOTO 93

          WRITE(lout,*)line
          WRITE(6,*)line
! call akspnum to check if  emitted species is known
          CALL akspnum(line,chrsp,numsp,isp)

          IF(isp.eq.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            species is unknown in line :'
            GOTO 90
          ENDIF

          surface_emi(isurf)%emission(i)%activefg = .TRUE.
          surface_emi(isurf)%emission(i)%name =         &
                 trim(line(1:index(line, ' ')-1))
          surface_emi(isurf)%emission(i)%index = isp
          ! molec by default
          surface_emi(isurf)%emission(i)%unit = 'molec cm-2'

          ntem = i_val(line,1,nchar,2,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword,'(2)'
          ENDIF

          DO j = 1,ntem
            write(6,*) j
            READ (lin,'(a76)') line
            surface_emi(isurf)%emission(i)%table(j,1) =  &
                              r_val(line,1,nchar,1,ierr)
            surface_emi(isurf)%emission(i)%table(j,2) =  &
                              r_val(line,1,nchar,2,ierr)
            ! catch input errors
            IF(j.GT.1)THEN
              IF(surface_emi(isurf)%emission(i)%table(j,1).LE. &
                 surface_emi(isurf)%emission(i)%table(j-1,1)) ierr = 1
            ENDIF

            IF(ierr.NE.0)THEN
              lokerr=.true.
              WRITE(lout,*)' --error--  while reading ',keyword,'(3)'
            ENDIF
          ENDDO
          surface_emi(isurf)%emission(i)%npoints = ntem
        ENDDO
        GOTO 90

! initial concentration
      ELSE IF(keyword.eq.'REAC')THEN
        IF (nbox.eq.0) THEN
          lokerr = .true.
          WRITE(lout,*) '--error--, REAC was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF

        CALL akspnum(line,chrsp,numsp,isp)
        IF(isp.eq.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            species is unknown in line :'
          WRITE(lout,*) line
          GOTO 90
        ENDIF
        conc(isp,1)=r_val(line,1,nchar,2,ierr)
        IF (nbox.eq.1) THEN
          cbg(isp)=r_val(line,1,nchar,3,ierr)
        ENDIF
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*) line
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'         first value was not read correctly'
        ENDIF
        IF (nbox.EQ.2) THEN
          conc(isp,2)=r_val(line,1,nchar,3,ierr)
          cbg(isp)=r_val(line,1,nchar,4,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'       second value was not read correctly'
          ENDIF
        ENDIF
        GOTO 90

! initial concentration pour les especes des equilibres thermo
      ELSE IF(keyword.eq.'C0EQ')THEN
        IF (iscape.le.0) THEN
            lokerr=.true.
          WRITE(lout,*) '--error--, C0EQ requires SCAP = 1'
        ENDIF
        IF (nbox.eq.0) THEN
            lokerr=.true.
          WRITE(lout,*) '--error--, C0EQ was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        IF (iscape.eq.1) THEN
          call akspnum(line,chrsp,numsp,isp)
          IF(isp.eq.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            species is unknown in line :'
            WRITE(lout,*) line
            GOTO 90
          ENDIF
          conc(isp,1)=r_val(line,1,nchar,2,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'         first value was not read correctly'
          ENDIF
          IF (nbox.eq.2) THEN
            conc(isp,2)=r_val(line,1,nchar,3,ierr)
            IF(ierr.NE.0)THEN
              lokerr=.true.
              WRITE(lout,*)' --error--  while reading ',keyword
              WRITE(lout,*)'       second value was not read correctly'
            ENDIF
          ENDIF
        ENDIF
        GOTO 90

! chamber injections
      ELSE IF(keyword .eq. 'INJC') THEN
          ninject = i_val(line, 1, nchar,1, ierr)
        call init_injections(ninject)
        do i = 1,ninject
          read(lin,'(a76)') line
          inject_code = line(1:index(line,' ')-1)
          inject_tstart = r_val(line, 1, nchar, 2, ierr)
          inject_tstop = r_val(line, 1, nchar, 3, ierr)
          inject_conc = r_val(line, 1, nchar, 4, ierr)
          call set_injection(inject_code, inject_tstart, inject_tstop, &
                             inject_conc)
        enddo
        goto 90
! lights constraints (chamber)
      ELSE IF (keyword .eq. 'LGHT') THEN
        nlights = i_val(line, 1, nchar, 1, ierr)
        call init_lights(nlights)
        do i = 1, nlights
          read(lin, '(a76)') line
          light_ton = r_val(line, 1, nchar, 1, ierr)
          light_toff = r_val(line, 1, nchar, 2, ierr)
          call set_light(light_ton, light_toff)
        enddo
        goto 90
! (A/V) ratio for chamber experiment
! useful for wall loss parameterization
      ELSE IF (keyword .eq. 'AVRT') THEN
        A_V_ratio = r_val(line, 1, nchar, 1, ierr)
        if (ierr.ne.0) then
          lokerr = .true.
          write(lout, *) ' --error-- while reading ', keyword
        endif
        goto 90

! constrained concentrations
      ELSE IF(keyword.eq.'CONS')THEN

        if(index(line, 'TOP ') .ne. 0) then
          ibox = 2
        else if (index(line, 'BOT ') .ne. 0) then
          ibox = 1
        else
          write(lout,*) ' --error--'
          write(lout,*) ' must specify which box the CONS'
          write(lout,*) ' keyword is applied to'
          write(lout,*) ' e.g. CONS BOT 8 or CONS TOP 6'
          stop
        endif
        ncons = i_val(line,1,nchar,2,ierr)
        WRITE(6,*) ncons," CONSTRAINED CONCENTRATIONS"
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword,'(1)'
        ENDIF
        IF(ncons.GT.maxem)THEN
          lokerr=.true.
          WRITE(lout,*)' --error-- ncons > maxem: update akparameter'
        ENDIF

! read individual species headers
        DO i = 1,ncons
91        READ (lin,'(a76)') line

! check if the line is a comment
          IF(line(1:1).eq.'.'.or.line(1:1).eq.'/'.or. &
             line(1:1).eq.'!') GOTO 91

          WRITE(lout,*) line
          WRITE(6,*)line
! call akspnum to check if constrained species is known
          IF (line(1:4) .ne. 'NOx ') then
            CALL akspnum(line,chrsp,numsp,isp)
          else
            isp = -1
            noxfix(ibox) = 1
          endif
          IF(isp.eq.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            species is unknown in line :'
            GOTO 90
          ENDIF
          cons_spec(ibox, i)%activefg = .TRUE.
          cons_spec(ibox, i)%name = trim(line(1:index(line, ' ')-1))
          cons_spec(ibox, i)%index = isp
          cons_spec(ibox, i)%unit = 'molec cm-3'
          ntcon = i_val(line,1,nchar,2,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword,'(2)'
          ENDIF
          DO j = 1,ntcon
            READ (lin,'(a76)') line
            cons_spec(ibox, i)%table(j,1) = r_val(line,1,nchar,1,ierr)
            cons_spec(ibox, i)%table(j,2) = r_val(line,1,nchar,2,ierr)
            ! catch input errors
            IF(j.GT.1)THEN
              IF(cons_spec(ibox, i)%table(j,1).LE. &
                 cons_spec(ibox, i)%table(j-1,1)) ierr = 1
            ENDIF
          !WRITE(6,*)i,j,cons_spec(i)%table(j,1),cons_spec(i)%table(j,2)
            IF(ierr.NE.0)THEN
              lokerr=.true.
              WRITE(lout,*)' --error--  while reading ',keyword,'(3)'
            ENDIF
          ENDDO
          cons_spec(ibox, i)%npoints = ntcon
        ENDDO

        GOTO 90

! emissions
      ELSE IF(keyword.eq.'EMIS')THEN
        nemi = i_val(line,1,nchar,1,ierr)
        IF(nemi.GT.0) emis_fg=1
        WRITE(6,*)nemi," EMITTED SPECIES"
        WRITE(lout,*)nemi," EMITTED SPECIES"
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword,'(1)'
        ENDIF
        IF(nemi.GT.maxem)THEN
          lokerr=.true.
          WRITE(lout,*)' --error-- nemi > maxem: update akparameter'
        ENDIF

! read individual species headers
        DO i = 1,nemi
92        READ (lin,'(a76)') line

! check if the line is a comment
          IF(line(1:1).eq.'.'.or.line(1:1).eq.'/'.or. &
             line(1:1).eq.'!') GOTO 92

          WRITE(lout,*)line
          WRITE(6,*)line
! call akspnum to check if  emitted species is known
          IF (line(1:4) .ne. 'NOx ') then
            CALL akspnum(line,chrsp,numsp,isp)
          else
            isp = -1
          endif
          IF(isp.eq.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword
            WRITE(lout,*)'            species is unknown in line :'
            GOTO 90
          ENDIF
          emi_spec(i)%activefg = .TRUE.
          emi_spec(i)%name = trim(line(1:index(line, ' ')-1))
          emi_spec(i)%index = isp
          ! molec by default
          emi_spec(i)%unit = 'molec cm-2 s-1'
          ntem = i_val(line,1,nchar,2,ierr)
          IF(ierr.NE.0)THEN
            lokerr=.true.
            WRITE(lout,*)' --error--  while reading ',keyword,'(2)'
          ENDIF
          DO j = 1,ntem
            READ (lin,'(a76)') line
            emi_spec(i)%table(j,1) = r_val(line,1,nchar,1,ierr)
            emi_spec(i)%table(j,2) = r_val(line,1,nchar,2,ierr)
            !WRITE(6,*)i,j,emi_spec(i)%table(j,1),emi_spec(i)%table(j,2)
            IF(ierr.NE.0)THEN
              lokerr=.true.
              WRITE(lout,*)' --error--  while reading ',keyword,'(3)'
            ENDIF
          ENDDO
          emi_spec(i)%npoints = ntem
        ENDDO
        GOTO 90
      ELSE IF(keyword.eq.'NOEM')THEN
        noem_fg=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90
      ELSE IF(keyword.eq.'ISOP') THEN
        isop_fac = r_val(line, 1, nchar, 1, ierr)
        IF(ierr.ne.0) THEN
          lokerr=.true.
          WRITE(lout,*)' --error-- while reading ', keyword
        ENDIF
        GOTO 90
      ELSE IF(keyword.eq.'MTER') THEN
        mterp_fac = r_val(line, 1, nchar, 1, ierr)
        IF(ierr.ne.0) THEN
          lokerr=.true.
          WRITE(lout,*)' --error-- while reading ', keyword
        ENDIF
        GOTO 90
! starting time
      ELSE IF(keyword.eq.'TSTR')THEN
        tstart=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! stopping time
      ELSE IF(keyword.eq.'TSTP')THEN
        tstop=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! ntstep = number of model computational timesteps
! EITHER set ntstep using number (NPAS)
      ELSE IF(keyword.eq.'NPAS')THEN
        ntstep=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! OR set ntstep using timestep length (TLEN)
      ELSE IF(keyword.eq.'TLEN')THEN
        tlen = r_val(line,1,nchar,1,ierr)
        ntstep=INT((tstop-tstart)/tlen)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! SKIP input = # of timesteps SKIPPED between output points 
! nskip as used by code = # of timesteps between consecutive output points INCLUDING the output point
! EITHER: set nskip using number (SKIP)
      ELSE IF(keyword.eq.'SKIP')THEN
        nskip = i_val(line,1,nchar,1,ierr)
        nskip = nskip + 1
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! OR: set nskip using desired # of output values (NPRT)
      ELSE IF(keyword.eq.'NPRT')THEN
        ntprint=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! relative tolerance
      ELSE IF(keyword.eq.'RTOL')THEN
        nrtol=nrtol+1
        IF (nrtol.GT.1) THEN
          lokerr=.true.
          write(6,*) '--error--, >1 x RTOL keyword '
        ENDIF
        rtol=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! absolute tolerance
      ELSE IF(keyword.eq.'ATOL')THEN
        natol=natol+1
        IF (natol.GT.1) THEN
          lokerr=.true.
          write(6,*) '--error--, >1 x ATOL keyword '
        ENDIF
        atol=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! minimum solver timestep input value
      ELSE IF(keyword.eq.'DTMN')THEN
        ndtmin=ndtmin+1
        IF (ndtmin.GT.1) THEN
          lokerr=.true.
          write(6,*) '--error--, >1 x DTMN keyword '
        ENDIF
        dtmin=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! aerosol surface (cm2/cm3)
      ELSE IF(keyword.eq.'AERO')THEN
        saero=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! aerosol bulk activity coefficient, gamma
      ELSE IF(keyword.eq.'GAMM')THEN
        gamm=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
       GOTO 90

! nonvolatile aerosol molecular mass (g/mole)
      ELSE IF(keyword.eq.'NVMW')THEN
        Mp=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
       GOTO 90

! nonvolatile aerosol particle radius (cm)
      ELSE IF(keyword.eq.'NVRO')THEN
        Rpo=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
       GOTO 90

! aerosol particle radius (tabular input)
      ELSE IF(keyword.eq.'RPTB')THEN
        IF (nbox.eq.0) THEN
          lokerr=.true.
          WRITE(lout,*) '--error--, RHTB was read before'
          WRITE(lout,*) '           setting NBOX'
        ENDIF
        rpfix_fg = 1
        nrp=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword RPTB'
        ENDIF
        IF(nrp.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword RPTB'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,nrp
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         rptim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of NRPO for line :',i
         ENDIF
         DO ibox=1,nbox
           rpval(ibox,i)=r_val(line,1,nchar,ibox+1,ierr)
         ENDDO
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of NRPO for line :',i
         ENDIF
       ENDDO
       GOTO 90

! nonvolatile cseed aerosol concentration (molec/cm3)
      ELSE IF(keyword.eq.'SEED')THEN
        cnv=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
       GOTO 90

! time-varying nonvolatile cseed aerosol concentration (molec/cm3)
      ELSE IF(keyword.eq.'SEET')THEN
       nseed=i_val(line,1,nchar,1,ierr)
       IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword SEET'
       ENDIF
       IF(nseed.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword SEET'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,nseed
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         tseed(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of SEET for line :',i
         ENDIF
         cseed(i)=r_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading height parameter'
           WRITE(lout,*)'            of SEET for line :',i
         ENDIF
       ENDDO
       GOTO 90

! photolysis (localisation of the box on earth)
      ELSE IF(keyword.eq.'PPHO')THEN
        !PRINT*, line
        sla=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        slo=r_val(line,1,nchar,2,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        tz=r_val(line,1,nchar,3,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        iy=i_val(line,1,nchar,4,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        im=i_val(line,1,nchar,5,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        id=i_val(line,1,nchar,6,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          WRITE(lout,*)'            value was not read correctly'
        ENDIF
        GOTO 90

! photolysis adjustment factor (to approximate measurements)
      ELSE IF(keyword.eq.'JFAC')THEN
       !PRINT*, line
       njf=i_val(line,1,nchar,1,ierr)
       IF(ierr.NE.0)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword JFAC'
       ENDIF
       IF(njf.GT.mtim)THEN
         WRITE(lout,*)' --error--  while reading number of data'
         WRITE(lout,*)'            linked to keyword JFAC'
         WRITE(lout,*)'            maximum size is', mtim
       ENDIF
       DO i=1,njf
         READ (lin,'(a76)') line
         WRITE (lout,'(a76)') line
         jftim(i)=r_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading time parameter'
           WRITE(lout,*)'            of JFAC for line :',i
         ENDIF
         jfval(i)=r_val(line,1,nchar,2,ierr)
         IF(ierr.NE.0)THEN
           lokerr=.true.
           WRITE(lout,*)' --error--  while reading value parameter'
           WRITE(lout,*)'            of JFAC for line :',i
         ENDIF
       ENDDO
       GOTO 90

! flag and value for fixed solar zenith angle (SZA)
      ELSE IF(keyword.eq.'SZAF')THEN
       szafix=1
       szaval=r_val(line,1,nchar,1,ierr)
       IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
       ENDIF
       GOTO 90

!  mixing height (in cm) data
       ELSE IF(keyword.eq.'HBOX')THEN
         nhd=i_val(line,1,nchar,1,ierr)
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword HBOX'
         ENDIF
         IF(nhd.GT.mhd)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword HBOX'
           WRITE(lout,*)'            maximum size is', mhd
         ENDIF
         DO i=1,nhd
           READ (lin,'(a76)') line
           WRITE (lout,'(a76)') line
           boxt(i)=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of HBOX for line :',i
           ENDIF
           boxh(i)=r_val(line,1,nchar,2,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading height parameter'
             WRITE(lout,*)'            of HBOX for line :',i
           ENDIF
         ENDDO
         GOTO 90

! height of the top of the boxes
      ELSE IF(keyword.eq.'HTOP')THEN
        htop=r_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! tropospheric subsidence velocity (cm.s-1), exchanges between top box & background.
      ELSE IF (keyword.eq.'SUBS')THEN
        vs = r_val(line, 1, nchar, 1, ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

! cross-boundary transport velocity (cm.s-1)
       ELSE IF(keyword.eq.'MIX')THEN
         nmx=i_val(line,1,nchar,1,ierr)
         IF(nmx.GT.0) mixing_fg=1
         IF(ierr.NE.0)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword MIX'
         ENDIF
         IF(nmx.GT.mhd)THEN
           WRITE(lout,*)' --error--  while reading number of data'
           WRITE(lout,*)'            linked to keyword MIX'
           WRITE(lout,*)'            maximum size is', mhd
         ENDIF
         DO i=1,nmx
           READ (lin,'(a76)') line
           WRITE (lout,'(a76)') line
           mixt(i)=r_val(line,1,nchar,1,ierr)
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading time parameter'
             WRITE(lout,*)'            of MIX for line :',i
           ENDIF
           DO ibox=1,nbox
             mixv(ibox,i)=r_val(line,1,nchar,ibox+1,ierr)
           ENDDO
           IF(ierr.NE.0)THEN
             lokerr=.true.
             WRITE(lout,*)' --error--  while reading rate parameter'
             WRITE(lout,*)'            of MIX for line :',i
           ENDIF
         ENDDO
         GOTO 90

! keyword unknown
      ELSE
        WRITE(lout,*)' --error--  while reading ',keyword
        WRITE(lout,*)'            keyword unknown'
        lokerr=.true.
        GOTO 90

      ENDIF ! (checking the keyword)

!***************************
! done reading input       *
!***************************

! keyword END not found
99    WRITE(lout,*)' --error--  keyword END not found'
      lokerr=.true.

! keyworkd END is found
100   CONTINUE

!****************************
! check for necessary input *
!****************************

! check time stop
      IF(tstop.le.0.0)THEN
        WRITE(lout,*)' --error--  TSTP not properly specified'
        lokerr=.true.
      ENDIF

! check # of boxes
      IF(nbox.eq.0)THEN
        WRITE(lout,*)' --error--  NBOX not specified'
        lokerr=.true.
      ENDIF

! check input for temperature
      IF(ntk.EQ.0)THEN
      DO i=1,nbox
        IF(tempm(i).le.0.0)THEN
          WRITE(lout,*)' --error--  TEMP not properly specified'
          WRITE(lout,*)'             in box:',i
          lokerr=.true.
        ENDIF
      ENDDO
      ENDIF

! check input for relative humidity
      IF(nrh.EQ.0)THEN
      DO i=1,nbox
        IF (rhm(i)-rha(i).lt.0.0) THEN
          WRITE(lout,*)' --error--  rhum not properly specified'
          WRITE(lout,*)'            (cannot be less than 0)'
          WRITE(lout,*)'            in box:',ibox
          lokerr=.true.
        ENDIF
        IF (rhm(i)+rha(i).gt.100.0) THEN
          WRITE(lout,*)' --error--  rhum not properly specified'
          WRITE(lout,*)'            (cannot be greater than 100)'
          WRITE(lout,*)'            in box:',ibox
          lokerr=.true.
        ENDIF
      ENDDO
      ENDIF

! check that dtmin is not larger than tlen
      dtmin = MIN(dtmin,tlen)

! check that mixing height do not exceed HTOP
      IF (nbox.eq.2) THEN
        DO i=1,nhd
          IF (boxh(i).GE.htop) THEN
            WRITE(lout,*)' --error--  mixing height exceed HTOP'
            lokerr=.true.
          ENDIF
        ENDDO
      ENDIF

! SEAS and SURF MUST be specified together
      IF((iseas.le.0.AND.nsd.gt.0).OR.(iseas.gt.0.AND.nsd.le.0))THEN
        WRITE(lout,*)' --error--  SEAS and SURF BOTH required'
        lokerr=.true.
      ENDIF

! check that time for mixing height variation is set for length of run
! (note that this will not discriminate between daily repeating values
!  and values set for the length of the run)
      IF (boxt(1).ne.0. .and. boxt(1) .gt. tstart) THEN
         WRITE(lout,*)' --error--  mixing height does not start at'
         WRITE(lout,*)'            at time=0'
         lokerr=.true.
      ENDIF
      IF (boxt(nhd).lt.tstop.AND.boxt(nhd).ne.86400.) THEN
         WRITE(lout,*) '--error-- '// &
          'mixing height time must be >= tstop or = 86400'
         lokerr=.true.
      ENDIF

! stop IF error
      IF (lokerr) STOP

      CLOSE(lin)
! reset ibox to default
      ibox = 0

      END SUBROUTINE readkey6
