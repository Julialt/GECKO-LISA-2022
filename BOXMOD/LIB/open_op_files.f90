      SUBROUTINE open_op_files

! open output files (if not already open)

      USE flags_module,ONLY: iofmt_fg,soa_fg,jall_fg,printphoto_fg, &
                             iscape
      USE io_units_module,ONLY: lppf,lppa,lpff,lpaa,ljval, &
                                ljall,lsoa,lpbl,lro2
      USE printphoto_module,ONLY: n_printphoto,photoreac_names
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main

      IMPLICIT NONE

      LOGICAL :: exist
!---------------------------------------------------
! open ascii file for selected output jvals
      IF (printphoto_fg.AND.n_printphoto .GT. 0) THEN
      INQUIRE(file='outdat.jvals',exist=exist)
      IF(exist)THEN
        OPEN(ljval,file='outdat.jvals',form='formatted',action='write',&
                      status='old',position='append')
      ELSE
        OPEN(ljval,file='outdat.jvals',form='formatted',action='write',&
                      status='new')
        WRITE(ljval,*) "time(s), ", &
                       (TRIM(photoreac_names(k)),',', &
                       k=1,n_printphoto-1), &
                       TRIM(photoreac_names(n_printphoto))
      ENDIF
      ENDIF !(printphoto_fg.AND.n_printphoto .GT. 0) THEN

!---------------------------------------------------
! open binary files for box 1 and box 2 results

      IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN
        PRINT*,"opening binary output files"

        OPEN(lppf,file='outdat.ppf',form='unformatted',status='unknown')
        IF (soa_fg.EQ.1) THEN
          OPEN(lppa,file='outdat.ppa',form='unformatted', &
             status='unknown')
        ENDIF

        IF (nbox == 2) THEN
          OPEN(lpff,file='outdat.pff',form='unformatted', &
               status='unknown')
          IF (soa_fg.GE.1) THEN
            OPEN(lpaa,file='outdat.paa',form='unformatted', &
                 status='unknown')
          ENDIF
        ENDIF

! write binary output headers
        CALL write_op_hdrs

! open ascii file to output ALL jvals
        IF (jall_fg.NE.0) THEN
        INQUIRE(file='outdat.jvals',exist=exist)
        IF(exist)THEN
        OPEN(ljall,file='outdat.jvals',form='formatted',action='write',&
                      status='old',position='append')
        ELSE
        OPEN(ljall,file='outdat.jvals',form='formatted',action='write',&
                      status='new')
        ENDIF !(exist)
        ENDIF !(jall_fg)

      ENDIF !(iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2)

!---------------------------------------------------
! open NetCDF file for results output

      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        CALL setup_ncdf_op
      ENDIF !(iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) 

!---------------------------------------------------
! open files for results linked to thermodynamical calculation

        IF (iscape.eq.1) THEN
          OPEN(20,file='outdat.distri',status='unknown')
          OPEN(67,file='outdat.aerodep',status='unknown')
          OPEN(71,file='outdat.1rhdat',status='unknown')
          OPEN(72,file='outdat.1res',status='unknown')
          OPEN(73,file='outdat.p1dat',status='unknown')
          OPEN(74,file='outdat.1brhdat',status='unknown')
          OPEN(75,file='outdat.p1bdat',status='unknown')
          IF (nbox == 2) THEN
            OPEN(81,file='outdat.2rhdat',status='unknown')
            OPEN(82,file='outdat.2res',status='unknown')
            OPEN(83,file='outdat.p2dat',status='unknown')
            OPEN(84,file='outdat.2brhdat',status='unknown')
            OPEN(85,file='outdat.p2bdat',status='unknown')
          ENDIF
        ENDIF

!---------------------------------------------------
! open some other ascii files, write ascii hdrs

! PBL height
      INQUIRE(file='outdat.pbl',exist=exist)
      IF(exist)THEN
        OPEN(lpbl,file='outdat.pbl',form='formatted',action='write', &
                      status='old',position='append')
      ELSE
        OPEN(lpbl,file='outdat.pbl',form='formatted',action='write', &
                      status='new')
        WRITE(lpbl,*) "Time [s], PBL Height [cm]"
      ENDIF

! RO2 concentration
      INQUIRE(file='outdat.ro2',exist=exist)
      IF(exist)THEN
        OPEN(lro2,file='outdat.ro2',form='formatted',action='write', &
                      status='old',position='append')
      ELSE
        OPEN(lro2,file='outdat.ro2',form='formatted',action='write', &
                      status='new')
        WRITE(lro2,*)"t1(s), t2(s), PERO1, PERO2, PERO3, PERO4,"//&
                     "PERO5, PERO6, PERO7, PERO8, PERO9 (molec/cc @ t1)"
      ENDIF

! SOA diagnostics
      IF(soa_fg.NE.0)THEN
      INQUIRE(file='outdat.soa',exist=exist)
      IF(exist)THEN
        OPEN(lsoa,file='outdat.soa',form='formatted',action='write', &
                      status='old',position='append')
      ELSE
        OPEN(lsoa,file='outdat.soa',form='formatted',action='write', &
                      status='new')
        WRITE(lsoa,*) "Time [s], ctotaer, [molec/cc], Rp [cm],"// &
                     " mwaer(mean) [g/mole] "
      ENDIF
      ENDIF

!---------------------------------------------------
      END SUBROUTINE open_op_files
!---------------------------------------------------
