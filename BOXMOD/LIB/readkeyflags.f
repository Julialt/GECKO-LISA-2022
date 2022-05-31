      SUBROUTINE readkeyflags
!----------------------------------------------------
! PURPOSE: read indat.key for simulation flags ONLY:
!          returns before reading input data 
!----------------------------------------------------

      USE flags_module,ONLY: print_steadystate_fg,iofmt_fg ,prevflag, 
     &                       lagflag,depos_fg ,mixing_fg ,emis_fg, 
     &                       soa_fg,pvap_fg ,dimer_fg ,ro2_fg, 
     &                       OFR_fg,dyn_fg ,wall_fg ,icham, jall_fg,
     &                       seedtyp_fg,isopsoa_fg, tracer_fg, 
     &                       printphoto_fg,seedtyp_fg,icham,iscape,
     &                       vbs_fg, vbs_aging_fg, reacrate_fg
      USE io_units_module,ONLY: lout
      USE forcing_params_module,ONLY: nbox
      USE akparameter_module

      IMPLICIT NONE

* local
      LOGICAL lokerr
      CHARACTER(76) line
      CHARACTER(4)  keyword
      INTEGER  i, j, ibox, isp, lin, i_val, ierr

      INTEGER,PARAMETER :: nchar=76

****************************
* initialise               *
****************************
! LOGICAL
      lokerr=.false.

****************************
*  reading input           *
****************************
      lin=12
      WRITE(6,*) "read indat.key for flags ......"
      OPEN(lin, file='indat.key', form='formatted',  status='old')

      WRITE(lout,*)
      WRITE(lout,*)'    Keyword input:'
      WRITE(lout,*)

* read next input line
90    continue
      READ(lin,'(a4,(a))',end=99)keyword,line
      !WRITE(6,'(a4,(a))')keyword,line

* check if the line is a comment
      IF(keyword(1:1).eq.'.'.or.keyword(1:1).eq.'/'.or.
     &   keyword(1:1).eq.'!')GOTO 90

* check if the line is a number
      IF(keyword(1:1).eq.'0'.or.
     &   keyword(1:1).eq.'1'.or.
     &   keyword(1:1).eq.'2'.or.
     &   keyword(1:1).eq.'3'.or.
     &   keyword(1:1).eq.'4'.or.
     &   keyword(1:1).eq.'5'.or.
     &   keyword(1:1).eq.'6'.or.
     &   keyword(1:1).eq.'7'.or.
     &   keyword(1:1).eq.'8'.or.
     &   keyword(1:1).eq.'9'.or.
     &   keyword(1:1).eq.' ')GOTO 90

* check for end of flags section - skip to end & return
      IF(keyword.eq.'DATA')THEN
        WRITE(lout,*)"END OF FLAGS SECTION OF INDAT.KEY"
        GOTO 100

* check for end of file - skip to end & return
      ELSE IF(keyword(1:3).eq.'END')THEN
        WRITE(lout,*)"END OF INDAT.KEY ENCOUNTERED"
        GOTO 100

* check end of keyword section
!      ELSE IF(keyword(1:4).eq.'DATA')THEN
!        WRITE(lout,*)"END OF FLAG SECTION ENCOUNTERED"
!        GOTO 90

! FLAGS FOR SIMULATION OPTIONS !
* fresh or continuing simulation?
      ELSE IF(keyword.eq.'PREV')THEN
        prevflag=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* Binary/Netcdf i/o?
      ELSE IF(keyword.eq.'IFMT')THEN
        iofmt_fg=i_val(line,1,nchar,1,ierr)

        PRINT*,"iofmt_fg =",iofmt_fg

        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* # of boxes to solve?
      ELSE IF(keyword.eq.'NBOX')THEN
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

* Eulerian/Lagrangan?
      ELSE IF(keyword.eq.'LAGR')THEN
        lagflag=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* Deposition?
      ELSE IF(keyword.eq.'DEPO')THEN
        depos_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* calcul des equilibre thermo par scape
      ELSE IF(keyword.eq.'SCAP')THEN
        iscape=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
        ENDIF
        GOTO 90

* SOA module? Mass transfer?
      ELSE IF(keyword.eq.'SOAF')THEN
        soa_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* vapor pressure scheme?
      ELSE IF(keyword.eq.'PVAP')THEN
        pvap_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* print final timestep concentrations?
      ELSE IF(keyword.eq.'PSSF')THEN
        print_steadystate_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* output instantaneous reaction rates?
      ELSE IF(keyword.eq.'RRAT')THEN
        reacrate_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* dimerisation?
      ELSE IF(keyword.eq.'DIMR')THEN
        dimer_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* Dynamic gas -> aero transfer?
      ELSE IF(keyword.eq.'DYNF')THEN
        dyn_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* wall partitioning?
      ELSE IF(keyword.eq.'WALL')THEN
        wall_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* chamber identity?
      ELSE IF(keyword.eq.'CHAM')THEN
        icham=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* identity of non-volatile seed aerosol?
      ELSE IF(keyword.eq.'NVID')THEN
        seedtyp_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90
        
! isoprene soa formation ?
      ELSE IF(keyword.eq.'ISOA')THEN
        isopsoa_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90
! vbs calculation
      ELSE IF(keyword.eq.'VBSF')THEN
        vbs_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90
      ELSE IF (keyword.eq.'VBAG')THEN
        vbs_aging_fg =i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* OFR mode?
      ELSE IF(keyword(1:3).eq.'OFR')THEN
        OFR_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* output j-values for all reference cheomophores?
      ELSE IF(keyword.eq.'CROM')THEN
        jall_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* RO2+RO2 allowed?
      ELSE IF(keyword.eq.'RO2F')THEN
        ro2_fg=i_val(line,1,nchar,1,ierr)
        IF(ierr.NE.0)THEN
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading ',keyword
          STOP
        ENDIF
        GOTO 90

* keyword unknown
      ELSE
        WRITE(lout,*)' --error--  while reading ',keyword
        WRITE(lout,*)'            keyword unknown'
        lokerr=.true.
          STOP
        GOTO 90
      ENDIF

****************************
* done reading input       *
****************************

      CLOSE(lin)
      RETURN

* stop IF error
100   IF (lokerr) THEN
        WRITE(lout,*)"100,lokerr"
        STOP
      ENDIF
      CLOSE(lin)
      RETURN

* keyword END not found
99    WRITE(lout,*)' --error--  keyword END not found'
      lokerr=.true.
          STOP

      END SUBROUTINE readkeyflags

!==============================================================

