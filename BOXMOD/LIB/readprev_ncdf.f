* subroutine to read existing model output and extract all results    *
* Interpolates to input time & returns values for single time point   *
! NetCDF version, Julia Lee-Taylor, NCAR, Jan 2018.                   !
***********************************************************************
      SUBROUTINE readprev_ncdf(ncidp,ncido,
     &                    numsp,nsat,ndim,
     &                    iotime,oconc,oconcaer,ibox)

      USE flags_module
      USE akparameter_module
      USE time_mgmt_module, ONLY: idat
      USE forcing_params_module, ONLY: nbox
      IMPLICIT NONE

! input
      INTEGER:: ncidp,ncido !previous and output ids
      REAL   :: iotime
      INTEGER:: numsp,nsat,ndim,ibox

! output
      REAL,DIMENSION(maxsp):: oconc
      REAL,DIMENSION(mxsat):: oconcaer

* maximum number of data in the initialisation file
      INTEGER,PARAMETER :: mdat=80000

! arrays to read in from initialisation file
      INTEGER:: numspp,numsat,numdim,ndat,nboxp
      REAL,DIMENSION(mdat):: ptime
      !REAL,DIMENSION(mdat,mxsat)::  pconc,pconcaer
      REAL,DIMENSION(maxsp,2):: pconc
      REAL,DIMENSION(mxsat,2):: pconcaer
      CHARACTER(len=40):: line1,line2

! other internals
      INTEGER:: i,j
      REAL   :: incr

* -----------------------------
* initialize
* -----------------------------
      ptime = 0.
      pconc = 0.
      pconcaer = 0.
      oconc = 0.
      oconcaer = 0.

* -------------------------------------------------
! (file is already open)
* check mechanisms match and 
* check mech dimensions of the previous result file
* -------------------------------------------------
      print*,ncido,ncidp
      CALL eznc_get_globalatt(ncido,"mech_name",line1)
      CALL eznc_get_globalatt(ncidp,"mech_name",line2)
      IF (line1.ne.line2) THEN
        WRITE(6,*) '--error--, mechanism names do not match'
        STOP
      ENDIF

      CALL eznc_get_globalatt(ncido,"mech_date",line1)
      CALL eznc_get_globalatt(ncidp,"mech_date",line2)
      IF (line1.ne.line2) THEN
        WRITE(6,*) '--error--, mechanism dates do not match'
        STOP
      ENDIF

      CALL eznc_get_0Dint(ncidp,"numsp",numspp)
      IF (numspp.ne.numsp) THEN
        WRITE(6,*) '--error--, number of species differs between files'
        STOP
      ENDIF

      CALL eznc_get_0Dint(ncidp,"nsat",numsat)
      IF (numsat.ne.nsat.AND.soa_fg.NE.0) THEN
        WRITE(6,*) '--error--, number of pvap spp differs between files'
        STOP
      ENDIF

      CALL eznc_get_0Dint(ncidp,"ndim",numdim)
      IF (numdim.ne.ndim) THEN
        WRITE(6,*) '--error--, number of dimers differs between files'
        STOP
      ENDIF

      CALL eznc_get_dimension(ncidp,"nbox",nboxp)
      IF (nboxp .ne. nbox) then
        WRITE(6,*) '--error--, number of boxes differs from prev run'
        STOP
      ENDIF

* ----------------------------
* read data in the previous result file
* ----------------------------

! find number of times in previous file
      CALL eznc_get_dimension(ncidp,"ntout",ndat)
      IF (ndat.gt.mdat) THEN
        WRITE(6,*) '--error--, number of data exceed mdat'
        STOP
      ENDIF

! read time data
      CALL eznc_get_1Dreal(ncidp,"time",ndat,ptime(1:ndat),1,ndat)
      IF (ptime(ndat).LT.iotime) THEN
        WRITE(6,*)
     &  '--error--, prev data does not encompass desired time'
        WRITE(6,*) ptime(ndat),iotime
        STOP
      ENDIF

      DO idat = 1, ndat
        if(ptime(idat) .lt. iotime) cycle
        exit
      ENDDO
! idat is the timestep we're looking for
! get the corresponding concentration
      print*,'nboxp=', nboxp, "ibox-", ibox
      print*,'ndat=', ndat, 'idat=', idat
      call eznc_get_3Dreal(ncidp,"conc",maxsp,nboxp,ndat,
     &                     pconc(1:maxsp,1:2),
     &                     1,maxsp,
     &                     ibox,ibox,
     &                     idat-1,idat)
      if (soa_fg .eq. 1) then
        call eznc_get_3Dreal(ncidp,"caer",maxsp,nboxp,ndat,
     &                     pconcaer(1:mxsat,1:2),
     &                     1,mxsat,
     &                     ibox,ibox,
     &                     idat-1,idat)
      endif

* interpolate values for desired time
      incr = (iotime-ptime(idat-1))/(ptime(idat)-ptime(idat-1))

      oconc(1:numsp) = pconc(1:numsp,1) +
     &           incr*(pconc(1:numsp,2) -
     &                 pconc(1:numsp,1) )

      if (soa_fg .eq. 1) then
        oconcaer(1:nsat) = pconcaer(1:nsat,1) +
     &             incr*(pconcaer(1:nsat,2) -
     &                   pconcaer(1:nsat,1) )
      endif
      RETURN
      END
***********************************************************************
