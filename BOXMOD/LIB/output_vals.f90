      SUBROUTINE output_vals(cbox,ctotaerbox,maerbox)

      USE flags_module,ONLY: iofmt_fg,wall_fg,soa_fg,dimer_fg,dyn_fg
      USE io_units_module,ONLY: lppa,lppf,lpaa,lpff,lsoa
      USE time_mgmt_module
      USE printphoto_module
      USE NetCDF_vars_module
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE fundamental_consts_module,ONLY: multimass

      IMPLICIT NONE

! -------------------------------------
! INPUTS
      REAL    :: cbox(maxsp)  ! = conc(1:maxsp,ibox), retained
      REAL    :: ctotaerbox   ! = ctotaer(ibox), retained
      REAL    :: maerbox      ! = maer(ibox),retained

! LOCAL
      INTEGER :: lboxa      ! current box aerosol output unit : lppa or lpaa
      INTEGER :: lboxf      ! current box gas output unit : lppf or lpff

!      REAL::  tmpconc(maxsp)        ! only concs >= maskval
!      REAL::  tmpcro2(maxro2)       ! only concs >= maskval
      REAL::  tmpcaer(mxsat)        ! only concs >= maskval
      
!------------------------------------------------------------------
! set up output units to use for current box
      SELECT CASE (ibox)
        CASE (1)
          lboxa = lppa
          lboxf = lppf
        CASE (2)
          lboxa = lpaa
          lboxf = lpff
        CASE DEFAULT
          lboxa = lppa
          lboxf = lppa
      END SELECT

!--------------------------------------------------------
! Output gas phase values: NetCDF, binary
!--------------------------------------------------------
! NB: time,temp are calculated as scalars but written into 1D space
!     cbot etc are calculated as 1D arrays but written into 2D space

! NOTE : during development we encountered the following error:
! error in NetCDF routine: NF90_PUT_VAR cbox                  -60
! i.e. /* Math result not representable */
! => range of data is too great for data type.
! Potential solution: if cbox < maskvalue, set tmpconc = 0. for o/p
            !tmpconc(:)=cbox(:)
            !WHERE(cbox.LT.maskval) tmpconc = 0.
            !tmpcro2(:)=cro2(:)
            !WHERE(cro2.LT.maskval) tmpcro2 = 0.

! NetCDF output
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        !PRINT*,"starting gas write"
        CALL wrtopgas_ncdf(ncid_out,ibox,itout,cbox)
      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

! binary output
      IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

        WRITE(lboxf) time,temp,rh,(cbox(isp),isp=1,numsp)

      ENDIF ! (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

!--------------------------------------------------------
! write aerosol conditions 
!--------------------------------------------------------
      IF(soa_fg.GT.0)THEN
        !PRINT*,"starting soa write"
        !CALL getpvap_nan(nsat,temp,tb,dB,psat)
! NetCDF output
        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

          tmpcaer(:)=caer(:)
          WHERE(caer.LT.maskval) tmpcaer=0.
          CALL wrtopaer_ncdf(ncid_out,ibox,itout,nsat, &
                             maerbox,psat,caer,ctotaerbox)

        ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

! binary output
        IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN

! equilibrium approach
          IF (soa_fg.EQ.1) THEN
            WRITE(lboxa)time,temp,(cbox(idsat(isp)),isp=1,nsat), &
                 (caer(isp),isp=1,nsat),(psat(isp),isp=1,nsat), &
                maerbox
          ENDIF ! (soa_fg.EQ.1)

! mass transfer approach
          IF (soa_fg.EQ.2) THEN

            IF (wall_fg .eq. 1 .and. dimer_fg.EQ.1) then
              WRITE(lboxa)time,temp,(cbox(idgsat(isp)),isp=1,nsat), &
                         (cbox(idasat(isp)),isp=1,nsat), &
                         (cbox(iddsat(isp)),isp=1,ndim), &
                         (cbox(idwsat(isp)),isp=1,nsat), &
                         (psat(isp),isp=1,nsat), &
                         maerbox
            ELSE IF (dimer_fg.EQ.1) THEN
              WRITE(lboxa)time,temp,(cbox(idgsat(isp)),isp=1,nsat), &
                         (cbox(idasat(isp)),isp=1,nsat), &
                         (cbox(iddsat(isp)),isp=1,ndim), &
                         (0.0,isp=1,nsat), &
                         (psat(isp),isp=1,nsat), &
                         maer
            ELSE IF (wall_fg .eq. 1) then
              WRITE(lboxa)time,temp,(cbox(idgsat(isp)),isp=1,nsat), &
                     (cbox(idasat(isp)),isp=1,nsat), &
                     (cbox(idwsat(isp)),isp=1,nsat), &
                     (psat(isp),isp=1,nsat), &
                     maer
            ELSE 
              WRITE(lboxa)time,temp,(cbox(idgsat(isp)),isp=1,nsat), &
                     (cbox(idasat(isp)),isp=1,nsat), &
                     (0.0,isp=1,nsat), &
                     (psat(isp),isp=1,nsat), &
                     maer
            ENDIF ! (wall, dimer combinations)
          ENDIF ! (soa_fg.EQ.2)

        ENDIF ! (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2)

! ascii mass transfer output (for checking)

        IF (dyn_fg.EQ.1)THEN
          IF( time.GT.0 .AND. ctotaerbox.GT.small) THEN
            WRITE(lsoa,*)time,ctotaerbox,Rp, &
                 (maer/multimass+cnv*Mp)/ctotaerbox
          ELSE
            WRITE(lsoa,*)time,ctotaerbox,Rp,0.
          ENDIF
        ENDIF

      ENDIF ! (soa_fg.GT.0)
      !PRINT*,"after soa write"
! -------------------------------------
      END SUBROUTINE output_vals
! -------------------------------------

