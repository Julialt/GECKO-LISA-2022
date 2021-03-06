      SUBROUTINE get_envinp

      USE flags_module,ONLY: iofmt_fg, depos_fg, OFR_fg, iscape
      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: tstart,tstop,ntout,nskip, &
                                 ntstep,ntprint,delt
      USE akparameter_module,ONLY: maxsp
      USE NetCDF_vars_module,ONLY: ncid_in,max_outvar_size
      USE solver_params_module,ONLY: dtmin,dtmax,numit
      USE module_data_gecko_main

      IMPLICIT NONE

      REAL deltprint
      REAL outvar_size
!----------------------------------------------------------------------
! Read the input file *.key (ascii format)
!----------------------------------------------------------------------
      WRITE(6,*) 'reading the input file ...'

      CALL readkey6

      WRITE(6,*) 'done input (key) file read... nbox =',nbox

! itout = 1 corresponds to initial time (t=0.)
! => output index ntout = (ntstep/nskip) + 1
! OR = ntprint (if >0)
! NB: nskip defined within code as input variable "SKIP" + 1 (see spreadkey6)
! IE: model-defined nskip INCLUDES the output timestep at the END of the "skip"

! Timestep parameters
      delt=(tstop-tstart)/REAL(ntstep)
      dtmax=delt
      PRINT*,"delt",delt

      IF(ntprint.EQ.0)THEN
! nskip is used in the model as "skipped times PLUS ONE"
! but nskip is input as "# of skipped times" which is more intuitive
! The input definition is set in spreadkey6.f90
        nskip = nskip 
      ELSE
        IF(ntprint.EQ.1)ntprint = ntprint + 1
        deltprint=(tstop-tstart)/REAL(ntprint-1)
        nskip = INT(deltprint/delt) 
        PRINT*,"deltprint",deltprint
      ENDIF

      ntout = INT(ntstep/nskip)+1
      PRINT*,"nskip,ntprint,ntout",nskip,ntprint,ntout

      WRITE(lout,'(a15,e9.2)') 'dtmin         =',dtmin
      WRITE(lout,'(a15,i3)')   'numit         =',numit

!-------------------------------------------------------
! THIS SECTION PREVENTS PRE-DEFINED CONC ARRAY FROM BEING TOO LARGE
! SKIP THIS SECTION TO USE UNLIMITED ARRAY DIMENSION (RECOMMENDED)

! NetCDF limits fixed-size output arrays to < 2^31 - 4 large.
!      outvar_size = maxsp*ntout*nbox
!      WRITE(6,'(A15,ES10.3)')"  OUTVAR_SIZE =",outvar_size
!
!      IF(outvar_size.GT.max_outvar_size)THEN
!        WRITE(6,*) '--ERROR--, output array is TOO LARGE'
!        WRITE(6,*) ' outvar_size = maxsp*ntout*nbox'
!        WRITE(6,*) ' maxsp = ',maxsp
!        WRITE(6,*) ' ntout = ',ntout
!        WRITE(6,*) ' nbox  = ',nbox
!        WRITE(6,'(A11,ES10.3)')"  must be <",max_outvar_size
!        WRITE(6,*) '=> REDUCE SIZE OF REQUESTED OUTPUT AND RESUBMIT <='
!        STOP
!      ENDIF
!
!---------------------------------------------------------------
! read photolysis frequencies (if appropriate)
! all reads are ascii: NetCDF routine saves more info for later o/p
!---------------------------------------------------------------

      IF (numhv.gt.0) THEN
        WRITE(6,*) 'reading the photolysis data ...'

        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
          CALL readj4_ncdf(lout, &  !ncid_out,                  &
               numsp,chrsp,nt1chromo,chromo1cf,            &
               chromomedcf,ntmedchromo,chromotopcf,        &
               numtet,xang,                                &
               rat1pho,coef1pho,ratmedpho,coefmedpho,      & 
               rattoppho,coeftoppho,                       &
               mxjtab,mxsza,llin,                          &
               njtab,nsza,idjtab,jsza,jvref,jnam,jreac)
        ELSE
          CALL readj4(lout,numsp,chrsp,nt1chromo,chromo1cf,&
               chromomedcf,ntmedchromo,chromotopcf,        &
               numtet,xang,                                &
               rat1pho,coef1pho,ratmedpho,coefmedpho,      &
               rattoppho,coeftoppho)
        ENDIF !(iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

        WRITE(6,*) 'DONE PHOT READ...',nbox,ibox
      ENDIF ! (numhv.gt.0) THEN

!---------------------------------------------------------------
! read data that are required to compute deposition (AFTER readkey)
!---------------------------------------------------------------

      IF (depos_fg.GT.0) THEN
        WRITE(6,*) 'reading chem params for ground deposition ...'

        IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
          CALL readdep3_ncdf(ncid_in,lout,chrsp,numsp,iscape, &
                             ndepspe,depnamspe,depdatspe,iddepspe)
        ELSE ! i.e. NOT (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)
          CALL readdep3(lout,chrsp,numsp,iscape,   &
                        ndepspe,depnamspe,depdatspe,iddepspe)
        ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2)

        WRITE(6,*) 'INITIALIZING SURFACE DEPOSITION DATA ...'
        CALL surfdep_wesely89(lout, iseas,                   &
                              rik,rlu0k,rack,rgsSk,rgsOk,    &
                              rclSk,rclOk,z0k)

      ENDIF

!---------------------------------------------------------------
! write some screen output for checking
!---------------------------------------------------------------

      WRITE(6,*)'sumcfix,waterfix,noxfix,presfix =', &
                 sumcfix,waterfix,noxfix,presfix
      WRITE(6,*)'sza =',szaval
      WRITE(6,*)'cnv =',cnv

      IF (OFR_fg.GT.0) THEN
        WRITE(6,*)'USING OFR MODE'
        WRITE(6,*)'jh2o, jo2 =',jh2o,jo2
      ENDIF

      WRITE(6,*) 'END OF ENVIRONMENTAL DATA READ'

!==================================
      RETURN
      END SUBROUTINE get_envinp
!==================================
