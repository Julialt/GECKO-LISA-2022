      SUBROUTINE get_prev(cbox)
! read results of previous run (if supplied)

      USE flags_module,ONLY: iofmt_fg,soa_fg
      USE io_units_module,ONLY: lout
      USE time_mgmt_module
      USE printphoto_module
      USE NetCDF_vars_module,ONLY: ncid_prev,ncid_out
      USE akparameter_module
      USE solver_params_module,ONLY: dtmin,dtmax,numit
      USE forcing_params_module
      USE module_data_gecko_main

      IMPLICIT NONE

      REAL,DIMENSION(maxsp),intent(inout) :: cbox ! = conc(1:numsp,ibox)
      REAL    :: maerbox    ! = maer(ibox)
      REAL    :: ctotaerbox    ! = ctotaer(ibox)
!--------------------------------------------------------
! read previous ncdf file and check that dictionaries are compatible:
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

        CALL readprev_ncdf(ncid_prev,ncid_out, &
                           numsp,nsat,ndim, &
                           time,cbox,caer,ibox)
        IF(soa_fg.EQ.2)THEN ! need previous maer, ctotaer 
          CALL soa_dyn_update(cbox,ctotaerbox,maerbox)
        ENDIF
      ELSE

!--open previous binary output file (not working now: idsat superseded)
        CALL readprev(time,chrsp,idsat,ndim,cbox,caer,ibox)

      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

! -------------------------------------
      END SUBROUTINE get_prev
! -------------------------------------
