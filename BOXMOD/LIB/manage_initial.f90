      SUBROUTINE manage_initial

      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: time,tstart,tout,itout,timemod,&
                                 nskip,ntprint,ntout
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main,ONLY: ibox

      IMPLICIT NONE

!------------------------------------------------------
      PRINT*,"nskip,ntprint,ntout",nskip,ntprint,ntout

! -------------------------------------
      time=tstart
      itout = 1

! time in modulo 24h
      timemod = MODULO(time, 86400.)

! calculate time-independent derived quantities
      !IF(soa_fg.EQ.2) CALL calc_Cstar298

! do initial calculations/forcings/output for each box

      DO ibox = 1,nbox
        CALL solve_box
      ENDDO

      PRINT*,' Initial Time:',itout, timemod
      WRITE(lout,*)' Initial Time:',itout, timemod

      PRINT*,time,tout,itout,"initial o/p"

!------------------------------------------------------
      END SUBROUTINE manage_initial
!------------------------------------------------------
