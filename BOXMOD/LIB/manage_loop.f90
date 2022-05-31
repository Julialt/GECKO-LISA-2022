      SUBROUTINE manage_loop

      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: delt,iskip,itout,nskip,tout, &
                                 time,timemod,tstop
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main,ONLY: ibox

      IMPLICIT NONE

      INTEGER :: j

!------------------------------------------------------
      WRITE(lout,*)' Time integration:'
      PRINT*,' Time integration:'

!----------------------------------------------------
! Timestep is set in the input file EITHER as:
! SKIP: # of internal timesteps between output points
! .OR.
! NPRT: total # of output points (including initialization point) 
! EACH option calculates nskip.

! start integration loop 
! the "delt/2" term prevents precision issues giving false test results
      DO WHILE (tout+delt/2..lt.tstop) 

! manage time, time in modulo 24h, time at end of timestep
        tout   = time + delt
        timemod = MODULO(time, 86400.)

! update output indicator and output time index
        iskip = iskip + 1

        IF (iskip.EQ.nskip) THEN
          itout = itout + 1
          WRITE(6,*) time,tout,itout
        ELSE
          WRITE(6,*) time,tout
        ENDIF

! do forcing/calculation/output for each box
        DO ibox = 1,nbox
          CALL solve_box
        ENDDO

! update time after solving BOTH boxes
        time = tout
        
! manage time
        IF (iskip.EQ.nskip) THEN
           iskip=0  ! SKIP case
        ENDIF

      ENDDO

! end integration loop (previously 300 CONTINUE)
          
!------------------------------------------------------
      END SUBROUTINE manage_loop
!------------------------------------------------------
