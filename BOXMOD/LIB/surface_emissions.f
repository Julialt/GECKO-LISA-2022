! This routine calculates emissions flux depENDing on surfaces proportions
      SUBROUTINE surface_emissions(lout,timemod, xsurf, surface_emi,
     &                             eflux)

      USE akparameter_module
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: lout
      REAL, INTENT(IN)   :: timemod, xsurf(msur)
      TYPE(surface_data), INTENT(IN)  :: surface_emi(msur)
      REAL, INTENT(OUT)  :: eflux(maxsp)

      INTEGER    :: isurf, iemi, j
      TYPE(species_data)  :: s
      REAL       :: tempval, m, p, t1, t2, v1, v2
      LOGICAL    :: found


      DO isurf = 1, msur
        IF (xsurf(isurf) .EQ. 0.) CYCLE
        IF (surface_emi(isurf)%nemis .eq. 0) CYCLE
        DO iemi = 1, surface_emi(isurf)%nemis
          s = surface_emi(isurf)%emission(iemi)
          IF (.NOT. s%activefg) CYCLE
! linearly interpolate between the two closest points
          tempval = 0.
          found = .FALSE.
          DO j=1, s%npoints-1
            IF (s%table(j,1).LE.timemod .AND.
     &          s%table(j+1,1).GT.timemod) THEN
              found = .TRUE.
              t1 = s%table(j,1)
              t2 = s%table(j+1,1)
              v1 = s%table(j,2)
              v2 = s%table(j+1,2)
            ! v = m*t+p
              m  = (v2-v1)/(t2-t1)
              p  = v2-m*t2
              tempval = m*timemod+p
            ENDIF
          ENDDO
          IF (.NOT. found) CYCLE ! we Don't constrain species IF we can't
                                 ! find the right points

          SELECT CASE (s%unit)
            CASE ('molec cm-2')
              tempval = tempval
            CASE DEFAULT
              WRITE(lout,*) ' error3 in apply_constraint'
              WRITE(lout,*) ' unit ', s%unit
              WRITE(lout,*) ' is unknown'
              STOP
          END SELECT
          eflux(s%index)= eflux(s%index) + tempval*xsurf(isurf)
        ENDDO
      ENDDO

      END SUBROUTINE surface_emissions
