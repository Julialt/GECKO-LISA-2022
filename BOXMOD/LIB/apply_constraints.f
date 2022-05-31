      SUBROUTINE apply_constraints(lout, time, cbot,temp, sumc,eflux,
     &                     cons_spec, emi_spec,idNO,idNO2, idNO3)

      USE akparameter_module
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: lout, idNO, idNO2, idNO3
      REAL,    INTENT(IN)   :: time, temp, sumc
      TYPE(species_data), INTENT(IN) :: cons_spec(maxconst)
      TYPE(species_data), INTENT(IN) :: emi_spec(maxem)
      REAL,    INTENT(INOUT)  :: cbot(maxsp), eflux(maxsp)

      !local
      INTEGER      :: i,j
      REAL         :: pres, ppbfact, pptfact
      REAL         :: tempval, v1, v2, t1, t2, m, p
      LOGICAL      :: found
      TYPE(species_data)   :: s

      !calculate pressure (Pa)
      pres = sumc*8.32*temp*1e5/(6.022E+22)

      !calculate ppb factor : C(molec cm-3) = C(ppb) * ppbfact
      ppbfact = 6.022E23*pres/(8.314*temp*1e15)

      !calculate ppt factor : C(molec cm-3) = C(ppt) * ppbtact
      pptfact = ppbfact*1e-3


      ! start with constant species

C       DO i=1, maxconst
C         if (const_species(i)%activefg) THEN
C           tempval = 0.
C           SELECT CASE (const_species(i)%unit)
C             CASE ('molec cm-3')
C               tempval = const_species(i)%value
C             CASE ('ppb')
C               tempval = const_species(i)%value * ppbfact
C             CASE ('ppt')
C               tempval = const_species(i)%value * pptfact
C             CASE DEFAULT
C               WRITE(lout,*) ' error1 in apply_constraint'
C               WRITE(lout,*) ' unit ', const_species(i)%unit
C               WRITE(lout,*) ' is unknown'
C               STOP
C           END SELECT
C           IF (const_species(i)%name .eq. 'NOx') THEN
C             CALL noxconstraint(cbot, tempval, idNO, idNO2, idNO3)
C           ELSE
C             cbot(const_species(i)%index) = tempval
C           ENDIF
C         ENDIF
C       ENDDO

      ! then constrained species
      DO i = 1, maxconst
        s = cons_spec(i)
        IF (s%activefg) THEN
          ! linearly interpolate between the two closest points
          tempval = 0.
          found = .FALSE.
          DO j=1, s%npoints-1
            IF (s%table(j,1).le.time .and. s%table(j+1,1).gt.time) THEN
              found = .TRUE.
              t1 = s%table(j,1)
              t2 = s%table(j+1,1)
              v1 = s%table(j,2)
              v2 = s%table(j+1,2)
            ! v = m*t+p
              m  = (v2-v1)/(t2-t1)
              p  = v2-m*t2
              tempval = m*time+p
            ENDIF
          ENDDO
          IF (.not. found) CYCLE ! we don't constrain species if we can't
                                 ! find the right points

          SELECT CASE (s%unit)
            CASE ('molec cm-3')
              tempval = tempval
            CASE ('ppb')
              tempval = tempval * ppbfact
            CASE ('ppt')
              tempval = tempval * pptfact
            CASE DEFAULT
              WRITE(lout,*) ' error2 in apply_constraint'
              WRITE(lout,*) ' unit ', s%unit
              WRITE(lout,*) ' is unknown'
              STOP
          END SELECT
          IF (s%name .eq. 'NOx') THEN
            CALL noxconstraint(cbot, tempval, idNO, idNO2, idNO3)
          ELSE
            cbot(s%index) = tempval
          ENDIF
        ENDIF
      ENDDO

      ! constrain emissions
      DO i = 1, maxem
        s = emi_spec(i)
        IF (s%activefg) then
          ! linearly interpolate between the two closest points
          tempval = 0.
          found = .FALSE.
          DO j=1, s%npoints-1
            IF (s%table(j,1).lt.time .and. s%table(j+1,1).ge.time) THEN
              found = .TRUE.
              t1 = s%table(j,1)
              t2 = s%table(j+1,1)
              v1 = s%table(j,2)
              v2 = s%table(j+1,2)
            ! v = m*t+p
              m  = (v2-v1)/(t2-t1)
              p  = v2-m*t2
              tempval = m*time+p
            ENDIF
          ENDDO
          IF (.not. found) CYCLE ! we don't constrain species if we can't
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
          eflux(s%index)= tempval

        endif
      enddo


      END SUBROUTINE apply_constraints

!====================================================================
      SUBROUTINE noxconstraint(cbot, value, idNO, idNO2, idNO3)

      USE akparameter_module
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: idNO, idNO2, idNO3
      REAL, INTENT(IN)    :: value ! Nox constraint in molec cm-3
      REAL, INTENT(INOUT) :: cbot(maxsp)

      REAL                :: sumnox, diff, xno2, xno, xno3

!! Original version includes no3 in sumnox but NOT in rest of accounting !!
C       sumnox=cbot(idno)+cbot(idno2)+cbot(idno3)

      sumnox=cbot(idno)+cbot(idno2)

      diff=value-sumnox
      IF((cbot(idno)+cbot(idno2)) .ne. 0.) THEN
        xno2=cbot(idno2)/(cbot(idno)+cbot(idno2))
        xno=1.-xno2

        cbot(idno)=cbot(idno)+diff*xno
        cbot(idno2)=cbot(idno2)+diff*xno2
      ENDIF

!! Version fully including no3 !!
      !sumnox=cbot(idno)+cbot(idno2)+cbot(idno3)
      !diff=value-sumnox
      !IF((cbot(idno)+cbot(idno2)+cbot(idno3)) .ne. 0.) THEN
      !  xno =cbot(idno) /sumnox
      !  xno2=cbot(idno2)/sumnox
      !  xno3=cbot(idno3)/sumnox

      !  cbot(idno)=cbot(idno)+diff*xno
      !  cbot(idno2)=cbot(idno2)+diff*xno2
      !  cbot(idno3)=cbot(idno3)+diff*xno3
      !ENDIF

      END SUBROUTINE noxconstraint
