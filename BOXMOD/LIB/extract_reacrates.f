! This subroutine outputs rates for the selected reactions
! in units of molec-(n-1) cm-3(n-1) s-1

      SUBROUTINE extract_reacrates(qfor, idreacs, rate_tab)
        IMPLICIT NONE
        REAL, INTENT(IN)                 :: qfor(:)
        INTEGER, ALLOCATABLE,INTENT(IN)  :: idreacs(:)
        REAL, ALLOCATABLE,INTENT(INOUT)  :: rate_tab(:)
        
        INTEGER  :: i, ire
        
        ! check that idreacs and rate_tab dimensions match
        IF (SIZE(idreacs) .NE. SIZE(rate_tab)) THEN
          WRITE(6,*) '--error-- in extract_reacrates'
          WRITE(6,*) 'idreacs and rate_tab sizes do not match'
          WRITE(6,*) SIZE(idreacs), ' ', SIZE(rate_tab)
          STOP
        ENDIF
        
        DO i = LBOUND(idreacs,1), UBOUND(idreacs,1)
          ire = idreacs(i)
          IF (ire < LBOUND(qfor,1) .OR. ire > UBOUND(qfor,1)) THEN
            WRITE(6,*) '--error-- in extract_reacrates'
            WRITE(6,*) 'index of reaction ',i
            WRITE(6,*) 'is out of bounds for qfor'
            WRITE(6,*) ' idreacs(',i,')=',idreacs(i)
            WRITE(6,*) 'qfor range:', LBOUND(qfor,1),'-',UBOUND(qfor,1)
            STOP
          ELSE
            rate_tab(i) = qfor(ire)
          ENDIF
        ENDDO        
      END SUBROUTINE extract_reacrates
