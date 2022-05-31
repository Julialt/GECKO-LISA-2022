! this subroutine prints the reaction at the given index
! in a character*(300)

      FUNCTION printreaction(ire, chrsp,numstoi, 
     &   idrestoi, idpdstoi, restoicf, pdstoicf) RESULT(res)
     
        USE akparameter_module
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: ire
        INTEGER, INTENT(IN) :: idrestoi(maxre,mxleft)
        INTEGER, INTENT(IN) :: idpdstoi(maxre,mxright)
        INTEGER, INTENT(IN) :: numstoi(maxre,2)
        REAL   , INTENT(IN) :: restoicf(maxre,mxleft)
        REAL   , INTENT(IN) :: pdstoicf(maxre,mxright)
        CHARACTER(maxlsp), INTENT(IN) :: chrsp(maxsp)

! variable declaration for function result
        CHARACTER(maxreac_char) :: res
       
        INTEGER  :: i, j, k
        CHARACTER(5) :: temp_num
        
        res = ""
        
        ! loop over reactants
        DO i = 1, numstoi(ire, 1)
          IF (i .GT. 1) res = TRIM(res) // ' +'
          IF (restoicf(ire,i) .NE. 1) THEN
            WRITE(temp_num, '(f5.2)') restoicf(ire,i)
            res = TRIM(res) // ' ' // TRIM(temp_num)
          ENDIF
          res = TRIM(res) // ' ' //  TRIM(chrsp(idrestoi(ire,i)))
        ENDDO
        
        res = TRIM(res) // ' => '
        ! loop over products
        DO i = 1, numstoi(ire, 2)
          IF (i .GT. 1) res = TRIM(res) // ' +'
          IF (pdstoicf(ire,i) .NE. 1) THEN
            WRITE(temp_num, '(f5.2)') pdstoicf(ire,i)
            res = TRIM(res) // ' ' // TRIM(temp_num)
          ENDIF
          res = TRIM(res) // ' ' // TRIM(chrsp(idpdstoi(ire,i)))
        ENDDO 
           
        RETURN  
      END FUNCTION printreaction
