      SUBROUTINE printsteadystate

      USE akparameter_module,ONLY: maxsp
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main,ONLY: chrsp,conc
      USE steadystate_vars_module

      IMPLICIT NONE

      INTEGER i,j,ibox
      CHARACTER(3):: cnbox
!-----------------------------------------------------------

      WRITE(cnbox,'(i1)')nbox

      non_zero_species_number = 0
! find all non-zero species ie more than 1 molec cm-3...
      non_zero_species_number = COUNT(conc(:,1) > 1.)
      ALLOCATE(initspecnames(non_zero_species_number), &
               initspecconc(non_zero_species_number,nbox))
      j = 0
      DO i=1,maxsp
        IF (conc(i,1) > 1.) THEN
          j = j+1
          initspecnames(j) = chrsp(i)
          DO ibox = 1,nbox
            initspecconc(j,ibox) = conc(i,ibox)
          ENDDO
        ENDIF
      ENDDO
! we print the steadystate concentrations for all non-zero species
      OPEN(87, file = "steadystate.key", status = "unknown")
      DO i = 1, non_zero_species_number
        WRITE(87, '(A5, A8, '//cnbox//'(1X, ES9.2), 1X, ES9.2)') &
                    'REAC ', initspecnames(i), &
                             initspecconc(i,:), 0.0
      ENDDO
      CLOSE(87)

! Deallocations
      IF(ALLOCATED(initspecnames)) DEALLOCATE(initspecnames)
      IF(ALLOCATED(initspecconc))  DEALLOCATE(initspecconc)

      RETURN
!-----------------------------------------------------------
      END SUBROUTINE printsteadystate
!-----------------------------------------------------------
