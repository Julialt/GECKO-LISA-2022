      SUBROUTINE setup_printphoto

! prepare table of photolysis reactions indices we want to print

      USE printphoto_module
      USE module_data_gecko_main,ONLY: chrsp,idrestoi,idpdstoi,numstoi,&
                                       restoicf,pdstoicf,ire,numhv,idhv
      IMPLICIT NONE

      INTEGER :: i,j
!=====================================

      ! find number of photolysis reactions of inorganics
      n_printphoto = 0
      DO i = 1, numhv
        ire = idhv(i)
        IF(chrsp(idrestoi(ire, 1))(1:4) == "GO3 "  .OR. &
           chrsp(idrestoi(ire, 1))(1:5) == "GNO2 " .OR. &
           chrsp(idrestoi(ire, 1))(1:5) == "GNO3 " .OR. &
           chrsp(idrestoi(ire, 1))(1:6) == "GH2O2 " .OR. &
           chrsp(idrestoi(ire, 1))(1:6) == "GHNO2 " .OR. &
           chrsp(idrestoi(ire, 1))(1:6) == "GHNO3 " .OR. &
           chrsp(idrestoi(ire, 1))(1:6) == "GHNO4 " .OR. &
           chrsp(idrestoi(ire, 1))(1:6) == "GCH2O " .OR. &
           chrsp(idrestoi(ire, 1))(1:7) == "GN01000" .OR. & ! CH3ONO
           chrsp(idrestoi(ire, 1))(1:7) == "GD02000" .OR. & ! CH3CHO
           chrsp(idrestoi(ire, 1))(1:6) == "GGLYOX" .OR. & ! CHO3CHO
           chrsp(idrestoi(ire, 1))(1:7) == "GMGLYOX" ) THEN ! CH3COCHO
           n_printphoto = n_printphoto + 1
         ENDIF
      ENDDO

!-----------------------------------------------------
      IF (n_printphoto .GT. 0) THEN

        ALLOCATE(idprintphoto(n_printphoto),  &
                 photorates(n_printphoto),    &
                 photoreac_names(n_printphoto))
! cnjv = character version of (n_printphoto -1), for output formatting

        IF(n_printphoto.LT.10) THEN
          ALLOCATE(character(1) :: cnjv)
          WRITE(cnjv,'(i1)') n_printphoto
        ELSEIF(n_printphoto.LT.100) THEN
          ALLOCATE(character(2) :: cnjv)
          WRITE(cnjv,'(i2)') n_printphoto
        ENDIF

        j = 0
        DO i = 1, numhv
          ire = idhv(i)
          IF(chrsp(idrestoi(ire, 1))(1:4) == "GO3 "  .OR. &
             chrsp(idrestoi(ire, 1))(1:5) == "GNO2 " .OR. &
             chrsp(idrestoi(ire, 1))(1:5) == "GNO3 " .OR. &
             chrsp(idrestoi(ire, 1))(1:6) == "GH2O2 " .OR. &
             chrsp(idrestoi(ire, 1))(1:6) == "GHNO2 " .OR. &
             chrsp(idrestoi(ire, 1))(1:6) == "GHNO3 " .OR. &
             chrsp(idrestoi(ire, 1))(1:6) == "GHNO4 " .OR. &
             chrsp(idrestoi(ire, 1))(1:6) == "GCH2O " .OR. &
             chrsp(idrestoi(ire, 1))(1:7) == "GN01000" .OR. & ! CH3ONO
             chrsp(idrestoi(ire, 1))(1:7) == "GD02000" .OR. & ! CH3CHO
             chrsp(idrestoi(ire, 1))(1:6) == "GGLYOX" .OR. & ! CHO3CHO
             chrsp(idrestoi(ire, 1))(1:7) == "GMGLYOX" ) THEN ! CH3COCHO
             j = j + 1
             idprintphoto(j) = ire
             photoreac_names(j) = printreaction(ire, chrsp, &
               numstoi,idrestoi, idpdstoi, restoicf, pdstoicf)
          ENDIF
        ENDDO

      ENDIF !(n_printphoto.GT.0) 

!------------------------------------------------------
      RETURN
      END SUBROUTINE setup_printphoto
!=====================================
