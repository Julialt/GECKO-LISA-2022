      SUBROUTINE map_sp_ids
! find species i.d. for common reactants

      USE flags_module,ONLY: tracer_fg
      USE akparameter_module,ONLY: maxlsp
      USE module_data_gecko_main
      USE akparameter_module

      IMPLICIT NONE

      CHARACTER(LEN=maxlsp):: tnam
!---------------------------------------------------
      PRINT*,"starting species map"

      tnam = 'GCH3O2  '
      CALL akspnum(tnam,chrsp,numsp,idch3o2)
      IF (idch3o2.eq.0) GOTO 99

      tnam = 'GHO     '
      CALL akspnum(tnam,chrsp,numsp,idho)
      IF (idho.eq.0) GOTO 99

      tnam = 'GHO2    '
      CALL akspnum(tnam,chrsp,numsp,idho2)
      IF (idho2.eq.0) GOTO 99

      tnam = 'GNO     '
      CALL akspnum(tnam,chrsp,numsp,idno)
      IF (idno.eq.0) GOTO 99

      tnam = 'GNO2    '
      CALL akspnum(tnam,chrsp,numsp,idno2)
      IF (idno2.eq.0) GOTO 99

      tnam = 'GNO3    '
      CALL akspnum(tnam,chrsp,numsp,idno3)
      IF (idno3.eq.0) GOTO 99

      tnam = 'GH2O    '
      CALL akspnum(tnam,chrsp,numsp,idh2o)
      IF (idh2o.eq.0) GOTO 99

      tnam = 'GO2     '
      CALL akspnum(tnam,chrsp,numsp,ido2dic)
      IF (ido2dic.eq.0) GOTO 99

      tnam = 'GO3     '
      CALL akspnum(tnam,chrsp,numsp,ido3)
      IF (ido3.eq.0) GOTO 99

      tnam = 'GISOPRN '
      CALL akspnum(tnam,chrsp,numsp,idisop)
! NB: isoprene is not always present - 
! i.e. "failure to find" is not necessarily an error

      tnam = 'GSumRO2 '
      CALL akspnum(tnam,chrsp,numsp,idSumRO2)

      tnam = 'GSuRCO3 '
      CALL akspnum(tnam,chrsp,numsp,idSuRCO3)

      IF (idSuRCO3.EQ.0)THEN ! try again in case MG/FM format applies)
        tnam = 'GSumRCO3'
        CALL akspnum(tnam,chrsp,numsp,idSuRCO3)
      ENDIF

! define array of ids for rate tracking
      IF (tracer_fg .eq. 1) THEN
        ntr = 3
        idtr(1) = idho
        idtr(2) = idho2
        idtr(3) = idh2o

        WRITE(30,*)'time, trprodHO, trlossHO, trprodHO2, trlossHO2'
        WRITE(30,'(5(1x,ES12.5))') 0.,0.,0.,0.,0.
      ENDIF

      PRINT*,"end of species map"
      RETURN

!---------------------------------------------------------- 
! error messages

99    WRITE(6,*) '--error--, in subroutine map_sp_ids'
      WRITE(6,*) '           species '//tnam//' not found'
      STOP
!---------------------------------------------------------- 

      END SUBROUTINE map_sp_ids
