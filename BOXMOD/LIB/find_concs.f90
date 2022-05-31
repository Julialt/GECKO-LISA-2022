      SUBROUTINE find_concs(cbox)
!!-----------------------------------------------------------------
! rates of :
! dilution and exchange between the 2 boxes (box-specific)
! emission & deposition (box-specific) 
!!-----------------------------------------------------------------

      USE flags_module,ONLY: iofmt_fg,ro2_fg,dimer_fg
      USE io_units_module,ONLY: lro2
      USE time_mgmt_module,ONLY: iskip,nskip,time,tout
      USE akparameter_module
      USE forcing_params_module,ONLY: water
      USE module_data_gecko_main

      IMPLICIT NONE

      REAL                  :: totcro2
      REAL,DIMENSION(maxsp) :: cbox
!------------------------------------------------------
      !PRINT*,"starting find_concentrations"

! compute RO2 concentration
      cro2(1:nclro2) = 0.

      IF (ro2_fg.ne.0) THEN
        DO i=1,nclro2
          DO k=1,numchemro2(i)
            cro2(i)=cro2(i)+cbox(idchemro2(k,i))
          ENDDO
        ENDDO

! set CH3O2 concentration (and check for 0 value)
      cmeo2 = cbox(idch3o2)
      IF (cmeo2.LT.1E-13) cmeo2 = 1E-13

!!JMLT for GECKO-vs-MECHGEN comparisons
!!       set constrained SumRO2,Su(m)CRO3
!!       SuCRO2 = PERO9   must be assigned species index numsp
!!       SumRO2 = sum of MEPERO + PERO1:8, species index numsp-1
      totcro2 = cmeo2
      DO i=1, nclro2-1
        totcro2 = totcro2 + cro2(i)
      ENDDO
      IF(idSumRO2.GT.0)THEN
        IF(cbox(idSumRO2).GT.small.AND.totcro2.GT.small) THEN
        cmeo2 = cbox(idSumRO2) * cmeo2/totcro2
        DO i=1, nclro2-1
          cro2(i) = cbox(idSumRO2) * cro2(i)/totcro2
        ENDDO
        ENDIF
      ENDIF
      IF(idSuRCO3.GT.0)THEN
        IF(cbox(idSuRCO3).GT.small) THEN
          cro2(nclro2) = cbox(idSuRCO3)
        ENDIF
      ENDIF
! END JMLT for GECKO-vs-MECHGEN

! write time & pre-solver [RO2] to ascii output file
! (netCDF output is done POST-solver)
!        IF (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2) THEN
        !IF (time.EQ.0.OR.iskip.eq.nskip) &
        IF (iskip.EQ.0.OR.iskip.eq.nskip) &
          WRITE(lro2,'(11(ES13.3))') &
                time,tout,(cro2(i),i=1,nclro2)
!        ENDIF ! (iofmt_fg.EQ.1.OR.iofmt_fg.EQ.2)

      ENDIF

! compute dimer concentration
      IF (dimer_fg.EQ.1) THEN
        cdimer(1:ncldimer)=0.
        DO i=1,ncldimer
          DO k=1,numchemdimer(i)
            cdimer(i)=cdimer(i)+cbox(idchemdimer(k,i))
          ENDDO
        ENDDO
      ENDIF ! (dimer_fg.EQ.1)

! set sum of NOx for this box
      IF (noxfix(ibox) == 1) sumnox = cbox(idno) + cbox(idno2)

! apply water concentration constraint from input file 
! and O2 concentration constraint from sumc
! (overwrite calculation-updated concentrations)
      cbox(idh2o) = water(ibox)
      cbox(ido2dic) = sumc(ibox)*0.2

      !PRINT*,"after find_concentrations"
! -------------------------------------
      END SUBROUTINE find_concs
! -------------------------------------

