      SUBROUTINE soa_equil(cbox,ctotaerbox,maerbox)

      USE time_mgmt_module,ONLY: delt
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE fundamental_consts_module,ONLY: multimass

      IMPLICIT NONE

      REAL    :: maerbox    ! = maer(ibox)
      REAL    :: ctotaerbox ! = ctotaer(ibox)
      REAL,DIMENSION(maxsp) :: cbox ! = conc(1:numsp,ibox)

! -------------------------------------
      !PRINT*,"starting equilibrium soa section"

! compute dilution for SOA species (use a simple euler equation)
      DO isp=1,nsat
        caer(isp)=caer(isp)*exp(-rdil*delt)
      ENDDO

      CALL soapartition(temp,nsat,idsat,cbox,caer,psat,cnv)

! compute the new concentration of SOA (in molec/cc)
! compute the new mass concentration of SOA (in ug/m3)
      ctotaerbox = cnv
      maerbox = 0.

!$OMP PARALLEL DO private(isp) reduction(+:maerbox,ctotaerbox)
        DO isp = 1,nsat
          ctotaerbox = ctotaerbox+caer(isp)
          maerbox    = maerbox   +caer(isp)*wmol(idsat(isp))*multimass
        ENDDO
!$OMP END PARALLEL DO

      !PRINT*,"after equilibrium soa section"

! -------------------------------------
      END SUBROUTINE soa_equil
! -------------------------------------

