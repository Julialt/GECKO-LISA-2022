      SUBROUTINE soa_dyn_update(cbox,ctotaerbox,maerbox)

      USE flags_module,ONLY: dimer_fg
      USE akparameter_module,ONLY: maxsp
      USE module_data_gecko_main,ONLY: wmol,nsat,idasat,isp,cnv, &
                                            ndim,iddsat
      USE fundamental_consts_module,ONLY: multimass

      !$ USE OMP_LIB

      IMPLICIT NONE

      REAL,DIMENSION(maxsp),INTENT(IN) :: cbox  ! = conc(:,ibox)
      REAL,INTENT(OUT)    :: maerbox    ! = maer(ibox)
      REAL,INTENT(OUT)    :: ctotaerbox ! = ctotaer(ibox)

! -------------------------------------
      !PRINT*,"updating dynamic soa"

! compute the new concentration of SOA (in molec/cc)
! compute the new mass concentration of SOA (in ug/m3)
      ctotaerbox = cnv
      maerbox = 0.

!$OMP PARALLEL DO private(isp) reduction(+:maerbox,ctotaerbox)
      DO isp = 1,nsat
        ctotaerbox = ctotaerbox + cbox(idasat(isp))
        maerbox    = maerbox    + cbox(idasat(isp)) &
                                * multimass*wmol(idasat(isp))
      ENDDO
!$OMP END PARALLEL DO

      IF (dimer_fg.EQ.1) THEN
!$OMP PARALLEL DO private(isp) reduction(+:maerbox,ctotaerbox)
      DO isp = 1,ndim
        ctotaerbox = ctotaerbox + cbox(iddsat(isp))
        maerbox    =    maerbox + cbox(iddsat(isp)) &
                                * multimass*wmol(iddsat(isp))
      ENDDO
!$OMP END PARALLEL DO
      ENDIF ! (dimer_fg)

! -------------------------------------
      END SUBROUTINE soa_dyn_update
! -------------------------------------

