! -------------------------------------
      SUBROUTINE calc_psat
! -------------------------------------
! calculate saturation vapor pressure (Psat) using pvap file info

      USE flags_module,ONLY: pvap_fg, soa_fg
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main

      !$ use OMP_LIB
      IMPLICIT NONE

      REAL :: log10Psat
      REAL :: Trb,Tbtmp,dBtmp
! -------------------------------------
      log10Psat = 0.
      psat = 0.

! -------------------------------------
! Myrdal & Yalkowski method
! requires inputs mech.pmy with columns namsat(i), Tb(i), HBN(i), tau(i)
! Only available for equilibrium configuration
! -------------------------------------
      IF (pvap_fg.EQ.1) THEN 
!$OMP PARALLEL DO private(i, log10Psat)
        DO i=1,nsat
          log10Psat= - (86. + 0.4*tau(i) + 1421*HBN(i)) &
                     * (Tb(i)-temp)/(19.1*temp)          &
                     + (-90.0-2.1*tau(i))/19.1            &
                     * ((Tb(i)-temp)/temp-log(Tb(i)/temp)) 
            psat(i)=exp(log10Psat*log(10.))
         ENDDO
!$OMP END PARALLEL DO

!------------------------
! Nannoolal method
! requires inputs mech.pnan with columns namsat(i), Tb(i), dB(i)
! Available for both equilibrium and dynamic configurations
!------------------------
        ELSE IF (pvap_fg.EQ.2) THEN 
          IF (soa_fg.eq.1)THEN
! equilibrium mode: psat is indexed in 1:nsat order (as in dict)
!$OMP PARALLEL DO private(i,Trb,log10Psat)
            DO i=1,nsat
              Trb       = temp/Tb(i)
              log10Psat = (4.1012+dB(i))*((Trb-1.)/(Trb-(1./8.)))
              psat(i)   = exp(log10Psat*log(10.))
            ENDDO
!$OMP END PARALLEL DO

          ELSE ! implies (soa_fg.EQ.2)

! dynamic mode: psat is indexed in 1:numsp order (as in mech)
!$OMP PARALLEL DO private(i,Trb,Tbtmp,dBtmp,log10Psat)
            DO i=1,numsp
              IF(satid(i).NE.0)THEN
                Tbtmp = Tb(satid(i))  
                dBtmp = dB(satid(i))
                Trb   = temp/Tbtmp
                log10Psat = (4.1012+dBtmp)*((Trb-1.)/(Trb-(1./8.)))
                psat(i)   = exp(log10Psat*log(10.))
              ENDIF
            ENDDO
!$OMP END PARALLEL DO

          ENDIF ! (soa_fg.EQ.1 or 2)

!------------------------
! Simpol method
! requires inputs simpol.dat (31 rows each with 4x bk(n) columns)
! and mech.pvap with columns: namsat(i), simpgroup(i,0:30)
! Available for equilibrium (now) and dynamic (soon) configurations
!------------------------
        ELSE IF (pvap_fg.EQ.3) THEN 
          IF (soa_fg.eq.1)THEN
!$OMP PARALLEL DO private(i,j,log10Psat)
            DO i=1,nsat
              log10psat = 0.
              DO j=1,31
                IF (simpgroup(i,j).NE.0) THEN
                  log10psat = log10psat + simpgroup(i,j)* &
                                         ((bk(j,1)/temp)+ &
                                           bk(j,2)+       &
                                          (bk(j,3)*temp)+ &
                                          (bk(j,4)*log(temp)))
                ENDIF
              ENDDO
              psat(i)=exp(log10Psat*log(10.))
            ENDDO
!$OMP END PARALLEL DO

          ELSE ! implies (soa_fg.EQ.2)
! dynamic mode: psat is indexed in 1:numsp order (as in mech)

!$OMP PARALLEL DO private(i,j,log10Psat)
            DO i=1,numsp
              IF (satid(i).NE.0) THEN
                log10psat = 0.
                DO j=1,31
                  log10psat = log10psat + simpgroup(satid(i),j)* &
                                          ((bk(j,1)/temp)+ &
                                            bk(j,2)      + &
                                           (bk(j,3)*temp)+ &
                                           (bk(j,4)*log(temp)))
                ENDDO
                psat(i)   = exp(log10Psat*log(10.))
              ENDIF
            ENDDO
!$OMP END PARALLEL DO

          ENDIF ! (soa_fg.EQ.1 or 2)
        ENDIF ! (pvap_fg.EQ._)

! -------------------------------------
      END SUBROUTINE calc_psat
! -------------------------------------
