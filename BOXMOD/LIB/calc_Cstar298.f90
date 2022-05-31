      SUBROUTINE calc_Cstar298
! calculate Cstar at 298K for dynamic mode, using either
! Nannoolal or SIMPOL pvap scheme

      USE flags_module,ONLY: iofmt_fg,pvap_fg
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE NetCDF_vars_module,ONLY: ncid_out

      !$ use OMP_LIB
      IMPLICIT NONE

      REAL :: log10Psat
      REAL :: Trb,psat298
      REAL,PARAMETER :: T298 = 298. ! Temperature (K)
      REAL,PARAMETER :: Ratm = 0.0820578 ! atm L K-1 mol-1
! -------------------------------------
      log10Psat = 0.
      psat = 0.

! Myrdal & Yalkowski method is not used in Dynamic mode
! (The code is included` here just in case it's ever wanted)
!      IF (pvap_fg.EQ.1) THEN 
!!$OMP PARALLEL DO private(i, log10Psat)
!        DO i=1,numsp
!          log10Psat= - (86. + 0.4*tau(myid(i)) + 1421*HBN(myid(i))) &
!                     * (Tb(myid(i))-T298)/(19.1*T298)          &
!                     + (-90.0-2.1*tau(myid(i)))/19.1            &
!                     * ((Tb(myid(i))-T298)/T298-log(Tb(myid(i))/T298)) 
!            psat298=exp(log10Psat*log(10.))
!             Cstar298(i) = psat298*wmol(i)/(Ratm*T298) * 1.E+9
!         ENDDO
!!$OMP END PARALLEL DO

!------------------------
! Nannoolal method
!------------------------
!      ELSE IF (pvap_fg.EQ.2) THEN 

      IF (pvap_fg.EQ.2) THEN 
! dynamic mode: pvap info is indexed (1:numsp), mechanism order
! Cstar = Psat(T)(atm)*MW(g/mol)/[R(atm.L/K/mol)*T(K)]*1e3(L/m3)*1e6(ug/g)

!$OMP PARALLEL DO private(i,Trb,log10Psat,psat298)
         DO i=1,numsp
           IF(satid(i).NE.0)THEN
             Trb = T298/Tb(satid(i))
             log10Psat  = (4.1012+dB(satid(i)))*((Trb-1.)/(Trb-(1./8.)))
             psat298    = exp(log10Psat*log(10.))
             Cstar298(i)= psat298*wmol(i)/(Ratm*T298) * 1.E+9
           ENDIF
         ENDDO
!$OMP END PARALLEL DO

!------------------------
! Simpol method
!------------------------
       ELSE IF (pvap_fg.EQ.3) THEN 
!$OMP PARALLEL DO private(i,j,log10Psat,psat298)
!! This OMP has a problem: it restricts the size of 'i'
         DO i=1,numsp
           log10psat = 0
           IF(satid(i).NE.0)THEN
             DO j=1,31
               IF (simpgroup(satid(i),j).NE.0) THEN
                 log10psat = log10psat + simpgroup(satid(i),j)* &
                                         ((bk(j,1)/T298)+       &
                                           bk(j,2)+             &
                                          (bk(j,3)*T298)+       &
                                          (bk(j,4)*log(T298)))
               ENDIF
             ENDDO
             psat298     = exp(log10Psat*log(10.))
             Cstar298(i) = psat298*wmol(i)/(Ratm*T298) * 1.E+9
           ENDIF
         ENDDO
!$OMP END PARALLEL DO

       ENDIF ! (pvap_fg.EQ._)
! -------------------------------------
! write to Netcdf output file, if appropriate flag
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        CALL eznc_put_1Dreal(ncid_out,"Cstar298", &
                                       Cstar298(1:numsp), &
                                                1,numsp)
      ENDIF

! -------------------------------------
      END SUBROUTINE calc_Cstar298
! -------------------------------------

