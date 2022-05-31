!=======================================================================
      SUBROUTINE wall_test
! PURPOSE: calculate rates of gas-wall & wall-gas trasnfer reactions
!          for chamber parameterization 2: Krechmer et al (2016)
!          characterizing Ziemann chanber
! CAUTION: for use with  SIMPOL pvap parameterization
! REF: Krechmer, Pagonis, Ziemann & Jimene (2016), ES&T, 50, 11, 5757-5765.
! DOI: 10.1021/acs.est.6b00606
!=======================================================================

      USE flags_module,ONLY: pvap_fg
      USE io_units_module,ONLY: lout
      !USE constraints_module,ONLY: nan_to_sim
      USE forcing_params_module,ONLY: temp
      USE module_data_gecko_main,ONLY: numwin,numwou,idwin,idwou, &
                                       idrestoi,wmol,psat,qfor
      USE fundamental_consts_module,ONLY: Ratm

      IMPLICIT NONE

! internal variables
      INTEGER :: ire,i
      REAL :: Cw,Cstar,Cstar_eff

      REAL,PARAMETER :: ratewin = 2.0E-05 ! s-1, for all species
      REAL,PARAMETER :: mw_wall = 250.
! published values, Krechmer et al 2016
      !REAL,PARAMETER :: Cwmin = 16. ! derived wall "concentration", ug/m3 
      !REAL,PARAMETER :: Cwmax = 1.0E+4 ! derived wall "concentration", ug/m3 (1e4:3e4)
      !REAL,PARAMETER :: woucfa = 16.0 ! a for Cw = a * Cstar^b
      !REAL,PARAMETER :: woucfb = 0.6  ! b for Cw = a * Cstar^b
! Updated values, Jimenez, pers. comm. 2018/2019
      REAL,PARAMETER :: Cwmin = 10. ! derived wall "concentration", ug/m3
      REAL,PARAMETER :: Cwmax = 1.0E+4 ! derived wall "concentration", ug/m3
      REAL,PARAMETER :: woucfa = 19.0 ! a for Cw = a * Cstar^b
      REAL,PARAMETER :: woucfb = 0.6  ! b for Cw = a * Cstar^b

!-----------------------------------------------------------------------
! ONLY use this subroutine with SIMPOL pvap scheme

      IF(pvap_fg.NE.3)THEN
        WRITE (lout,*) '--error--, in subroutine wall_Krech2016_Z.f90'
        WRITE (lout,*) 'this wall parameterization requires SIMPOL'
        STOP
      ENDIF

!-----------------------------------------------------------------------
! setup variables....

!-------------
! Gas -> Wall 
!-------------
! The rate constant comes from the equilibrium constant
! with the wall and the definition of time constant to equilibration
! assuming 2 first order reactions for gas/wall partitioning.

      DO i=1,numwin
        ire = idwin(i)
        qfor(ire) = ratewin
      ENDDO ! i-1,numwim

!-------------
! Wall -> Gas
!-------------
! Use arrhcf(ire,1) as inverse of time constant to equilibration.
! The rate constant comes from the equilibrium constant
! with the wall and the definition of time constant to equilibration
! assuming 2 first order reactions for gas/wall partitioning.
! The time constant to equilibrium tau = 1/kgw + 1/kwg.
! Parameterization is for Cstar calculated at 298K
!   Keq = Cw (ug/m3) / Cstar (ug/m3)
!   where Cw = a * (Cstar ^ b), for various Cstar ranges
! => kwg = kgw * Cstar/Cw

      DO i=1,numwou
        ire = idwou(i)

        !Cstar = psat(idrestoi(ire,1)) & 
        !      * wmol(idrestoi(ire,1)) &
        !      /  (Ratm*temp) * 1.E+9
        Cstar_eff = psat(idrestoi(ire,1)) & 
              * mw_wall &
              /  (Ratm*temp) * 1.E+9

!! OLD COMMENT BELOW
! use simple Cstar correlation derived from alkanol products
! to convert NANNOOLAL c* to SIMPOL c* (local)
! because wall parameterization relies on SIMPOL C*
        !IF(pvap_fg.EQ.2)THEN
        !  Cstar = nan_to_sim(Cstar)
        !ENDIF !(pvap_fg)

! Cw : derived wall "concentration" in units of ug/m3
        Cw = MAX(Cwmin,(MIN(Cwmax,woucfa * Cstar_eff ** woucfb)))

! For Cstar/Cw, units cancel; ratewin units are s-1
        qfor(ire) = ratewin * Cstar_eff / Cw
        IF(i.EQ.1) WRITE(52,*) idwou(i),Cstar_eff,Cw,ratewin,qfor(ire)

      ENDDO ! i=1,numwou

!=======================================================================
      END SUBROUTINE wall_test
!=======================================================================
