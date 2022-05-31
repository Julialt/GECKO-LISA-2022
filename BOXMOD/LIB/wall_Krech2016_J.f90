!=======================================================================
      SUBROUTINE wall_Krech2016_J
! PURPOSE: calculate rates of gas-wall & wall-gas trasnfer reactions
!          for chamber parameterization 3: Krechmer et al (2016)
!          characterizing Jimenez chanber
! CAUTION: for use with  SIMPOL pvap parameterization
! REF1: Krechmer, Pagonis, Ziemann & Jimenez (2016),ES&T,50,11,5757-5765.
! DOI: 10.1021/acs.est.6b00606
! see subroutine wall_Krech2016_Z for original published coeffs in
! Ziemann chamber
! REF2: Xiaoxi Liu, submitted manuscript, 2019
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
      REAL :: Cw,Cstar_eff

      REAL,PARAMETER :: ratewin = 1.000E-03 ! s-1, for all species
      REAL,PARAMETER :: mw_wall = 250.
! min & max derived wall "concentrations", ug/m3
      REAL,PARAMETER :: Cwmin = 0.1 ! 
      REAL,PARAMETER :: Cwmax = 1.0E+4 
! coeffs a & b for eqn Cw = a * Cstar^b, Cstar > woucfc
      REAL,PARAMETER :: woucfa1 = 16.0 ! a (Krechmer 2019)
      REAL,PARAMETER :: woucfb1 = 0.6  ! b (Krechmer 2016)
! coeffs a & b for eqn Cw = a * Cstar^b, Cstar < woucfc
      REAL,PARAMETER :: woucfa2 = 6.49  ! a, Xiaoxi Liu ms, Cstar < 1e3
      REAL,PARAMETER :: woucfb2 = 0.75  ! b, Xiaoxi Liu ms, Cstar < 1e3
! woucfc = transition point between 2 different wall equations
! for woucfa1 = 19
      REAL,PARAMETER :: woucfc  = 4e2 ! transition btw 2 sets of coeffs
      REAL,PARAMETER :: chamfac = 0.75 ! conversion from PJZ(1.5) to
!                                      ! JLJ(2) chamber surface-to-vol ratio

!-----------------------------------------------------------------------
! ONLY use this subroutine with SIMPOL pvap scheme
! NB: a * b coeffs were derived by Xiaoxi Liu valid using SIMPOL pvap(i)

      IF(pvap_fg.NE.3)THEN
        WRITE (lout,*) '--error--, in subroutine wall_Krech2016_J.f90'
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
! Parameterization is for Cstar calculated per SIMPOL at actual T
! using effective molecular wt of CHAMBER WALL material
!   Keq = Cw (ug/m3) / Cstar (ug/m3)
!   where Cw = a * (Cstar ^ b), for various Cstar ranges
! => kwg = kgw * Cstar/Cw

      DO i=1,numwou
        ire = idwou(i)

        Cstar_eff = psat(idrestoi(ire,1)) * mw_wall & 
              /  (Ratm*temp) * 1.E+9

! use simple Cstar correlation derived from alkanol products
! to convert NANNOOLAL c* to SIMPOL c* (local)
! because wall parameterization relies on SIMPOL C*
!        IF(pvap_fg.EQ.2)THEN 
!          Cstar = nan_to_sim(Cstar)
!        ENDIF !(pvap_fg)

! Cw : derived wall "concentration" in units of ug/m3
        IF(Cstar_eff.GT.woucfc)THEN
          Cw = MAX(Cwmin,(MIN(Cwmax,woucfa1 * Cstar_eff ** woucfb1)))
        ELSE
          Cw = MAX(Cwmin,(MIN(Cwmax,woucfa2 * Cstar_eff ** woucfb2)))
        ENDIF
! adjust Cw for chamber surface-to-volume ratio
        Cw = Cw*chamfac

! For Cstar/Cw, units cancel; ratewin units are s-1
        qfor(ire) = ratewin * Cstar_eff / Cw

! DEBUG/DIAGNOSTIC O/P
        IF(i.EQ.63) WRITE(52,*) idwou(i),ratewin,Cstar_eff,Cw

      ENDDO ! i=1,numwou

!=======================================================================
      END SUBROUTINE wall_Krech2016_J
!=======================================================================
