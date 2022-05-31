!=======================================================================
      SUBROUTINE wall_MatZie2010
! PURPOSE: calculate rates of gas-wall & wall-gas trasnfer reactions
!          for chamber parameterization 1: Matsunaga & Ziemann (2010)
!=======================================================================

      USE forcing_params_module,ONLY: temp
      USE module_data_gecko_main,ONLY: numwin,numwou,idwin,idwou, &
                                       winfac,weqfac,woucf, &
                                       idrestoi,wmol,psat,arrhcf,qfor
      USE fundamental_consts_module,ONLY: Ratm

      IMPLICIT NONE

! internal variables
      INTEGER :: ire,i
      REAL :: lnwin,lnweq
      REAL :: lnrate,lnpvap,arrhcfwin
      REAL :: lnmultifacgw,cmw

!-----------------------------------------------------------------------
! setup variables....

      lnwin = LOG(winfac)
      lnweq = LOG(weqfac)

      lnmultifacgw = LOG(Ratm*1E-3*temp) ! RT in atm.m3.mol-1

!-------------
! Gas -> Wall 
!-------------
! Use arrhcf(ire,1) as inverse of time constant to equilibration.
! The rate constant comes from the equilibrium constant
! with the wall and the definition of time constant to equilibration
! assuming 2 first order reactions for gas/wall partitioning.

      DO i=1,numwin
        ire = idwin(i)
        !arrhcfwin = arrhcf(ire,3)
        arrhcfwin = arrhcf(ire,1)
        lnrate    = arrhcfwin+lnwin
        qfor(ire) = EXP(lnrate)
      ENDDO ! i-1,numwim

!-------------
! Wall -> Gas
!-------------
! Use arrhcf(ire,1) as inverse of time constant to equilibration.
! The rate constant comes from the equilibrium constant
! with the wall and the definition of time constant to equilibration
! assuming 2 first order reactions for gas/wall partitioning.
! The time constant to equilibrium tau = 1/kgw + 1/kwg.
! The equilibrium constant
! Keq=kgw/kwg=(RT/Pvap)*(cmw)=(RT/Pvap)*(Cw/Mw*gamma)
! see paper, Matsunaga & Ziemann, AST, 2010, p887.
! R=0.0820578 (atm L K-1 mol-1),and 1E-3 m3/L
! Here, cmw is given as the third auxilliary information "woucf(3,*)"
! (input parameter)
! vapor pressure computed with the Nanoonal SAR.
! woucf(1,*)=Bp, woucf(2,*)=dB : but now we calculate psat elsewhere.

      DO i=1,numwou
        ire = idwou(i)

! vapor pressure, natural LOG
        lnpvap = LOG(psat(idrestoi(ire,1)))

        cmw = LOG(woucf(3,i))+lnweq

        !lnrate = arrhcf(ire,3)+lnpvap-lnmultifacgw-cmw+lnwin
        lnrate = arrhcf(ire,1)+lnpvap-lnmultifacgw-cmw+lnwin
        qfor(ire) = EXP(lnrate)
        IF(i.EQ.1) WRITE(52,*) idwin(i),qfor(idwin(i)), &
                              idwou(i),qfor(idwou(i))
      ENDDO ! i=1,numwou

!=======================================================================
      END SUBROUTINE wall_MatZie2010
!=======================================================================
