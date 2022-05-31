! --------------------------------------------------
! compute rate constant for mass transfer "reaction"
! --------------------------------------------------
      SUBROUTINE mtrat(ratfac,winfac,weqfac,
     &                 numain,numaou,numwin,numwou,
     &                 idain,idaou,idwin,idwou,
     &                 temp,ctotaerbox,cnv,sumcbox,
     &                 gamm,Mp,Rpo,maerbox, psat,
     &                 wmol,imtr,idrestoi,
     &                 Rp,qfor)

      USE akparameter_module
      USE flags_module,ONLY: wall_fg,seedtyp_fg,icham
      USE solver_params_module,ONLY: dtmin
      USE module_data_gecko_main,ONLY: Cstar298,dvsp,chrsp,small
      USE fundamental_consts_module,ONLY: Ratm,Navo
      IMPLICIT NONE

! INPUT
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      REAL     temp
      REAL     ratfac,winfac,weqfac
      REAL     lnfac
      REAL     lnwin,lnweq
      REAL     wmol(maxsp)
      INTEGER  imtr
      INTEGER  idrestoi(maxre,mxleft)
      REAL     ctotaerbox,cnv,Mp,Rpo
      REAL     sumcbox  ! local value for sumcbox(ibox)
      REAL     gamm  ! bulk activity coeff for SOA 
      REAL,INTENT(IN)::  maerbox  
      REAL,INTENT(IN)::  psat(maxsp)

! OUTPUT
      REAL     Rp
      REAL     qfor(maxre)

! LOCAL
      INTEGER  ire, i
      REAL     qln 
      REAL     lnmultifacgp,lngamma,dvi
      REAL     lnpvap,lncaer,lnkeq,lnrate
      REAL     Cstar
      REAL     teq_approx  ! factor to allow equilibrium to be reached
                           ! on the timescale of dtmin. (Too-fast
                           ! equilibration leads to overshooting and
                           ! anomalous mass creation). Mathematical
                           ! solution from Andrew Conley, NCAR, Dec 2020.

! flag selector for specific chamber parameterization:
! 1 = Matsunaga & Ziemann 2010
! 2 = Krechmer et al 2016 (Ziemann chamber at CU)
! 3 = Krechmer et al 2016 / Li pers comm 2018 (Jimenez chamber at CU)
! 4 = Huang et al 2018 (universal parameterization)
!===============================================
! return if no mass transfer equations
      IF (numain+numwin.EQ.0) RETURN

* initialize
      !small= TINY(1.0)

      teq_approx=1.
      lnfac=LOG(ratfac)
      lnwin=LOG(winfac)
      lnweq=LOG(weqfac)
      write(44,*) ctotaerbox, cnv

! do gas/aerosol equilbirum only if aerosol already exists
      IF(ctotaerbox.GT.small) THEN

! --------------
! gas -> aerosol
! --------------

      DO i = 1,numain
        ire = idain(i)
        dvi = dvsp(idrestoi(ire,1))

       CALL mtdyn_rp
     &         (idrestoi,wmol,imtr,
!     &          ctotaerbox,cnv,
     &          cnv,
     $          sumcbox,Mp,Rpo,maerbox,dvi,
     &          ire,temp,Rp,qln)

        lnrate = qln + lnfac
        qfor(ire) = EXP(lnrate)
      ENDDO
      
! --------------
! aerosol -> gaz
! --------------

! reverse reaction comes from the equilibrium constant Keq=k_in/k_out
! where k_in is the forward reaction (see gas -> aerosol )
! The equilibrium constant Keq=Caer/Cgas (concentration in molecule/cm3 of air)
! is given by : Keq=(RT*Caer)/(N*Pvap) with R the gas constant,
!              N the Avogadro number and Pvap the vapor pressure (atm).
! R=0.0820578 (atm L K-1 mol-1), N=6.02214E23 (molec/mol) and 1000 cm3/L
! vapor pressure computed with the Nannoolal SAR.
! Aoucf(1,*) is the Bp, Aoucf(2,*) is dB for the species

      !lnmultifacgp = LOG(0.0820578*1000*temp/6.02214E23)
      lnmultifacgp = LOG(Ratm*1000.*temp/Navo)

      IF(seedtyp_fg.EQ.1)then ! inorganic seed
        lncaer=LOG(AMAX1((ctotaerbox-cnv),1.))
      ELSE                    ! organic seed
        lncaer=LOG(ctotaerbox)
      ENDIF

      DO i=1,numaou

        ire = idaou(i)
        dvi = dvsp(idrestoi(ire,1))

! compute the equilibrium constant by Pankow(1994)
        if (psat(idrestoi(ire,1)) == 0.) then
!          print *, 'error--'
!          print *, chrsp(idrestoi(ire,1)), ' has no psat'
          lnpvap = 30.
        else
          lnpvap = LOG(psat(idrestoi(ire,1)))
        endif
        lngamma = LOG(gamm)
        lnkeq = lnmultifacgp+lncaer-lnpvap-lngamma

! compute the rate constant for the transformation
        CALL mtdyn_rp
     &       (idrestoi,wmol,imtr,
!     &        ctotaerbox,cnv,
     &        cnv,
     &        sumcbox,Mp,Rpo,maerbox,dvi,
     &        ire,temp,Rp,qln)

        lnrate = qln+lnfac-lnkeq
        qfor(ire) = EXP(lnrate)

      ENDDO

! --------------
! apply fix for rates << dtmin to ensure equilibrium (not overshooting)
! ASSUME that for every AIN rxn, there is also an AOU rxn.
!      PRINT*,"mtrat gas-particle rates (raw):", qfor(110),qfor(111) 
!      DO i=1,numain
!        IF(1/(qfor(idain(i))+qfor(idaou(i))).LT.dtmin)THEN
!          teq_approx = 1/(dtmin*(qfor(idain(i))+qfor(idaou(i))))
!          qfor(idain(i)) = qfor(idain(i))*teq_approx
!          qfor(idaou(i)) = qfor(idaou(i))*teq_approx
!        ENDIF
!      ENDDO
!      PRINT*,"mtrat gas-part rates (revised):", qfor(110),qfor(111) 

! --------------
      ELSE ! no seed mass, no aerosol partitioning
        DO i=1,numain
          qfor(idain(i)) = 0.
        ENDDO
        DO i=1,numaou
          qfor(idaou(i)) = 0.
        ENDDO
      ENDIF

! --------------
! gas -> wall
!     and
! wall -> gas
! Are dealt with in chamber-specific subroutines, 
! that calculate qfor(ire) for all idwin and idwou reactions
! --------------

      IF(wall_fg .eq. 1) THEN
        SELECT CASE(icham) 
          CASE(1)  ! Matsunaga & Ziemann, 2010
            CALL wall_MatZie2010 
          CASE(2) ! Krechmer et al 2016, Ziemann chamber
            CALL wall_Krech2016_Z 
          CASE(3) ! Krechmer et al 2016, Jimenez chamber
            CALL wall_Krech2016_J 
          CASE(4) ! Huang et al. 2018, universal parameterization
            CALL wall_Huang2018
          CASE(99)  ! test routine
            CALL wall_test
          CASE DEFAULT ! no partitioning if chamber not specified
           DO i=1,numwin
             qfor(idwin(i)) = 0.
           ENDDO
           DO i=1,numwou
             qfor(idwou(i)) = 0.
           ENDDO
 
        END SELECT

! --------------
! apply fix for rates << dtmin to ensure equilibrium (not overshooting)
! ASSUME that for every AIN rxn, there is also an AOU rxn.
!        PRINT*,"mtrat gas-wall rates (raw):", qfor(112),qfor(113) 
!        DO i=1,numwin
!          IF(1/(qfor(idwin(i))+qfor(idwou(i))).LT.dtmin)THEN
!            teq_approx = 1/(dtmin*(qfor(idwin(i))+qfor(idwou(i))))
!            qfor(idwin(i)) = qfor(idwin(i))*teq_approx
!            qfor(idwou(i)) = qfor(idwou(i))*teq_approx
!          ENDIF
!        ENDDO
!        PRINT*,"mtrat gas-wall rates (revised):", qfor(112),qfor(113) 
! --------------

      ELSE ! no wall partitioning

        DO i=1,numwin
          qfor(idwin(i)) = 0.
        ENDDO
        DO i=1,numwou
          qfor(idwou(i)) = 0.
        ENDDO

      ENDIF

! comment this statement out to get diagnostic output
      RETURN

!------------------------------------------------
! DIAGNOSTICS !
! generate diagnostic output for reaction:

!      i=4 !  the ith "iain" reaction : K01000
!        ire = idain(i)

!        lnpvap = LOG(psat(idrestoi(ire,1)))
!        lnkeq = lnmultifacgp+lncaer-lnpvap
!        Cstar =  psat(idrestoi(ire,1))*wmol(idrestoi(ire,1))
!     &           /(Ratm*298) * 1.E+9

!        WRITE(51,*) 
!     & i,ire,qfor(idain(i)),qfor(idaou(i)),EXP(lnkeq),EXP(lnpvap),Cstar!,
!     &  qfor(idwin(i)),qfor(idwou(i)),ctotaerbox,Cstar
!        WRITE(51,*) 
!     &  i,ire,idain(i),idaou(i)


!      i=9 !  the ith "iain" reaction : N02008
!        ire = idain(i)

!        lnpvap = LOG(psat(idrestoi(ire,1)))
!        lnkeq = lnmultifacgp+lncaer-lnpvap
!        Cstar =  psat(idrestoi(ire,1))*wmol(idrestoi(ire,1))
!     &           /(Ratm*298) * 1.E+9

!        WRITE(56,*) 
!     &  i,ire,qfor(idain(i)),qfor(idaou(i)),EXP(lnkeq),EXP(lnpvap),Cstar!,
!     &  qfor(idwin(i)),qfor(idwou(i)),ctotaerbox,Cstar
! END DIAGNOSTICS !
!------------------------------------------------

      END
