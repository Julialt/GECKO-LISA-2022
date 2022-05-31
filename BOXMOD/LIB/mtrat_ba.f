! compute rate constant for mass transfer "reaction"
! (not dynamic)
! --------------------------------------------------
      SUBROUTINE mtrat_ba(ratfac,winfac,weqfac,
     &                 numain,numaou,numwin,numwou,
     &                 idain,idaou,idwin,idwou,aoucf,woucf,
     &                 arrhcf,temp,ctotaer,
     &                 qfor, wall_fg)
      USE akparameter_module
      IMPLICIT NONE

! INPUT
      INTEGER  wall_fg
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      REAL     aoucf(2,maxt),woucf(3,maxt)
      REAL     arrhcf(maxre,3)
      REAL     temp
      REAL     ctotaer
      REAL     ratfac,winfac,weqfac
      REAL     lnfac
      REAL     lnwin,lnweq
!      REAL     wmol(maxsp)

! OUTPUT
      REAL     qfor(maxre)

! LOCAL
      INTEGER  ire, i
      REAL     multifac,loga10,Trb,logpvap,lncaer,lnpvap,lnkeq,rate
      REAL     cmw


! return if no mass transfer equations
      IF (numain+numwin.EQ.0) RETURN

      lnfac=log(ratfac)
      lnwin=log(winfac)
      lnweq=log(weqfac)

! --------------
! gas -> aerosol
! --------------
! make the program later, use arrhcf(ire,1) as the forward rate constant
      DO i=1,numain
        ire=idain(i)
        rate=arrhcf(ire,1)+lnfac
        qfor(ire)=exp(rate)
      ENDDO

! --------------
! aerosol -> gaz
! --------------
! reverse reaction comes from the equilibrium constant Keq=k_in/k_out
! where k_in is the forward reaction (see gas -> aerosol )
! The equilibrium constant Keq=Caer/Cgas (concentration in molecule/cm3 of air)
! is given by : Keq=(RT*Caer)/(N*Pvap) with R the gas constant, N the Avogadro
! number and Pvap the vapor pressure.
! R=0.0820578 (atm L K-1 mol-1), N=6.02214E23 and 1000 cm3/L
! vapor pressure computed with the Nanoonal SAR.
! Aoucf(1,*) is the Bp, Aoucf(2,*) is dB for the species
      multifac = log(0.0820578*1000*temp/6.02214E23)
      loga10=log(10.)
      lncaer=log(ctotaer)
      DO i=1,numaou
        ire=idaou(i)
!        Trb=temp/Tb(i)
        Trb=temp/aoucf(1,i)
! compute the vapor pressure, log base 10 and switch to natural log
        logpvap = (4.1012+aoucf(2,i))*((Trb-1.)/(Trb-(1./8.)))
        lnpvap=logpvap*loga10

! compute the equilibrium constant
! statment to be used if aerosol concentration is computed OUTSIDE 2step - see iter3
        lnkeq=multifac+lncaer-lnpvap
! statment to be used if aerosol concentration is computed INSIDE 2step -see iter3
c        lnkeq=multifac-lnpvap

! compute the rate constant for the transformation
! (note that arrhcf(ire,1) is the rate constant for the forward reaction)
        rate=arrhcf(ire,1)-lnkeq +lnfac
        qfor(ire)=exp(rate)
      ENDDO
      
!      write(6,*) qfor(idain(1)),"<->",qfor(idaou(1))


! --------------
! gas -> wall
! --------------
! make the program later, use arrhcf(ire,1) as the forward rate constant
      if (wall_fg .eq. 0) return

      DO i=1,numwin
        ire=idwin(i)
        rate=arrhcf(ire,1)+lnwin
        qfor(ire)=exp(rate)
      ENDDO

! --------------
! wall -> gaz
! --------------
! reverse reaction comes from the equilibrium constant with the wall
! The equilibrium constant Keq=Cwall/Cgas=(RT/Pvap)*(cmw)=(RT/Pvap)*(Cw/Mw*gamma)
! see Ziemann paper, ast, 2010, p887).
! Here, cmw is given as the third auxilliary information "woucf(3,*)" (input parameter)
! R=0.0820578 (atm L K-1 mol-1),and 1E-3 m3/L
! vapor pressure computed with the Nanoonal SAR.
! woucf(1,*) is the Bp, woucf(2,*) is dB for the species
      multifac = log(0.0820578*1E-3*temp) ! RT in atm.m3.mol-1
      loga10=log(10.)
      DO i=1,numwou
        ire=idwou(i)
!        Trb=temp/Tb(i)
        Trb=temp/woucf(1,i)
!        logPvap(i) = (4.1012+dB(i))*((Trb-1.)/(Trb-(1./8.)))
! compute the vapor pressure, log base 10 and switch to natural log
        logpvap = (4.1012+woucf(2,i))*((Trb-1.)/(Trb-(1./8.)))
        lnpvap=logpvap*loga10

! compute the rate constant for the transformation
! (note that arrhcf(ire,1) is the rate constant for the forward reaction)
        cmw=log(woucf(3,i))+lnweq
        rate=arrhcf(ire,1)+lnpvap-multifac-cmw+lnwin
        qfor(ire)=exp(rate)
!        write(*,*) i,' ',ire,' ',qfor(ire)
      ENDDO

      END
