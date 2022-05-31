!**********************************************************************
! calcul de quelques forcages (hauteur de la boite, temperature,      *
! concentration de M, humidite relative et concentration en eau)      *
! new in version 4: calculation of emissions                          *
!                   and photolysis adjustment factor jfac             *
! new in version 5: - interpolation of constrained concentrations     *
!                   - ability to fix sumc,  in keyfile           *
!                                                                     *
! HAUTEUR DE BOITE : le calcul suppose l-existence de 2 boites au     *
!                    maximum. La hauteur du sommet de la boite 1      *
!                    varie en fonction du temps, le sommet de la      *
!                    boite 2 a une hauteur fixe                       *
! TEMPERATURE : Le calcul suppose un profil sinusoidal en             *
!   fonction du temps OR may be entered as a table of values          *
! HUMIDITE REALATIVE : Le calcul suppose un profil sinusoidal en      *
!   fonction du temps OR may be entered as a table of values          *
!     >>>>>>Tabular values OVERRIDE sine calculations<<<<<<           *
! Attention : le temps en entree est le temps modulo (24*3600)        *
! input variables entrees :                                           *
!     => mbox  : nombre maximum de boite dans le modele               *
!     => ibox : numero de la boite traitee                           *
!     => timemod  : temps local (modulo 24*3600)                         *
!     => emid : indices of emitted species                            *
!     => emtim : times for specified emissions                        *
!     => emval : values for specified emissions                       *
!     => conid : indices of constrained species                       *
!     => contim : times for constrained concentrations                *
!     => conval : values for constrained concentrations               *
!     => tempm(i): valeur moyenne de la temperature                   *
!                dans la boite i (en K)                               *
!     => tempa(i): amplitude de variation de la temperature           *
!                  autour de la valeur moyennne dans la boite (en K)  *
!     => temptm(i): fixe le temps pour le maximum de temp (en s)      *
!     => ntk   : number of tabular temperatures                       *
!     => tktim : time for tabular temperatures (s)                    *
!     => tkval : tabular temperatures (K)                             *
!     => rhm(i): valeur moyenne de l-humidite relative                *
!                dans la boite i (en %)                               *
!     => rha(i): amplitude de variation de l-humidite relative        *
!                autour de la valeur moyennne dans la boite (en %)    *
!     => rhtm(i): fixe le temps pour le maximum d-humidite (en s)     *
!     => nrh   : number of tabular RH values                          *
!     => rhtim : time for tabular RH values (s)                       *
!     => rhval : tabular RH values (%)                                *
!     => diltim : time for tabular dilution values (s)                *
!     => dilval : tabular dilution values (s-1)                       *
!     => sumcfix : flag for fixed concentration of M                  *
!     => waterfix : flag for fixed concentration of H2O               *
! input/output variables
!     => sumc   : conencentration de M                                *
!     => water   : conencentration de H2O                             *
!     => jfac   : factor to multiply photolysis rates to match obsns  *
!     => conc   : concentrations (constraints apply to box 1 only)    *
! output variables sorties                                            *
!      height, dhdt are now constrained in separate subroutine        *
!      pbl_height_forcage.f                                           *
!     => vmix   : mixing velocity across box boundary at time         *
!     => cnv    : background aerosol concentration at time            *
!     => temp   : temperature local                                   *
!     => rh     : humidite relative (en %)                            *
!     => eflux(i)  : emissions (molecules/cm2)                        *
! parametre de la routine                                             *
!     => fconv => (2*pi)/(24*3600)                                    *
!     => pi_2  => pi/2                                                *
!***********************************************************************
      SUBROUTINE forcage6(timein,conc,organic_mass)

      USE flags_module,ONLY: isopsoa_fg
      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: timemod
      USE akparameter_module
      USE forcing_params_module
      USE constraints_module,ONLY: constrained_thing, update_constraint
      USE inorganic_aer_module,ONLY: aer_population, init_aer_population
      USE module_data_gecko_main,ONLY: ibox,idno,idno2,idno3, &
                                       vmix,inorg_aer

      IMPLICIT NONE

! input
      REAL,INTENT(IN) :: timein
      REAL,INTENT(IN) :: organic_mass

! in/out
      REAL,INTENT(INOUT),DIMENSION(maxsp) :: conc

! local
      INTEGER i,j
      REAL    b,hmean,dxdt
      REAL    c1,c2,c3,TC,pvap_H2O,psat_H2O

      REAL,PARAMETER :: fconv=7.2722E-5
      REAL,PARAMETER :: pi_2=1.5708
      REAL,PARAMETER :: po=1.01325  ! bar
      REAL,PARAMETER :: ho=0.     ! m

!-------------------------------------------------------------------
! use real time ("timein") for emissions
! use time modulo 24h ("timemod") for physical contraints

      !PRINT*,"timein,timemod",timein,timemod

      !PRINT*,"! calculate the advective / diffusive mixing rate"
      IF(nmx.GT.0)THEN
! check time interval is valid
        IF(mixt(nmx).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -advection-'
            WRITE (lout,*) 'upper limit for time not found'
            STOP
        ENDIF
        DO i=1,nmx-1
          IF (mixt(i+1).GT.timemod) THEN
            dxdt=(mixv(ibox,i+1)-mixv(ibox,i))/(mixt(i+1)-mixt(i))
            vmix(ibox)=ABS(mixv(ibox,i) + dxdt*(timemod-mixt(i)))
            EXIT
          ENDIF
        ENDDO
      ENDIF

      !PRINT*,"! calculate the background aerosol"
      IF(nseed.GT.0)THEN
! check time interval is valid
        IF(tseed(nseed).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -seed-aerosol-'
            WRITE (lout,*) 'upper limit for time not found'
            STOP
        ENDIF
        DO i=1,nseed-1
          IF (tseed(i+1).GT.timemod) THEN
            dxdt=(cseed(i+1)-cseed(i))/(tseed(i+1)-tseed(i))
            cnv = cseed(i) + dxdt*(timemod-tseed(i))
            EXIT
          ENDIF
        ENDDO
      ENDIF

      !print*,"! calculate the dilution"
      IF(ndil.gt.0) then
        dilfix = 1
        if(diltim(ndil).lt.timemod) then
          write(lout,*) '--error--, in forcage6 -dilution-'
          write(lout,*) 'upper limit for time not found'
          stop
        endif
        do i=1,ndil
          if (diltim(i+1).gt.timemod) then
            dxdt=(dilval(i+1)-dilval(i))/(diltim(i+1)-diltim(i))
            dilconst=dilval(i) + dxdt*(timemod-diltim(i))
            exit
          endif
        enddo
      endif

!-JFAC------------------------------

      !PRINT*,"! calculate the photolysis adjustment factor"
      IF(njf.GT.0)THEN
! check time interval is valid
        IF(jftim(njf).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -photolysis-factor-'
            WRITE (lout,*) 'upper limit for time not found'
            WRITE (lout,*) '=> edit JFAC in user input file'
            STOP
        ENDIF
        DO i=1,njf-1
          IF (jftim(i+1).GE.timemod) THEN
            dxdt=(jfval(i+1)-jfval(i))/(jftim(i+1)-jftim(i))
            jfac=jfval(i) + dxdt*(timemod-jftim(i))
            EXIT
          ENDIF
        ENDDO
      ENDIF

!-TEMPERATURE------------------------------

      !PRINT*,"! calculate or interpolate the temperature"
      !PRINT*,"ibox = ",ibox
      !PRINT*,"ntk = ",ntk
      IF(ntk.GT.0)THEN
! check time interval is valid
        IF(tktim(ntk).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -temperature-'
            WRITE (lout,*) 'upper limit for time not found'
            STOP
        ENDIF
        DO i=1,ntk-1
          IF (tktim(i+1).GT.timemod) THEN
            dxdt=(tkval(ibox,i+1)-tkval(ibox,i))/(tktim(i+1)-tktim(i))
            temp=tkval(ibox,i) + dxdt*(timemod-tktim(i))
            EXIT
          ENDIF
        ENDDO
      ELSE
      !PRINT*,"temptm(ibox) = ",temptm(ibox)
      !PRINT*,"tempa(ibox) = ",tempa(ibox)
        b=pi_2-fconv*temptm(ibox)
        temp=tempm(ibox)+tempa(ibox)*SIN((fconv*timemod)+b)
      ENDIF

!-PRESSURE & SUMC------------------------------

      !PRINT*,"! calculate pressure (for sumc estimate)"
      IF (presfix .EQ. 1) THEN 

        pres = prconst(ibox)

      ELSE
        IF(npr.NE.0)THEN ! interpolate tabular pressure
! check time interval is valid
          IF(prtim(npr).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -PRTB-'
            WRITE (lout,*) 'upper limit for time not found'
            STOP
          ENDIF
          DO i=1,npr-1
            IF (prtim(i+1).GT.timemod) THEN
              dxdt=(prval(ibox,i+1)-prval(ibox,i))/(prtim(i+1)-prtim(i))
              pres=prval(ibox,i) + dxdt*(timemod-prtim(i))
              EXIT
            ENDIF
          ENDDO
        ELSE ! calculate pressure (bar) at box mid ht
! comment JMLT 2020
! barometric formula p ~ po.EXP(-((h-h0).g.M)/(Tb.R*)) where:
! p = pressure(bar) at altitude h(m)
! po = pressure at sealevel (h0) = 1.01325 (bar used here)
! g = Earth-surface gravitational acceleration =  9.80665 (m.s-2)
! M = Molar mass of dry air = 0.02896 (kg.mol-1)
! To = standard temperature = 288.16 (K)
! R* = Unversal gas constant = 8.31446 (J.mol-1.K-1) or (kg.m2.K-1.mol-1.s-2)

          IF(ibox == 1) THEN 
            hmean = height/2
          ELSE IF(ibox == 2) THEN
            hmean = (htop + height)/2
          ENDIF
! original equation uses hmean in km.
! and converts constant divisor to km units also
          hmean=hmean*1E-5
          pres=po*EXP(-(hmean-ho)/7.)
! JMLT 2020:calculate pressure in bar using factor 1000/((g.M)/(To.Ro))
!          pres=po*EXP(-hmean/8.436)
        ENDIF

      ENDIF

      !PRINT*,"! calculate [M] ( = sumc), unless NDEN specified in key file"
! units: sumc(molec.cm-3)=pres(Pa/1e5)*Navo(molec.cm-3/10)/
!                                      (R0(m3.pa.K-1.mol-1).T(K))
! conversion 1/(1e5*10) compensates for the cm-3/m3 conversion
      IF(sumcfix.NE.1) sumc(ibox)=(pres*6.022E+22)/(8.31446*temp)

!-RH-------------------------------

      !PRINT*,"! interpolate or calculate (sine function) RH"
      IF(nrh.GT.0)THEN
! check time interval is valid
        IF(rhtim(nrh).lt.timemod)THEN
            WRITE (lout,*) '--error--, in forcage6 -RH-'
            WRITE (lout,*) 'upper limit for time not found'
            STOP
        ENDIF
        DO i=1,nrh-1
          IF (rhtim(i+1).GT.timemod) THEN
            dxdt=(rhval(ibox,i+1)-rhval(ibox,i))/(rhtim(i+1)-rhtim(i))
            rh=rhval(ibox,i) + dxdt*(timemod-rhtim(i))
            EXIT
          ENDIF
        ENDDO
      ELSE
        b=pi_2-fconv*rhtm(ibox)
        rh=rhm(ibox)+rha(ibox)*SIN((fconv*timemod)+b)
      ENDIF

      ! constrain Rp particles radiu
      IF (nrp .GT. 0) THEN
        IF (rptim(nrp) .LT. timemod) THEN
          write(lout,*) '--error--, in forcage6 -RP-'
          write(lout,*) 'upper limit for time not found'
          STOP
        ENDIF
        DO i=1,nrp-1
          IF (rptim(i+1).GT.timemod) THEN
            dxdt=(rpval(ibox,i+1)-rpval(ibox,i))/(rptim(i+1)-rptim(i))
            rpfix=rpval(ibox,i) + dxdt*(timemod - rptim(i))
            EXIT
          ENDIF
        ENDDO
      ENDIF

      !PRINT*,"! calculate [H2O], unless specified in input.key file"
! RH = 100 x pvap(H2O)/psat(H2O)
! "Magnus formula": psat(H2O) = c1 * EXP(c2*T/(c3+T)) (units Pa)
!  where: T is expressed in degrees ! : T(K)-273.16
!         c1 = 610.94, c2 = 17.625, c3 = 243.04
! Ref: Alduchov & Eskridge 1996, J.Appl.Meteorol. 35 601-609.
! quoted by: Lawrence, 2005, BAMS DOI:10.1175/BAMS-86-2-225
!  pvap(Pa) = RT(K) ! n/V (units: n = moles, V = 1m3)
!  R = 8.314 J.K-1.mol-1
!  "water" unites = molec/cc -> conversion 1e6/6.022e23 to moles/m3
      c1 = 610.94
      c2 = 17.625
      c3 = 243.04
      TC = temp-273.16
      psat_H2O = c1*EXP(c2*TC/(c3+TC))

      !PRINT*,"waterfix = ",waterfix
      IF(waterfix.EQ.1)THEN
        !PRINT*,psat_H2O
        pvap_H2O = water(ibox) * (8.314*temp*1.e6/6.02214e23)
        rh = 100. * pvap_H2O / psat_H2O
      ELSE
        !PRINT*,"temp = ",temp
        pvap_H2O = psat_H2O * rh/100.
        water(ibox) = pvap_H2O / (8.314*temp*1.e6/6.02214e23)
      ENDIF ! (waterfix.EQ.1)

! ORIGINAL FORMULAE (use grammes, calories - rather opaque)
!      IF(waterfix.EQ.1)THEN
!        rh = water(ibox) /
!     &      ( 6.1078*EXP( -1.*(597.3-0.57*(temp-273.16)) * 18./1.986 *
!     &        (1./ temp-1./273.16)) * 10./(1.38E-16*(temp)) )
!      ELSE
!        water(ibox) = rh *
!     &      ( 6.1078*EXP(-1.*(597.3-0.57*(temp-273.16))*18./1.986*
!     &        (1./ temp-1./273.16))*10./(1.38E-16*(temp)) )
!!       water(ibox)=1.55E17  ! default !
!      ENDIF

      IF (isopsoa_fg .EQ. 1) then
        !PRINT*,"! interpolate aerosol things"
        CALL update_constraint(ph_const(ibox), timemod)
        CALL update_constraint(sulfate_const(ibox), timemod)
        CALL update_constraint(nitrate_const(ibox), timemod)
        CALL update_constraint(kappa_const(ibox), timemod)
        CALL update_constraint(naer_const(ibox), timemod)

        !PRINT*,"! update aerosol properties"
        CALL init_aer_population(inorg_aer(ibox),   &
                           1e3, & ! hardcoded aerosol density (kg/m-3)
                           kappa_const(ibox)%val,   &
                           naer_const(ibox)%val,    &
                           ph_const(ibox)%val,      &
                           sulfate_const(ibox)%val, &
                           nitrate_const(ibox)%val, &
                           organic_mass, & ! from model chemistry
                           rh ) 
      ENDIF

      CALL apply_constraints(lout, timemod, conc, temp, sumc(ibox), &
                       eflux, cons_spec, emi_spec, idno, idno2, idno3)


      END
