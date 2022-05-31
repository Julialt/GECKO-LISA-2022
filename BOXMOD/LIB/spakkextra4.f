      SUBROUTINE akkextra4
     &           (chrsp,numsp,idrestoi,extracf,wmol,Vd,
     &            ire,iex,temp,sumc,ibox,water,height,saero,qln,
     &            inorg_aer)

      USE flags_module,ONLY: isopsoa_fg,isopsoa_fac,vbs_fg,vbs_aging_fg
      USE time_mgmt_module, ONLY: daytime_fg
      USE akparameter_module
      USE inorganic_aer_module, ONLY : aer_population
      USE time_mgmt_module, ONLY: daytime_fg
      USE vbs_module, ONLY: kvbs, kvbs3, kvbs4, kvbs5
      USE module_data_gecko_main, ONLY: arrhcf, small
      IMPLICIT NONE

* INPUT
      CHARACTER(maxlsp) chrsp(maxsp)
      INTEGER numsp
      INTEGER ire,iex,ibox
      INTEGER idrestoi(maxre,mxleft)
      REAL    temp,height,saero
      REAL,DIMENSION(mbox) :: sumc,water
      REAL    extracf(maxaux,maxextra)
! USE REAL UNDER COMPILATION OPTION real-8
!      DOUBLE PRECISION    extracf(maxaux,maxextra)
      REAL    wmol(maxsp)
      REAL    Vd(maxsp)
      TYPE(aer_population), intent(in) :: inorg_aer(mbox)

* OUTPUT
      REAL    qln

* LOCAL
      INTEGER nsp,ino2
      !REAL    small
      REAL    xk,xk0,xk2,xk3,factor,gamma,xv,seff,q1,q2
      REAL    water_mono, water_dimer, kd, h2o, sumcb
      REAL Dg, alpha, Heff, khp, knuc, khso4
      REAL hp, so4mm, no3m, hso4m, kaq, R, wet_rad
      REAL y_lonox, y_hinox, y, A, B, C
      REAL tmp
      !small=1E-32
     
      sumcb = sumc(ibox)

* calculate water monomers vs dimers concentration
* dimerization constant [atm-1] according to Scribano et al., 2006
      kd = 4.7856e-4*exp(1851/temp-5.10485e-3*temp) ! atm^(-1)
!convert atm-1 to molec-1 cm3
      kd = kd*8.314*temp*1e6/(1.01325e5*6.02E23)

      h2o = water(ibox)

      water_dimer = (1+2*kd*h2o-sqrt(1+4*kd*h2o))/(2*kd)
      water_mono = h2o-water_dimer

**********************************************************
* identificateur=400 --> reaction multiplier par RO2     *
* qln en entree est la bonne valeur de k pour la         *
* reaction d'ordre 2. A multiplier par RO2 (common block)*
**********************************************************
c      IF (nint(extracf(1,iex)).EQ.400) THEN
c        qln=qln+log(cro2(1))
c      ELSE IF (nint(extracf(1,iex)).EQ.401) THEN
c        qln=qln+log(cro2(1))
c      ELSE IF (nint(extracf(1,iex)).EQ.402) THEN
c        qln=qln+log(cro2(2))
c      ELSE IF (nint(extracf(1,iex)).EQ.403) THEN
c        qln=qln+log(cro2(3))
c      ELSE IF (nint(extracf(1,iex)).EQ.404) THEN
c        qln=qln+log(cro2(4))
c      ELSE IF (nint(extracf(1,iex)).EQ.405) THEN
c        qln=qln+log(cro2(5))
c      ELSE IF (nint(extracf(1,iex)).EQ.406) THEN
c        qln=qln+log(cro2(6))
c      ELSE IF (nint(extracf(1,iex)).EQ.407) THEN
c        qln=qln+log(cro2(7))
c      ELSE IF (nint(extracf(1,iex)).EQ.408) THEN
c        qln=qln+log(cro2(8))
c      ELSE IF (nint(extracf(1,iex)).EQ.409) THEN
c        qln=qln+log(cro2(9))

**********************************************************
* identificateur=300 --> reaction avec O2                *
* qln en entree est la bonne valeur de k pour la         *
* reaction d'ordre 2. A multiplier par O2 (0.2*M)        *
**********************************************************
c      ELSE IF (nint(extracf(1,iex)).EQ.300) THEN
c        qln=qln+log(sumcb*0.2)

**********************************************************
* identificateur=100 --> 0 + O2 -> O3                    *
* qln en entree est la bonne valeur de k pour la         *
* recation d'ordre 3. A multiplier par M et O2 (0.2*M)   *
**********************************************************
      IF (nint(extracf(1,iex)).EQ.100) THEN
        qln=qln+log(sumcb)+log(sumcb*0.2)

**********************************************************
* identificateur=500 --> reaction avec H2O               *
* identificateur=501 --> 2 reaction :                    *
*                          1 avec H2O (q1)               *
*                          2 avec H2O et M (q2)          *
*            les parametres de 2 sont dans extracf       *
**********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.500) THEN
        IF(h2o.GT.0.0)THEN
          qln=qln+log(water_mono)
        ELSE
          qln=-9999
        ENDIF
      ELSE IF (nint(extracf(1,iex)).EQ.501) THEN
        IF(h2o.GT.0.0)THEN
          q1=exp(qln)*water_mono
          q2=extracf(2,iex)*EXP(-extracf(4,iex)/temp)*water_mono*sumcb
          qln=log(q1+q2)
        ELSE
          qln=-9999
        ENDIF
**********************************************************
* identificateur=502 --> reaction avec H2O dimer          *
      ELSE IF (nint(extracf(1,iex)).EQ.502) THEN
        IF(h2o.GT.0.0)THEN
          qln = qln + log(water_dimer)
        ELSE
          qln=-9999
        ENDIF

**********************************************************
* identificateur=550 --> reaction OH + HNO3              *
* k=k0 + k3M/(1+K3M/K2)                                  *
*    k0 : donnees par les parametres de Arrhenius (qln)  *
*    k2 : parametres 2,3,4 de extracf                    *
*    k3 : parametres 5,6,7 de extracf                    *
**********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.550) THEN
        xk0=exp(qln)
        xk2=extracf(2,iex)*(temp**extracf(3,iex))*
     &     EXP(-extracf(4,iex)/temp)
        xk3=extracf(5,iex)*(temp**extracf(6,iex))*
     &     EXP(-extracf(7,iex)/temp)
        xk=xk0 + xk3*sumcb / (1.+ ((xk3*sumcb)/xk2) )
        qln=log(xk)

*********************************************************
* reaction de formation heterogene de HONO              *
* identificateur=600                                    *
*    parametrisation selon Harisson et al., 1996        *
*    cf. these I.Bey pour plus de detail                *
*********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.600) THEN
        IF (ibox.EQ.1) THEN
          qln=qln-log(height)
        ELSE
          qln=log(small)
        ENDIF

*********************************************************
* reaction heterogene a la surface du sol.              *
*                       suppose pas de limitation       *
*                       par transport de masse en       *
*                       en phase gazeuse                *
*                       k=0.25*gamma*v*Aeff/h           *
*                       gamma est en log dans qln       *
*                       h est donne en entree de        *
*                       la routine                      *
*                       Aeff (cm2 effectif / cm2 au sol)*
*                       et donnee juste apres           *
*                       l'ID 650                        *
*                                                       *
* identificateur=650                                    *
*********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.650) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(11,*) '--error--, heterogeneous reaction with'
          write(11,*) 'code 650 is not a first order reaction'
          stop
        ENDIF
        IF (ibox.EQ.1) THEN
          CALL akspnum('NO2 ',chrsp,numsp,ino2)
          IF (ino2.eq.0) THEN
            WRITE(6,*) '--error--, in akkextra'
            WRITE(6,*) '           NO2 not found'
            STOP
          ENDIF
          nsp=idrestoi(ire,1)
          write(11,*) 'nsp=',nsp
          xv=( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
* passe de m/s en cm/s
          xv=xv*100.
          gamma=exp(qln)
          seff=extracf(2,iex)
          xk=0.25*xv*gamma*seff
          xk=xk/height
          write(11,'(A3,E10.3)') 'xk1=',xk
          qln=log(xk)
        ELSE
          qln=log(small)
        ENDIF

*********************************************************
* reaction heterogene a la surface du sol.              *
* UTILISATION DES VITESSES DE DEPOT                     *
*                       gamma est en log dans qln       *
*                       h est donne en entree de        *
*                       la routine                      *
*                                                       *
* identificateur=620                                    *
*********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.620) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(11,*) '--error--, heterogeneous reaction with'
          write(11,*) 'code 620 is not a first order reaction'
          stop
        ENDIF
        IF (ibox.EQ.1) THEN
          nsp=idrestoi(ire,1)
          write(11,*) 'nsp=',nsp
          IF (vd(nsp).gt.0) THEN
            write(11,*) 'Vd(nsp)=',Vd(nsp)
            factor=exp(qln)
            xk=(vd(nsp)/height)*factor
            write(11,'(A3,E10.3)') 'xk=',xk
            qln=log(xk)
          ELSE
            qln=log(small)
          ENDIF
        ELSE
          qln=log(small)
        ENDIF

*********************************************************
* reaction heterogene - suppose pas de limitation       *
*                       par transport de masse en       *
*                       en phase gazeuse                *
*                       k=0.25*gamma*v*S                *
*                       gamma est en log dans qln       *
*                       S est donne en entree           *
* identificateur=700                                    *
*********************************************************
      ELSE IF (nint(extracf(1,iex)).EQ.700) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 700 is not a first order reaction'
          stop
        ENDIF
        nsp=idrestoi(ire,1)
        xv=( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
* passe de m/s en cm/s
        xv=xv*100.
        gamma=exp(qln)
        xk=0.25*xv*gamma*saero
        qln=log(xk)
c        write(11,*) 'nsp=',nsp
c        write(11,*) 'xv=',xv
c        write(11,*) 'gamma=',gamma
c        write(11,*) 'xk=',xk
*********************************************************
* reactive uptake of isoprene products                  *
*  k = saero/(r/Dg + 4/(xv*gamma))                      *
*      saero = wet aerosol surface area (cm2/cm3)       *
*      r = wet aerosol radius (cm)                      *
*     Dg = gas diffusion coefficient (cm2/s)            *
*     xv = mean molecular speed (cm/s)                  *
*  gamma = uptake coefficient (dimensionless)           *
*  gamma = 1/(1/alpha + 3*xv/(4*r*R*T*Heff*kaq))        *
*   alpha = accomodation coefficient (dimensionless)    *
*       R = 0.08206 L atm / K / mol                     *
*       T = temperature (K)                             *
*    Heff = effective Henrys law coefficient (M/atm)    *
*     kaq = pseudo first order aqueous phase reaction   *
*  rate constant for conversion to non-volatile products*
*  (1/s)                                                *
*  kaq is calculated depending on product (Marais et al. 2016)*
*********************************************************

* identifier = 711 : IEPOX, MEPOX parameterisation from *
* make tetrol
* Marais et al. 2016
      ELSE IF (nint(extracf(1,iex)).EQ.711) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 711 is not a first order reaction'
          stop
        ENDIF

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
          nsp = idrestoi(ire,1)
! mean molecular speed
          xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
! passe de m/s en cm/s
          xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
         Dg=0.22*(18./wmol(nsp))**0.5
! alpha is known from second extra parameter 
          alpha = extracf(2, iex)
! Heff is known from third extra parameter
          Heff = extracf(3, iex)
! to calculate kaq, kH+ is known from parameters 4
! all in units M^-n-1 s-1
          khp = extracf(4, iex)
! hp = [H+], nuc = [SO4--] + [NO3-], hso4 = [HSO4-]
! all in M
          hp = inorg_aer(ibox)%hp
          kaq = khp*hp
          
          R = 0.08206
          wet_rad = inorg_aer(ibox)%wet_radius
          saero = inorg_aer(ibox)%saero
          gamma = 1/(1/alpha + 3*xv/(4*wet_rad*R*temp*Heff*kaq))
! for more acidic aerosol, gamma(MEPOX) assumed to be 30 times lower
! than IEPOX
! use scaling factor extracf(7) to account for this (1 for iepox, 1/30 for mepox)
          if (inorg_aer(ibox)%ph < 4) then
            gamma = gamma*extracf(5, iex)
          endif
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)
        endif
* identifier = 712 : IEPOX, MEPOX parameterisation from *
* make organosulfate reacting with HSO4-
* Marais et al. 2016
      ELSE IF (nint(extracf(1,iex)).EQ.712) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 712 is not a first order reaction'
          stop
        ENDIF

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else        
          nsp = idrestoi(ire,1)
! mean molecular speed
          xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
! passe de m/s en cm/s
          xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
         Dg=0.22*(18./wmol(nsp))**0.5
! alpha is known from second extra parameter 
          alpha = extracf(2, iex)
! Heff is known from third extra parameter
          Heff = extracf(3, iex)
! to calculate kaq, kkhso4- is known from parameters 4
! all in units M^-n-1 s-1
          khso4 = extracf(4, iex)
! hp = [H+], nuc = [SO4--] + [NO3-], hso4 = [HSO4-]
! all in M
          hp = inorg_aer(ibox)%hp
          hso4m = inorg_aer(ibox)%hso4m
          kaq = khso4*hso4m
          
          R = 0.08206
          wet_rad = inorg_aer(ibox)%wet_radius
          saero = inorg_aer(ibox)%saero
          gamma = 1/(1/alpha + 3*xv/(4*wet_rad*R*temp*Heff*kaq))
! for more acidic aerosol, gamma(MEPOX) assumed to be 30 times lower
! than IEPOX
! use scaling factor extracf(7) to account for this (1 for iepox, 1/30 for mepox)
          if (inorg_aer(ibox)%ph < 4) then
            gamma = gamma*extracf(5, iex)
          endif
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)        
        endif

* identifier = 713 : IEPOX, MEPOX parameterisation from *
* make organosulfate SO4--
* Marais et al. 2016
      ELSE IF (nint(extracf(1,iex)).EQ.713) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 713 is not a first order reaction'
          stop
        ENDIF

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
         nsp = idrestoi(ire,1)
! mean molecular speed
         xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
! passe de m/s en cm/s
          xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
          Dg=0.22*(18./wmol(nsp))**0.5
! alpha is known from second extra parameter 
         alpha = extracf(2, iex)
! Heff is known from third extra parameter
          Heff = extracf(3, iex)
! to calculate kaq, knuc, is known from parameters 4
! all in units M^-n-1 s-1
          knuc = extracf(4, iex)
! hp = [H+], nuc = [SO4--], hso4 = [HSO4-]
! all in M
          hp = inorg_aer(ibox)%hp
          so4mm = inorg_aer(ibox)%so4mm
          kaq = knuc*so4mm*hp
          
          R = 0.08206
          wet_rad = inorg_aer(ibox)%wet_radius
          saero = inorg_aer(ibox)%saero
          gamma = 1/(1/alpha + 3*xv/(4*wet_rad*R*temp*Heff*kaq))
! for more acidic aerosol, gamma(MEPOX) assumed to be 30 times lower
! than IEPOX
! use scaling factor extracf(7) to account for this (1 for iepox, 1/30 for mepox)
          if (inorg_aer(ibox)%ph < 4) then
            gamma = gamma*extracf(5, iex)
          endif
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)
        endif
* identifier = 714 : IEPOX, MEPOX parameterisation from *
* make organonitrates with NO3-
* Marais et al. 2016
      ELSE IF (nint(extracf(1,iex)).EQ.714) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 714 is not a first order reaction'
          stop
        ENDIF

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
         nsp = idrestoi(ire,1)
! mean molecular speed
         xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
! passe de m/s en cm/s
         xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
         Dg=0.22*(18./wmol(nsp))**0.5
! alpha is known from second extra parameter 
         alpha = extracf(2, iex)
! Heff is known from third extra parameter
          Heff = extracf(3, iex)
! to calculate kaq, knuc, is known from parameters 4
! all in units M^-n-1 s-1
         knuc = extracf(4, iex)
! hp = [H+], nuc = [NO3-], hso4 = [HSO4-]
! all in M
          hp = inorg_aer(ibox)%hp
          no3m = inorg_aer(ibox)%no3m
          kaq = knuc*no3m*hp
          
          R = 0.08206
          wet_rad = inorg_aer(ibox)%wet_radius
          saero = inorg_aer(ibox)%saero
          gamma = 1/(1/alpha + 3*xv/(4*wet_rad*R*temp*Heff*kaq))
! for more acidic aerosol, gamma(MEPOX) assumed to be 30 times lower
! than IEPOX
! use scaling factor extracf(7) to account for this (1 for iepox, 1/30 for mepox)
          if (inorg_aer(ibox)%ph < 4) then
            gamma = gamma*extracf(5, iex)
          endif
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)
        endif

* identifier = 702 : case when kaq is direcly constrained from mechanism
* used for ISOPN for instance
* Marais et al. 2016        
      ELSE IF (nint(extracf(1,iex)).EQ.702) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 702 is not a first order reaction'
          stop
        ENDIF        

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
         nsp = idrestoi(ire,1)
! mean molecular speed
          xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
* passe de m/s en cm/s
         xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
         Dg=0.22*(18./wmol(nsp))**0.5
! alpha is known from second extra parameter 
          alpha = extracf(2, iex)
! Heff is known from third extra parameter
          Heff = extracf(3, iex)
!  kaq is known from parameter 4
! in units s-1
          kaq = extracf(4, iex)
          R = 0.08206
          wet_rad = inorg_aer(ibox)%wet_radius
          saero = inorg_aer(ibox)%saero
          gamma = 1/(1/alpha + 3*xv/(4*wet_rad*R*temp*Heff*kaq))
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)
        endif

* identifier = 703 : case when gamma is direcly constrained from mechanism
* used for ISOPOOH for instance
* Marais et al. 2016        
      ELSE IF (nint(extracf(1,iex)).EQ.703) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 703 is not a first order reaction'
          stop
        ENDIF        

        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
          nsp = idrestoi(ire,1)
! mean molecular speed
          xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
* passe de m/s en cm/s
          xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2  
          Dg=0.22*(18./wmol(nsp))**0.5
! gamma is known from second parameter 
          gamma = extracf(2, iex)
          xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
          qln = log(xk)       
        endif
* identifier = 704 : case when gamma is direcly constrained from mechanism
* but is different between day and night
* used GLY and MGLY
* Marais et al. 2016        
      ELSE IF (nint(extracf(1,iex)).EQ.704) THEN
* check that the reaction is a first order reaction
        IF (idrestoi(ire,2).ne.0) THEN
          write(*,*) '--error--, heterogeneous reaction with'
          write(*,*) 'code 704 is not a first order reaction'
          stop
        ENDIF
        if (isopsoa_fg .eq. 0) then
          qln = -9999.
        else
          nsp = idrestoi(ire,1)
! mean molecular speed
          xv = ( (8.*8.32*temp)/(3.14159*wmol(nsp)*0.001) )**0.5
! passe de m/s en cm/s
          xv = xv*100.  
! gas diffusion, scaled from water
! using Dg(H2O) as a reference : Dg(i)=Dg(H2O)*(M(H2O)/M(i))^1/2
          Dg=0.22*(18./wmol(nsp))**0.5
! day gamma is second parameter, night gamma is third parameter
         IF (.not. daytime_fg) then
           gamma = extracf(3, iex)
         ELSE
           gamma = extracf(2, iex)
         ENDIF
         saero = inorg_aer(ibox)%saero
         xk = isopsoa_fac*saero/(wet_rad/Dg + 4/(xv*gamma))
         qln = log(xk)      
       endif
!treat special reactions 900 for VBS counter
      ELSE IF (nint(extracf(1,iex)).EQ.900) THEN
         if (vbs_fg .gt. 0) then
           A = exp(arrhcf(ire, 1))
           B = arrhcf(ire, 3)
           C = arrhcf(ire, 2)
           y_lonox = extracf(3, iex)
           y_hinox = extracf(2, iex)
           tmp = kvbs(A, B, C, temp, y_lonox, y_hinox)
           if (tmp .ne. 0.) then
             qln = log(tmp)
           else
             qln = -9999
           endif
         else
           qln = -9999
         endif
        
      ELSE IF (nint(extracf(1,iex)).eq.903) THEN
         if (vbs_fg .gt. 0) then
           A = exp(arrhcf(ire, 1))
           B = arrhcf(ire, 3)
           C = arrhcf(ire, 2)
           y = extracf(2, iex)
           tmp = kvbs3(A, B, C, temp, y) 
           if (tmp .ne. 0.) then
             qln = log(tmp)
           else
             qln = -9999.
           endif  
         else
           qln = -9999
         endif     
      
      ELSE IF (nint(extracf(1,iex)).eq.904) THEN
         if (vbs_fg .gt. 0 .and. vbs_aging_fg .gt. 0) then
           A = exp(arrhcf(ire, 1))
           tmp = kvbs4(A)
           if (tmp .ne. 0.) then
             qln = log(tmp)
           else
             qln = -9999
           endif
         else
           qln = -9999
         endif
      
      ELSE IF (nint(extracf(1,iex)).eq.905) THEN
         if (vbs_fg .gt. 0 .and. vbs_aging_fg .gt. 0) then
           A = exp(arrhcf(ire, 1))
           tmp = kvbs5(A)
           if (tmp .ne. 0.) then
             qln = log(tmp)
           else
             qln = -9999
           endif  
         else
           qln = -9999
         endif   
      
*********************************************************
* identificateur inconnu                                *
*********************************************************
      ELSE
        WRITE (*,'(a21)')'--error-- in akkextra'
        WRITE (*,'(a24,i4,a44)')
     &      '          the extra case',nint(extracf(1,iex)),
     &      'could not be treated'
        STOP
      ENDIF
      END
