      SUBROUTINE mtdyn_rp
     &           (idrestoi,wmol,imtr,
!     &            ctotaerbox,cnv,
     &            cnv,
     $            sumcbox,Mp,Rpo,maerbox,dvi,
     &            ire,temp,Rp,qln)
!--------------------------------------------
! TEST VERSION FOR NEW ktd PARAMETERIZATION
!--------------------------------------------
      USE akparameter_module
      USE flags_module,ONLY: seedtyp_fg
      USE forcing_params_module,ONLY: rpfix_fg, rpfix
      USE module_data_gecko_main,ONLY: small
      IMPLICIT NONE

* INPUT
      INTEGER ire,imtr
      INTEGER idrestoi(maxre,mxleft)
      REAL    wmol(maxsp)
      REAL    cnv     ! conc nonv-volatile aerosol molecules (molec/cc)
!      REAL    ctotaerbox ! conc nvseed + soa (molec/cc)
      REAL    dvi     ! dimensionless diffusion volume of molecule
      REAL,INTENT(IN):: maerbox    ! total mass of soa (ug/m3)
      REAL    Mp      ! molecular mass of seed particle material (g/mol)
      REAL    Rpo     ! initial particle radius of seed (cm)
      REAL    sumcbox ! number density of atmosphere (=sumc(ibox))
      REAL    temp    ! temperature (K)

* OUTPUT
      REAL    Rp ! particle radius, for output
      REAL    qln

* LOCAL
      INTEGER nsp

      REAL    ci,Cp !,Rp
      REAL    Dg,Dg2,ktd,kti,ktot,invktot
      REAL    rhop ! seed density, g/cm3
!      REAL    Vol  ! particle volume, 
!      REAL    small
      REAL    FS, Kn, mfp, mAB

! AIR MOLECULAR CHARACTERISTICS
      REAL,PARAMETER :: Mair = 28.56   ! average molwt of air
      REAL,PARAMETER :: dvair = 18.02  ! dimensionless diffn vol: air molecules
! calculated as SUM(n(i)V(i)) where V(N2)=18.5, V(O2)=16.3, V(ar) =
! 16.2, V(CO2)=26.9. Abundances = 78.09%,20.95%,0.93%,0.04%.
! Reference: Polin, Prausnitz & O'Connell (2001), Table 11-1: 5th edition of 
! Reid, Prausnitz & Poling (1987), The properties of gases and liquids,
! McGraw Hill, New York 
      REAL,PARAMETER :: Sair = 4.0E-19 ! Collision cross-section, m2
! e.g. Benyoucef & Tahri, Can J. Phys. 95, 346-352 (2017)

! AEROSOL SPECIES CHARACTERISTICS
      REAL,PARAMETER :: alpha = 1.     ! accommodation coefficient
      REAL,PARAMETER :: rhoaer = 1.3   ! SOA density, g/cm3  assumption !
!      REAL,PARAMETER :: Rpo  = 1.25E-5 ! initial particle radius (cm) = 125 nm
!                                       ! Rpo is now passed in as argument

! FUNDAMENTAL CONSTANTS
      REAL,PARAMETER :: Natm = 2.46E+19 ! molec cm-3 at 1 atm, 298K
      REAL,PARAMETER :: Navo = 6.02214E+23 ! molec.mol-1
      REAL,PARAMETER :: Pi   = 3.14159
      REAL,PARAMETER :: Rkg  = 8.31446  ! kg.m2.s-2.K-1.mol-1

! initialization
      !small=TINY(1.0)
      ci=0.
      Cp=0.
      Dg=0.
      Rp=0.
      kti=0.
      ktd=0.
      ktot=0.
      invktot=0.

! if cnv = 0, we would expect NO condensation.
      IF(cnv.LE.small) RETURN  ! no particles

      nsp = idrestoi(ire,1)

! density, g/cm3 
      IF(seedtyp_fg.EQ.1)THEN ! inorganic seed
        rhop = 1.77
      ELSE  ! organic seed 
        rhop = 0.9
      ENDIF ! mixed seed not considered

!----------------------------------------------------
!k_in from diffusion limited mass transfer equation : 
!----------------------------------------------------
!k_in = (12*Mp*Dg(i)*Caer/(N*rhop*Dp**2) where:
!     Mp = particle mean molecular mass, 
!     Dg = the diffusion coefficient,
!     N = Avogadro number (6.02214E23)
!     rhop = particle density, 
!     rpo = initial particle radius,
!     rp = particle radius (calculated below),
!     Rkg = 8.314 (gas constant, units of Kg m2 s-2 K-1 mol-1),
!           factor of 0.001 converts g/kg

!     c(i) = average velocity (cm.s-1):
!     c(i)=(8*R*T/(Pi*wmol(i)))**0.5
! 100 converts ci from m s-1 to cm s-1 
! 0.001 converts wmol from g to kg

      ci = 100. * ((8.*Rkg*temp)/(Pi*wmol(nsp)*0.001))**0.5

! ====Dg: original equation==============================
!     Dg = diffusion coefficient, depends on wmol(i) : cm2 molec-1 s-1
!Dg(i) = (1/sumcbox)*(1/Sair)*c(i)*(1/(1+(wmol(i)/Mair)))**0.5 where:
!     Sair = cross section in air (m2): 1e4 converts Sair from m2 to cm2

      Dg = (1./3.)*(1./sumcbox)*(ci/(Sair*1.e4))*
     &     (1./(1.+(wmol(nsp)/Mair)))**0.5

      !IF(nsp.EQ.1821)print*,"Dg1 = ",Dg

! ====Dg: equation of Fuller(1966)=============================
! for diffusion of trace gas A (=species) in bath gas B (= air)

! Reduced mass m(A,B) = 1/(1/mw(A)+1/mw(B) : g.mol-1
! Tang et al (2015) eqn 8 , doi:10.5194/acp-15-5585-2015 has '2' as
! numerator BUT Seinfeld & Pandis show '1', as does the Encyclopaedia
! of Physics (Lerner & Trigg, 1991).
      mAB = 1. / (1./wmol(nsp) + 1./Mair)  

! dvi (aka Va) is calculated by generator, output as ascii list.
! dimensionless diffusion volume of air dvair (aka Vb) declared above
! Reid et al 1987, cited by Tang et al. (2015) eqn 9

! Diffusion coefficient Dg = Dp(X).P  [Tang, eqn 6]
! Dg = pressure-independent diffusion coeff (diffusivity): Torr cm2 s-1
!    = 1.0868 T^1.75 /[sqrt(mAB)*(Va^1/3 + Vb^1/3)^2] : Torr.cm2.s-1
! (comment updated 201207 JMLT to correctly include denominator square)
! Torr-to-atm factor = 1/760 .........................> cm2.s-1
! Fuller at al (1966) cited by Tang et al (2015) eqn 7 

      Dg2 = 1.0868*temp**(7./4.) / (760.*sumcbox/Natm) /
     &  ( mAB**0.5 *(dvi**(1./3.) + dvair**(1./3.))**2.) 
      
! =============================================================
! compute particle number Cp (this version changes with ctotaerbox)
c      Cp= (Mp*ctotaerbox)/(Navo*(4/3)*Pi*rhop*(Rpo**3))

! compute particle number Cp (does not change)
        Cp = (Mp*cnv)/(Navo*(4./3.)*Pi*rhop*(Rpo**3.))

! compute particle radius Rp (varies with maerbox from previous timestep)
! JMLT, 180831: 

! Radius Rp = (3/4*vol/pi)**(1/3) units = cm
! Vol       = particle volume = Vol(seed) + Vol(SOA) (cm3/particle)
! Vol(seed) =    Mp    *     cnv     /     Navo     /  rhop   /   Cp
!   units   = (g/mole) * (molec/cm3) / (molec/mole) / (g/cm3) / (particles/cm3)
! Cp        = # aerosol particles (particles/cm3)
! Mp        = seed molwt =>
! Mp*cnv/Navo = mass of seed (g/cm3)
! rhop      = density of seed (g/cm3)
! Vol(SOA)  =   maerbox  *     1e-12       / rhoaer  /     Cp
!   units   = (ug/m3) * (g/cm3/(ug/m3)) / (g/cm3) / (particles/cm3)
! maerbox      = mass of SOA (ug/m3)
! 1e-12     = conversion from ug/m3 to g/cm3
! rhoaer    = density of SOA (g/cm3)
        IF (rpfix_fg == 0) THEN
          Rp = ( 3./4. *(maerbox*1.e-12/rhoaer + Mp*cnv/Navo/rhop)
     &            / Cp / Pi )**(1./3.)
        ELSE
          Rp = rpfix
        ENDIF

!      PRINT*,"Rp,Rpo,maerbox",Rp,Rpo,maerbox

      IF(Rp.lt.small) RETURN ! zero radius -> transfer rates = 0.
                               ! (can occur in 1st timestep)

! approximate ktd

      ktd = 4.*Pi*Dg*Cp*Rp

!----------------------------------------------------
!k_in from gas/particle collision (used for approximate method only): 
!----------------------------------------------------
!k_in = (1/4)*gamma*c(i)*Caer where: (That's not what's written!)

!     c(i) = average velocity (m.s-1)(calculated above)

       kti = (alpha*ci*1.e2*Pi*(Rp**2)*Cp)

!----------------------------------------------------
! qln depends on whether at a limit, or a combination of both processes

      IF (imtr.EQ.1) THEN
c        WRITE(6,*) 'mass transfer limited by diffusion'
        qln = log(ktd)
      ELSE IF (imtr.EQ.2) THEN
c        WRITE(6,*) 'mass transfer described by gas/particles collision'
        qln = log(kti)
      ELSE
!        WRITE(6,*) 'mass transfer approach'
        invktot = (1./ktd)+(1./kti)
        ktot = 1./invktot
        qln = log(ktot) 
      ENDIF

!----------------------------------------------------------
! To calculate ktot including FS (Seinfeld & Pandis, Ch12),
! Comment out following line
!      GOTO 400
!----------------------------------------------------------
! mfp = mean free path : cm.molec-1
! Kn = Knudson number : molec-1
! Fuchs-Sutugin factor (1971) : no units

      mfp = 3.*Dg2/ci
      Kn = mfp / Rp
      FS = (1.+Kn)/
     & (1. + 0.377*Kn + (4./3.)*(Kn + Kn**2.)/alpha )

! Gas-aerosol transfer rate ktot:
! Kulmala M, Wagner PE (2001) Mass accommodation and uptake coefficients
! — a quantitative comparison. Journal of Aerosol Science 32(7):833-841.
! Seinfeld JH, Pandis SN (2006) Atmospheric chemistry and physics:
! from air pollution to climate change.  (John Wiley & Sons, Inc.,
! Hoboken, New Jersey), 2nd Ed, pp 537-577.
      ktot = 4.*Pi*Dg2*Cp*Rp*FS
      qln = log(ktot) 

! 400  CONTINUE

!------------------------------------------------
! comment this statement out to get diagnostic output
      RETURN

!------------------------------------------------
! DIAGNOSTICS !

!      PRINT*,"Cp,maerbox,Rp,ktot:",Cp,maerbox,Rp,ktot

      !IF(ire.EQ.98.OR.ire.eq.97)THEN ! C12t dodecanol
!      WRITE(45,*) ire,nsp,wmol(nsp),dvi
!      WRITE(43,*) Rpo,Rp,Dg,Dg2,Cp,FS,qln
!      WRITE(41,*) temp,sumcbox,Natm,mAB,dvi,dvair
!      ENDIF

! K01000 in 10-11-diemthyleicosane_3g
!      IF(ire.EQ.1408) THEN
!        WRITE(43,*) nsp,wmol(nsp),ktd,kti,FS,ktot,Rp,maerbox
!      ENDIF
! N02008 in 10-11-diemthyleicosane_3g
!      IF(ire.EQ.2206) THEN
!        WRITE(45,*) nsp,wmol(nsp),ktd,kti,FS,ktot,Rp,maerbox
!      ENDIF

! END DIAGNOSTICS !
      END
