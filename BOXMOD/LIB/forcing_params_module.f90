      MODULE forcing_params_module

! module for environmental parameters as used in readkey and forcage

      USE constraints_module,ONLY: constrained_thing
      USE akparameter_module,ONLY: mhd,mbox,mtim,msur,maxem,maxsp, &
                                   species_data,surface_data, maxconst

      IMPLICIT NONE

! parameters for photolysis calculation (date, coordinate, ...)
      REAL     sla,slo,tz
      INTEGER  iy,im,id

! forcing parameters (mixing height, temperature, humidity, ...)
      INTEGER  nbox
      INTEGER  nhd,nmx,ntk,nrh,nws,npr
      INTEGER  ndil
      REAL     htop
      REAL     boxt(mhd),boxh(mhd)
      REAL     mixt(mhd),mixv(mbox,mhd)
      REAL     tempm(mbox),tempa(mbox),temptm(mbox)
      REAL     tktim(mtim),tkval(mbox,mtim),temp
      REAL     prtim(mtim),prval(mbox,mtim),prconst(mbox),pres
      REAL     rhm(mbox),rha(mbox),rhtm(mbox)
      REAL     rhtim(mtim),rhval(mbox,mtim),rh
      REAL     windm,winda,windtm
      REAL     wstim(mtim),wsval(mtim)
      REAL     vs  ! tropospheric subsidence velocity
      REAL     height,dhdt,sumc(mbox),water(mbox)
!      REAL     dhdt2 !,temptop
      REAL     dilconst
      REAL     diltim(mtim), dilval(mtim)
      REAL     sza

! flags for fixed value inputs
      INTEGER  waterfix,sumcfix, noxfix(mbox), presfix, dilfix
! parameter of the ground surface
      INTEGER  nsd, isurf
      REAL     psurf(mhd,msur),surft(mhd)
      REAL     xsurf(msur)
      REAL     ResAerdep ! surface resistances relevant for aerosol deposition
! parameter for ground deposition
      INTEGER  iseas
! thermodynamic
!      INTEGER  iscape

! constrained concentrations
      TYPE(species_data) :: cons_spec(mbox, maxconst)
! species emissions
      TYPE(species_data) :: emi_spec(maxem)
      REAL     eflux(maxsp)

! photolysis adjustment factors
      INTEGER  njf
      REAL     jftim(mtim),jfval(mtim),jfac
! flag and value for fixed SZA option
      INTEGER  szafix
      REAL     szaval
! aerosol surface aera
      REAL     saero
! non-volatile seed aerosol
      INTEGER  nseed
      REAL     cseed(mtim),tseed(mtim)
      REAL     cnv ! seed conc: enter in *.key under SEED
      REAL     gamm ! OA bulk activity coefficient: enter in *.key under SEED
      REAL     Mp  ! seed molwt: enter in *.key under SEED
! constrain aerosol radius
      integer  nrp, rpfix_fg
      REAL     Rpo! seed particle radius: enter in *.key under SEED
      REAL     Rp  ! aerosol particle radius (output)
      REAL     rptim(mtim), rpval(mbox, mtim)
      REAL     Rpfix 

! time constraints
      TYPE(constrained_thing) :: ph_const(mbox)
      TYPE(constrained_thing) :: sulfate_const(mbox)
      TYPE(constrained_thing) :: nitrate_const(mbox)
      TYPE(constrained_thing) :: kappa_const(mbox)
      TYPE(constrained_thing) :: naer_const(mbox)

      END MODULE forcing_params_module
!========================================================

      MODULE OFR_params_module
! forcing params for OFR Simulations

      IMPLICIT NONE

      REAL     jo2,jh2o,jo3,f185,f254
      REAL,PARAMETER :: sigmo2_185 = 1.1e-20
      REAL,PARAMETER :: sigmo3_254 = 1.03e-17
      REAL,PARAMETER :: sigmh2o_185 = 6.78e-20

      END MODULE OFR_params_module
