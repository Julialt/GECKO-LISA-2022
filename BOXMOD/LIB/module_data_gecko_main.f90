      MODULE module_data_gecko_main

      USE akparameter_module 
      USE inorganic_aer_module, only : aer_population, update_deposition
      USE prodloss_module
      USE forcing_params_module
      USE OFR_params_module

      IMPLICIT NONE

! species concentrations
      REAL     conc(maxsp,mbox)
      REAL     cbot_sav(maxsp) ! bottom box saved for top box exchange
      REAL     cbg(maxsp)      ! background

! variables linked to the chemical scheme
      CHARACTER(maxlsp) chrsp(maxsp)

      INTEGER  numsp, numre, num_n
      INTEGER  num_m, numfo, numhv, numcvar, numextra
      INTEGER  mx12stoi
      INTEGER  nauxpar(maxaux)
      INTEGER  numstoi(maxre,2)
      INTEGER  idrestoi(maxre,mxleft)
      INTEGER  idpdstoi(maxre,mxright)
      INTEGER  id_m(max_m)
      INTEGER  idfo(maxfo,3)
      INTEGER  idhv(maxhv), idcvar(maxcvar)
      INTEGER  idextra(maxextra)
      INTEGER  nself,idselfreac(mself,2)
      INTEGER  numo2,ido2(maxo2)
      INTEGER  numiso,idiso(maxiso)
      INTEGER  nummeo2,idmeo2(mxrpero)
      INTEGER  numreacro2(maxro2),idreacro2(mxrpero,maxro2)
      INTEGER  numreacdimer(maxdimer),idreacdimer(mxrdimer,maxdimer)

      INTEGER  chromotopcf(mtopchromo),numtopchromo(mtopchromo)
      INTEGER  idtopchromo(mtopchromo,msptopchromo)
      INTEGER  chromomedcf(mmedchromo),nummedchromo(mmedchromo)
      INTEGER  idmedchromo(mmedchromo,mspmedchromo)
      INTEGER  ntmedchromo
      INTEGER  nt1chromo,chromo1cf(mchromo),id1chromo(mchromo)
      INTEGER  nchrom,idchrom(mchromo)

      REAL     restoicf(maxre,mxleft),pdstoicf(maxre,mxright)
      REAL     arrhcf(maxre,3)
      REAL     focf(maxaux+3,maxfo)
      REAL     hvcf(maxhv),hvfact(maxhv),cvarcf(maxcvar)
      DOUBLE PRECISION     extracf(maxaux,maxextra)
      REAL     isocf(maxaux,maxiso)
      REAL     wmol(maxsp)

      INTEGER  nclro2, numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
      REAL     cro2(maxro2)
      INTEGER  ncldimer, numchemdimer(maxdimer)
      INTEGER  idchemdimer(mxrdimer,maxdimer)
      REAL     cdimer(maxdimer)

      REAL     qfor(maxre),rate(maxre)

! photolysis input data, for re-output
! NB: hardwiring max sizes of nsza and njtab
      INTEGER,PARAMETER:: mxjtab=150,mxsza=15,llin=76
      INTEGER  njtab,nsza,idjtab(mxjtab)
      REAL     jsza(mxsza),jvref(mxsza,mxjtab)
      CHARACTER(maxcoe)  jnam(mxjtab)
      CHARACTER(llin)  jreac(mxjtab)

! photolysis frequencies data
      INTEGER  numtet
      REAL     xang(maxang)
      REAL     rat1pho(mchromo,maxang),coef1pho(mchromo,nlo)
      REAL     ratmedpho(mmedchromo,maxang),coefmedpho(mmedchromo,nlo)
      REAL     rattoppho(mtopchromo,maxang),coeftoppho(mtopchromo,nlo)
! photolysis rates at specific times
      REAL     temprat1pho(mchromo,maxang)
      REAL     tempratmedpho(mmedchromo,maxang)
      REAL     temprattoppho(mtopchromo,maxang)

! stoichiometric coefficient as a function of T
      INTEGER  ntype(maxcvar), numcoe(maxcvar), nopc(maxcvar)
      INTEGER  ndatopc(maxcvar,mopc), nposopc(maxcvar,mopc,mpos)
      REAL     valcoe(nset,maxcoe,maxcvar)

! physical process rates
      REAL     rem(maxsp),rdep(maxsp)
      REAL     rex(maxsp),rdil,vmix(mbox)

! emission and deposition data relative to land use (urban, forest ...)
      INTEGER  nshcur, nenott, idnh3, idisop, idapin,idbpin
      INTEGER  idehcur(maxem),idenott(maxem)
      REAL     ecour,enoxur,ehcur
      REAL     ehcdatur(maxem,2),senoxur(3)
      REAL     enott(maxem,msur)
      REAL     enh3(msur),eisop(msur),eapin(msur),ebpin(msur)
      REAL     cscoef3(4,25)
      REAL     rik(msur),rlu0k(msur),rack(msur),rgsSk(msur)
      REAL     rgsOk(msur),rclSk(msur),rclOk(msur),z0k(msur)
      REAL     Vd(mxdep)
      TYPE(surface_data)  :: surface_emi(msur)
      REAL     isop_fac, mterp_fac

! SOA 
      INTEGER   nsat,ndim    
! pvap species indices:
! idsat = index in mech file of pvap species (i)
! satid = index in pvap file of mech species (i) 
! difid = index in diffusion volume file of mech species (i) 
      INTEGER,DIMENSION(mxsat):: idsat
      INTEGER,DIMENSION(maxsp):: satid,difid
      CHARACTER(maxlsp),DIMENSION(mxsat) :: namsat,namdif
! aerosol input params
      REAL      simpgroup(mxsat,31),bk(31,4)
      REAL      Tb(mxsat),HBN(mxsat),tau(mxsat),dB(mxsat)
      REAL      difvol(mxsat)   ! dyn readin: (1:nsat)
! aerosol derived quantities 
      REAL      psat(maxsp)     ! equilib uses (1:nsat), dyn needs (1:numsp)
      REAL      caer(mxsat)     ! equilib derived: (1:nsat)
      REAL      Cstar298(maxsp) ! dyn derived: (1:numsp)
      REAL      dvsp(maxsp)     ! dyn rearranged: (1:numsp)
      REAL      maer(mbox)
      REAL      ctotaer(mbox)
      TYPE(aer_population) :: inorg_aer(mbox)
! aerosol constants
      REAL      multifac

! chemical parameters for ground deposition
      CHARACTER(maxlsp) depnamspe(mxdep)
      INTEGER  ndepspe, iddepspe(mxdep)
      REAL     depdatspe(mxdep,3)

! mass transfer parameters
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      REAL     aoucf(2,maxt),woucf(3,maxt),wincf(3,maxt)
      REAL     ratfac, winfac,weqfac
      INTEGER  idgsat(mxsat),idasat(mxsat),idwsat(mxsat),iddsat(mxsat)
      INTEGER  imtr

      INTEGER  i,j,k,ibox,isp
      INTEGER  idno,idno2,idno3,idch3o2, ire
      INTEGER  idh2o,ido2dic,ido3,idho,idho2
      INTEGER  idSumRO2,idSuRCO3
      REAL     sumnox

      REAL     cmeo2
!      REAL     dilu
!      REAL(KIND=4)  :: tdat(2)

! aqueous phase

! deals with gas<->liq mass transfer
!      REAL tralphain(maxtr),tralphaout(maxtr)
!      REAL trdeltahin(maxtr),trdeltahout(maxtr)
!      REAL trhenryin(maxtr),trhenryout(maxtr)
!      INTEGER numtrin,numtrout
!      INTEGER idtrin(maxtr),idtrout(maxtr)

! hydration
!      REAL hydcf(maxhyd)
!      INTEGER numhyd,idhyd(maxhyd)
! acid dissociation
!      REAL ka(maxacid)
!      INTEGER numacid,idacid(maxacid)

! Aq ; aqueous phase oxidation
!      REAL ohaq(mxohaq,3)
!      INTEGER  numohaq,idohaq(mxohaq)
!      INTEGER idreacro2_aq(mxrpero,maxro2)
!      INTEGER nrpero_aq(mxrpero)

! arrays for tracer production/loss rate tracking
! total production/loss rate of each tracer for current timestep
      INTEGER  ntr
      INTEGER  idtr(mtr)
      REAL     trprod(maxsp),trloss(maxsp)
      REAL     puits(45),source(45),kicovi
      INTEGER  nbsource,nbpuits
      REAL     sum_source,sum_puits

      CHARACTER(18) filename

! array for indexing production and loss reaction of each species
! type spec_reac_map is defined in prodloss.h
      TYPE(spec_reac_map),TARGET   :: lpmap(maxsp) ! loss_production map
      INTEGER               :: ispec, nprod, nloss

! parameter for smallest number
! version previously used in:
! mtdyn_rp, mtrat, spakkrat6, spdeosition3, wrtopaer_ncdf
! DO NOT USE UNIVERSALLY: IT CRASHES iter4.f
      !REAL,PARAMETER :: small = TINY(1.0)
! version previously used in: iter4
      REAL,PARAMETER :: small = 1e-30
! version previously used in: spakkextra4, spakkrat6, spinterp5 
      !REAL,PARAMETER :: small = 1e-32
! version previously used in: vbs_module.f90
      !REAL,PARAMETER :: small = 1e-35

      END MODULE module_data_gecko_main
!======================================================================
      MODULE io_units_module
!     output file numbers.
      IMPLICIT NONE

! NB: units otherwise in use = {7,12,19,21,29,32,39,43,44,50,99)
      INTEGER,PARAMETER :: lpbl=10, lout=11, lppf=13, lpff=14, &
                           lsoa=15, ljval=16,lppa=17, lpaa=18, &
                           lread=20,lro2=38, ljall=22

      END MODULE io_units_module
!======================================================================
      MODULE NetCDF_vars_module
!     variables for NetCDF file handling
      USE akparameter_module,ONLY: maxsp,maxro2,mxsat
      IMPLICIT NONE

!--File IDs (NetCDF)
      INTEGER :: ncid_in,ncid_out,ncid_prev
!--Attribute text
      CHARACTER(len= 8) :: attname
      CHARACTER(len=40) :: line1,line2
! minimum concentration value to be output for Netcdf o/p file
      REAL,PARAMETER :: maskval=0.0
!      REAL::  tmpconc(maxsp)        ! only concs >= maskval
!      REAL::  tmpcro2(maxro2)       ! only concs >= maskval
!      REAL::  tmpcaer(mxsat)        ! only concs >= maskval
! maximum array size value to be output for Netcdf o/p file
      REAL,PARAMETER :: max_outvar_size = 5.252E+8

      END MODULE NetCDF_vars_module
!======================================================================
      MODULE printphoto_module
! arrays for photolysis rates printing
      USE akparameter_module,ONLY: maxreac_char
      IMPLICIT NONE

      INTEGER :: n_printphoto
      INTEGER,DIMENSION(:),ALLOCATABLE :: idprintphoto
      REAL,DIMENSION(:),ALLOCATABLE    :: photorates
      CHARACTER(LEN=maxreac_char),ALLOCATABLE :: photoreac_names(:)
      CHARACTER(LEN=maxreac_char),EXTERNAL :: printreaction
      CHARACTER(LEN=:),ALLOCATABLE :: cnjv ! string = # of j-values

      END MODULE printphoto_module
!======================================================================
      MODULE fundamental_consts_module
      IMPLICIT NONE

      REAL,PARAMETER :: Navo = 6.02214E+23 ! mol-1
      !REAL,PARAMETER :: multimass = 1.66E-12 ! = 1E12/6.02214E23
      REAL,PARAMETER :: multimass = 1.E+12 / Navo
      REAL,PARAMETER :: Ratm = 0.0820578 ! atm L K-1 mol-1
      REAL,PARAMETER :: Natm = 2.46E19 ! molec cm-3 at  atm, 298K
      REAL,PARAMETER :: PI = 3.14159265

      END MODULE fundamental_consts_module
!======================================================================
      MODULE solver_params_module
!     parameters for the solver
      IMPLICIT NONE

      INTEGER :: numit
      REAL    :: rtol, atol
      REAL    :: dtmin,dtmax

      END MODULE solver_params_module
!======================================================================
      MODULE steadystate_vars_module
! for printing end concentrations in namelist readable by other simulations
      USE akparameter_module,ONLY: maxlsp
      IMPLICIT NONE

      CHARACTER(maxlsp), ALLOCATABLE :: initspecnames(:)
      REAL, ALLOCATABLE        :: initspecconc(:,:)
      INTEGER                  :: non_zero_species_number

      END MODULE steadystate_vars_module
!======================================================================
      MODULE time_mgmt_module
!     parameters for the simulation (start, end time; timestep)
      IMPLICIT NONE

      INTEGER :: ntstep,ntprint
      INTEGER :: ntout,itout    
! datapoint in previous output file where new calculation starts
      INTEGER :: idat
! itout replaces nosav 
      INTEGER :: iskip, nskip
      REAL    :: tstart,tstop, delt
      REAL    :: time, time_save, tout
      REAL    :: timemod,timemod2
      REAL    :: tlocal
      REAL    :: etime, duree
      REAL(KIND=4)  :: tdat(2)
! flag to know if we're in daytime
      LOGICAL :: daytime_fg

      END MODULE time_mgmt_module
!======================================================================
