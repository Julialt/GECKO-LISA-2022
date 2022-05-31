      MODULE module_data_gecko_main

      USE module_data_gecko_kind, ONLY:  r8
      USE netcdf
      USE sorting, only:sort_species, srtid_chrsp
      USE constraint_module, only : constrained_thing
      USE inorganic_aer_module, only : aer_population
      !$ use OMP_LIB


      IMPLICIT NONE
      INCLUDE 'akparameter.h'
      INCLUDE 'flags.h'
      INCLUDE 'prodloss.h'

* species concentrations
      REAL     cbot(maxsp), cbot_sav(maxsp) ! bottom box
      REAL     ctop(maxsp) ! top box
      REAL     cbg(maxsp)       ! background

* variables linked to the chemical scheme
      CHARACTER*(maxlsp) chrsp(maxsp)

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
NTEGER  numiso,idiso(maxiso)
      INTEGER  nummeo2,idmeo2(mxrpero)
      INTEGER  numreacro2(maxro2),idreacro2(mxrpero,maxro2)
      INTEGER  numreacdimer(maxdimer),idreacdimer(mxrdimer,maxdimer)

      INTEGER  chromotopcf(mtopchromo),numtopchromo(mtopchromo)
      INTEGER  idtopchromo(mtopchromo,msptopchromo)

      INTEGER  chromomedcf(mmedchromo),nummedchromo(mmedchromo)
      INTEGER  idmedchromo(mmedchromo,mspmedchromo)
      INTEGER  ntmedchromo

      INTEGER  nt1chromo,chromo1cf(mchromo),id1chromo(mchromo)


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

* photolysis input data, for re-output
! NB: hardwiring max sizes of nsza and njtab
      INTEGER,PARAMETER:: mxjtab=150,mxsza=15,llin=76
      INTEGER  njtab,nsza,idjtab(mxjtab)
      REAL     jsza(mxsza),jvref(mxsza,mxjtab)
      CHARACTER*(maxcoe)  jnam(mxjtab)
      CHARACTER*(llin)  jreac(mxjtab)

* photolysis frequencies data
      INTEGER  numtet
      REAL     xang(maxang)
      REAL     rat1pho(mchromo,maxang),coef1pho(mchromo,nlo)
      REAL     ratmedpho(mmedchromo,maxang),coefmedpho(mmedchromo,nlo)
      REAL     rattoppho(mtopchromo,maxang),coeftoppho(mtopchromo,nlo)
* photolysis rates at specific times
      REAL     temprat1pho(mchromo,maxang)
      REAL     tempratmedpho(mmedchromo,maxang)
      REAL     temprattoppho(mtopchromo,maxang)

* stoichiometric coefficient as a function of T
      INTEGER  ntype(maxcvar), numcoe(maxcvar), nopc(maxcvar)
      INTEGER  ndatopc(maxcvar,mopc), nposopc(maxcvar,mopc,mpos)
      REAL     valcoe(nset,maxcoe,maxcvar)

* parameters for the simulation (starting time, date, coordinate, ...)
      INTEGER  nbox,ntstep,ntprint
      REAL     tstart,tstop
      REAL     sla,slo,tz
      INTEGER  iy,im,id


* solver parameters
      INTEGER  numit
      REAL     time, time2, time_save
      REAL     delt,deltprint,dtmin,dtmax
      REAL     tprint, tout
      REAL     rtol, atol
      REAL     rem(maxsp),rdep(maxsp)
      REAL     rex(maxsp),rdil,rmix

* forcing parameters (mixing height, temperature, humidity, ...)
      INTEGER  nflag
      INTEGER  nhd,nmx,ntk,nrh,nws
      REAL     timemod,timemod2
      REAL     htop
      REAL     boxt(mhd),boxh(mhd)
      REAL     mixt(mhd),mixr(mhd)
      REAL     tempm(mbox),tempa(mbox),temptm(mbox)
      REAL     tktim(mtim),tkval(mbox,mtim),temp
      REAL     rhm(mbox),rha(mbox),rhtm(mbox)
      REAL     rhtim(mtim),rhval(mbox,mtim),rh,rh2
      REAL     windm,winda,windtm
      REAL     wstim(mtim),wsval(mtim)
      REAL     vs  ! tropospheric subsidence velocity
      REAL     height,dhdt,tempbot,sumc(mbox),water(mbox)
      REAL     dhdt2,temptop
      REAL     pres, dilconst
! time constraints
      TYPE(constrained_thing) :: ph_const(mbox)
      TYPE(constrained_thing) :: sulfate_const(mbox)
      TYPE(constrained_thing) :: nitrate_const(mbox)
      TYPE(constrained_thing) :: kappa_const(mbox)
      TYPE(constrained_thing) :: naer_const(mbox)
* flags for fixed value inputs
      INTEGER  waterfix,sumcfix, noxfix(mbox), presfix, dilfix
* forcing params for OFR
      REAL     jo2,jh2o,jo3,f185,f254

* parameter of the ground surface
      INTEGER  nsd
      REAL     psurf(mhd,msur),surft(mhd)
      REAL     xsurf(msur)

* emission and deposition data relative to land use (urban, forest ...)
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
      REAL     isop_fac

* constrained concentrations
      TYPE(species_data) :: cons_spec(mbox, maxem)
* species emissions
      TYPE(species_data) :: emi_spec(maxem)
      REAL     eflux(maxsp)
* species emissions
!      INTEGER  emid(maxem)
!      REAL emtim(maxem,mtim),emval(maxem,mtim)
* constrained concentrations
!      INTEGER  conid(maxem),ncons
!      REAL     contim(maxem,mtim),conval(maxem,mtim)
* photolysis adjustment factors
      INTEGER  njf
      REAL     jftim(mtim),jfval(mtim),jfac
* flag and value for fixed SZA option
      INTEGER  szafix
      REAL     szaval

      REAL  sza
* thermodynamic
      INTEGER  iscape

* aerosol surface aera
      REAL     saero

* ouput file numbers. 
! NB: units otherwise in use = {7,12,19,21,29,32,39,43,44,50,99)
      INTEGER,PARAMETER :: lpbl=10, lout=11, lppf=13, lpff=14,
     &                     lsoa=15, ljval=16,lppa=17, lpaa=18

* SOA thermodynamic
      REAL     psat(mxsat),psattop(mxsat)
      CHARACTER*(maxlsp) namsat(mxsat)
      REAL               Tb(mxsat),HBN(mxsat),tau(mxsat),logPvap(mxsat)
      INTEGER            nsat,ndim,idsat(mxsat)
      REAL               caer(mxsat),cupt(mxsat),kup
      REAL               caertop(mxsat)
      REAL               maerbot, maertop
      REAL               multifac
      REAL,PARAMETER :: multimass = 1.66E-12 ! = 1E12/6.02214E23
      REAL               cnv ! eed conc: enter in *.key under SEED
      REAL               Mp  ! seed molwt: enter in *.key under SEED
      REAL               Rpo ! seed particle radius: enter in *.key
under SEED
      REAL               Rp  ! aerosol particle radius (output)
      REAL               dB(mxsat)
      REAL               simpgroup(mxsat,0:30),bk(0:30,4)
      TYPE(aer_population) :: inorg_aer(mbox)

* parameter for ground deposition
      INTEGER  iseas

* chemical parameter for ground deposition
      CHARACTER*(maxlsp) depnamspe(mxdep)
      INTEGER  ndepspe, iddepspe(mxdep)
      REAL     depdatspe(mxdep,3)

* non-volatile seed aerosol
      INTEGER            nseed
      REAL               cseed(mtim),tseed(mtim)
! mass transfer parameters
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      REAL     aincf(3,maxt),aoucf(3,maxt),woucf(3,maxt),wincf(3,maxt)
      REAL     ctotaer
      REAL     ratfac, winfac,weqfac
      INTEGER  idgsat(mxsat),idasat(mxsat),idwsat(mxsat),iddsat(mxsat)
      INTEGER  imtr

! for printing end concentrations in namelist readable by other
! simulations
      CHARACTER*(maxlsp), ALLOCATABLE :: initspecnames(:)
      REAL, ALLOCATABLE        :: initspecconcbot(:), initspecconctop(:)
      INTEGER                  :: non_zero_species_number

* PARAMETERS FOR OFR SIMULATION
      REAL,PARAMETER :: sigmo2_185 = 1.1e-20
      REAL,PARAMETER :: sigmo3_254 = 1.03e-17
      REAL,PARAMETER :: sigmh2o_185 = 6.78e-20

* local
      INTEGER  nosav
      INTEGER  i,j,k,isp,idno,idno2,idno3,idch3o2, ire
      INTEGER  idh2o,ido2dic,ido3,idho,idho2
      INTEGER  iskip, nskip
      REAL     cmeo2
      REAL     dilu
      REAL     etime, duree
      REAL(KIND=4)  :: tdat(2)
      REAL     sumnox,diff,rap,xno2,xno

* arrays for tracer production/loss rate tracking
* total production/loss rate of each tracer for current timestep
      INTEGER  ntr
      INTEGER  idtr(mtr)
      REAL     trprod(maxsp),trloss(maxsp)
      REAL     puits(45),source(45),kicovi
      INTEGER  nbsource,nbpuits
      REAL     sum_source,sum_puits

      CHARACTER*18 filename

! array for indexing production and loss reaction of each species
! type spec_reac_map is defined in prodloss.h
      TYPE(spec_reac_map),TARGET   :: lpmap(maxsp) ! loss_production map
      INTEGER               :: ispec, nprod, nloss

! photolysis rates printing
      INTEGER, DIMENSION(:), ALLOCATABLE :: idprintphoto
      REAL, DIMENSION(:), ALLOCATABLE    :: photorates
      CHARACTER(LEN=maxreac_char), ALLOCATABLE :: photoreac_names(:)
      CHARACTER(LEN=maxreac_char) :: printreaction
      CHARACTER(LEN=:),ALLOCATABLE :: cnjv ! string = # of j-values 
      INTEGER                      :: n_printphoto
!==variables for NetCDF file handling
!--File IDs (NetCDF)
      INTEGER :: ncid_in,ncid_out,ncid_prev
!--Attribute text
      CHARACTER(len=40) :: line1,line2
!==time output variables (for NetCDF output)
      INTEGER :: ntout,itout
! minimum concentration value to be output for Netcdf o/p file
      REAL,PARAMETER :: maskval=0.0
      REAL::  tmpconc(maxsp)        ! only concs >= maskval
      REAL::  tmpcro2(maxro2)       ! only concs >= maskval
      REAL::  tmpcaer(mxsat)        ! only concs >= maskval

! mandatory interface for extract_reacrates    
      INTERFACE
        SUBROUTINE extract_reacrates(qfor, idreacs, rate_tab )
          REAL, INTENT(IN)                 :: qfor(:)
          INTEGER, ALLOCATABLE,INTENT(IN)  :: idreacs(:)
          REAL, ALLOCATABLE,INTENT(INOUT)  :: rate_tab(:)
        END SUBROUTINE extract_reacrates
      END INTERFACE

!-----------------------------------------------------------
      END MODULE module_data_gecko_main
!-----------------------------------------------------------
