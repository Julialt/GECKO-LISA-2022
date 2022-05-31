
module link_ncdf
  implicit none
  contains
      SUBROUTINE wrtlinkncdf(lout,numsp,numre,num_n,num_m, &
                  numfo,numhv,numcvar,numextra,numo2,nummeo2, &
                  numain,numaou,numwin,numwou, nauxpar, &
                  numstoi,numiso,mx12stoi,mx1stoi,mx2stoi, &
                  chrsp,id_n,id_m, &
                  idfo,idhv,idcvar,idextra,ido2,idmeo2, &
                  idiso,idain,idaou, idwin,idwou, &
                  idrestoi,restoicf,idpdstoi,pdstoicf, &
                  arrhcf,focf,hvcf,hvfact,cvarcf,extracf, &
                  isocf,aoucf,woucf,wincf, &
                  nrpero,idreacro2,nrdimer,idreacdimer,wmol)

!===================================================================
! PURPOSE: create a NetCDF version of a GECKO-A mechanism
!          (aka a link file) to be read by the box model.
!          File outdat.nc is analogous to binary file outdat.akli
!          WITH THE ADDITION OF INFORMATION FROM FILES:
!          akparameter_module
!          {mech}.dict
!          {mech}.pvap
!          {mech}.Henry
!          The latter 3 files must be made available as links (indat.*)
!          and are directly read from within this subroutine.
! AUTHOR: Julia Lee-Taylor, NCAR, 18 Oct 2017
! ADAPTED FROM: original routine inca.f
!===================================================================
      USE sorting, ONLY : find_species_index
      USE akparameter_module
      USE general_module
      
      IMPLICIT NONE
      INCLUDE 'common.h'

      !INTEGER llink, lout
      INTEGER lout
      INTEGER numsp, numre, num_n, numo2, nummeo2,numiso
      INTEGER num_m, numfo, numhv, numcvar, numextra
! these dimensions are no longer used
      INTEGER mx12stoi, mx1stoi, mx2stoi, mx3stoi

      INTEGER numain, numaou, numwin, numwou
      INTEGER nauxpar(maxaux)
      INTEGER numstoi(maxre,2)
      INTEGER idrestoi(maxre,mxleft)
      INTEGER idpdstoi(maxre,mxright)
      INTEGER id_n(maxre), id_m(max_m)
      INTEGER idfo(maxfo,3)
      INTEGER idhv(maxhv), idcvar(maxcvar)
      INTEGER ido2(maxo2)
      INTEGER idiso(maxiso)
      INTEGER idmeo2(mxrpero)
      INTEGER idextra(maxextra)
      INTEGER idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)
      !INTEGER nrpero(mxrpero)
      INTEGER nrpero(maxro2)
      INTEGER idreacro2(mxrpero,maxro2)
      !INTEGER nrdimer(mxrdimer)
      INTEGER nrdimer(maxdimer)
      INTEGER idreacdimer(mxrdimer,maxdimer)

      REAL focf(maxaux+3,maxfo)
      REAL extracf(maxaux,maxextra)
! NB: must USE compilation option real-8
!      DOUBLE PRECISION extracf(maxaux,maxextra)
      REAL isocf(maxaux,maxiso)
      REAL restoicf(maxre,mxleft)
      REAL pdstoicf(maxre,mxright)
      REAL arrhcf(maxre,3)
      REAL hvcf(maxhv), hvfact(maxhv), cvarcf(maxcvar)
      REAL wmol(maxsp)
      REAL aoucf(2,maxt),woucf(3,maxt),wincf(3,maxt)

      CHARACTER(LEN=maxlsp)  chrsp(maxsp)

      INTEGER i,ire,j,k

      REAL tmp1dreal(3)

!==variables for NetCDF file
      INTEGER ncid

      INTEGER,PARAMETER:: maxaux3 = maxaux+3
      INTEGER,PARAMETER:: dim3 = 3, dim4 = 4

!--dictionary variables
! length of string & functionalities list in the dictionary
! dicnam, formula, functionalities list lengths of a given species
! parameters may be loaded in general.h
!      INTEGER,PARAMETER :: ldi=146,lco=6,lfo=(ldi-26),lfl=15

      INTEGER :: mxdic, max_n
      INTEGER :: iddic(maxsp)
      INTEGER :: igen(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo), target :: chem(maxsp)
      CHARACTER(LEN=lfl) :: code(maxsp)

!-- kOH,O3,NO3 variables
      INTEGER :: nkOH,nkO3,nkNO3
      INTEGER,DIMENSION(maxsp) :: idkOH,idkO3,idkNO3
      INTEGER,DIMENSION(maxsp) :: kOHid,kO3id,kNO3id
      REAL,DIMENSION(maxsp,3) :: kOHdat,kO3dat,kNO3dat
      CHARACTER(LEN=maxlsp),DIMENSION(maxsp) :: namkOH,namkO3,namkNO3

!--pvap variables
      INTEGER :: nsat
      INTEGER,DIMENSION(mxsat) :: idsat,idnan,idsim
      INTEGER,DIMENSION(maxsp) :: satid,nanid,simid,difid
      INTEGER,PARAMETER :: mxsim = 31
      REAL,DIMENSION(mxsat) :: Tb,dB,difvol
      REAL,DIMENSION(mxsat,2) :: nandat
      REAL,DIMENSION(mxsat,31) :: simpgroup
      REAL,DIMENSION(31,4) :: bk
      CHARACTER(LEN=maxlsp),DIMENSION(mxsat) :: namnan,namsim,namdif

!--Henry variables

! dimension mxdep now included in akparameter_module
!      INTEGER,PARAMETER :: mxinorgdep=10
!      INTEGER,PARAMETER :: mxdep=mxsat+mxinorgdep
!NB: dimension for dep in BOXMOD is given as maxsp
!    however this is FAR larger than needed!
! => MXDEP IS A NEW DIMENSION IN THIS VERSION.
      INTEGER :: ndep
      INTEGER :: iddep(mxdep)
      CHARACTER(LEN=maxlsp) :: namdep(mxdep)
      REAL,DIMENSION(mxdep,3) :: depdat

!--RO2 variables
      INTEGER nclro2
!      INTEGER,PARAMETER :: nclro2 = 9
      INTEGER numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
      INTEGER ipos

!--settings input variables: 1) declarations as in GECKO-A/RUN/main.f
      REAL :: critvp
      REAL :: cutoff_default
      REAL :: cutoff_OH
      REAL :: cutoff_O3
      REAL :: cutoff_NO3
      REAL :: cutoff_PAN
      REAL :: cutoff_HV
      REAL :: cutoff_RO
      REAL :: cutoff_RO2
      REAL :: cutoff_RCOO2
      REAL :: temp
      INTEGER :: maxgen

!--settings input variables: 2) declarations in common.h
      NAMELIST  /userparams/  critvp, cutoff_default, cutoff_OH, &
        cutoff_O3, cutoff_NO3, cutoff_PAN, cutoff_HV, cutoff_RO, &
        cutoff_RO2, cutoff_RCOO2, temp, maxgen, & !nrewind, &
        z_conformer_prop
      NAMELIST /flags/ wtflag,wtopeflag,per3_typefg,high_NOxfg, &
       zero_NOXfg,iflag,pvap_fg,kdissfg,kisomfg,hoadd_c1_fg,flgdhf, &
       autoox_fg,dimer_fg,masstransfg, spinup_fg, spinup_mech, &
       aeropart_fg, wallpart_fg, criegee_fg, stab_criegee_fg, &
       ro2ho2_fg, isopsoa_fg, aerophot_fg, OFR_fg, VBS_fg, kill_fg, &
       ro2no2_fg, pvapcd3_fg

!--test variables
      CHARACTER(LEN=ldi) :: line
      CHARACTER*20    :: filnam

!---------------------------
! CHECK THE INPUT ARGUMENTS
! --------------------------

      PRINT*,"numsp = ",numsp
      PRINT*,"numre = ",numre
      PRINT*,"num_m = ",num_m
      PRINT*,"num_n = ",num_n
      PRINT*,"numfo = ",numfo," : maxfo = ",maxfo
      PRINT*,"numhv = ",numhv
      PRINT*,"numcvar = ",numcvar
      PRINT*,"numextra = ",numextra
      PRINT*,"numo2 = ",numo2
      PRINT*,"nummeo2 = ",nummeo2
      PRINT*,"numain = ",numain
      PRINT*,"numaou = ",numaou
      PRINT*,"numwin = ",numwin
      PRINT*,"numwou = ",numwou
      PRINT*,"numiso = ",numiso

! ----------------------
!++CREATE AND OPEN NETCDF FILE IN DEFINE MODE
! ----------------------

      WRITE(6,*) 'in subroutine wrtncdf and writing the output ...'

      CALL open_ncfile_new("outdat.nc",ncid)
      PRINT*,"ncid = ",ncid

! ----------------------
!++DEFINE GLOBAL ATTRIBUTES
! ----------------------
!+++++CALL eznc_def_sysglobatt(ncid,AttName,syscall)
!AND>>CALL eznc_def_globalatt(ncid,AttName,your_text_here)

      CALL eznc_def_globalatt(ncid,"file_description", &
           "GECKO-A mechanism packaged for boxmodel input")
!++mechanism name ASSUMING mech is in a self-named directory
!++...and that this program is being run in that directory.
      CALL eznc_def_sysglobatt(ncid,"mechanism", &
            "`pwd | awk -F/ '{print $NF}'`")

      CALL eznc_def_sysglobatt(ncid,"date","`date`")
      CALL eznc_def_sysglobatt(ncid,"user","$USER")
      CALL eznc_def_sysglobatt(ncid,"platform","$HOSTNAME")

!++{Git_branch, Git_commit:
!-- if git calls do not work on your platform, substitute the following:
!++ reads ascii file generated by script run_interp.bash
      filnam="indat.gitinfo"
      OPEN(20,FILE=filnam,STATUS='OLD')
      READ(20,'(a)') line
      CALL eznc_def_globalatt(ncid,"Git_branch", line)
      READ(20,'(a)') line
      CALL eznc_def_globalatt(ncid,"Git_commit",line)
      CLOSE(20)

!++{Git_branch, Git_commit}
!        CALL eznc_def_sysglobatt(ncid,"Git_branch",
!     &       "`git rev-parse --abbrev-ref HEAD`")
!        CALL eznc_def_sysglobatt(ncid,"Git_commit",
!     &       "`git describe`")

! ----------------------
!++DEFINE MAX DIMENSIONS OF ARRAYS (per akparameter_module)
! ----------------------
!+++++CALL eznc_def_dim(ncid,DimName,DimSize)

      CALL eznc_def_dim(ncid,"maxlsp",maxlsp)
      CALL eznc_def_dim(ncid,"maxsp",maxsp)
      CALL eznc_def_dim(ncid,"mxleft",mxleft)
      CALL eznc_def_dim(ncid,"dim3",dim3)
      CALL eznc_def_dim(ncid,"dim4",dim4)
      CALL eznc_def_dim(ncid,"mxright",mxright)
!      CALL eznc_def_dim(ncid,"mx1stoi",mx1stoi)
!      CALL eznc_def_dim(ncid,"mx2stoi",mx2stoi)
!      CALL eznc_def_dim(ncid,"mx12stoi",mx12stoi)
      CALL eznc_def_dim(ncid,"maxre",maxre)
      CALL eznc_def_dim(ncid,"max_m",max_m)
      CALL eznc_def_dim(ncid,"maxfo",maxfo)
      CALL eznc_def_dim(ncid,"maxhv",maxhv)
      CALL eznc_def_dim(ncid,"maxcvar",maxcvar)
      CALL eznc_def_dim(ncid,"maxextra",maxextra)
      CALL eznc_def_dim(ncid,"maxo2",maxo2)
      CALL eznc_def_dim(ncid,"maxiso",maxiso)
      CALL eznc_def_dim(ncid,"maxt",maxt)
      CALL eznc_def_dim(ncid,"maxaux",maxaux)
      CALL eznc_def_dim(ncid,"maxaux3",maxaux3)
      CALL eznc_def_dim(ncid,"maxro2",maxro2)
      CALL eznc_def_dim(ncid,"mxro2cl",mxro2cl)
      CALL eznc_def_dim(ncid,"mxrpero",mxrpero)
      CALL eznc_def_dim(ncid,"maxdimer",maxdimer)
      CALL eznc_def_dim(ncid,"mxrdimer",mxrdimer)
      CALL eznc_def_dim(ncid,"mxsat",mxsat)
      CALL eznc_def_dim(ncid,"mxsim",mxsim)
      CALL eznc_def_dim(ncid,"mxdep",mxdep)
      CALL eznc_def_dim(ncid,"lco",lco)
      CALL eznc_def_dim(ncid,"lfo",lfo)
      CALL eznc_def_dim(ncid,"lfl",lfl)
      CALL eznc_def_dim(ncid,"maxcoe",maxcoe)
      CALL eznc_def_dim(ncid,"nset",nset)
      CALL eznc_def_dim(ncid,"maxang",maxang)
      CALL eznc_def_dim(ncid,"mchromo",mchromo)
      CALL eznc_def_dim(ncid,"mtopchromo",mtopchromo)
      CALL eznc_def_dim(ncid,"msptopchromo",msptopchromo)
      CALL eznc_def_dim(ncid,"mmedchromo",mmedchromo)
      CALL eznc_def_dim(ncid,"mspmedchromo",mspmedchromo)
      CALL eznc_def_dim(ncid,"nlo",nlo)
      CALL eznc_def_dim(ncid,"mbox",mbox)
      CALL eznc_def_dim(ncid,"mhd",mhd)
      CALL eznc_def_dim(ncid,"msur",msur)
      CALL eznc_def_dim(ncid,"mopc",mopc)
      CALL eznc_def_dim(ncid,"mpos",mpos)
      CALL eznc_def_dim(ncid,"mself",mself)
      CALL eznc_def_dim(ncid,"mes",mes)
      CALL eznc_def_dim(ncid,"maxem",maxem)
      CALL eznc_def_dim(ncid,"mtim",mtim)
      CALL eznc_def_dim(ncid,"mtr",mtr)
      CALL eznc_def_dim(ncid,"maxconst",maxconst)
      CALL eznc_def_dim(ncid,"maxinput",maxinput)

! -------------------------------------------
!++DEFINE USER FLAGS

      CALL eznc_def_0Dreal(ncid,"critvp")
      CALL eznc_def_localatt(ncid,"critvp","description", &
         "vapor pressure exponent below which "// &
         "species partition completely to aerosol")

      CALL eznc_def_0Dreal(ncid,"cutoff_default")
      CALL eznc_def_localatt(ncid,"cutoff_default","description", &
         "branching ratio below which reaction pathways not considered")

      CALL eznc_def_0Dreal(ncid,"cutoff_OH")
      CALL eznc_def_0Dreal(ncid,"cutoff_O3")
      CALL eznc_def_0Dreal(ncid,"cutoff_NO3")
      CALL eznc_def_0Dreal(ncid,"mech_temp")

      CALL eznc_def_0Dint(ncid,"maxgen")
      CALL eznc_def_localatt(ncid,"maxgen","description", &
         "maximum number of generations allowed")

      CALL eznc_def_0Dint(ncid,"z_conformer_prop")
      CALL eznc_def_localatt(ncid,"z_conformer_prop","description", &
         "temporary value for Z-E conformation of -CH=CH- alkenes")

      CALL eznc_def_0Dint(ncid,"per3_typefg")
      CALL eznc_def_localatt(ncid,"per3_typefg","description", &
        "flag to consider only 3 types of peroxys "// &
        "for recombination reactions")

      CALL eznc_def_0Dint(ncid,"high_NOxfg")
      CALL eznc_def_localatt(ncid,"high_NOxfg","description", &
         "flag to consider ONLY reactions with NOx "// &
         "for peroxys and acylperoxys ")

      CALL eznc_def_0Dint(ncid,"zero_NOxfg")
      CALL eznc_def_localatt(ncid,"zero_NOxfg","description", &
         "flag to NOT consider reactions with NOx "// &
         "for peroxys and acylperoxys ")

      CALL eznc_def_0Dint(ncid,"isomer_fg")
      CALL eznc_def_localatt(ncid,"isomer_fg","description", &
         "iflag allowing substitution by isomers (1) or not (0)")

      CALL eznc_def_0Dint(ncid,"flgdhf")
      CALL eznc_def_localatt(ncid,"flgdhf","description", &
         "flag to activate DHF formation : 1=DHF running ")

      CALL eznc_def_0Dint(ncid,"autoox_fg")
      CALL eznc_def_localatt(ncid,"autoox_fg","description", &
         "flag to activate auto-oxidation of peroxy radicals")

      CALL eznc_def_0Dint(ncid,"dimer_fg")
      CALL eznc_def_localatt(ncid,"dimer_fg","description", &
         "flag to activate dimerisation")

      CALL eznc_def_0Dint(ncid,"aeropart_fg")
      CALL eznc_def_localatt(ncid,"aeropart_fg","description", &
         "flag to activate aerosol partitioning")

      CALL eznc_def_0Dint(ncid,"wallpart_fg")
      CALL eznc_def_localatt(ncid,"wallpart_fg","description", &
         "flag to activate wall partitioning")

      CALL eznc_def_0Dint(ncid,"pvap_fg")
      CALL eznc_def_localatt(ncid,"pvap_fg","description", &
         "flag to select vapor pressure estimation method: "// &
         "1 = JR-MY; 2 = nannoolal; 3 = SIMPOL-1")

      CALL eznc_def_0Dint(ncid,"kdissfg")
      CALL eznc_def_localatt(ncid,"kdissfg","description", &
         "flag to select alkoxy decomposition estimation method: "// &
         "1 = Atkinson 2007; 2 = Vereecken 2009")

      CALL eznc_def_0Dint(ncid,"kisomfg")

      CALL eznc_def_0Dint(ncid,"hoadd_c1_fg")
      CALL eznc_def_localatt(ncid,"hoadd_c1_fg","description", &
         "flag for OH add-on alkene constant rate estimation method: "// &
         "1 = Peeters 1997; 2 = Ziemann 2009; "// &
         "3 = Magnify : Mike Jenkin SAR")

      CALL eznc_def_0Dint(ncid,"criegee_fg")
      CALL eznc_def_localatt(ncid,"criegee_fg","description", &
         "flag for criegee intermediates decomposition estimation "// &
         "method: 1 = Old version; 2 = New Version, CMV 2016")

      CALL eznc_def_0Dint(ncid,"stab_criegee_fg")
      CALL eznc_def_localatt(ncid,"criegee_fg","description", &
         "flag for stable criegee intermediates bimolecular reactions"// &
         " : 1 = Old version; 2 = New Version, CMV 2016")

      CALL eznc_def_0Dint(ncid,"masstransfg")
      CALL eznc_def_localatt(ncid,"masstransfg","description", &
         "flag to select the dynamic or thermodynamic mass transfer "// &
         "method : 1 = thermo; 2 = dynamic")

      CALL eznc_def_0Dint(ncid,"ro2ho2_fg")
      CALL eznc_def_localatt(ncid,"ro2ho2_fg","description", &
         "flag to select RO2 + HO2 reaction type ")
     
      CALL eznc_def_0Dint(ncid,"ro2no2_fg")
      CALL eznc_def_localatt(ncid,"ro2no2_fg","description", &
         "flag to select RO2 + NO2 reaction ")
     
      CALL eznc_def_0Dint(ncid,"ro2dep_fg")
      CALL eznc_def_localatt(ncid,"ro2dep_fg","description", &
         "flag to activate RO2 deposition ")

      CALL eznc_def_0Dint(ncid,"ro2cond_fg")
      CALL eznc_def_localatt(ncid,"ro2cond_fg","description", &
         "flag to activate RO2 condensation ")

      CALL eznc_def_0Dint(ncid,"isopsoa_fg")
      CALL eznc_def_localatt(ncid,"isopsoa_fg","description", &
         "flag to activate isoprene soa formation ")
     
      CALL eznc_def_0Dint(ncid,"aerophot_fg")
      CALL eznc_def_localatt(ncid,"aerophot_fg","description", &
         "flag to select photolysis of aerosol species")
     
      CALL eznc_def_0Dint(ncid,"OFR_fg")
      CALL eznc_def_localatt(ncid,"OFR_fg","description", &
         "flag to include OFR-specific reactions") 
     
      CALL eznc_def_0Dint(ncid,"pvapcd3_fg")
      CALL eznc_def_localatt(ncid,"pvapcd3_fg","description", &
         "flag to estimate pvap params for species with >2 = bonds") 
     
!--OLD CODE
!--    WRITE(llink)maxlsp,numsp,numre,num_n,
!     &            num_m,numfo,numhv,numcvar,numextra,numo2,nummeo2,
!     &            numain,numaou,numwin,numwou,
!     &            numiso,mx12stoi,mx1stoi,mx2stoi,
!     &            maxaux,
!     &            (nauxpar(i),i=1,maxaux)

! ----------------------
!++DEFINE ACTUAL SIZES OF ARRAYS AS INDEPENDENT SCALARS
! (VALUES ARE WRITTEN IN A LATER SECTION, IN "WRITE" MODE)
!!!CAUTION!!! THESE VALUES MAY BE ZERO
!->           DO NOT ATTEMPT TO DEFINE THEM AS DIMENSIONS!
! ----------------------
!+++++CALL eznc_def_0Dint(ncid,VarName)
!+++++CALL eznc_def_localatt(ncid,VarName,AttName,att_text)

      CALL eznc_def_0Dint(ncid,"numsp")
      CALL eznc_def_localatt(ncid,"numsp","title","actual number of species")

      CALL eznc_def_0Dint(ncid,"numre")
      CALL eznc_def_localatt(ncid,"numre","title","actual number of reactions")

      CALL eznc_def_0Dint(ncid,"num_n")
      CALL eznc_def_localatt(ncid,"num_n","title","actual number of simple thermal reactions")

      CALL eznc_def_0Dint(ncid,"num_m")
      CALL eznc_def_localatt(ncid,"num_m","title","actual number of +M reactions")

      CALL eznc_def_0Dint(ncid,"numfo")
      CALL eznc_def_localatt(ncid,"numfo","title","actual number of fo reactions")

      CALL eznc_def_0Dint(ncid,"numhv")
      CALL eznc_def_localatt(ncid,"numhv","title","actual number of hv reactions")

      CALL eznc_def_0Dint(ncid,"numcvar")
      CALL eznc_def_localatt(ncid,"numcvar","title","actual number of cvar reactions")

      CALL eznc_def_0Dint(ncid,"numextra")
      CALL eznc_def_localatt(ncid,"numextra","title","actual number of EXTRA reactions")

      CALL eznc_def_0Dint(ncid,"numo2")
      CALL eznc_def_localatt(ncid,"numo2","title","actual number of reactions with OXYGEN")

      CALL eznc_def_0Dint(ncid,"nummeo2")
      CALL eznc_def_localatt(ncid,"nummeo2","title","actual number of MeO2 reactions")

      CALL eznc_def_0Dint(ncid,"numain")
      CALL eznc_def_localatt(ncid,"numain","title","actual number of species undergoing phase equilibrium")

      CALL eznc_def_0Dint(ncid,"numaou")
      CALL eznc_def_localatt(ncid,"numaou","title","actual number of species undergoing particle-gas transfer")

      CALL eznc_def_0Dint(ncid,"numwin")
      CALL eznc_def_localatt(ncid,"numwin","title","actual number of species undergoing gas-wall transfer")

      CALL eznc_def_0Dint(ncid,"numwou")
      CALL eznc_def_localatt(ncid,"numwou","title","actual number of species undergoing wall-gas transfer")

      CALL eznc_def_0Dint(ncid,"ndim")
      CALL eznc_def_localatt(ncid,"ndim","title","actual number of species undergoing dimerization")

      CALL eznc_def_0Dint(ncid,"numiso")
      CALL eznc_def_localatt(ncid,"numiso","title","actual number of isomerizations")

!++OTHER PARAMETERS FROM AKPARAMETER.H
      CALL eznc_def_0Dint(ncid,"nsat")
      CALL eznc_def_localatt(ncid,"nsat","title","actual number of species for which Psat is evaluated")

      CALL eznc_def_0Dint(ncid,"ndep")
      CALL eznc_def_localatt(ncid,"ndep","title","actual number of species for which Vdep is evaluated")

      CALL eznc_def_0Dint(ncid,"nclro2")
      CALL eznc_def_localatt(ncid,"nclro2","title","actual number of classes of RO2")

!NOW USED AS A DIMENSION
!      CALL eznc_def_0Dint(ncid,"mxdic")
!      CALL eznc_def_localatt(ncid,"mxdic","title",
!     $    "actual number of species in dictionary")

! ----------------------
!++DEFINE VARIABLE ARRAYS
! ----------------------
!++INTEGER/REAL/DOUBLE/CHARACTER ARRAYS (1/2-D)
!++NB: DIMENSIONS ARE DEFINED HERE IN SAME ORDER AS FORTRAN ORDERING
!      EXCEPT FOR CHAR(dim,nchars)
!  !!!BUT!!! listed in the OTHER direction in the NetCDF file header !!!

!+++++CALL eznc_def_1Dint(ncid,VarName,DimName)
! (AND LATER...)
!+++++CALL eznc_put_1Dint(ncid,VarName,var,start,end)

      CALL eznc_def_1Dint(ncid,"nauxpar","maxaux")

!+++++CALL eznc_def_1Dchar(ncid,VarName,nchars,nvars)
      CALL eznc_def_1Dchar(ncid,"chrsp","maxlsp","maxsp")

      CALL eznc_def_1Dint(ncid,"id_n","maxre")
      CALL eznc_def_localatt(ncid,"id_n","title","rxn index of the ith thermal reaction")
      CALL eznc_def_localatt(ncid,"id_n","actual_size","(num_n)")

      CALL eznc_def_1Dint(ncid,"id_m","max_m")
      CALL eznc_def_localatt(ncid,"id_m","title","rxn index of the ith +M reaction")
      CALL eznc_def_localatt(ncid,"id_m","actual_size","(num_m)")

!+++++CALL eznc_def_2Dint(ncid,VarName,DimName1,DimName2)

      CALL eznc_def_2Dint(ncid,"idfo","maxfo","dim3")
      CALL eznc_def_localatt(ncid,"idfo","title","rxn index of the ith fall-off reaction")
      CALL eznc_def_localatt(ncid,"idfo","note","2nd number = 2, 3rd number = 0")
      CALL eznc_def_localatt(ncid,"idfo","actual_size","(numfo,3)")

      CALL eznc_def_1Dint(ncid,"idhv","maxhv")
      CALL eznc_def_localatt(ncid,"idhv","title","rxn index of the ith photolysis reaction")
      CALL eznc_def_localatt(ncid,"idhv","actual_size","(numhv)")

      CALL eznc_def_1Dint(ncid,"idcvar","maxcvar")
      CALL eznc_def_localatt(ncid,"idcvar","title","rxn index of the ith cvar reaction")
      CALL eznc_def_localatt(ncid,"idcvar","actual_size","(numcvar)")

      CALL eznc_def_1Dint(ncid,"idextra","maxextra")
      CALL eznc_def_localatt(ncid,"idextra","title","rxn index of the ith extra reaction")
      CALL eznc_def_localatt(ncid,"idextra","actual_size","(numextra)")

      CALL eznc_def_1Dint(ncid,"ido2","maxo2")
      CALL eznc_def_localatt(ncid,"ido2","title","rxn index of the ith reaction with 'OXYGEN'")
      CALL eznc_def_localatt(ncid,"ido2","actual_size","(numo2)")

      CALL eznc_def_1Dint(ncid,"idmeo2","mxrpero")
      CALL eznc_def_localatt(ncid,"ido2","title","rxn index of the ith reaction with CH3O2")
      CALL eznc_def_localatt(ncid,"idmeo2","actual_size","(nummeo2)")

      CALL eznc_def_1Dint(ncid,"idiso","maxiso")
      CALL eznc_def_localatt(ncid,"idiso","title","rxn index of the ith isomerization reaction")
      CALL eznc_def_localatt(ncid,"idiso","actual_size","(numiso)")

      CALL eznc_def_1Dint(ncid,"idain","maxt")
      CALL eznc_def_localatt(ncid,"idain","actual_size","(numain)")

      CALL eznc_def_1Dint(ncid,"idaou","maxt")
      CALL eznc_def_localatt(ncid,"idaou","actual_size","(numaou)")

      CALL eznc_def_1Dint(ncid,"idwin","maxt")
      CALL eznc_def_localatt(ncid,"idwin","actual_size","(numwin)")

      CALL eznc_def_1Dint(ncid,"idwou","maxt")
      CALL eznc_def_localatt(ncid,"idwou","actual_size","(numwou)")


!+++++CALL eznc_def_2Dreal(ncid,VarName,DimName1,DimName2)

      CALL eznc_def_2Dreal(ncid,"arrhcf","maxre","dim3")
      CALL eznc_def_localatt(ncid,"arrhcf","actual_size","(numre,3)")
      CALL eznc_def_localatt(ncid,"arrhcf","definition","for k=A*(T)^n*exp(-E/RT): 1 => ln(A); 2 =>  n; 3 => E/R")
      CALL eznc_def_localatt(ncid,"arrhcf","units","1) ln(molec.cm3.s) ;2) (none) ;3) Kelvin")

      CALL eznc_def_2Dint(ncid,"numstoi","maxre","mxleft")
      CALL eznc_def_localatt(ncid,"numstoi","title","numstoi(ire,1:2) = # stoi.coeffs. for rxn ire, both sides")

      CALL eznc_def_2Dint(ncid,"idrestoi","maxre","mxleft")
      CALL eznc_def_localatt(ncid,"idrestoi","title","idrestoi(i,k) = chrsp index of kth reactant in rxn i")

      CALL eznc_def_2Dreal(ncid,"restoicf","maxre","mxleft")
      CALL eznc_def_localatt(ncid,"restoicf","title","restoicf(i,k) = stoi.coeff. of kth reactant in rxn i")

      CALL eznc_def_2Dint(ncid,"idpdstoi","maxre","mxright")
      CALL eznc_def_localatt(ncid,"idpdstoi","title","idpdstoi(i,k) = chrsp index of kth product in rxn i")

      CALL eznc_def_2Dreal(ncid,"pdstoicf","maxre","mxright")
      CALL eznc_def_localatt(ncid,"pdstoicf","title","pdstocf(i,k) = stoi.coeff. of kth product in rxn i")

      CALL eznc_def_2Dreal(ncid,"focf","maxaux3","maxfo")
      CALL eznc_def_localatt(ncid,"focf","title","focf(j,i) = jth coefficient for the ith fo reaction")

!+++++CALL eznc_def_1Dreal(ncid,VarName,DimName)

      CALL eznc_def_1Dreal(ncid,"hvcf","maxhv")
      CALL eznc_def_localatt(ncid,"hvcf","title","lookup table label for the ith hv reaction")

      CALL eznc_def_1Dreal(ncid,"hvfact","maxhv")
      CALL eznc_def_localatt(ncid,"hvfact","title","coefficient for the ith hv reaction")

      CALL eznc_def_1Dreal(ncid,"cvarcf","maxcvar")
      CALL eznc_def_localatt(ncid,"cvarcf","title","label for the ith cvar reaction")

!+++++CALL eznc_def_2Ddbl(ncid,VarName,DimName1,DimName2)
! now using compilation option real-8 so extracf is a real.
      CALL eznc_def_2Ddbl(ncid,"extracf","maxaux","maxextra")
      CALL eznc_def_localatt(ncid,"extracf","title","extracf(j,i) = jth coefficient for the ith EXTRA rxn")

      CALL eznc_def_2Dreal(ncid,"isocf","maxaux","maxiso")
      CALL eznc_def_2Dreal(ncid,"aoucf","mxleft","maxt")
      CALL eznc_def_2Dreal(ncid,"woucf","dim3","maxt")
      CALL eznc_def_2Dreal(ncid,"wincf","dim3","maxt")

      CALL eznc_def_1Dint(ncid,"nrpero","maxro2")
      CALL eznc_def_localatt(ncid,"nrpero","title","actual number of reactions involving ith PERO")

      CALL eznc_def_2Dint(ncid,"idreacro2","mxrpero","maxro2")
      CALL eznc_def_localatt(ncid,"idreacro2","title","idreacro2(j,i) = rxn index of jth rxn involving PEROi")

      CALL eznc_def_1Dint(ncid,"nrdimer","maxdimer")
      CALL eznc_def_localatt(ncid,"nrdimer","title","actual number of reactions involving ith dimer")

      CALL eznc_def_2Dint(ncid,"idreacdimer","mxrdimer","maxdimer")
      CALL eznc_def_localatt(ncid,"idreacdimer","title","idreacdimer(j,i) = rxn index of jth rxn involving DIM_i")

      CALL eznc_def_1Dreal(ncid,"wmol","maxsp")
      CALL eznc_def_localatt(ncid,"wmol","congruence","chrsp")

!-----------------------------------------------------------------
      CALL eznc_def_0Dint(ncid,"nkOH")
!      CALL eznc_def_1Dchar(ncid,"namkOH","maxlsp","mxsat")
      CALL eznc_def_1Dint(ncid,"idkOH","mxsat")
      CALL eznc_def_localatt(ncid,"idkOH","title","chrsp index of ith kOH species")
      CALL eznc_def_1Dint(ncid,"kOHid","maxsp")
      CALL eznc_def_localatt(ncid,"kOHid","title","kOH index of ith species in chrsp list")
      CALL eznc_def_2Dreal(ncid,"kOHdat","mxsat","dim3")
      CALL eznc_def_localatt(ncid,"kOHdat","title","kOH arrhenius parameters")
!      CALL eznc_def_localatt(ncid,"kOHdat","congruence","namkOH")

      CALL eznc_def_0Dint(ncid,"nkNO3")
!      CALL eznc_def_1Dchar(ncid,"namkNO3","maxlsp","mxsat")
      CALL eznc_def_1Dint(ncid,"idkNO3","mxsat")
      CALL eznc_def_localatt(ncid,"idkNO3","title","chrsp index of ith kNO3 species")
      CALL eznc_def_1Dint(ncid,"kNO3id","maxsp")
      CALL eznc_def_localatt(ncid,"kNO3id","title","kNO3 index of ith species in chrsp list")
      CALL eznc_def_2Dreal(ncid,"kNO3dat","mxsat","dim3")
      CALL eznc_def_localatt(ncid,"kNO3dat","title","kNO3 arrhenius parameters")
!      CALL eznc_def_localatt(ncid,"kNO3dat","congruence","namkNO3")

      CALL eznc_def_0Dint(ncid,"nkO3")
!      CALL eznc_def_1Dchar(ncid,"namkO3","maxlsp","mxsat")
      CALL eznc_def_1Dint(ncid,"idkO3","mxsat")
      CALL eznc_def_localatt(ncid,"idkO3","title","chrsp index of ith kO3 species")
      CALL eznc_def_1Dint(ncid,"kO3id","maxsp")
      CALL eznc_def_localatt(ncid,"kO3id","title","kO3 index of ith species in chrsp list")
      CALL eznc_def_2Dreal(ncid,"kO3dat","mxsat","dim3")
      CALL eznc_def_localatt(ncid,"kO3dat","title","kO3 arrhenius parameters")
!      CALL eznc_def_localatt(ncid,"kO3dat","congruence","namkO3")

!-----------------------------------------------------------------
!      CALL eznc_def_1Dchar(ncid,"namnan","maxlsp","mxsat")
!      CALL eznc_def_1Dint(ncid,"idnan","mxsat")
!      CALL eznc_def_localatt(ncid,"idnan","title",
!     $    "chrsp index of ith species in Nannoolal file")
!      CALL eznc_def_1Dint(ncid,"nanid","maxsp")
!      CALL eznc_def_localatt(ncid,"nanid","title",
!     $    "Nannoolal file index of ith species in chrsp list")
!-----------------------------------------------------------------
      CALL eznc_def_1Dchar(ncid,"namsat","maxlsp","mxsat")
      CALL eznc_def_1Dint(ncid,"idsat","mxsat")
      CALL eznc_def_localatt(ncid,"idsat","title","chrsp index of ith saturation vapor pressure species")
      CALL eznc_def_1Dint(ncid,"satid","maxsp")
      CALL eznc_def_localatt(ncid,"satid","title","saturation vapor pressure index of ith species in chrsp list")
!-----------------------------------------------------------------

      CALL eznc_def_2Dreal(ncid,"nandat","mxsat","mxleft")
      CALL eznc_def_localatt(ncid,"nandat","congruence","namnan")
      CALL eznc_def_localatt(ncid,"nandat","actual_size","(nsat,2)")
      CALL eznc_def_localatt(ncid,"nandat","title","Parameters for Nannoolal vapor pressure calculation.")
      CALL eznc_def_localatt(ncid,"nandat","reference","Nannoolal 2004, 2008")
      CALL eznc_def_localatt(ncid,"nandat","title_1","Boiling point, Tb")
      CALL eznc_def_localatt(ncid,"nandat","units_1","Kelvin")
      CALL eznc_def_localatt(ncid,"nandat","title_2","Slope parameter, dB")
      CALL eznc_def_localatt(ncid,"nandat","definition_2","sum_i(N_i.C_i+GI)-0.176055;GI=sum group interaction contribs Ci")
      CALL eznc_def_localatt(ncid,"nandat","units_2","none")
      CALL eznc_def_localatt(ncid,"nandat","equation","logPsat = (4.1012+dB)*((Trb-1.)/(Trb-(1./8.))) ; Trb=T/Tb")
! IN CASE YOU WANT TO DEFINE NANNOOLAL PVAP PARAMS SEPARATELY....
!      CALL eznc_def_1Dreal(ncid,"Tb","mxsat")
!      CALL eznc_def_localatt(ncid,"Tb","congruence","namnan")
!      CALL eznc_def_localatt(ncid,"Tb","actual_size","(nsat)")
!      CALL eznc_def_localatt(ncid,"Tb","reference","Nannoolal 2004")
!      CALL eznc_def_localatt(ncid,"Tb","title","Boiling point")
!      CALL eznc_def_localatt(ncid,"Tb","units","Kelvin")
!      CALL eznc_def_1Dreal(ncid,"dB","mxsat")
!      CALL eznc_def_localatt(ncid,"dB","actual_size","(nsat)")
!      CALL eznc_def_localatt(ncid,"dB","reference","Nannoolal 2004")
!      CALL eznc_def_localatt(ncid,"dB","definition"",")
!      CALL eznc_def_localatt(ncid,"dB","units","Kelvin")
!      CALL eznc_def_localatt(ncid,"dB","congruence","namnan")
!-----------------------------------------------------------------
!      CALL eznc_def_1Dchar  (ncid,"namsim","maxlsp","mxsat")
!      CALL eznc_def_1Dint   (ncid,"idsim","mxsat")
!      CALL eznc_def_localatt(ncid,"idsim","title",
!     $    "chrsp index of ith species in SIMPOL file")
!      CALL eznc_def_1Dint   (ncid,"simid","maxsp")
!      CALL eznc_def_localatt(ncid,"simid","title",
!     $    "SIMPOL file index of ith species in chrsp list")
!-----------------------------------------------------------------
      CALL eznc_def_2Dreal  (ncid,"simpgroup","mxsat","mxsim")
      CALL eznc_def_localatt(ncid,"simpgroup","congruence","namsim")
      CALL eznc_def_localatt(ncid,"simpgroup","actual_size","(nsat,31)")
      CALL eznc_def_localatt(ncid,"simpgroup","title","Group contributions for SIMPOL vapor pressure calculation.")
      CALL eznc_def_localatt(ncid,"simpgroup","units","none")
      CALL eznc_def_localatt(ncid,"simpgroup","equation","log10Psat = sum over 31 groups(j) of [simpgroup(simid(i),j)* " &
                          // "((bk(j,1)/T) + bk(j,2) + (bk(j,3)*T) + (bk(j,4)*alog(T)))")
      CALL eznc_def_2Dreal  (ncid,"bk","mxsim","dim4")
      CALL eznc_def_localatt(ncid,"bk","congruence","simpgroup")
      CALL eznc_def_localatt(ncid,"bk","actual_size","(31,4)")
      CALL eznc_def_localatt(ncid,"bk","title","Group coefficients for SIMPOL vapor pressure calculation.")
      CALL eznc_def_localatt(ncid,"bk","units","none")
      CALL eznc_def_1Dchar  (ncid,"namdif","maxlsp","mxsat")
      CALL eznc_def_1Dreal  (ncid,"difvol","mxsat")
      CALL eznc_def_localatt(ncid,"difvol" ,"title","dimensionless diffusion volume of species")
      CALL eznc_def_1Dint   (ncid,"difid" ,"maxsp")
      CALL eznc_def_localatt(ncid,"difid" ,"title","difvol file index of ith species in chrsp list")

      CALL eznc_def_1Dchar(ncid,"namdep","maxlsp","mxdep")
      CALL eznc_def_1Dint(ncid,"iddep","mxdep")
      CALL eznc_def_localatt(ncid,"iddep","title","chrsp index of ith species in namdep array")

      CALL eznc_def_2Dreal(ncid,"depdat","mxdep","dim3")
      CALL eznc_def_localatt(ncid,"depdat","congruence","namdep")
      CALL eznc_def_localatt(ncid,"depdat","actual_size","(ndep,3)")
      CALL eznc_def_localatt(ncid,"depdat","title","Parameters for deposition velocity calculation.")
      CALL eznc_def_localatt(ncid,"depdat","reference","Wesely, 1989")
      CALL eznc_def_localatt(ncid,"depdat","definition_1","Diffusion coefficient ratio D_H2O/D_x")
      CALL eznc_def_localatt(ncid,"depdat","definition_2","Henry constant effective at pH=7")
      CALL eznc_def_localatt(ncid,"depdat","definition_3","Reactivity factor (0, 0.1 or 1)")

      CALL eznc_def_1Dint(ncid,"numchemro2","maxro2")
      CALL eznc_def_localatt(ncid,"numchemro2","title","actual number of species comprising each PEROi")

      CALL eznc_def_2Dint(ncid,"idchemro2","mxro2cl","maxro2")
      CALL eznc_def_localatt(ncid,"idchemro2","title","idchemro2(i,j) = chrsp index of ith PEROj")

! ------------------------------------
!++SWITCH TO WRITE (i.e. 'DATA') MODE
! ------------------------------------
      CALL switch_ncfile_to_data_mode(ncid)

! ----------------------
!++WRITE ("PUT") VARIABLES
!++VARIABLES ARE WRITTEN IN ORDER AS ALREADY DEFINED
! ----------------------
!++WRITE USER FLAGS

      CALL read_flags()
      print*,"done read_flags"

      CALL read_userparams( critvp, cutoff_default , cutoff_OH , & 
        cutoff_O3, cutoff_NO3, cutoff_PAN, cutoff_HV , &
        cutoff_RO , cutoff_RO2, cutoff_RCOO2, temp, &
        maxgen)
      print*,"done read_userparams"

      CALL eznc_put_0Dint(ncid,"z_conformer_prop",z_conformer_prop)
      CALL eznc_put_0Dint(ncid,"per3_typefg",per3_typefg)
      CALL eznc_put_0Dint(ncid,"high_NOxfg",high_NOxfg)
      CALL eznc_put_0Dint(ncid,"zero_NOxfg",zero_NOxfg)
      CALL eznc_put_0Dint(ncid,"isomer_fg",iflag)
      CALL eznc_put_0Dint(ncid,"flgdhf",flgdhf)
      CALL eznc_put_0Dint(ncid,"autoox_fg",autoox_fg)
      CALL eznc_put_0Dint(ncid,"dimer_fg",dimer_fg)
      CALL eznc_put_0Dint(ncid,"aeropart_fg",aeropart_fg)
      CALL eznc_put_0Dint(ncid,"wallpart_fg",wallpart_fg)
      CALL eznc_put_0Dint(ncid,"pvap_fg",pvap_fg)
      CALL eznc_put_0Dint(ncid,"kdissfg",kdissfg)
      CALL eznc_put_0Dint(ncid,"kisomfg",kisomfg)
      CALL eznc_put_0Dint(ncid,"hoadd_c1_fg",hoadd_c1_fg)
      CALL eznc_put_0Dint(ncid,"criegee_fg",criegee_fg)
      CALL eznc_put_0Dint(ncid,"stab_criegee_fg",stab_criegee_fg)
      CALL eznc_put_0Dint(ncid,"masstransfg",masstransfg)
      CALL eznc_put_0Dint(ncid,"maxgen",maxgen)
      CALL eznc_put_0Dint(ncid,"ro2ho2_fg",ro2ho2_fg)
      CALL eznc_put_0Dint(ncid,"ro2no2_fg",ro2no2_fg)
      CALL eznc_put_0Dint(ncid,"ro2dep_fg",ro2dep_fg)
      CALL eznc_put_0Dint(ncid,"ro2cond_fg",ro2cond_fg)
      CALL eznc_put_0Dint(ncid,"isopsoa_fg",isopsoa_fg)
      CALL eznc_put_0Dint(ncid,"aerophot_fg",aerophot_fg)
      CALL eznc_put_0Dint(ncid,"OFR_fg",OFR_fg)
      CALL eznc_put_0Dint(ncid,"pvapcd3_fg",pvapcd3_fg)

      CALL eznc_put_0Dreal(ncid,"critvp",critvp)
      CALL eznc_put_0Dreal(ncid,"cutoff_default",cutoff_default)
      CALL eznc_put_0Dreal(ncid,"cutoff_OH",cutoff_OH)
      CALL eznc_put_0Dreal(ncid,"cutoff_O3",cutoff_O3)
      CALL eznc_put_0Dreal(ncid,"cutoff_NO3",cutoff_NO3)
      CALL eznc_put_0Dreal(ncid,"mech_temp",temp)

!++INTEGER VALUES
!+++++CALL eznc_put_0Dint(ncid,VarName,var)

      CALL eznc_put_0Dint(ncid,"numsp",numsp)
      CALL eznc_put_0Dint(ncid,"numre",numre)
      CALL eznc_put_0Dint(ncid,"num_m",num_m)
      CALL eznc_put_0Dint(ncid,"num_n",num_n)
      CALL eznc_put_0Dint(ncid,"numhv",numhv)
      CALL eznc_put_0Dint(ncid,"numcvar",numcvar)
      CALL eznc_put_0Dint(ncid,"numextra",numextra)
      CALL eznc_put_0Dint(ncid,"numfo",numfo)
      CALL eznc_put_0Dint(ncid,"numo2",numo2)
      CALL eznc_put_0Dint(ncid,"nummeo2",nummeo2)
      CALL eznc_put_0Dint(ncid,"numain",numain)
      CALL eznc_put_0Dint(ncid,"numaou",numaou)
      CALL eznc_put_0Dint(ncid,"numwin",numwin)
      CALL eznc_put_0Dint(ncid,"numwou",numwou)
      CALL eznc_put_0Dint(ncid,"numiso",numiso)

!+++++CALL eznc_put_1Dint(ncid,VarName,var,start,end)

      CALL eznc_put_1Dint(ncid,"nauxpar",nauxpar,1,maxaux)

!++INTEGER ARRAYS (1-D)
!+++++CALL eznc_put_1Dint(ncid,VarName,var,start,end)

!--WRITE(llink) (id_n(i),i=1,num_n),
!--     &             (id_m(i),i=1,num_m),
!--     &             (idhv(i),i=1,numhv),
!--     &             (idcvar(i),i=1,numcvar),
!--     &             (idextra(i),i=1,numextra)
!--WRITE(llink) (ido2(i),i=1,numo2)
!--WRITE(llink) (idmeo2(i),i=1,nummeo2)
!--WRITE(llink) (idiso(i),i=1,numiso)
!--WRITE(llink) (idain(i),i=1,numain)
!--WRITE(llink) (idaou(i),i=1,numaou)
!--WRITE(llink) (idwin(i),i=1,numwin)
!--WRITE(llink) (idwou(i),i=1,numwou)
!--WRITE(llink) (nrpero(k),k=1,maxro2)
!--WRITE(llink) (nrdimer(k),k=1,maxdimer)

      CALL eznc_put_1Dint(ncid,"id_n",id_n,1,num_n)
      CALL eznc_put_1Dint(ncid,"id_m",id_m,1,num_m)
      CALL eznc_put_1Dint(ncid,"idhv",idhv,1,numhv)
      CALL eznc_put_1Dint(ncid,"idcvar",idcvar,1,numcvar)
      CALL eznc_put_1Dint(ncid,"idextra",idextra,1,numextra)
      CALL eznc_put_1Dint(ncid,"ido2",ido2,1,numo2)
      CALL eznc_put_1Dint(ncid,"idmeo2",idmeo2,1,nummeo2)
      CALL eznc_put_1Dint(ncid,"idiso",idiso,1,numiso)
      CALL eznc_put_1Dint(ncid,"idain",idain,1,numain)
      CALL eznc_put_1Dint(ncid,"idaou",idaou,1,numaou)
      CALL eznc_put_1Dint(ncid,"idwin",idwin,1,numwin)
      CALL eznc_put_1Dint(ncid,"idwou",idwou,1,numwou)
      CALL eznc_put_1Dint(ncid,"nrpero",nrpero,1,maxro2)
      CALL eznc_put_1Dint(ncid,"nrdimer",nrdimer,1,maxdimer)

! ----------------------
!++INTEGER ARRAYS (2-D: SPECIFY ARRAY SUB-SECTION)
! ----------------------

!--   ((idfo(i,k),k=1,3),i=1,numfo),
      CALL eznc_put_2Dint(ncid,"idfo",idfo(1:numfo,1:3),1,numfo,1,3)

!--WRITE(llink) ((idreacro2(i,k),i=1,nrpero(k)),k=1,maxro2)
!      DO k=1,maxro2
!        CALL eznc_put_2Dint(ncid,"idreacro2",
!     $                            idreacro2(1:nrpero(k),k),
!     $                                      1,nrpero(k),k,k)
!      ENDDO
        CALL eznc_put_2Dint(ncid,"idreacro2", &
                                  idreacro2(1:mxrpero,1:maxro2), &
                                            1,mxrpero,1,maxro2)

!--WRITE(llink) ((idchemdimer(i,k),i=1,nrdimer(k)),k=1,maxdimer)
!      DO k=1,maxdimer
!        CALL eznc_put_2Dint(ncid,"idreacdimer",
!     $                            idreacdimer(1:nrdimer(k),k),
!     $                                        1,nrdimer(k),k,k)
!      ENDDO
        CALL eznc_put_2Dint(ncid,"idreacdimer",idreacdimer(1:mxrdimer,1:maxdimer),1,mxrdimer,1,maxdimer)

!--WRITE(llink) ((numstoi(ire,k),k=1,2),ire=1,numre)
      CALL eznc_put_2Dint(ncid,"numstoi",numstoi(1:numre,1:2),1,numre,1,2)

!++INTEGER ARRAYS (2-D, IRREGULAR:
!++        USE LOOP & "START" TO SPECIFY ARRAY SUB-SECTIONS)

!--WRITE(llink) ((idrestoi(ire,i),i=1,numstoi(ire,1)),ire=1,numre)
!--             ((idpdstoi(ire,i),i=1,numstoi(ire,2)),ire=1,numre)
!      print*,"eznc_put_2Dreal:","idrestoi & idpdstoi"
!      DO ire=1,maxre
!        CALL eznc_put_2Dint(ncid,"idrestoi",
!     $                            idrestoi(ire,1:numstoi(ire,1)),
!     $                                     ire,ire,1,numstoi(ire,1))
!        CALL eznc_put_2Dint(ncid,"idpdstoi",
!     $                            idpdstoi(ire,1:numstoi(ire,2)),
!     $                                     ire,ire,1,numstoi(ire,2))
!      ENDDO

        CALL eznc_put_2Dint(ncid,"idrestoi",idrestoi(1:maxre,1:mxleft),1,maxre,1,mxleft)
        CALL eznc_put_2Dint(ncid,"idpdstoi",idpdstoi(1:maxre,1:mxright),1,maxre,1,mxright)

! ----------------------
!++REAL ARRAYS (1-D)
! ----------------------
!+++++CALL eznc_put_1Dreal(ncid,VarName,var,start,end)

!--WRITE(llink) (hvcf(k),k=1,numhv)
!--WRITE(llink) (wmol(i),i=1,numsp)
!--WRITE(llink) (hvfact(k),k=1,numhv)
!--WRITE(llink) (cvarcf(k),k=1,numcvar)

      CALL eznc_put_1Dreal(ncid,"hvcf",hvcf,1,numhv)
      CALL eznc_put_1Dreal(ncid,"wmol",wmol,1,numsp)
      CALL eznc_put_1Dreal(ncid,"hvfact",hvfact,1,numhv)
      CALL eznc_put_1Dreal(ncid,"cvarcf",cvarcf,1,numcvar)
      !CALL eznc_put_1Dreal(ncid,"Tb",Tb,1,nsat)
      !CALL eznc_put_1Dreal(ncid,"dB",dB,1,nsat)

! ----------------------
!++REAL ARRAYS (2-D; USE "START" & "END" TO SPECIFY ARRAY SUB-SECTION)
! ----------------------
!+++++>SUBROUTINE eznc_put_2Dreal(ncid,VarName,
!+++++$                                var(start1:end1,start2:end2),
!+++++$                                    start1,end1,start2,end2)
!NB: variables start*, count* describe location of array subsection in
!    completed array.

!--WRITE(llink)((arrhcf(ire,k),k=1,3),ire=1,numre),

      CALL eznc_put_2Dreal(ncid,"arrhcf",arrhcf(1:numre,1:3),1,numre,1,3)

!--WRITE(llink)((wincf(i,k),i=1,3),k=1,numwin)
      CALL eznc_put_2Dreal(ncid,"wincf",wincf(1:3,1:numwin),1,3,1,numwin)

!--WRITE(llink)((woucf(i,k),i=1,3),k=1,numwou)
      CALL eznc_put_2Dreal(ncid,"woucf",woucf(1:3,1:numwou),1,3,1,numwou)

!--WRITE(llink)((aoucf(i,k),i=1,2),k=1,numaou)
      CALL eznc_put_2Dreal(ncid,"aoucf",aoucf(1:2,1:numaou),1,2,1,numaou)

!--WRITE(llink)((focf(i,k),i=1,maxaux+3),k=1,numfo)
      CALL eznc_put_2Dreal(ncid,"focf",focf(1:maxaux3,1:numfo),1,maxaux3,1,numfo)

!--WRITE(llink)((isocf(i,k),i=1,maxaux),k=1,numiso)
      CALL eznc_put_2Dreal(ncid,"isocf",isocf(1:maxaux,1:numiso),1,maxaux,1,numiso)

!!!(THESE 2 ARE IRREGULAR: USE A LOOP)
!--WRITE(llink)((restoicf(ire,i),i=1,numstoi(ire,1)),ire=1,numre)
!--            ((pdstoicf(ire,i),i=1,numstoi(ire,2)),ire=1,numre)

!      DO ire=1,numre
!        CALL eznc_put_2Dreal(ncid,"restoicf",
!     $                             restoicf(ire,1:numstoi(ire,1)),
!     $                                      ire,ire,1,numstoi(ire,1))
!        CALL eznc_put_2Dreal(ncid,"pdstoicf",
!     $                             pdstoicf(ire,1:numstoi(ire,2)),
!     $                                      ire,ire,1,numstoi(ire,2))
!      ENDDO

        CALL eznc_put_2Dreal(ncid,"restoicf",restoicf(1:maxre,1:mxleft),1,maxre,1,mxleft)
        CALL eznc_put_2Dreal(ncid,"pdstoicf",pdstoicf(1:maxre,1:mxright),1,maxre,1,mxright)

! ----------------------
!++CHARACTER ARRAYS
! ----------------------

!--WRITE(llink) (chrsp(i),i=1,numsp)
      CALL eznc_put_1Dchar(ncid,"chrsp",chrsp,maxlsp,1,numsp)

!** DOUBLE PRECISION VARIABLE!
!** NOW AS A REAL WITH COMPILATION OPTION real-8
!--WRITE(llink)((extracf(i,k),i=1,maxaux),k=1,numextra)
!      CALL eznc_put_2Ddbl(ncid,"extracf",extracf,
      CALL eznc_put_2Dreal(ncid,"extracf",extracf(1:maxaux,1:numextra),1,maxaux,1,numextra)

! --------------------------------------------------
! OPEN AND READ DICTIONARY, PVAP, HENRY, RO2, SETTINGS FILES
! (doing this part LAST because it s time-consuming)
! --------------------------------------------------
      CALL read_lists(chrsp,mxdic,  dicnam, iddic, chem,code,igen, &
                            nkOH,   namkOH, idkOH, kOHid, kOHdat, &
                            nkNO3,  namkNO3,idkNO3,kNO3id,kNO3dat, &
                            nkO3,   namkO3, idkO3, kO3id, kO3dat, &
                            nsat,   idsat, satid, &
                            namnan, idnan,  nanid, nandat, &
                            namsim, idsim,  simid, simpgroup, bk, &
                            namdif,         difid, difvol, &
                            ndep,  namdep,    iddep,   depdat, &
                            nclro2,numchemro2,idchemro2)
      CALL eznc_def_dim(ncid,"mxdic",mxdic)     
      
      CALL eznc_def_1Dchar(ncid,"dicnam","lco","mxdic")

      CALL eznc_def_1Dint(ncid,"iddic","mxdic")
      CALL eznc_def_localatt(ncid,"iddic","title","FIRST chrsp index corresponding to ith species in dicnam")

      CALL eznc_def_1Dchar(ncid,"chem","lfo","mxdic")
      CALL eznc_def_localatt(ncid,"chem","title","full chemical formula of species")
      CALL eznc_def_localatt(ncid,"chem","congruence","dicnam")

      CALL eznc_def_1Dchar(ncid,"code","lfl","mxdic")
      CALL eznc_def_localatt(ncid,"code","title","codes for species functional groups")
      CALL eznc_def_localatt(ncid,"code","congruence","dicnam")
      
      CALL eznc_def_1Dint(ncid,"igen","mxdic")
      CALL eznc_def_localatt(ncid,"igen","title", &
          "number of stable generations to get to this species")
      CALL eznc_def_localatt(ncid,"igen","congruence","dicnam")
      
      CALL eznc_put_1Dint(ncid,"iddic",iddic,1,mxdic)
      CALL eznc_put_1Dint(ncid,"igen",igen,1,mxdic)
      CALL eznc_put_1Dchar(ncid,"dicnam",dicnam,lco,1,mxdic)
      CALL eznc_put_1Dchar(ncid,"chem",chem,lfo,1,mxdic)
      CALL eznc_put_1Dchar(ncid,"code",code,lfl,1,mxdic)
      print*,"done put_dict"

! integrated reaction-rate parameters
      CALL eznc_put_0Dint(ncid,"nkOH",nkOH)
      CALL eznc_put_1Dint(ncid,"idkOH",idkOH,1,nkOH)
      CALL eznc_put_1Dint(ncid,"kOHid",kOHid,1,numsp)
!      CALL eznc_put_1Dchar(ncid,"namkOH",namkOH,maxlsp,1,nkOH)
      CALL eznc_put_2Dreal(ncid,"kOHdat",kOHdat(1:nkOH,1:3),1,nkOH,1,3)

      CALL eznc_put_0Dint(ncid,"nkNO3",nkNO3)
      CALL eznc_put_1Dint(ncid,"idkNO3",idkNO3,1,nkNO3)
      CALL eznc_put_1Dint(ncid,"kNO3id",kNO3id,1,numsp)
!      CALL eznc_put_1Dchar(ncid,"namkNO3",namkNO3,maxlsp,1,nkNO3)
      CALL eznc_put_2Dreal(ncid,"kNO3dat",kNO3dat(1:nkNO3,1:3),1,nkNO3,1,3)

      CALL eznc_put_0Dint(ncid,"nkO3",nkO3)
      CALL eznc_put_1Dint(ncid,"idkO3",idkO3,1,nkO3)
      CALL eznc_put_1Dint(ncid,"kO3id",kO3id,1,numsp)
!      CALL eznc_put_1Dchar(ncid,"namkO3",namkO3,maxlsp,1,nkO3)
      CALL eznc_put_2Dreal(ncid,"kO3dat",kO3dat(1:nkO3,1:3),1,nkO3,1,3)

! pvap-related parameters
      CALL eznc_put_0Dint(ncid,"nsat",nsat)
! diffusion volume
      CALL eznc_put_1Dint(ncid,"difid",difid,1,numsp)
      CALL eznc_put_1Dchar(ncid,"namdif",namdif,maxlsp,1,nsat)
      CALL eznc_put_1Dreal(ncid,"difvol", &
                                 difvol(1:nsat), &
                                        1,nsat)
! pvap parameters: Universal
      CALL eznc_put_1Dint(ncid,"idsat",idsat,1,nsat)
      CALL eznc_put_1Dint(ncid,"satid",satid,1,numsp)
      CALL eznc_put_1Dchar(ncid,"namsat",namnan,maxlsp,1,nsat)

!! pvap parameters: Nannoolal
!      CALL eznc_put_1Dint(ncid,"idnan",idnan,1,nsat)
!      CALL eznc_put_1Dint(ncid,"nanid",nanid,1,numsp)
!      CALL eznc_put_1Dchar(ncid,"namnan",namnan,maxlsp,1,nsat)
      CALL eznc_put_2Dreal(ncid,"nandat",nandat(1:nsat,1:2),1,nsat,1,2)
! pvap parameters: SIMPOL
!      CALL eznc_put_1Dint(ncid,"idsim",idsim,1,nsat)
!      CALL eznc_put_1Dint(ncid,"simid",simid,1,numsp)
!      CALL eznc_put_1Dchar(ncid,"namsim",namsim,maxlsp,1,nsat)
      CALL eznc_put_2Dreal(ncid,"simpgroup",simpgroup(1:nsat,1:mxsim),1,nsat,1,mxsim)
      CALL eznc_put_2Dreal(ncid,"bk",bk(1:mxsim,1:dim4),1,mxsim,1,dim4)
      print*,"done put_pvap"

!      CALL read_dep(chrsp,ndep,namdep,iddep,depdat)
      CALL eznc_put_0Dint(ncid,"ndep",ndep)
      CALL eznc_put_1Dint(ncid,"iddep",iddep,1,ndep)
      CALL eznc_put_1Dchar(ncid,"namdep",namdep,maxlsp,1,ndep)
      CALL eznc_put_2Dreal(ncid,"depdat",depdat(1:ndep,1:3),1,ndep,1,3)
      print*,"done put_dep"



      CALL eznc_put_0Dint(ncid,"nclro2",nclro2)
      CALL eznc_put_1Dint(ncid,"numchemro2",numchemro2,1,maxro2)
!--WRITE(llink) ((idchemro2(i,k),i=1,nrpero(k)),k=1,maxro2)
      print*,"eznc_put_2Dint: idchemro2"
! BUG !
!      DO k=1,nclro2
!        CALL eznc_put_2Dint(ncid,"idchemro2",
!     $                            idchemro2(:,k),
!     $                                      1,numchemro2(k),k,k)
!      ENDDO
        CALL eznc_put_2Dint(ncid,"idchemro2",idchemro2(1:mxro2cl,1:maxro2),1,mxro2cl,1,maxro2)
      print*,"done put_ro2"

1000  CONTINUE
! ----------------------
!++CLOSE NetCDF Dataset
! ----------------------

      CALL close_ncfile(ncid)

      RETURN
      END SUBROUTINE wrtlinkncdf

! ===================================================================

      SUBROUTINE read_lists(chrsp,ndic,dicnam,iddic,chem,code,igen, &
                            nkOH,  namkOH, idkOH, kOHid, kOHdat, &
                            nkNO3, namkNO3,idkNO3,kNO3id,kNO3dat, &
                            nkO3,  namkO3, idkO3, kO3id, kO3dat, &
                            nsat,  idsat,satid, &
                            namnan,idnan,nanid,nandat, &
                            namsim,idsim,simid,simpgroup,bk, &
                            namdif,      difid,difdat, &
                            ndep,namdep,iddep,      depdat, &
                            nclro2,numchemro2,idchemro2)

!PURPOSE: read ascii dictionary file  for writing to NCDF.
!PURPOSE: read ascii file of Nannoolal pvap parameters for writing to NCDF.
!PURPOSE: read ascii file of SIMPOL pvap parameters for writing to NCDF.
!PURPOSE: read ascii file of molecular diffusion volumes for writing to NCDF.
!PURPOSE: read ascii file of Henry parameters for writing to NCDF.
!PURPOSE: read ascii files ('X-files') of RO2 indices for writing to NCDF.
!AUTHOR: julial, NCAR, 2018-03-15.

!....................................

      USE akparameter_module
      USE general_module, 
      USE sorting, only: sort_species,switch_sorting_to_chem,srtid_chem
      IMPLICIT NONE

! length of string & functionalities list in the dictionary
! dicnam, formula, functionalities list lengths of a given species
! following parameters are loaded from general.h
!     INTEGER,PARAMETER :: ldi=146,lco=6,lfo=(ldi-26),lfl=15
! define dimension mxdep
! WHY? : dimension for dep in BOXMOD is given as maxsp
!    mxsat is defined in akparameter_module
!    however this is FAR larger than needed!
! => MXDEP IS ALSO NOW DEFINED IN akparameter_module
      INTEGER,PARAMETER :: ninorgdep=8

!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
      INTEGER nclro2

!OUTPUT VARIABLES
      INTEGER :: ndic
      INTEGER :: iddic(maxsp)
      INTEGER :: igen(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo) :: chem(maxsp)
      CHARACTER(LEN=lfl) :: code(maxsp)
      INTEGER :: nsat
      INTEGER :: idsat(mxsat),satid(maxsp)
      INTEGER :: idnan(mxsat),nanid(maxsp)
      INTEGER :: idsim(mxsat),simid(maxsp)
      INTEGER ::              difid(maxsp)
      CHARACTER(LEN=maxlsp),DIMENSION(mxsat) :: namnan,namsim,namdif
      REAL,DIMENSION(mxsat) :: difdat
      REAL,DIMENSION(mxsat,2) :: nandat
      REAL,DIMENSION(mxsat,31) :: simpgroup
      REAL,DIMENSION(31,4) :: bk
      INTEGER :: ndep
      INTEGER :: iddep(mxdep)
      CHARACTER(LEN=maxlsp) :: namdep(mxdep)
      REAL,DIMENSION(mxdep,3) :: depdat
      INTEGER numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
!-- kOH,O3,NO3 variables
      INTEGER :: nkOH,nkO3,nkNO3
      INTEGER,DIMENSION(maxsp) :: idkOH,idkO3,idkNO3
      INTEGER,DIMENSION(maxsp) :: kOHid,kO3id,kNO3id
      REAL,DIMENSION(maxsp,3) :: kOHdat,kO3dat,kNO3dat
      CHARACTER(LEN=maxlsp),DIMENSION(maxsp) :: namkOH,namkO3,namkNO3

!INTERNAL VARIABLES
      INTEGER :: i,j,k,ipos,ch1,ch2
      REAL    :: Tb,dB
      INTEGER,PARAMETER :: lout = 11
      INTEGER,PARAMETER :: nhdr_pvap=1
      INTEGER,PARAMETER :: nhdr_Henry=10
      CHARACTER*20      :: filnam
      CHARACTER(LEN=ldi)   :: line
      CHARACTER(LEN=maxlsp) ro2sp(mxrpero,maxro2)

! initialize

      idchemro2 = 0
      simid = 0
      nanid = 0
      idsat = 0
      idsim = 0
      idnan = 0
      chem = ""

! read dictionary info
      CALL read_dict(chrsp,ndic,dicnam,iddic,chem,code,igen)

! read Nannoolal dicnams, pvap params
      CALL read_pvap_nan(chrsp,nsat,namnan,idnan,nanid,nandat)

      idsat(:)=idnan(:)
      satid(:)=nanid(:)

! read ascii files of SIMPOL pvap parameters for writing to NCDF.
      CALL read_pvap_sim(chrsp,nsat,namsim,idsim,simid,simpgroup,bk)

! read diffusion volume dicnams, diffusion volumes
      CALL read_difvol(chrsp,nsat,namdif,difid,difdat)

! read info for deposition
      CALL read_dep(chrsp,ndep,namdep,iddep,depdat)

! read RO2 information
      CALL read_ro2(chrsp,nclro2,numchemro2,idchemro2)

! read kOH parameters
! switch the sorting module to chem mode to allow looking up formulas
      call switch_sorting_to_chem(chem, ndic)
      do i = lbound(chem,1), ubound(chem,1)
        write(96,*) trim(chem(i))
      enddo
! sort chem
      write(6,*) 'sort the list of formulas ...'
      call sort_species()
      write(6,*) '      end of sort'
      do i = lbound(srtid_chem,1), ubound(srtid_chem,1)
        write(97,*) trim(chem(srtid_chem(i)))
      enddo
      
      write(6,*) ' -- start reading kOH ...'
      CALL read_kOH(chrsp,ndic,chem,dicnam,iddic,nkOH,namkOH,idkOH,kOHid,kOHdat)
      write(6,*) ' ... done'

      write(6,*) ' -- start reading kNO3 ...'
      CALL read_kNO3(chrsp,ndic,chem,dicnam,iddic,nkNO3,namkNO3,idkNO3,kNO3id,kNO3dat)
      write(6,*) ' ... done'

      write(6,*) ' -- start reading kO3 ...'
      CALL read_kO3(chrsp,ndic,chem,dicnam,iddic,nkO3,namkO3,idkO3,kO3id,kO3dat)
      write(6,*) ' ... done'
      END SUBROUTINE read_lists

! ===================================================================
      SUBROUTINE read_dict(chrsp,ndic,dicnam,iddic,chem,code,igen)

!PURPOSE: read ascii dictionary file  for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-10-20.
      USE sorting, ONLY : find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
! length of string & functionalities list in the dictionary
! dicnam, formula, functionalities list lengths of a given species
! following parameters are loaded from general.h
!     INTEGER,PARAMETER :: ldi=149,lco=6,lfo=(ldi-26),lfl=15
!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
! OUTPUT
      INTEGER :: ndic
      INTEGER :: iddic(maxsp)
      INTEGER :: igen(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo) :: chem(maxsp)
      CHARACTER(LEN=lfl) :: code(maxsp)
! INTERNAL VARIABLES
      INTEGER i,j,ch1,ch2,ipos
      CHARACTER(LEN=ldi)::  line
      CHARACTER(LEN=20) :: filnam
      CHARACTER(LEN=3)  :: cgen
!................................................

      ndic=0
      igen = 0
      filnam="indat.dict"
      OPEN(20,FILE=filnam,STATUS='OLD')

      DO i=1,maxsp
         READ(20,'(a)') line
         IF (line(1:4).EQ.'****') EXIT
         !IF (line(1:5).EQ.'     ') CYCLE
         !IF (line(1:5).EQ.'HV   ') CYCLE
         !IF (line(1:5).EQ.'M    ') CYCLE
         !IF (line(1:5).EQ.'(M)  ') CYCLE
         ndic=ndic+1
         ch1 = 1
         ch2 = lco
         READ(line(ch1:ch2),'(a)') dicnam(ndic)
         ch1 = 10
         ch2 = ch1+lfo-1
         READ(line(ch1:ch2),'(a)') chem(ndic)
         ch1 = ch2+2+1
         ch2 = ch1+lfl-1
         READ(line(ch1:ch2),'(a)') code(ndic)
         ch1 = ch2+1
         ch2 = ch1+3-1
         READ(line(ch1:ch2),'(a)') cgen
         IF(cgen.NE."   ") READ(cgen,*) igen(ndic)
      ENDDO
      WRITE(6,*) 'number of species in the dictionary=',ndic
      WRITE(6,*) 'maxsp =',maxsp
      IF(maxsp.LT.ndic) STOP
      CLOSE (20)

! set up indexing
      DO i=1,ndic
        ipos = find_species_index('G'//dicnam(i)(1:lco), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, dic species unidentified in readdic'
          WRITE(6,'(a)') dicnam(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF
        iddic(i) = ipos
      ENDDO

      END SUBROUTINE read_dict
!===================================================================
      SUBROUTINE read_difvol(chrsp,nsat,namdif,difid,difdat)
!PURPOSE: read ascii file of diffusion volumes for writing to NCDF.
!AUTHOR: julial, NCAR, 2019-02-05.
      USE sorting, ONLY : find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
      INTEGER :: nsat
!OUTPUT VARIABLES
      INTEGER :: difid(maxsp)
      REAL,DIMENSION(mxsat) :: difdat
      CHARACTER(LEN=maxlsp) :: namdif(mxsat)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER*20    :: filnam
!....................................

      filnam="indat.difv"
      OPEN(20,FILE=filnam,STATUS='OLD')

! read dicnams, data
      DO j=1,mxsat
        READ(20,'(a)') line
        IF(line(1:3).EQ.'END')EXIT
        READ(line,*,err=9) namdif(j), difdat(j)
        CYCLE ! bypass error handling if read successful
9       WRITE(lout,*) '-error--, while reading indat.dif'
        WRITE(lout,*) '        , at line : ',j
        STOP
      ENDDO
      CLOSE(20)
! set up indexing: 
! difid = dif index of chrsp species
      DO i=1,nsat
        ipos = find_species_index(namdif(i), .false.)
        if (ipos .eq. 0) then
          WRITE(6,*) '--error--, sat species unidentified in readdif'
          WRITE(6,'(a)') namdif(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        endif
! identify GAS species
        difid(ipos) = i          
! loop to catch AER species
        DO j=ipos+1,maxsp
          IF(namdif(i)(2:7).EQ.chrsp(j)(2:7)) THEN
            difid(j)=i
            EXIT
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE read_difvol
!===================================================================
      SUBROUTINE read_pvap_nan(chrsp,nsat,namnan,idnan,nanid,nandat)
!PURPOSE: read ascii file of Nannoolal pvap parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-10-20.
      USE sorting, ONLY: find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
!OUTPUT VARIABLES
      INTEGER :: nsat
      INTEGER :: idnan(mxsat)
      INTEGER :: nanid(maxsp)
      REAL,DIMENSION(mxsat,2) :: nandat
      CHARACTER(LEN=maxlsp) :: namnan(mxsat)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      REAL    :: Tb,dB
      INTEGER,PARAMETER :: nhdr_pvap=1
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER*20    :: filnam
!....................................

      filnam="indat.pnan"
      OPEN(20,FILE=filnam,STATUS='OLD')
      DO j=1,nhdr_pvap
        READ(20,'(a)') line
      ENDDO

! read Nannoolal dicnams, pvap params
      IF(line(1:6).EQ.'NANNOO')THEN
        DO j=1,mxsat
          nsat = j-1
          READ(20,'(a)') line
          IF(line(1:3).EQ.'END')EXIT
          READ(line,*,err=9) namnan(j), Tb, dB
          nandat(j,1)=Tb
          nandat(j,2)=dB
          CYCLE ! bypass error handling if read successful
9         WRITE(lout,*) '-error--, while reading indat.pnan'
          WRITE(lout,*) '        , at line : ',j
          STOP
        ENDDO
      ELSE
        WRITE(lout,*) 'pvap data is not Nannoolal format'
        WRITE(lout,*) '=> must add appropriate read/write to wrtncdf.f'
        STOP
      ENDIF
      CLOSE(20)

! set up indexing: 
! idnan = chrsp index of nan species
! nanid = nan index of chrsp species

      DO i = 1, nsat

! find 'G' species
        ipos = find_species_index(namnan(i), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, sat species unidentified in readnan'
          WRITE(6,'(a)') namnan(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF
        idnan(i) = ipos
        nanid(ipos) = i        

! find 'A' species
        ipos = find_species_index('A'//namnan(i)(2:len(namnan(i))),.TRUE.)
        IF(ipos.EQ.0) THEN
          WRITE(6,*) '--error--, sat species unidentified in readnan'
          WRITE(6,'(a)') 'A'//namnan(i)(2:len(namnan(i)))
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF
        nanid(ipos) = i        

! find 'W' species (if any)
        ipos = find_species_index('W'//namnan(i)(2:len(namnan(i))),.TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'W' species
        nanid(ipos) = i        

      ENDDO

      END SUBROUTINE read_pvap_nan
!===================================================================
      SUBROUTINE read_pvap_sim(chrsp,nsat,namsim,idsim,simid,simpgroup,bk)
!PURPOSE: read ascii files of SIMPOL pvap parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2019-02-04
      USE sorting, ONLY: find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
!OUTPUT VARIABLES
      INTEGER :: nsat
      INTEGER :: idsim(mxsat)
      INTEGER :: simid(maxsp)
      REAL,DIMENSION(mxsat,31) :: simpgroup
      REAL,DIMENSION(31,4) :: bk
      CHARACTER(LEN=maxlsp) :: namsim(mxsat)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      INTEGER,PARAMETER :: nhdr_pvap=1
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER*20    :: filnam
!....................................

      filnam="indat.psim"
      OPEN(20,FILE=filnam,STATUS='OLD')
      DO j=1,nhdr_pvap
        READ(20,'(a)') line
      ENDDO

! read SIMPOL dicnams, group contributions
      IF(line(1:6).EQ.'SIMPOL')THEN
        DO i=1,mxsat
          nsat = i-1
          READ(20,'(a)') line
          IF(line(1:3).EQ.'END')EXIT
          READ(line,'(A7,2x,31(f3.0))',err=890) namsim(i),(simpgroup(i,j),j=1,31)
          CYCLE ! bypass error handling if read successful
890       WRITE(lout,*) '-error--, while reading indat.psim'
          WRITE(lout,*) '        , at line : ',i
          STOP
        ENDDO
      ELSE
        WRITE(lout,*) 'pvap data is not SIMPOL format'
        WRITE(lout,*) '=> must add appropriate read/write to wrtncdf.f'
        STOP
      ENDIF
      CLOSE(20)

! set up indexing: 
! idsim = chrsp index of sim species
! simid = sim index of chrsp species
      DO i = 1, nsat
        ipos = find_species_index(namsim(i), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, sat species unidentified in readsat'
          WRITE(6,'(a)') namsim(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF
        idsim(i) = ipos
        simid(ipos) = i        
      ENDDO

! read SIMPOL group coefficients
      filnam="simpol.dat"
      OPEN(20,FILE=filnam,STATUS='OLD')
      DO i=1,31
        READ(20,'(f10.4,2x,f9.6,2x,f11.8,2x,f11.8)') bk(i,1),bk(i,2), bk(i,3),bk(i,4)
      ENDDO
      CLOSE(20)
      END SUBROUTINE read_pvap_sim
! ===================================================================
      SUBROUTINE read_dep(chrsp,ndep,namdep,iddep,depdat)
!PURPOSE: read ascii file of Henry parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-10-20.
      USE sorting, ONLY : find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
! define dimension mxdep
! WHY? : dimension for dep in BOXMOD is given as maxsp
!    mxsat is defined in akparameter_module
!    however this is FAR larger than needed!
! => MXDEP IS ALSO NOW DEFINED IN akparameter_module
      INTEGER,PARAMETER :: ninorgdep=8
      !INTEGER,PARAMETER :: mxdep=mxsat+ninorgdep
!INPUT
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
!OUTPUT
      INTEGER :: ndep
      INTEGER :: iddep(mxdep)
      CHARACTER(LEN=maxlsp) :: namdep(mxdep)
      REAL,DIMENSION(mxdep,3) :: depdat
!OTHER INTERNAL VALUES
      INTEGER i,j,k,ipos
      INTEGER,PARAMETER :: nhdr_Henry=10
      INTEGER,PARAMETER :: lout = 11
      CHARACTER*20      :: filnam
      CHARACTER(LEN=ldi)   :: line
!.................................................

      filnam="indat.Henry"
      OPEN(20,FILE=filnam,STATUS='OLD')
      DO j=1,nhdr_Henry
        READ(20,'(a)') line
      ENDDO

! read inorganic and organic species dicnams & Henry params,
      ndep=0
      DO j=1,mxdep
         READ(20,'(a)') line
         IF (line(1:3).eq.'END') EXIT
         ndep=ndep+1
         !READ(line,'(1pe7.1,2x,1pe7.1,2x,1pe7.1,2x,a16)',err=10)
         READ(line,*,err=10) (depdat(j,k),k=1,3),namdep(j)
         CYCLE ! bypass error handling if read successful
10       WRITE(lout,*) '-error--, while reading indat.Henry'
         WRITE(lout,*) '        , at line : ',j
         STOP
      ENDDO
      CLOSE (20)

! set up indexing
      DO i=1,ndep
        ipos = find_species_index(namdep(i), .true.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, dep species unidentified in readdep'
          WRITE(6,'(a)') namdep(i)
          WRITE(6,*) 'in file ', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
          STOP
        ENDIF
        iddep(i) = ipos
      ENDDO

      END SUBROUTINE read_dep
! ===================================================================
      SUBROUTINE read_ro2(chrsp,nclro2,numchemro2,idchemro2)
!PURPOSE: read ascii files ('X-files') of RO2 indices for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-10-20.
      USE sorting, ONLY : find_species_index
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      INTEGER nclro2
      CHARACTER(LEN=maxlsp)  chrsp(maxsp)
!OUTPUT
      INTEGER numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
!INTERNAL
      INTEGER i,ipos,j,k
      INTEGER,PARAMETER :: lout = 11
      CHARACTER*20      :: filnam
      CHARACTER(LEN=ldi)   :: line
      CHARACTER(LEN=maxlsp) ro2sp(mxrpero,maxro2)
! initialize
      numchemro2 = 0
      nclro2 = 9
      DO k=1,nclro2
! open file
        filnam(1:)='indat_.ro2'
        WRITE(filnam(6:6),'(i1)') k
        WRITE(6,*) '    opening ',filnam
        OPEN (20, file=filnam,STATUS='OLD')
! read RO2 in the class
        DO i=1,mxrpero+1
          READ(20,'(a)',END=59) line
          IF (line(1:3).EQ.'***') EXIT
          IF (i.GT.mxrpero) THEN
            WRITE(6,*) '--error-- number of RO2 in file > mxrpero'
            WRITE(6,*)   filnam
            STOP
          ENDIF
          IF (INDEX(line,' ').LE.1) THEN
            WRITE(6,*) '--error--, lecture RO2 in file ', filnam
            STOP
          ENDIF
          ro2sp(i,k)=line(1:7)
          numchemro2(k)=i
59      ENDDO
        CLOSE(20)
! set up indexing
        DO i=1,numchemro2(k)
          ipos = find_species_index(ro2sp(i,k), .true.)
          IF(ipos.EQ.0)THEN
            WRITE(6,*) '--error--, RO2 species unidentified in readro2'
            WRITE(6,*) i,numchemro2(k),ro2sp(i,k)
            WRITE(6,*) 'in file ', filnam
            WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE'
            WRITE(6,*) 'SPECIES MAY NOT BE SORTED'
            STOP
          ENDIF
          idchemro2(i,k) = ipos
        ENDDO
      ENDDO
      END SUBROUTINE read_ro2

! ===================================================================
      SUBROUTINE read_flags()
!PURPOSE: read user input file "userparams.input" ("GECKO_INPUTS/settings_*")
!AUTHOR: julial, NCAR, 2018-03-07.

      USE akparameter_module
      USE general_module
      IMPLICIT NONE
      INCLUDE 'common.h'


!--settings input variables: 2) declarations in common.h

      !Declare namelist to be used for user flags
      NAMELIST /flags/ wtflag,wtopeflag,per3_typefg,high_NOxfg, &
       zero_NOXfg,iflag,pvap_fg,kdissfg,kisomfg,hoadd_c1_fg,flgdhf, &
       autoox_fg,dimer_fg,masstransfg, spinup_fg, spinup_mech, &
       aeropart_fg, wallpart_fg, criegee_fg, stab_criegee_fg, &
       ro2ho2_fg, isopsoa_fg, aerophot_fg, OFR_fg , VBS_fg, kill_fg, &
       ro2no2_fg, pvapcd3_fg


! open & read file

      OPEN (20, file="userparams.input",STATUS='OLD')
      READ(20, nml  = flags )
      CLOSE(20)

      END SUBROUTINE read_flags

! ===================================================================
      SUBROUTINE read_userparams( critvp, cutoff_default , cutoff_OH , &
        cutoff_O3, cutoff_NO3, cutoff_PAN, cutoff_HV , &
        cutoff_RO , cutoff_RO2, cutoff_RCOO2, temp, &
        maxgen)
!PURPOSE: read user input file "userparams.input" ("GECKO_INPUTS/settings_*")
!AUTHOR: julial, NCAR, 2018-03-07.

      USE akparameter_module
      USE general_module
      IMPLICIT NONE
      INCLUDE 'common.h'

!--settings input variables: 1) declarations as in GECKO-A/RUN/main.f
      REAL :: critvp
      REAL :: cutoff_default
      REAL :: cutoff_OH
      REAL :: cutoff_O3
      REAL :: cutoff_NO3
      REAL :: cutoff_PAN
      REAL :: cutoff_HV
      REAL :: cutoff_RO
      REAL :: cutoff_RO2
      REAL :: cutoff_RCOO2
      REAL :: temp
      INTEGER :: maxgen

!--settings input variables: 2) declarations in common.h
      !Declare namelist to be used for user parameters file
      NAMELIST  /userparams/  critvp, cutoff_default, cutoff_OH, &
         cutoff_O3, cutoff_NO3, cutoff_PAN, cutoff_HV, cutoff_RO, &
         cutoff_RO2, cutoff_RCOO2, temp, maxgen, &!nrewind, &
         z_conformer_prop

! open & read file
      OPEN (20, file="userparams.input",STATUS='OLD')
      READ(20, nml  = userparams )
      CLOSE(20)

      END SUBROUTINE read_userparams
!===================================================================
      SUBROUTINE read_kOH(chrsp,ndic,chem,dicnam,iddic, &
                           nkOH,namkOH,idkOH,kOHid,kOHdat)
!PURPOSE: read ascii file of kOH parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-11-04

! DBUG: dicnam desn't appear to be what we need!

      USE sorting, ONLY: find_species_index, find_formula_index, &
                switch_sorting_to_chrsp,switch_sorting_to_chem
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      INTEGER :: ndic,iddic(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo), target :: chem(maxsp)
      CHARACTER(LEN=maxlsp), target :: chrsp(maxsp)
!OUTPUT VARIABLES
      INTEGER :: nkOH
      INTEGER :: idkOH(maxsp)
      INTEGER :: kOHid(maxsp)
      REAL,DIMENSION(maxsp,3) :: kOHdat
      CHARACTER(LEN=maxlsp) :: namkOH(maxsp)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      REAL    :: tarrhcf(3)
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER(LEN=20)  :: filnam
      CHARACTER(LEN=lfo) :: tchem
!....................................

      nkOH = 0
      filnam="indat.kOH"
      OPEN(20,FILE=filnam,STATUS='OLD')

! read kOH Arrhenius coeffs, formulae
      DO i=1,maxsp
        READ(20,'(a)') line
        IF(line(1:3).EQ.'END')EXIT
        READ(line,*,err=9) tarrhcf(:),tchem
        nkOH = nkOH+1
        kOHdat(i,1:3)=tarrhcf(:)

! set up indexing: 
! idkOH = chrsp index of kOH species
! kOHid = kOH index of chrsp species

!        WRITE(10,*)i,kOHdat(i,1:3),tchem(1:LEN_TRIM(tchem))

! find 'chrsp' from knowledge of chem
        call switch_sorting_to_chem(chem, ndic)    
        j = find_formula_index(tchem)
        
        if (j .le. 0) then
          write(6,*) ' error in read_kOH, '
          write(6,*) ' could not find ', trim(tchem), ' in chem '
          cycle
        endif
        
!        WRITE(11,*)i,j,tchem(1:LEN_TRIM(tchem)),' ',chem(j)(1:LEN_TRIM(chem(j)))

        namkOH(i)=dicnam(j)(1:lco)

!        WRITE(12,*)i,j,iddic(j),namkOH(i),tchem(1:LEN_TRIM(tchem))
        

! find 'G' species indices
        call switch_sorting_to_chrsp(chrsp)
        ipos = find_species_index('G'//namkOH(i), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, species unidentified in read_kOH'
          WRITE(6,'(a)') chem
          WRITE(6,*) 'in file ', filnam
          STOP
        ENDIF
        IF(ipos.NE.iddic(j))THEN
          WRITE(6,*) '--error--, species unidentified in read_kOH'
          WRITE(6,'(a)') ipos,'<>',iddic(j)
          STOP
        ENDIF

        idkOH(i) = ipos
        kOHid(ipos) = i        

! find 'A' species
        ipos = find_species_index('A'//namkOH(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'A' species matches
        kOHid(ipos) = i        

! find 'W' species (if any)
        ipos = find_species_index('W'//namkOH(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'W' species matches
        kOHid(ipos) = i        
      enddo      

      CLOSE(20)

      RETURN

9     WRITE(lout,*) '-error--, while reading indat.pkOH'
      WRITE(lout,*) '        , at line : ',j
      STOP

      END SUBROUTINE read_kOH
!===================================================================
      SUBROUTINE read_kO3(chrsp,ndic,chem,dicnam,iddic, &
                           nkO3,namkO3,idkO3,kO3id,kO3dat)
!PURPOSE: read ascii file of kO3 parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-11-04

! DBUG: dicnam desn't appear to be what we need!

      USE sorting, ONLY: find_species_index, find_formula_index, &
                switch_sorting_to_chrsp,switch_sorting_to_chem
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      INTEGER :: ndic,iddic(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo) :: chem(maxsp)
      CHARACTER(LEN=maxlsp) :: chrsp(maxsp)
!OUTPUT VARIABLES
      INTEGER :: nkO3
      INTEGER :: idkO3(maxsp)
      INTEGER :: kO3id(maxsp)
      REAL,DIMENSION(maxsp,3) :: kO3dat
      CHARACTER(LEN=maxlsp) :: namkO3(maxsp)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      REAL    :: tarrhcf(3)
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER(LEN=20)  :: filnam
      CHARACTER(LEN=lfo) :: tchem
!....................................

      nkO3 = 0
      filnam="indat.kO3"
      OPEN(20,FILE=filnam,STATUS='OLD')

! read kO3 Arrhenius coeffs, formulae
      DO i=1,maxsp
        READ(20,'(a)') line
        IF(line(1:3).EQ.'END')EXIT
        READ(line,*,err=9) tarrhcf(:),tchem
        nkO3 = nkO3+1
        kO3dat(nkO3,1:3)=tarrhcf(:)

! set up indexing: 
! idkO3 = chrsp index of kO3 species
! kO3id = kO3 index of chrsp species

!        WRITE(10,*)i,kO3dat(i,1:3),tchem(1:LEN_TRIM(tchem))


! find 'chrsp' from knowledge of chem
        call switch_sorting_to_chem(chem, ndic)    
        j = find_formula_index(tchem)
        
        if (j .le. 0) then
          write(6,*) ' error in read_kO3, '
          write(6,*) ' could not find ', trim(tchem), ' in chem '
          cycle
        endif
                
        namkO3(i)=dicnam(j)(1:lco)

        call switch_sorting_to_chrsp(chrsp)
        ipos = find_species_index('G'//namkO3(i), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, species unidentified in read_kO3'
          WRITE(6,'(a)') chem
          WRITE(6,*) 'in file ', filnam
          STOP
        ENDIF
        IF(ipos.NE.iddic(j))THEN
          WRITE(6,*) '--error--, species unidentified in read_kO3'
          WRITE(6,'(a)') ipos,'<>',iddic(j)
          STOP
        ENDIF

        idkO3(i) = ipos
        kO3id(ipos) = i        

! find 'A' species
        ipos = find_species_index('A'//namkO3(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'A' species matches
        kO3id(ipos) = i        

! find 'W' species (if any)
        ipos = find_species_index('W'//namkO3(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'W' species matches
        kO3id(ipos) = i       

      ENDDO
      CLOSE(20)

      RETURN

9       WRITE(lout,*) '-error--, while reading indat.pkO3'
        WRITE(lout,*) '        , at line : ',j
        STOP

      END SUBROUTINE read_kO3
!===================================================================
      SUBROUTINE read_kNO3(chrsp,ndic,chem,dicnam,iddic, &
                           nkNO3,namkNO3,idkNO3,kNO3id,kNO3dat)
!PURPOSE: read ascii file of kNO3 parameters for writing to NCDF.
!AUTHOR: julial, NCAR, 2017-11-04

! DBUG: dicnam desn't appear to be what we need!

      USE sorting, ONLY: find_species_index, find_formula_index, &
                switch_sorting_to_chrsp,switch_sorting_to_chem
      USE akparameter_module
      USE general_module
      IMPLICIT NONE
!INPUT
      INTEGER :: ndic,iddic(maxsp)
      CHARACTER(LEN=lco) :: dicnam(maxsp)
      CHARACTER(LEN=lfo) :: chem(maxsp)
      CHARACTER(LEN=maxlsp) :: chrsp(maxsp)
!OUTPUT VARIABLES
      INTEGER :: nkNO3
      INTEGER :: idkNO3(maxsp)
      INTEGER :: kNO3id(maxsp)
      REAL,DIMENSION(maxsp,3) :: kNO3dat
      CHARACTER(LEN=maxlsp) :: namkNO3(maxsp)
!INTERNAL VARIABLES
      INTEGER :: i,j,j1,ipos,k,k2
      REAL    :: tarrhcf(3)
      INTEGER,PARAMETER :: lout = 11
      CHARACTER(LEN=ldi) :: line
      CHARACTER(LEN=20)  :: filnam
      CHARACTER(LEN=lfo) :: tchem
!....................................

      nkNO3 = 0
      filnam="indat.kNO3"
      OPEN(20,FILE=filnam,STATUS='OLD')

! read kNO3 Arrhenius coeffs, formulae
      DO i=1,maxsp
        READ(20,'(a)') line
        IF(line(1:3).EQ.'END')EXIT
        READ(line,*,err=9) tarrhcf(:),tchem
        nkNO3 = nkNO3+1
        kNO3dat(nkNO3,1:3)=tarrhcf(:)

! set up indexing: 
! idkNO3 = chrsp index of kNO3 species
! kNO3id = kNO3 index of chrsp species

 !       WRITE(10,*)i,knO3dat(i,1:3),tchem(1:LEN_TRIM(tchem))

! find 'chrsp' from knowledge of chem
        call switch_sorting_to_chem(chem, ndic)    
        j = find_formula_index(tchem)
        
        if (j .le. 0) then
          write(6,*) ' error in read_kO3, '
          write(6,*) ' could not find ', trim(tchem), ' in chem '
          cycle
        endif
               
        namkNO3(i)=dicnam(j)(1:lco)

! find 'G' species indices
        call switch_sorting_to_chrsp(chrsp)
        ipos = find_species_index('G'//namkNO3(i), .false.)
        IF(ipos.EQ.0)THEN
          WRITE(6,*) '--error--, species unidentified in read_kNO3'
          WRITE(6,'(a)') chem
          WRITE(6,*) 'in file ', filnam
          STOP
        ENDIF
        IF(ipos.NE.iddic(j))THEN
          WRITE(6,*) '--error--, species unidentified in read_kNO3'
          WRITE(6,'(a)') ipos,'<>',iddic(j)
          STOP
        ENDIF

        idkNO3(i) = ipos
        kNO3id(ipos) = i        

! find 'A' species
        ipos = find_species_index('A'//namkNO3(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'A' species matches
        kNO3id(ipos) = i        

! find 'W' species (if any)
        ipos = find_species_index('W'//namkNO3(i), .TRUE.)
        IF(ipos.EQ.0) CYCLE ! no 'W' species matches
        kNO3id(ipos) = i          

      ENDDO
      CLOSE(20)

      RETURN

9       WRITE(lout,*) '-error--, while reading indat.pkNO3'
        WRITE(lout,*) '        , at line : ',j
        STOP

      END SUBROUTINE read_kNO3
end module link_ncdf
