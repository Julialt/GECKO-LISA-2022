***********************************************************************
* SUBROUTINE spakinit9_ncdf : this routine reads the chemical scheme
* produced by the interpreter "inca". 
* NetCDF VERSION, julial, NCAR, Oct 2017
* Analogous to spakinit9.f
* INPUT ARGUMENTS : 
* OUTPUT : working array sizes and values from mechanism file
!! NOTE !! Array size values (dimensions) are in BOTH akparameter_mod AND
!!         NetCDF file. We are using akparameter_mod here to avoid needing
!!         to use allocatable arrays. Future model versions could do
!!         that and avoid keeping track of multiple ascii input files.
!!         Commented calls to "get_dimension" are placeholders for that
!!         future work. They are not required if akparameter_mod is included.
*                  lout = Fortran unit # for diagnostic output (errors)
*              (currently hardwired, not passed, in SUBROUTINE handle_err)
************************************************************************
      SUBROUTINE spakinit9_ncdf(ncid,
     1                     numsp,numre,num_n,num_m,numfo,numhv,
     2                     numcvar,numextra,numo2,nummeo2,nrpero,
     3                     numiso,idiso,
     3                     mx12stoi,nauxpar,numstoi,idrestoi,idpdstoi,
     4                     id_m,idfo,
     5                     idhv,idcvar,idextra,ido2,idmeo2,idreacro2,
     5                     idreacdimer,nrdimer,
     6                     nself,idselfreac,
     7                nt1chromo,chromo1cf,id1chromo,
     7                ntmedchromo,chromomedcf,nummedchromo,idmedchromo,
     7                chromotopcf,numtopchromo,idtopchromo,
     8                     restoicf,pdstoicf,arrhcf,focf,
     9                     hvcf,hvfact,cvarcf,extracf,isocf,
     1                     numain, numaou, numwin, numwou,
     1                     idain,idaou,idwin,idwou,
     1                     !aoucf,woucf,wincf,
     1                     woucf,wincf,
     o                     wmol,chrsp)

      USE netcdf
      USE akparameter_module
      USE flags_module,ONLY: OFR_fg
      USE module_data_gecko_main,ONLY: difvol,nchrom,idchrom

      IMPLICIT NONE

      CHARACTER(maxlsp) chrsp(maxsp)

      INTEGER  numsp, numre, num_n
      INTEGER  num_m, numfo, numhv, numcvar, numextra
      INTEGER  mx12stoi
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  nauxpar(maxaux)
      INTEGER  numstoi(maxre,2)
      INTEGER  idrestoi(maxre,mxleft),idpdstoi(maxre,mxright)
      INTEGER  id_n(maxre)
      INTEGER  id_m(max_m)
      INTEGER  idfo(maxfo,3)
      INTEGER  idhv(maxhv), idcvar(maxcvar)
      INTEGER  idextra(maxextra)
      INTEGER  nself,idselfreac(mself,2)
      INTEGER  numo2,ido2(maxo2)
      INTEGER  numiso,idiso(maxiso)
      INTEGER  nummeo2,idmeo2(mxrpero)
      INTEGER  nrpero(maxro2),idreacro2(mxrpero,maxro2)
      INTEGER  nrdimer(maxdimer),idreacdimer(mxrdimer,maxdimer)
      INTEGER  idain(maxt),idaou(maxt),idwin(maxt),idwou(maxt)

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
      REAL     extracf(maxaux,maxextra)
      !DOUBLE PRECISION     extracf(maxaux,maxextra)
      REAL     isocf(maxaux,maxiso)
      REAL     wmol(maxsp)
      !REAL     aoucf(2,maxt),woucf(3,maxt),wincf(3,maxt)
      REAL     woucf(3,maxt),wincf(3,maxt)


* internal
      INTEGER  ntchromo,chromocf(mchromo),numchromo(mchromo)
      INTEGER  label, mxlab, cphvcf(maxhv)
      LOGICAL  lostop
      INTEGER  ierr, i, j, mxlsp, mxsp, k, ire, mxaux, index
      INTEGER  idcf
      INTEGER  mx1stoi, mx2stoi
      INTEGER  tmp_fg

* NetCDF 
      INTEGER         ncid,status,nDims,nVars, nGlobalAtts, unlimdimid
      INTEGER         lout
      CHARACTER(20)  filename, text

* -----------------------------
      lostop=.false.
      ierr=0

      print*,"in spakinit9_ncdf"

* -----------------------------
* initialization of all values
* -----------------------------

* initialization of the strings
      DO i=1,maxsp
        chrsp(i)=' '
      ENDDO

* initialization of the integer variables
      numsp=0
      numre=0
      num_n=0
      num_m=0
      numfo=0
      numhv=0
      numcvar=0
      numextra=0
      numo2=0
      numiso=0
      nummeo2=0
! defined as DIMENSIONS in netcdf link file: read from
! akparameter_module
!      mx12stoi=0
!      mx1stoi=0
!      mx2stoi=0

      nauxpar(:)=0
      numstoi(:,:)=0
      idrestoi(:,:)=0
      idpdstoi(:,:)=0
      id_n(:)=0
      id_m(:)=0
      idfo(:,:)=0
      idhv(:)=0
      idcvar(:)=0
      idextra(:)=0
      ido2(:)=0
      idmeo2(:)=0

* initialization of real variables
      arrhcf(:,:)=0.
      restoicf(:,:)=0.
      pdstoicf(:,:)=0.
      focf(:,:)=0.
      hvcf(:)=0.
      hvfact(:)=0.
      cvarcf(:)=0.
      extracf(:,:)=0.
      isocf(:,:)=0
      wmol(:)=0.

* ------------------------------------------
* reading the chemical scheme (NetCDF file)
* checking some dimensions against akparameter_module versions
* ------------------------------------------
! check various things about the dataset
!      status = NF90_INQUIRE(ncid, nDims, nVars, nGlobalAtts, unlimdimid)
!      text = "NF90_INQUIRE, "//filename
!      if (status /= nf90_noerr) call eznc_handle_err(status,text)
!      print*,"nDims =",nDims
!      print*,"nVars =",nVars

! switch to DATA mode for read
      CALL switch_ncfile_to_data_mode(ncid)

* mxlsp = NetCDF dimension maxlsp, for diagnostic comparison with akparameter
* species length
      CALL eznc_get_dimension(ncid,"maxlsp",mxlsp)
        IF (mxlsp.NE.maxlsp) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, length of species names given by the '
          WRITE(6,*) '         chemical interpreter does not correspond'
          WRITE(6,*) '         to the length used in this  program'
          WRITE(6,*) '         see parameter MAXLSP (in akparameter)'
          WRITE(6,*) 'mxlsp= ',mxlsp," maxlsp= ",maxlsp
          lostop=.true.
        ENDIF

* mxsp = NetCDF dimension maxsp, for diagnostic comparison with akparameter
* species number
      CALL eznc_get_dimension(ncid,"maxsp",mxsp)
        IF (mxsp.NE.maxsp) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, max number of  species given by the '
          WRITE(6,*) '         chemical interpreter does not correspond'
          WRITE(6,*) '         to the value used in this  program'
          WRITE(6,*) '         see parameter MAXSP (in akparameter)'
          WRITE(6,*) 'mxsp= ',mxsp," maxsp= ",maxsp
          lostop=.true.
        ENDIF

* numsp
      CALL eznc_get_0Dint(ncid,"numsp",numsp)
        IF (numsp.GT.maxsp) THEN
          WRITE(6,*)
          WRITE(6,*) '--error--, number of species is larger '
          WRITE(6,*) '           than the parameter MAXSP'
          WRITE(6,*) '           see parameter MAXSP (in akparameter)'
          WRITE(6,*) 'numsp= ',numsp," maxsp= ",maxsp
          lostop=.true.
        ENDIF

* numre
* find and check the number of reactions
!      CALL eznc_get_dimension(ncid,"maxre",maxre)
      CALL eznc_get_0Dint(ncid,"numre",numre)
        IF (numre.GT.maxre) THEN
          WRITE(6,*)
          WRITE(6,*) '--error--, number of reactions is larger '
          WRITE(6,*) '           than the parameter MAXRE'
          WRITE(6,*) '           see parameter MAXRE (in akparameter)'
          lostop=.true.
        ENDIF

* num_n
* find and check the number of "normal" reactions
      CALL eznc_get_0Dint(ncid,"num_n",num_n)
        IF (num_n.GT.maxre) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of simple thermal reactions > '
          WRITE(6,*) '         than the parameter MAXRE'
          WRITE(6,*) '         see parameter MAXRE (in akparameter)'
          lostop=.true.
        ENDIF

* num_m
* find and check the number of reactions with "M"
!      CALL eznc_get_dimension(ncid,"max_m",max_m)
      CALL eznc_get_0Dint(ncid,"num_m",num_m)
        IF (num_m.GT.max_m) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of third body reactions > '
          WRITE(6,*) '         than the parameter MAX_M'
          WRITE(6,*) '         see parameter MAX_M (in akparameter)'
          lostop=.true.
        ENDIF

* numfo
* find and check the number of fall_off reactions
!      CALL eznc_get_dimension(ncid,"maxfo",maxfo)
      CALL eznc_get_0Dint(ncid,"numfo",numfo)
        IF (numfo.GT.maxfo) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of falloff reactions is larger '
          WRITE(6,*) '         than the parameter MAXFO'
          WRITE(6,*) '         see parameter MAXFO (in akparameter)'
          lostop=.true.
        ENDIF

* numhv
* find and check the number of photolysis reactions
!      CALL eznc_get_dimension(ncid,"maxhv",maxhv)
      CALL eznc_get_0Dint(ncid,"numhv",numhv)
        IF (numhv.GT.maxhv) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of hv reactions is larger '
          WRITE(6,*) '         than the parameter MAXHV'
          WRITE(6,*) '         see parameter MAXHV (in akparameter)'
          lostop=.true.
        ENDIF

* numcvar
* find and check the number of cvar reactions
!      CALL eznc_get_dimension(ncid,"maxcvar",maxcvar)
      CALL eznc_get_0Dint(ncid,"numcvar",numcvar)
        IF (numcvar.GT.maxcvar) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of cvar reactions is larger '
          WRITE(6,*) '         than the parameter MAXCVAR'
          WRITE(6,*) '         see parameter MAXCVAR (in akparameter)'
          lostop=.true.
        ENDIF

* numextra
* find and check the number of extra reactions
!      CALL eznc_get_dimension(ncid,"maxextra",maxextra)
      CALL eznc_get_0Dint(ncid,"numextra",numextra)
        IF (numextra.GT.maxextra) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of extra reactions '
          WRITE(6,*) '         than the parameter MAXEXTRA'
          WRITE(6,*) '        see parameter MAXEXTRA (in akparameter)'
          lostop=.true.
        ENDIF

* numo2
* find and check the number of o2 reactions
!      CALL eznc_get_dimension(ncid,"maxo2",maxo2)
      CALL eznc_get_0Dint(ncid,"numo2",numo2)
        IF (numo2.GT.maxo2) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of o2 reactions '
          WRITE(6,*) '        greater than the parameter maxo2'
          WRITE(6,*) '        see parameter maxo2 (in akparameter)'
          lostop=.true.
        ENDIF

* numiso
* find and check the number of isomerisation reactions
!      CALL eznc_get_dimension(ncid,"maxiso",maxiso)
      CALL eznc_get_0Dint(ncid,"numiso",numiso)
        IF (numiso.GT.maxiso) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of isomerisation reactions '
          WRITE(6,*) '        greater than the parameter maxiso'
          WRITE(6,*) '        see parameter maxiso (in akparameter)'
          lostop=.true.
        ENDIF

* nummeo2
* find and check the number of ch3o2 reactions (RO2+CH3O2)
!      CALL eznc_get_dimension(ncid,"mxrpero",mxrpero)
      CALL eznc_get_0Dint(ncid,"nummeo2",nummeo2)
        IF (nummeo2.GT.mxrpero) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of ch3o2 reactions '
          WRITE(6,*) '        greater than the parameter mxrpero'
          WRITE(6,*) '        see parameter mxrpero (in akparameter)'
          WRITE(6,*) 'nummeo2=',nummeo2
          lostop=.true.
        ENDIF

* mxaux = NetCDF version of maxaux
* find the number of auxiliary information given
      CALL eznc_get_dimension(ncid,"maxaux",mxaux)
        IF (mxaux.NE.maxaux) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of different types '
          WRITE(6,*) '         of auxiliary information given by'
          WRITE(6,*) '         the chemical interpreter does not '
          WRITE(6,*) '         correspond to the number used in this'
          WRITE(6,*) '         program'
          WRITE(6,*) '         see parameter MAXAUX (in akparameter)'
          lostop=.true.
        ENDIF

!>>>> CALL eznc_get_1Dint(ncid,VarName,dim,values,start1,end1)
* (nauxpar(i),i=1:maxaux)
      CALL eznc_get_1Dint(ncid,"nauxpar",maxaux,
     $                     nauxpar(1:maxaux),1,maxaux)

* numain
! find and check the number of transfer reactions
!      CALL eznc_get_dimension(ncid,"maxt",maxt)
      CALL eznc_get_0Dint(ncid,"numain",numain)
        IF (numain.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of gas -> aero'
          WRITE(6,*) '         reactions is greater than maxt'
          WRITE(6,*) '         numain=',numain
          lostop=.true.
        ENDIF
* numaou
      CALL eznc_get_0Dint(ncid,"numaou",numaou)
        IF (numaou.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of  aero -> gas'
          WRITE(6,*) '         reactions is greater than maxt'
          WRITE(6,*) '         numaou=',numaou
          lostop=.true.
        ENDIF
* numwin
      CALL eznc_get_0Dint(ncid,"numwin",numwin)
        IF (numwin.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of gas -> wall'
          WRITE(6,*) '         reactions is greater than maxt'
          WRITE(6,*) '         numwin=',numwin
          lostop=.true.
        ENDIF
* numwou
      CALL eznc_get_0Dint(ncid,"numwou",numwou)
        IF (numwou.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of  wall -> gas'
          WRITE(6,*) '         reactions is greater than maxt'
          WRITE(6,*) '         numwou=',numwou
          lostop=.true.
        ENDIF

! Defined as DIMENSIONS in NetCDF link file.
! No need to interrogate since not used in this program unit.
* mx12stoi
* mx1stoi
* mx2stoi

* stop on error
        IF (lostop) THEN
          WRITE(6,*)
          WRITE(6,*) '--error exit--,   in subroutine SPAKINIT3'
          WRITE(6,*) '  make sure that the interpreter "inca"'
          WRITE(6,*) '  and your application program'
          WRITE(6,*) '  have been compiled using the same version'
          WRITE(6,*) '  of the parameter file (akparemeter.h)'
          WRITE(6,*)
          STOP
        ENDIF

* reading rest of link-file

* 1-D INTEGERS
!>>>> CALL eznc_get_1Dint(ncid,VarName,dim,values,start1,end1)
* (nrpero(k),k=1,maxro2), with size check
      CALL eznc_get_1Dint(ncid,"nrpero",maxro2,
     $                     nrpero(1:maxro2),1,maxro2)

      j=0
      DO k=1,maxro2
        IF (nrpero(k).gt.j) j=nrpero(k)
      ENDDO
      IF (j.gt.mxrpero) THEN
        WRITE(6,*) '--err--, number of RO2+RO2 reaction for a given '
        WRITE(6,*) '         type is greater than mxrpero'
        WRITE(6,*) '         see parameter MXRPERO (in akparameter)'
        WRITE(6,*) '         number found=',j
        WRITE(6,*) '         mxrpero=',mxrpero
      ENDIF
* (nrdimer(k),k=1,maxdimer), with size check
* (same as "numreacdimer" in spakinit9.f)
      CALL eznc_get_1Dint(ncid,"nrdimer",maxdimer,
     $                     nrdimer(1:maxdimer),1,maxdimer)

      j=0
      DO k=1,maxdimer
        IF (nrdimer(k).gt.j) j=nrdimer(k)
      ENDDO
      IF (j.gt.mxrdimer) THEN
        WRITE(6,*) '--err--, number of dimer reaction for a given '
        WRITE(6,*) '         type is greater that mxrdimer'
        WRITE(6,*) '         see parameter MXRDIMER (in akparameter)'
        WRITE(6,*) '         number found=',j
        WRITE(6,*) '         mxrdimer=',mxrdimer
      ENDIF
* (id_n(i),i=1,num_n)
      CALL eznc_get_1Dint(ncid,"id_n",maxre,
     $                     id_n(1:num_n),1,num_n)
* (id_m(i),i=1,num_m),
      CALL eznc_get_1Dint(ncid,"id_m",max_m,
     $                     id_m(1:num_m),1,num_m)
* (idhv(i),i=1,numhv)
      CALL eznc_get_1Dint(ncid,"idhv",maxhv,
     $                     idhv(1:numhv),1,numhv)
* (idcvar(i),i=1,numcvar),
      CALL eznc_get_1Dint(ncid,"idcvar",maxcvar,
     $                     idcvar(1:numcvar),1,numcvar)
* (idextra(i),i=1,numextra)
      CALL eznc_get_1Dint(ncid,"idextra",maxextra,
     $                     idextra(1:numextra),1,numextra)
* (ido2(i),i=1,numo2)
      CALL eznc_get_1Dint(ncid,"ido2",maxo2,
     $                     ido2(1:numo2),1,numo2)
* (idmeo2(i),i=1,nummeo2)
      CALL eznc_get_1Dint(ncid,"idmeo2",mxrpero,
     $                     idmeo2(1:nummeo2),1,nummeo2)
* (idiso(i),i=1,numiso)
      CALL eznc_get_1Dint(ncid,"idiso",maxiso,
     $                     idiso(1:numiso),1,numiso)
! read id number for mass transfer equation
* (idain(i),i=1,numain)
      CALL eznc_get_1Dint(ncid,"idain",maxt,
     $                     idain(1:numain),1,numain)
* (idaou(i),i=1,numaou)
      CALL eznc_get_1Dint(ncid,"idaou",maxt,
     $                     idaou(1:numaou),1,numaou)
* (idwin(i),i=1,numwin)
      CALL eznc_get_1Dint(ncid,"idwin",maxt,
     $                     idwin(1:numwin),1,numwin)
* (idwou(i),i=1,numwou)
      CALL eznc_get_1Dint(ncid,"idwou",maxt,
     $                     idwou(1:numwou),1,numwou)

* 2-D INTEGERS
* ((idfo(i,k),k=1,3),i=1,numfo),
      CALL eznc_get_2Dint(ncid,"idfo",maxfo,3,
     $                     idfo(1:numfo,1:3),
     $                          1,numfo,1,3)

* (numstoi(ire,k),k=1,2,ire=1,numre)
      CALL eznc_get_2Dint(ncid,"numstoi",maxre,2,
     $                     numstoi(1:numre,1:2),
     $                             1,numre,1,2)

* (idrestoi(ire,i),i=1,numstoi(ire,1)),ire=1,numre)
      idrestoi = 0
      call eznc_get_2Dint(ncid, "idrestoi",maxre, mxleft,
     $                     idrestoi(1:numre,1:mxleft),
     $                      1,numre,1,mxleft)

* (idpdstoi(ire,i),i=1,numstoi(ire,2)),ire=1,numre)
      idpdstoi = 0
       CALL eznc_get_2Dint(ncid,"idpdstoi",maxre,mxright,
     $                    idpdstoi(1:numre,1:mxright),
     $                    1,numre,1,mxright) 

* ((idreacdimer(i,k),i=1,nrdimer(k)),k=1,maxdimer)
      CALL eznc_get_2Dint(ncid,"idreacdimer",mxrdimer,maxdimer,
     $                   idreacdimer(1:mxrdimer,1:maxdimer),
     $                               1,mxrdimer,1,maxdimer)

* ((idreacro2(i,k),i=1,nrpero(k)),k=1,maxro2)
      CALL eznc_get_2Dint(ncid,"idreacro2",mxrpero,maxro2,
     $                   idreacro2(1:mxrpero,1:maxro2),
     $                             1,mxrpero,1,maxro2)

* 1-D REALS
* (hvfact(k),k=1,numhv)
      CALL eznc_get_1Dreal(ncid,"hvfact",maxhv,
     $                      hvfact(1:numhv),1,numhv)
* (hvcf(k),k=1,numhv)
      CALL eznc_get_1Dreal(ncid,"hvcf",maxhv,
     $                      hvcf(1:numhv),1,numhv)
* (cvarcf(k),k=1,numcvar)
      CALL eznc_get_1Dreal(ncid,"cvarcf",maxcvar,
     $                      cvarcf(1:numcvar),1,numcvar)
* (wmol(i),i=1,numsp)
      CALL eznc_get_1Dreal(ncid,"wmol",maxsp,
     $                      wmol(1:numsp),1,numsp)

* 2-D REALS
* ((arrhcf(ire,k),k=1,3),ire=1,numre)
!      PRINT*,"eznc_get_2Dreal: arrhcf"
      CALL eznc_get_2Dreal(ncid,"arrhcf", maxre,3,
     $                     arrhcf(1:numre,1:3),
     $                            1,numre,1,3)

* (restoicf(ire,i),i=1,numstoi(ire,1)),ire=1,numre)
      PRINT*,"eznc_get_2Dreal: restoicf"
      CALL eznc_get_2Dreal(ncid,"restoicf",maxre,mxleft,
     $                    restoicf(1:numre,1:mxleft),
     $                       1,numre,1,mxleft)
!      DO ire=1,numre
!        CALL eznc_get_2Dreal(ncid,"restoicf", maxre,mxleft,
!     $                        restoicf(ire,1:numstoi(ire,1)),
!     $                             ire,ire,1,numstoi(ire,1))
!      ENDDO

* (pdstoicf(ire,i),i=1,numstoi(ire,2)),ire=1,numre)
      PRINT*,"eznc_get_2Dreal: pdstoicf"
      CALL eznc_get_2Dreal(ncid,"pdstoicf",maxre,mxright,
     $                    pdstoicf(1:numre,1:mxright),
     $                       1,numre,1,mxright)
      
!      DO ire=1,numre
!        CALL eznc_get_2Dreal(ncid,"pdstoicf",maxre,mxright,
!     $                        pdstoicf(ire,1:numstoi(ire,2)),
!     $                             ire,ire,1,numstoi(ire,2))
!      ENDDO

* ((focf(i,k),i=1,maxaux+3),k=1,numfo)
      PRINT*,"eznc_get_2Dreal: focf"
      CALL eznc_get_2Dreal(ncid,"focf",maxaux+3,maxfo,
     $                      focf(1:maxaux+3,1:numfo),
     $                           1,maxaux+3,1,numfo)

* ((isocf(i,k),i=1,maxaux),k=1,numiso)
      PRINT*,"eznc_get_2Dreal: isocf"
      CALL eznc_get_2Dreal(ncid,"isocf",maxaux,maxiso,
     $                      isocf(1:maxaux,1:numiso),
     $                            1,maxaux,1,numiso)

! coefficients for mass transfer equations
*((aoucf(i,k),i=1,2),k=1,numaou)
!      PRINT*,"eznc_get_2Dreal: aoucf"
!      CALL eznc_get_2Dreal(ncid,"aoucf",2,maxt,
!     $                      aoucf(1:2,1:numaou),
!     $                            1,2,1,numaou)

*((woucf(i,k),i=1,3),k=1,numwou)
      PRINT*,"eznc_get_2Dreal: woucf"
      CALL eznc_get_2Dreal(ncid,"woucf",3,maxt,
     $                      woucf(1:3,1:numwou),
     $                            1,3,1,numwou)

*((wincf(i,k),i=1,3),k=1,numwin)
      PRINT*,"eznc_get_2Dreal: wincf"
      CALL eznc_get_2Dreal(ncid,"wincf",3,maxt,
     $                      wincf(1:3,1:numwin),
     $                            1,3,1,numwin)

* 2-D NF90_DOUBLE but written as "real"
* ((extracf(i,k),i=1,maxaux),k=1,numextra)
      PRINT*,"eznc_get_2Dreal: extracf"
      CALL eznc_get_2Dreal(ncid,"extracf",maxaux,maxextra,
     $                      extracf(1:maxaux,1:numextra),
     $                              1,maxaux,1,numextra)


* 1-D CHARACTER
* (chrsp(i),i=1,numsp)
      PRINT*,"eznc_get_1Dchar: chrsp"
      CALL eznc_get_1Dchar(ncid,"chrsp",maxlsp,maxsp,
     $                      chrsp,1,numsp)

* (dummyid_n,i=1,num_n)
!      ??? not in NetCDF file, not passed to main routine.
!      is apparently a placeholder for reading id_n, which is not 
!      passed either

* ----------------------
* get OFR_fg from mech, check compatible with keyfile

      PRINT*,"eznc_get_0Dint: OFR_fg",OFR_fg
      CALL eznc_get_0Dint(ncid,"OFR_fg",tmp_fg)

      IF(tmp_fg.NE.OFR_fg)THEN
        PRINT*,"eznc_get_0Dint: OFR_fg",tmp_fg
        WRITE(lout,*)'--error-- OFR_fg differs between'
        WRITE(lout,*)'mechanism and keyfile. STOPPING.'
        STOP
      ENDIF

* -------------------------------------------
* set or calculate same often used variables
* -------------------------------------------

* check that all falloff expressions follow a TROE expression
      DO i=1,numfo
        IF (idfo(i,2).ne.2) THEN
          WRITE(6,'(a38)') '--error-- this falloff case was not written'
          STOP
        ENDIF
      ENDDO

* Check that all stoichiometric coefficients in the reactant
* side of the reaction are only 1 or 2. If equal to 2, then
* store the info into a table. Purpose is to win time in the routine
* involved in time integration
      nself=0
      DO 100 ire=1,numre
        ierr=0
        DO j=1,numstoi(ire,1)
          IF (restoicf(ire,j).eq.1.) GOTO 100
          IF (restoicf(ire,j).eq.2.) THEN
            nself=nself+1
            ierr=ierr+1
            IF (nself.GT.mself) THEN
              WRITE(6,*) '--error--, in spakinit9.'
              WRITE(6,*) 'The number of self reaction exceed the'
              WRITE(6,*) 'maximum value allowed'
              WRITE(6,*) '(see mself in akparameter) '
              STOP
            ENDIF
            IF (ierr.gt.1) THEN
              WRITE(6,*) '--error--, in spakinit9.'
              WRITE(6,*) 'two reactant are reacting with a stoe=2'
              WRITE(6,*) '(i.e. the reaction is at least 4 orders)'
              WRITE(6,*) 'such a case is not allowed '
              STOP
            ENDIF
            idselfreac(nself,1)=ire
            idselfreac(nself,2)=j
          ELSE
            WRITE(6,*) '--error--, in spakinit9.'
            WRITE(6,*) 'stoichiometric coefficient at the reactant'
            WRITE(6,*) 'side is not 1 or 2 for the reaction :',ire
            WRITE(6,*) 'this is not allowed '
            STOP
          ENDIF
        ENDDO
100   CONTINUE

* ----------------------
* CHECK EXTRA REACTIONS
* ----------------------

!      WRITE(6,*) '     checking extra reactions ....'
c      cpnumextra=0
c      DO i=1,maxro2
c        nrpero(i)=0
c      ENDDO
c      DO i=1,maxextra
c        cpidextra(i)=0
c        DO j=1,maxaux
c          cpextracf(j,i)=0.
c        ENDDO
c      ENDDO


c      DO 201 i=1,numextra
c        ire=idextra(i)
c
c
* RO2+RO2 reactions
c        IF (nint(extracf(1,i)).EQ.401) THEN
c          nrpero(1)=nrpero(1)+1
c          idreacro2(nrpero(1),1)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.402) THEN
c          nrpero(2)=nrpero(2)+1
c          idreacro2(nrpero(2),2)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.403) THEN
c          nrpero(3)=nrpero(3)+1
c          idreacro2(nrpero(3),3)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.404) THEN
c          nrpero(4)=nrpero(4)+1
c          idreacro2(nrpero(4),4)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.405) THEN
c          nrpero(5)=nrpero(5)+1
c          idreacro2(nrpero(5),5)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.406) THEN
c          nrpero(6)=nrpero(6)+1
c          idreacro2(nrpero(6),6)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.407) THEN
c          nrpero(7)=nrpero(7)+1
c          idreacro2(nrpero(7),7)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.408) THEN
c          nrpero(8)=nrpero(8)+1
c          idreacro2(nrpero(8),8)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.409) THEN
c          nrpero(9)=nrpero(9)+1
c          idreacro2(nrpero(9),9)=ire
c          GOTO 201
c
c* other type of extra reaction
c        ELSE
c          cpnumextra=cpnumextra+1
c          cpidextra(cpnumextra)=ire
c          DO j=1,maxaux
c            cpextracf(j,cpnumextra)=extracf(j,i)
c          ENDDO
c        ENDIF
c201   CONTINUE
c      DO i=1,maxro2
c        WRITE(6,*) 'nrpero ', i,' ',nrpero(i)
c      ENDDO
c      DO i=1,9
c        DO j=1,nrpero(i)
c         WRITE(90+i,'(i8)') idreacro2(j,i)
c        ENDDO
c      ENDDO
c      STOP
* restore the table for extra reaction
c      numextra=cpnumextra
c      DO i=1,numextra
c        idextra(i)=cpidextra(i)
c        DO j=1,maxaux
c          extracf(j,i)=cpextracf(j,i)
c        ENDDO
c      ENDDO


* ----------------------
* CHECK HV REACTIONS
* ----------------------

* check the reactions with hv. To decrease CPU time, the work is
* done on the chromophore only (instead of the reaction). Here, the
* various chromophores are defined based on the reaction dataset.
* Many chromophores are used only once (e.g. NO2 photolysis). Few
* chromophores are used for a large number of reactions
* (e.g. ROOH photolysis). To save memory space of chromophores,
* 3 tables of chromophores are considered.

* initialize the table
* ------------------------------
      ntchromo=0
      mxlab=0
      DO i=1,mchromo
        chromocf(i)=0
        numchromo(i)=0
      ENDDO
      DO i=1,numhv
        cphvcf(i)=nint(hvcf(i))
      ENDDO

      DO i=1,mtopchromo
        chromotopcf(i)=0
        numtopchromo(i)=0
        DO j=1,msptopchromo
          idtopchromo(i,j)=0
        ENDDO
      ENDDO

      DO i=1,mmedchromo
        chromomedcf(i)=0
        nummedchromo(i)=0
        DO j=1,mspmedchromo
          idmedchromo(i,j)=0
        ENDDO
        ntmedchromo=0
      ENDDO

      nt1chromo=0
      DO i=1,mchromo
        chromo1cf(i)=0
        id1chromo(i)=0
      ENDDO

* loop over the hv reaction to get the number of reaction hidden below
* each label (table numchromo), the total number of chromophore
* (ntchromo) and the label table (chromocf)
* ------------------------------
      DO 400 i=1,numhv
        idcf=cphvcf(i)

* check if chromo already exist. If yes goto next (after storing data)
        DO j=1,ntchromo
          IF (idcf.eq.chromocf(j)) THEN
            numchromo(j)=numchromo(j)+1
            GOTO 400
          ENDIF
        ENDDO

* if that point is reached, chromo does not exist => add to the list
        ntchromo=ntchromo+1
        IF (ntchromo.gt.mchromo) THEN
          WRITE(6,*) '--error--, number of chromophore (i.e. labels)'
          WRITE(6,*) '           exceed mchromo. Change akparameter'
          STOP
        ENDIF
        chromocf(ntchromo)=idcf
        numchromo(ntchromo)=1
400   CONTINUE


* get the most used label (mtop)
* ------------------------------
      DO j=1,mtopchromo

* identify the label
        mxlab=0
        DO i=1,ntchromo
          IF (numchromo(i).gt.mxlab) THEN
            mxlab=numchromo(i)
            label=chromocf(i)
            index=i
          ENDIF
        ENDDO
        IF (mxlab.gt.msptopchromo) THEN
          WRITE(6,*) '--error--, the number of species in the'
          WRITE(6,*) '           "top" table of chromophore'
          WRITE(6,*) '           exceed size nsptopchromo'
          WRITE(6,*) '           change akparameter accordingly'
          STOP
        ENDIF

* store the reaction number i into idtopchromo and the label into
* chromotopcf
        k=0
        DO i=1,numhv
          IF (cphvcf(i).eq.label) THEN
            k=k+1
            idtopchromo(j,k)=i
            cphvcf(i)=0
          ENDIF
        ENDDO
        chromotopcf(j)=label
        numtopchromo(j)=mxlab
        IF (mxlab.NE.k) STOP '-error1- in spakinit.f, numbers not ok'

* erase the label and get next label
        numchromo(index)=0
      ENDDO


* get the label used only once
* ------------------------------
      k=0
      DO 500 j=1,ntchromo
        IF (numchromo(j).eq.1) THEN
          numchromo(j)=0
          label=chromocf(j)
          k=k+1
          DO i=1,numhv
            IF (cphvcf(i).eq.label) THEN
              id1chromo(k)=i
              chromo1cf(k)=chromocf(j)
              cphvcf(i)=0
              GOTO 500
            ENDIF
          ENDDO
        ENDIF
500   CONTINUE
      nt1chromo=k


* store other labels
* ------------------------------
      ntmedchromo=0
      DO j=1,ntchromo
        IF (numchromo(j).ne.0) THEN
          mxlab=numchromo(j)
          IF (mxlab.gt.mspmedchromo) THEN
            WRITE(6,*) 'mxlab, mspmed chromo ', mxlab,mspmedchromo
            WRITE(6,*) '--error--, the number of species in the'
            WRITE(6,*) '           "med" table of chromophore'
            WRITE(6,*) '           exceed size mspmedchromo'
            WRITE(6,*) '           change akparameter accordingly'
            STOP
          ENDIF

          numchromo(j)=0
          label=chromocf(j)

          ntmedchromo=ntmedchromo+1
          IF (ntmedchromo.gt.mmedchromo) THEN
            WRITE(6,*) 'ntmedchromo,mmedchromo ',ntmedchromo,mmedchromo
            WRITE(6,*) '--error--, the number of ID for chomophore'
            WRITE(6,*) '           in "med" table of chromophore'
            WRITE(6,*) '           exceed size mspmedchromo'
            WRITE(6,*) '           change akparameter accordingly'
            STOP
          ENDIF
          chromomedcf(ntmedchromo)=label
          nummedchromo(ntmedchromo)=mxlab

          k=0
          DO i=1,numhv
            IF (cphvcf(i).eq.label) THEN
              k=k+1
              idmedchromo(ntmedchromo,k)=i
              cphvcf(i)=0
            ENDIF
          ENDDO
          IF (mxlab.NE.k) STOP '-err2- in spakinit.f, numbers not ok'
        ENDIF
      ENDDO

! set up parameters for all j-value output
* -----------------------------------------
      nchrom = nt1chromo+ntmedchromo+mtopchromo

      i = 1 
      k = mtopchromo
      idchrom(i:k) = chromotopcf(1:mtopchromo)
      i = k+1 
      k = mtopchromo+ntmedchromo
      idchrom(i:k) = chromomedcf(1:ntmedchromo)
      i = k+1 
      k = mtopchromo+ntmedchromo+nt1chromo
      idchrom(i:k) = chromo1cf(1:nt1chromo)

* check that every hv reaction has been treated
* ------------------------------
      DO i=1,numhv
        IF (cphvcf(i).ne.0) THEN
          WRITE(6,*) '--error--, the  HV reaction was not treated'
          WRITE(6,*) '           label=',cphvcf(i)
          STOP
        ENDIF
      ENDDO

      mxlab=mtopchromo+ntmedchromo+nt1chromo

! DEBUG
      print*,"mxlab = ",mxlab
      print*,"ntchromo = ",ntchromo
      IF (mxlab.NE.ntchromo) STOP '-error-, hv reactions are missing'

* ----------------------

      END

* -------------------------------
