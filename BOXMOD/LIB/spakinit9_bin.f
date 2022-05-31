******************************************************************
* SUBROUTINE spakinit9 : this routine reads the chemical scheme
* produced by the interpreter "inca".
*
* INPUT : nothing
* OUTPUT :
* (comment : all the variables that were initially in the akcommon.h
* are now given as output of this routine : need to make the
* comment in the future)
******************************************************************
      SUBROUTINE spakinit9_bin(numsp,numre,num_n,num_m,numfo,numhv,
     2                     numcvar,numextra,numo2,nummeo2,numreacro2,
     3                     numiso,idiso,
     3                     mx12stoi,nauxpar,numstoi,idrestoi,idpdstoi,
     4                     id_m,idfo,
     5                     idhv,idcvar,idextra,ido2,idmeo2,idreacro2,
     5                     idreacdimer,numreacdimer,
     6                     nself,idselfreac,
     7                nt1chromo,chromo1cf,id1chromo,
     7                ntmedchromo,chromomedcf,nummedchromo,idmedchromo,
     7                chromotopcf,numtopchromo,idtopchromo,
     8                     restoicf,pdstoicf,arrhcf,focf,
     9                     hvcf,hvfact,cvarcf,extracf,isocf,
     1                     numain, numaou, numwin, numwou,
     1                     idain,idaou,idwin,idwou,
     1                     aoucf,woucf,wincf,
     o                     wmol,chrsp) !,
!     &            tralphain,tralphaout,
!     &                trdeltahin,trdeltahout,
!     &            trhenryin,trhenryout,
!     &            numtrin,numtrout,
!     &            idtrin,idtrout,
!     &            hydcf,numhyd,idhyd,
!     &            ka,numacid,idacid,
!     &            ohaq,numohaq,idohaq,
!     &            idreacro2_aq,nrpero_aq)
      !$ use OMP_LIB
      USE akparameter_module
      USE module_data_gecko_main,ONLY: nchrom,idchrom
      IMPLICIT NONE

      CHARACTER(maxlsp) chrsp(maxsp)

      INTEGER  numsp, numre, num_n
      INTEGER  num_m, numfo, numhv, numcvar, numextra
      INTEGER  mx12stoi
      INTEGER  numain, numaou, numwin, numwou
      INTEGER  nauxpar(maxaux)
      INTEGER  numstoi(maxre,2)
      INTEGER  idrestoi(maxre,mxleft),idpdstoi(maxre,mxright)
      INTEGER  dummyid_n, id_m(max_m)
      INTEGER  idfo(maxfo,3)
      INTEGER  idhv(maxhv), idcvar(maxcvar)
      INTEGER  idextra(maxextra)
      INTEGER  nself,idselfreac(mself,2)
      INTEGER  numo2,ido2(maxo2)
      INTEGER  numiso,idiso(maxiso)
      INTEGER  nummeo2,idmeo2(mxrpero)
      INTEGER  numreacro2(maxro2),idreacro2(mxrpero,maxro2)
      INTEGER  numreacdimer(maxdimer),idreacdimer(mxrdimer,maxdimer)
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
! USE REAL UNDER COMPILATION OPTION real-8
      !DOUBLE PRECISION     extracf(maxaux,maxextra)
      REAL     isocf(maxaux,maxiso)
      REAL     wmol(maxsp)
      REAL     aoucf(2,maxt),woucf(3,maxt),wincf(3,maxt)

!      REAL tralphain(maxtr),tralphaout(maxtr)
!      REAL trdeltahin(maxtr),trdeltahout(maxtr)
!      REAL trhenryin(maxtr),trhenryout(maxtr)
!      INTEGER numtrin,numtrout
!      INTEGER idtrin(maxtr),idtrout(maxtr)

!      REAL hydcf(maxhyd)
!      INTEGER numhyd,idhyd(maxhyd)
!      REAL ka(maxacid)
!      INTEGER numacid,idacid(maxacid)

!      REAL ohaq(mxohaq,3)
!      INTEGER  numohaq,idohaq(mxohaq)
!      INTEGER idreacro2_aq(mxrpero,maxro2)
!      INTEGER nrpero_aq(mxrpero)

* internal
      INTEGER  ntchromo,chromocf(mchromo),numchromo(mchromo)
      INTEGER  label, mxlab, cphvcf(maxhv)
      LOGICAL  lostop
      INTEGER  ierr, i, j, mxlsp, k, ire, mxaux, llink, index
      INTEGER  idcf
      INTEGER  mx1stoi, mx2stoi

c      INTEGER  cpnumextra, cpidextra(maxextra)
c      REAL     cpextracf(maxaux,maxextra)

      lostop=.false.
      ierr=0

      print*,"in spakinit9"
*open the chemical file
      llink=12
      OPEN(llink,file='indat.li',  form='unformatted',status='old')

* -----------------------------
* initialization of all values
* -----------------------------

* initialization of the strings
      DO i=1,maxsp
        chrsp(i)=' '
      ENDDO

* initialization of the integer variable
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
      mx12stoi=0
      mx1stoi=0
      mx2stoi=0

      nauxpar(:)=0
      numstoi(:,:)=0
      idrestoi(:,:)=0
      idpdstoi(:,:)=0
      dummyid_n=0
      id_m(:)=0
      idfo(:,:)=0
      idhv(:)=0
      idcvar(:)=0
      idextra(:)=0
      ido2(:)=0
      idmeo2(:)=0

! initialization of real variables
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
* reading the chemical scheme (binary file)
* ------------------------------------------

* read some key input values and check that the sizes are ok
      READ(llink) mxlsp, numsp, numre, num_n,
     &            num_m,numfo,numhv,numcvar,numextra,numo2,nummeo2,
     &            numain,numaou,numwin,numwou,
     &            numiso,mx12stoi, mx1stoi, mx2stoi,
!     &            numtrin, numtrout, numhyd, numohaq, numacid,
     &            mxaux,
     &           (nauxpar(i),i=1,maxaux)

* check for the length of the species
        IF (mxlsp.GT.maxlsp) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, length of species names given by the '
          WRITE(6,*) '         chemical interpreter does not correspond'
          WRITE(6,*) '         to the length used in this  program'
          WRITE(6,*) '         see parameter MAXLSP (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of species
        IF (numsp.GT.maxsp) THEN
          WRITE(6,*)
          WRITE(6,*) '--error--, number of species is larger '
          WRITE(6,*) '           than the parameter MAXSP'
          WRITE(6,*) '           see parameter MAXSP (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of reaction
        IF (numre.GT.maxre) THEN
          WRITE(6,*)
          WRITE(6,*) '--error--, number of reaction is larger '
          WRITE(6,*) '           than the parameter MAXRE'
          WRITE(6,*) '           see parameter MAXRE (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of reaction with "M"
        IF (num_m.GT.max_m) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of third body reaction is larger '
          WRITE(6,*) '         than the parameter MAX_M'
          WRITE(6,*) '         see parameter MAX_M (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of fall_off reactions
        IF (numfo.GT.maxfo) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of falloff reaction is larger '
          WRITE(6,*) '         than the parameter MAXFO'
          WRITE(6,*) '         see parameter MAXFO (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of photolysis reactions
        IF (numhv.GT.maxhv) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of hv reaction is larger '
          WRITE(6,*) '         than the parameter MAXHV'
          WRITE(6,*) '         see parameter MAXHV (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of cvar reactions
        IF (numcvar.GT.maxcvar) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of cvar reaction is larger '
          WRITE(6,*) '         than the parameter MAXCVAR'
          WRITE(6,*) '         see parameter MAXCVAR (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of extra reactions
        IF (numextra.GT.maxextra) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of extra reaction '
          WRITE(6,*) '         than the parameter MAXEXTRA'
          WRITE(6,*) '        see parameter MAXEXTRA (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of o2 reactions
        IF (numo2.GT.maxo2) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of o2 reaction '
          WRITE(6,*) '        greater than the parameter maxo2'
          WRITE(6,*) '        see parameter maxo2 (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of isomerisation reactions
        IF (numiso.GT.maxiso) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of isomerisation reaction '
          WRITE(6,*) '        greater than the parameter maxiso'
          WRITE(6,*) '        see parameter maxiso (akparameter mod)'
          lostop=.true.
        ENDIF

* check the number of ch3o2 reactions (RO2+CH3O2)
        IF (nummeo2.GT.mxrpero) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, number of ch3o2 reaction '
          WRITE(6,*) '        greater than the parameter mxrpero'
          WRITE(6,*) '        see parameter mxrpero (akparameter mod)'
          WRITE(6,*) 'nummeo2=',nummeo2
          lostop=.true.
        ENDIF

* check the maximum number of species in a reaction
c        IF (mx3stoi.GT.maxstoi) THEN
c          WRITE(6,*)
c          WRITE(6,*) '--err--, number of reaction partners is '
c          WRITE(6,*) '         larger than the parameter MAXSTOI'
c          WRITE(6,*) '         see parameter MAXSTOI (akparameter mod)'
c          lostop=.true.
c        ENDIF

* check the number of auxiliary information given
        IF (mxaux.NE.maxaux) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of different types '
          WRITE(6,*) '         of auxiliary information given by'
          WRITE(6,*) '         the chemical interpreter does not '
          WRITE(6,*) '         correspond to the number used in this'
          WRITE(6,*) '         program'
          WRITE(6,*) '         see parameter MAXAUX (akparameter mod)'
          lostop=.true.
        ENDIF

! check the number of tranfer reactions
        IF (numain.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of gas -> aero'
          WRITE(6,*) '         reaction is greater than maxt'
          WRITE(6,*) '         numain=',numain
          lostop=.true.
        ENDIF
        IF (numaou.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of  aero -> gas'
          WRITE(6,*) '         reaction is greater than maxt'
          WRITE(6,*) '         numaou=',numaou
          lostop=.true.
        ENDIF
        IF (numwin.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of gas -> wall'
          WRITE(6,*) '         reaction is greater than maxt'
          WRITE(6,*) '         numwin=',numwin
          lostop=.true.
        ENDIF
        IF (numwou.gt.maxt) THEN
          WRITE(6,*)
          WRITE(6,*) '--err--, the number of  wall -> gas'
          WRITE(6,*) '         reaction is greater than maxt'
          WRITE(6,*) '         numwou=',numwou
          lostop=.true.
        ENDIF

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
      READ(llink) (chrsp(i),i=1,numsp)
      READ(llink) (dummyid_n,i=1,num_n),(id_m(i),i=1,num_m),
     &            ((idfo(i,k),k=1,3),i=1,numfo),
     &            (idhv(i),i=1,numhv),(idcvar(i),i=1,numcvar),
     &            (idextra(i),i=1,numextra)!,(idtrin(i),i=1,numtrin),
!     &            (idtrout(i),i=1,numtrout),(idhyd(i),i=1,numhyd),
!     &            (idohaq(i),i=1,numohaq),(idacid(i),i=1,numacid)


      READ(llink) (ido2(i),i=1,numo2)
      READ(llink) (idmeo2(i),i=1,nummeo2)
      READ(llink) (idiso(i),i=1,numiso)

! read id number for mass transfer equation
      READ(llink) (idain(i),i=1,numain)
      READ(llink) (idaou(i),i=1,numaou)
      READ(llink) (idwin(i),i=1,numwin)
      READ(llink) (idwou(i),i=1,numwou)

      READ(llink)((arrhcf(ire,k),k=1,3),ire=1,numre)
      READ(llink)((numstoi(ire,k),k=1,2),ire=1,numre)

      READ(llink)((idrestoi(ire,i),restoicf(ire,i),i=1,numstoi(ire,1)),
     &            (idpdstoi(ire,i),pdstoicf(ire,i),i=1,numstoi(ire,2)),
     &             ire=1,numre)

      READ(llink) ((focf(i,k),i=1,maxaux+3),k=1,numfo)
      READ(llink) (hvcf(k),k=1,numhv)
      READ(llink) (hvfact(k),k=1,numhv)
      READ(llink) (cvarcf(k),k=1,numcvar)
      READ(llink) ((extracf(i,k),i=1,maxaux),k=1,numextra)
      READ(llink) ((isocf(i,k),i=1,maxaux),k=1,numiso)

! read coefficients for mass transfer equations
      READ(llink)((aoucf(i,k),i=1,2),k=1,numaou)
      READ(llink)((woucf(i,k),i=1,3),k=1,numwou)
      READ(llink)((wincf(i,k),i=1,3),k=1,numwin)

* read data for RO2+RO2. Check that enougth size has been given
      READ(llink) (numreacro2(k),k=1,maxro2)
      j=0
      DO k=1,maxro2
        IF (numreacro2(k).gt.j) j=numreacro2(k)
      ENDDO
      IF (j.gt.mxrpero) THEN
        WRITE(6,*) '--err--, number of RO2+RO2 reaction for a given '
        WRITE(6,*) '         type is greater that mxrpero'
        WRITE(6,*) '         see parameter MAXAUX (akparameter mod)'
        WRITE(6,*) '         number found=',j
        WRITE(6,*) '         mxrpero=',mxrpero
      ENDIF
      READ(llink) ((idreacro2(i,k),i=1,numreacro2(k)),k=1,maxro2)

* read data for dimerisation. Check that enougth size has been given
      READ(llink) (numreacdimer(k),k=1,maxdimer)
      j=0
      DO k=1,maxdimer
        IF (numreacdimer(k).gt.j) j=numreacdimer(k)
      ENDDO
      IF (j.gt.mxrdimer) THEN
        WRITE(6,*) '--err--, number of dimer reaction for a given '
        WRITE(6,*) '         type is greater that mxrdimer'
        WRITE(6,*) '         see parameter MAXAUX (akparameter mod)'
        WRITE(6,*) '         number found=',j
        WRITE(6,*) '         mxrdimer=',mxrdimer
      ENDIF
      READ(llink) ((idreacdimer(i,k),i=1,numreacdimer(k)),k=1,maxdimer)
!      READ(llink) (trhenryin(i),i=1,numtrin),
!     &           (tralphain(i),i=1,numtrin),(trdeltahin(i),i=1,numtrin)
!      READ(llink)(trhenryout(i),i=1,numtrout),
!     &       (tralphaout(i),i=1,numtrout),(trdeltahout(i),i=1,numtrout)

!      READ(llink) (hydcf(i),i=1,numhyd)

!      READ(llink) ((ohaq(i,k),k=1,3),i=1,numohaq)
!      READ(llink) (nrpero_aq(k),k=1,3)
!      READ(llink) ((idreacro2_aq(i,k),i=1,nrpero_aq(k)),k=1,3)

!      READ(llink) (ka(i),i=1,numacid)

* reading the molecular weight
      READ(llink)(wmol(i),i=1,numsp)

* close the file
      CLOSE(llink)

* -------------------------------------------
* set or calculate same often used variables
* -------------------------------------------

* check that all falloff expression follow a TROE expression
      DO i=1,numfo
        IF (idfo(i,2).ne.2) THEN
          WRITE(6,'(a38)') '--error-- this falloff case was not written'
          STOP
        ENDIF
      ENDDO

* Check that all stochiometrique coefficient in the reactant
* side of the reaction are only 1 or 2. If equal 2, then
* store the info into a table. Purpose is to win time in the routine
* involded in time integration
      nself=0
!$OMP PARALLEL DO private(ire, j, ierr)
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
              WRITE(6,*) '(see mself akparameter mod) '
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
!$OMP END PARALLEL DO

* ----------------------
* CHECK EXTRA REACTION
* ----------------------

      WRITE(6,*) '     checking extra reaction ....'
c      cpnumextra=0
c      DO i=1,maxro2
c        numreacro2(i)=0
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
c          numreacro2(1)=numreacro2(1)+1
c          idreacro2(numreacro2(1),1)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.402) THEN
c          numreacro2(2)=numreacro2(2)+1
c          idreacro2(numreacro2(2),2)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.403) THEN
c          numreacro2(3)=numreacro2(3)+1
c          idreacro2(numreacro2(3),3)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.404) THEN
c          numreacro2(4)=numreacro2(4)+1
c          idreacro2(numreacro2(4),4)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.405) THEN
c          numreacro2(5)=numreacro2(5)+1
c          idreacro2(numreacro2(5),5)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.406) THEN
c          numreacro2(6)=numreacro2(6)+1
c          idreacro2(numreacro2(6),6)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.407) THEN
c          numreacro2(7)=numreacro2(7)+1
c          idreacro2(numreacro2(7),7)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.408) THEN
c          numreacro2(8)=numreacro2(8)+1
c          idreacro2(numreacro2(8),8)=ire
c          GOTO 201
c        ELSE IF (nint(extracf(1,i)).EQ.409) THEN
c          numreacro2(9)=numreacro2(9)+1
c          idreacro2(numreacro2(9),9)=ire
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
c        WRITE(6,*) 'numreacro2 ', i,' ',numreacro2(i)
c      ENDDO
c      DO i=1,9
c        DO j=1,numreacro2(i)
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
* CHECK HV REACTION
* ----------------------

* check the reaction with hv. To decrease CPU time, the work is
* made on the chormophore only (instead of the reaction). Here, the
* various chromophore are defined based on the reaction dataset.
* Many chromophore are used only once (e.g. NO2 photolysis). Few
* chromophores are used for a large number of reactions
* (e.g. ROOH photolysis). To save memory space of chromophore, 3
* tables of chromophore are considered.

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
          WRITE(6,*) '--error--, the  HV reaction was not traeted'
          WRITE(6,*) '           label=',cphvcf(i)
          STOP
        ENDIF
      ENDDO

      mxlab=mtopchromo+ntmedchromo+nt1chromo
! DEBUG
      print*,"mxlab = ",mxlab
      print*,"ntchromo = ",ntchromo
      IF (mxlab.NE.ntchromo) STOP '-error-, hv reaction are missing'
      
* -------------------------
* write some data
* -------------------------


c      write(6,*) 'numhv=',numhv
c      write(6,*) 'ntchromo=',ntchromo
c      j=0
c      do i=1,ntchromo
c         j=j+numchromo(i)
c      ENDDO
c      write(6,*) 'total stock=',j
c
c      DO i=1,ntchromo
c       write(47,*) 'chromocf= ',chromocf(i), ' numchromo= ',numchromo(i)
c        DO j=1,numchromo(i)
c          write(47,*) '    ',j,idchromo(i,j),idhv(idchromo(i,j))
c        ENDDO
c      ENDDO
c      STOP

c      write(6,*) 'numextra=',numextra
c      write(6,'(10(i7))') (numreacro2(i),i=1,maxaux)
c      do i=1,numreacro2(1)
c        write(47,*) i, idreacro2(i,1)
c      enddo
c      STOP

c      do i=1,numextra
c        write(48,'(i8,i8,10(E11.3))')
c     &      i, idextra(i),(extracf(k,i),k=1,maxaux)
c      enddo
c      stop

c      write(6,*) 'self=',nself
c      DO i=1,nself
c       write(6,*) 'i=',i
c       write(6,*) 'selfreac(i,1)=',idselfreac(i,1)
c       write(6,*) 'selfreac(i,2)=',idselfreac(i,2)
c      ENDDO
c      STOP

      END
