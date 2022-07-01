      PROGRAM inca

      USE sorting, ONLY : sort_species, srtid_chrsp
      USE akparameter_module
      IMPLICIT NONE
      INCLUDE 'general.h' 
      INTEGER linedim, lenread, lenline
      INTEGER maxall
      PARAMETER (linedim=2, lenread=160, lenline=linedim*lenread)
      PARAMETER (maxall = maxre+maxsp)

      INTEGER lin, lout
      INTEGER llink
      INTEGER isp, i, ii, ipos, k, j, is
      INTEGER icvar, nlines, ire
      INTEGER numsp, numre, num_n, numo2, nummeo2,numiso
      INTEGER num_m, numfo, numhv, numcvar, numextra
      INTEGER mx12stoi, mx1stoi, mx2stoi, mx3stoi
      INTEGER numain, numaou, numwin, numwou

      INTEGER nauxpar(maxaux)
      INTEGER numstoi(maxre,2)
c      INTEGER idstoi(maxre,maxstoi,2)
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
      INTEGER itype(maxre)
      !INTEGER nrpero(mxrpero)
      INTEGER nrpero(maxro2)
      INTEGER idreacro2(mxrpero,maxro2)
      !INTEGER nrdimer(mxrdimer)
      INTEGER nrdimer(maxdimer)
      INTEGER idreacdimer(mxrdimer,maxdimer)

      REAL focf(maxaux+3,maxfo)
      REAL extracf(maxaux,maxextra)
!      DOUBLE PRECISION extracf(maxaux,maxextra)
      REAL isocf(maxaux,maxiso)
      REAL xauxcf(0:maxaux,maxaux)
!      DOUBLE PRECISION xauxcf(0:maxaux,maxaux)
c      REAL stoicf(maxre,maxstoi,2)
      REAL restoicf(maxre,mxleft)
      REAL pdstoicf(maxre,mxright)
      REAL arrhcf(maxre,3)
      REAL hvcf(maxhv), hvfact(maxhv), cvarcf(maxcvar)
      REAL wmol(maxsp)
! note that AIN reactions do not have additional coeffs, so AINCF not needed.
! DITTO AOU
      REAL woucf(1,maxt),wincf(1,maxt)

      CHARACTER*(maxlsp)  chrsp(maxsp)
      CHARACTER*(lenline) line
      CHARACTER*(lenread) inline,outline(linedim)

      LOGICAL loreply, lostop, locheck(maxall)
      LOGICAL lo_m, lofo, lohv, loaux, locvar, loextra, lo_o2, lopero
      LOGICAL loain,loaou,lowin,lowou
      LOGICAL lo_meo2,lo_iso,lodimer
      INTEGER ipero,idimer

      EQUIVALENCE (line,outline)
      INTEGER lenstr

      REAL            etime, duree, tdat(2)
      CHARACTER*10    date,begtime,endtime
      REAL            telstd
      REAL            time_diff

      INTEGER  ntchromo,chromocf(mchromo),numchromo(mchromo),idcf
      INTEGER  mxlab
! flag for i/o format (0:NetCDF; 1:binary; 2:in=NetCDF, out=binary)
      INTEGER :: iofmt_fg = 0
* ----------
* INITIALIZE
* ----------

* check that mxleft is lesser than mxright (required for getreac)
      IF (mxright.lt.mxleft) THEN
        WRITE(6,*) '--error--, check akparameter_module'
        WRITE(6,*) 'mxright must be greater than mxleft'
        STOP
      ENDIF

* loreply write data interpreted in the *.akoi output file
* set loreply to true if writing needed
      loreply=.false.

      telstd=0.
      DO isp=1,maxsp
        chrsp(isp)=' '
      ENDDO

      lostop =.false.
      loaux  =.false.
      lo_m   =.false.
      lofo   =.false.
      lohv   =.false.
      locvar =.false.
      loextra=.false.
      lo_o2  =.false.
      lo_meo2=.false.
      lopero =.false.
      lodimer =.false.
      lo_iso =.false.
      loain=.false.
      loaou=.false.
      lowin=.false.
      lowou=.false.
      line=' '
      lin=15
      lout=16
      llink=11

* initialise numre,numfo, numcvar, num_m,numextra, numhv
      numre=0
      numfo=0
      numcvar=0
      num_m=0
      numextra=0
      numhv=0
      numo2=0
      numiso=0
      nummeo2=0
      numsp=0
      numain=0
      numaou=0
      numwin=0
      numwou=0
      mx12stoi=0
      mx1stoi=0
      mx2stoi=0
      mx3stoi=0
      nauxpar=0
      xauxcf=0.
      nrpero=0
      nrdimer=0
      idrestoi=0
      idpdstoi=0

* Special reaction are given using a keyword (EXTRA, HV, ...).
* For each keyword, a given number of data must be read (data are
* given after the each special reaction between slahes). The
* number of data that must be read for the various cases are given
* below. Old keyword "SRI" (iidaux=3) and "LT" (iidaux=4) are not
* used anymore (and number are free if additional keyword are required).

* keyword  FALLOFF (previously LOW and TROE)   (iidaux=1)
      nauxpar(1)=7

* keyword  (previously TROE used 2 lines)
      nauxpar(2)=0

* a virer apres debuggage
      nauxpar(3)=5
      nauxpar(4)=2

* keyword  HV    (iidaux=5)
      nauxpar(5)=2

* keyword  EXTRA (iidaux=6)
      nauxpar(6)=10

* keyword  CVAR  (iidaux=7)
      nauxpar(7)=1

* keyword  AOU  (iidaux=8)
      nauxpar(8)=0

* keyword  WOU  (iidaux=9)
      nauxpar(9)=1

* keyword  WIN  (iidaux=10)
      nauxpar(10)=1

* keyword  ISOM (iidaux=11)
      nauxpar(11)=5

* ---------------------------
* OPEN INPUT AND OUTPUT FILES
* ---------------------------

      OPEN(lout ,file='outdat.akoi')

* ----------------------------------
* READ THE LISTS OF SPECIES & MASSES
* ----------------------------------
      DO i=1,3
      line=' '

      !------------------------------
      SELECT CASE (i)

        CASE (1)
          OPEN(lin  ,file='indat.sp_gas',status='OLD')
        CASE (2)
          OPEN(lin  ,file='indat.sp_part',status='OLD')
        CASE (3)
          OPEN(lin  ,file='indat.sp_wall',status='OLD')
      END SELECT

* search for 'PHASE' - first occurrence
1000  CONTINUE
      READ(lin,'(a)',end=9100)inline
      IF (loreply) WRITE(lout,'(a)')inline
      CALL cleanline(inline,outline(1),lenread)

      IF (lenstr(line,lenline).EQ.0) GOTO 1000
      IF (index(line,"!").EQ.1) GOTO 1000
      IF (index(line,"SPECIES").NE.0) GOTO 1000
      IF (index(line,"PHASE").NE.1) GOTO 9210

* read species list => When all species read => jump to 1300
      SELECT CASE (i)
        CASE (1)
          WRITE(6,*) 'reading gas species ...'
        CASE (2)
          WRITE(6,*) 'reading particle species ...'
        CASE (3)
          WRITE(6,*) 'reading wall species ...'
      END SELECT

1110  CONTINUE
      line=' '
      READ(lin,'(a)',END=9220) inline
      IF (loreply) WRITE(lout,'(a)') inline
      CALL cleanline(inline,outline(1),lenread)

* search for 'PHASE' - 2nd occurrence
      !ipos=index(line,'END')
      ipos=index(line,'PHASE')
      IF (ipos.eq.1) GOTO 1300

* get species names
      CALL getspec (
     1   line, lenline, lout,
     2   maxsp, maxlsp, numsp,
     3   lostop,
     4   wmol, chrsp)
      GOTO 1110

1300  CONTINUE
      CLOSE(lin)
      ENDDO

      ! created sorted index of species array
      write(6,*) 'sort the list of species ...'
      call sort_species(chrsp)
      write(6,*) '      end of sort'
      WRITE(6,*) numsp

* check species names
      WRITE(6,*) 'checking the species ...'
      CALL chkspec (
     1   lout,
     2   maxsp, numsp, maxlsp,
     3   wmol, chrsp,
     4   loreply, lostop)

*  stop if error found
      IF (lostop) THEN
        WRITE(lout,*)
        WRITE(lout,*)' -- error exit --  before reading reaction data'
        WRITE(lout,*)
        CLOSE(lout)
        STOP ' ERROR : in the list of species => check file '
      ENDIF

* ---------------------------
* READ THE LIST OF REACTIONS
* ---------------------------
      OPEN(lin  ,file='indat.mech',status='OLD')

* set all locheck(i) to .false.
      DO i=1,maxre
        locheck(i)=.false.
      ENDDO

* search for keyword 'REACTIONS'
1999  CONTINUE
      line=' '
      READ(lin,'(a)',END=9400) inline
      IF (loreply) WRITE(lout,'(a)') inline
      CALL cleanline(inline,outline(1),lenread)
      IF (lenstr(line,lenline).eq.0) GOTO 1999
      IF (index(line,"!").EQ.1) GOTO 1999
      IF (index(line,"SPECIES").NE.0) GOTO 1999
      IF (index(line,'REACTIONS').NE.1) GOTO 9410
      WRITE(6,*) 'reading the reactions ...'

* read reactions. Label 2000 is the reentry point to read next reaction
* ---------------------------------------------------------------------
2000  CONTINUE
c      WRITE(6,*) numre
      line=' '
      nlines=1
      READ(lin,'(a)',end=9420) inline
      CALL cleanline(inline,outline(nlines),lenread)
! DEBUG
      !PRINT*,outline(nlines)

* if keyword END found => goto 5000 (check first the last reaction)
      ipos=index(line,'END')
      IF (ipos.eq.1) THEN
         CALL chkaux(
     1      maxaux, maxfo, maxre,
     2      maxextra, maxcvar,
     3      max_m, maxhv, maxo2, mxrpero, maxro2,maxiso,
     3      maxt,mxrdimer,maxdimer,
     4      numre, numfo, numcvar, num_m,
     5      numextra, numhv, numo2, nummeo2,numiso,
     5      numain, numaou, numwin, numwou,
     6      lout,
     7      idfo, idhv,
     9      idextra, idcvar, id_m, ido2, idmeo2,idiso,
     9      idain,idaou,idwin,idwou,
     9      cvarcf, extracf,isocf,
     8      focf, hvcf, hvfact,
     8      woucf,wincf,
     &      nrpero,idreacro2,
     &      nrdimer,idreacdimer,
     7      lostop, locheck, lo_m, lofo,lo_iso,
     6      lohv,loaux,locvar,loextra,lo_o2,lo_meo2,lopero,ipero,
     6      lodimer,idimer,
     7      loain,loaou,lowin,lowou,
     l      itype, xauxcf)
         GOTO 5000
      ENDIF

* find out if input line is a reaction or auxiliary information
* -------------------------------------------------------------

2100  CONTINUE

* read the input line, starting from the end
      DO 2200 i=lenline,1,-1

* blanks
        IF (line(i:i).eq.' ') THEN
          GOTO 2200

* auxiliary information
        ELSE IF(line(i:i).eq.'/')THEN
          IF(loreply)WRITE(lout,'(a)')inline
           CALL getaux (
     1       maxaux,
     1       line, lenline, lout,
     1       lostop, loaux,
     1       nauxpar,
     1       xauxcf)

          GOTO 2000

* continue line for reaction written using more than one line
        ELSE IF (line(i:i).eq.'+'.or.
! JMLT 2021: dashes are ALLOWED for MECHGEN names
     &           (line(i:i).eq.'-'.AND.lco.NE.8).or.
     &           line(i:i+1).eq.'=>') THEN

          IF (loreply.and.nlines.eq.1)
     &        WRITE(lout,'(a16,i7)') 'REACTION NUMBER ',numre+1
          IF (loreply) WRITE(lout,'(a)')inline

          nlines=nlines+1
          IF (nlines.gt.linedim) THEN
            WRITE(lout,*)
            WRITE(lout,*)'   --error--   too many input lines',
     &                   ' for a single reaction equation',
     &                   '(linedim=',linedim,')'
            WRITE(lout,*)
            lostop=.true.
            GOTO 2000
          ENDIF

          READ (lin,'(a)',END=9420) inline
          CALL cleanline(inline,outline(nlines),lenread)
          GOTO 2100

* otherwise : check the auxiliary information for the reaction
* numre (the reaction that was read previously) and read the
* reaction numre+1.
        ELSE

          CALL chkaux(
     1      maxaux, maxfo, maxre,
     2      maxextra, maxcvar,
     3      max_m, maxhv, maxo2, mxrpero, maxro2,maxiso,
     3      maxt,mxrdimer,maxdimer,
     4      numre, numfo, numcvar, num_m,
     5      numextra, numhv, numo2, nummeo2,numiso,
     5      numain, numaou, numwin, numwou,
     6      lout,
     7      idfo, idhv,
     9      idextra, idcvar, id_m, ido2, idmeo2,idiso,
     9      idain,idaou,idwin,idwou,
     9      cvarcf, extracf,isocf,
     8      focf, hvcf, hvfact,
     8      woucf,wincf,
     &      nrpero,idreacro2,
     &      nrdimer,idreacdimer,
     7      lostop, locheck, lo_m, lofo,lo_iso,
     6      lohv,loaux,locvar,loextra,lo_o2,lo_meo2,lopero,ipero,
     6      lodimer,idimer,
     7      loain,loaou,lowin,lowou,
     l      itype, xauxcf)

          IF (loreply.and.nlines.eq.1)
     &        WRITE(lout,'(a16,i7)') 'REACTION NUMBER ',numre+1
          IF (loreply) WRITE(lout,'(a)') inline

          CALL getreac(
     1       line, lenline, lout,
     3       numre, numsp, numstoi,
     4       idrestoi, idpdstoi, restoicf, pdstoicf, chrsp,
     5       mx12stoi, mx1stoi, mx2stoi, mx3stoi,
     6       lostop, locheck, lo_m, lofo, lo_o2, lo_meo2,lo_iso,
     7       lohv, locvar, loextra, lopero, ipero,lodimer,idimer,
     7       loain,loaou,lowin,lowou,
     8       itype,
     9       arrhcf
     1       )
          GOTO 2000
        ENDIF

2200  CONTINUE

* read next reaction => goto 2000
      IF (loreply) WRITE(lout,'(a)') inline
      GOTO 2000

* ------------------------------------
* FINAL CHECK OF THE CHEMICAL SCHEMES
* ------------------------------------

5000  CONTINUE

      IF (lostop) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   -- error exit -- before reac'
        WRITE(lout,*)
        CLOSE(lin)
        CLOSE(lout)
        STOP 'ERROR 2: in the list of reactions => check file '
      ENDIF

      CALL chkreac (
     1    lout,
     2    maxre, mxleft, mxright, mxrpero, maxro2,numiso,
     3    num_n, num_m, numfo, numhv, numre, numo2, nummeO2, nrpero,
     3    mxrdimer,nrdimer,maxdimer,
     4    numcvar, numextra, itype, id_n,
     4    numain, numaou, numwin, numwou,
     5    numstoi, idrestoi, idpdstoi,
     6    lostop, locheck)

* check hv reaction
* ----------------------
* check the reaction with hv. To decrease CPU
* time, the work is made on the chormophore only
* (instead of the reaction). Here, the various
* chromophore are defined based on the reaction
* dataset.
      WRITE(6,*) '     checking hv reaction ....'

      ntchromo=0
      mxlab=0
      DO i=1,mchromo
        chromocf(i)=0
        numchromo(i)=0
      ENDDO

      DO 400 i=1,numhv
        idcf=nint(hvcf(i))

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

* write info about hv reaction and check the size of the table
      DO i=1,ntchromo
        IF (mxlab.lt.numchromo(i)) mxlab=numchromo(i)
      ENDDO
      WRITE (lout,*) ' '
      WRITE (lout,*) 'total number of HV reaction :',numhv
      WRITE (lout,*) 'number of chromophore (label):',ntchromo
      WRITE (lout,*) 'max number of species in a given label:',mxlab
      WRITE (lout,*) ' '
      WRITE (lout,*) 'occurence (column 3) of each label (column 2):'
      DO i=1,ntchromo
        WRITE (lout,'(i3,2x,i6,2x,i6)') i,chromocf(i), numchromo(i)
      ENDDO
      DO i=1,ntchromo
        IF (ntchromo.gt.mchromo) THEN
          WRITE(6,*) '--error-- in size of the table. Parameter
     &                ntspchromo is underestimated. see *.akoi file'
          WRITE(lout,*) '--error--, number of species that belong to'
          WRITE(lout,*) '          the chromophore (label):',idcf
          WRITE(lout,*) '          exceed mspchromo. '
          WRITE(lout,*) '          Change akparameter accordingly'
c          STOP
        ENDIF
      ENDDO

* stop if error found
      IF (lostop) THEN
        WRITE(lout,*)
        WRITE(lout,*)'   -- error exit -- after reading reaction data'
        WRITE(lout,*)
        CLOSE(lin)
        CLOSE(lout)
        STOP 'ERROR 3: in the list of reactions => check file '
      ENDIF

* calculate arrhcf(i,1)=log(arrhcf(i,1))

      DO i=1,numre
        IF(arrhcf(i,1).GT.1.0E-50)THEN
          arrhcf(i,1)=log(arrhcf(i,1))
        ELSE
          arrhcf(i,1)=-999.
        ENDIF
      ENDDO

* stoichiometric coefficients read in CVAR reaction does not have
* anay sence for the reaction products (purpose of CVAR is to set
* these coefficients has a function of temperature). Put hugge value
* here (remove comment if wanted).
c      DO icvar = 1, numcvar
c        ire = idcvar(icvar)
c        DO is = 1, numstoi(ire,2)
c          stoicf(ire,is,2) = 9.e9
c        ENDDO
c      ENDDO

* ----------------------
* CREATING LINK-FILE
* ----------------------

      WRITE(6,*) 'writing the output ...'

* -------------------------------
* WRITE-TO-LINK USED TO BE HERE!
* -------------------------------
! DEBUG !
      GOTO 81
! END DEBUG !
        CALL wrtlinkbin(llink,numsp,numre,num_n,num_m,
     &            numfo,numhv,numcvar,numextra,numo2,nummeo2,
     &            numain,numaou,numwin,numwou, nauxpar,
     &            numstoi,numiso,mx12stoi,mx1stoi,mx2stoi,
     &            chrsp,id_n,id_m,
     &            idfo,idhv,idcvar,idextra,ido2,idmeo2,
     &            idiso,idain,idaou,idwin,idwou,
     &            idrestoi,restoicf,idpdstoi,pdstoicf,
     &            arrhcf,focf,hvcf,hvfact,cvarcf,extracf,
     &            isocf,woucf,wincf,
     &            nrpero,idreacro2,nrdimer,idreacdimer,wmol)

81    CONTINUE

! DEBUG !
! END DEBUG !
      IF (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN
        CALL wrtlinkncdf(lout,numsp,numre,num_n,num_m,
     &             numfo,numhv,numcvar,numextra,numo2,nummeo2,
     &             numain,numaou,numwin,numwou, nauxpar,
     &             numstoi,numiso,mx12stoi,mx1stoi,mx2stoi,
     &             chrsp,id_n,id_m,
     &             idfo,idhv,idcvar,idextra,ido2,idmeo2,
     &             idiso,idain,idaou,idwin,idwou,
     &             idrestoi,restoicf,idpdstoi,pdstoicf,
     &             arrhcf,focf,hvcf,hvfact,cvarcf,extracf,
     &             isocf,woucf,wincf,
     &             nrpero,idreacro2,nrdimer,idreacdimer,wmol)
      ENDIF ! (iofmt_fg.EQ.0.OR.iofmt_fg.EQ.2) THEN

82    CONTINUE

      WRITE(lout,*)
      WRITE(lout,*)'   -- NORMAL EXIT --  LINK-FILE HAS BEEN CREATED'
      WRITE(lout,*)

      do i = 1, numre
        do j = 1, mxright
          if (idpdstoi(i,j) > numsp) then
            stop 'wrong index in idrestoi'
          endif
        enddo
      enddo

      WRITE(6,*) 'number of reaction : ', numre
      WRITE(6,*) 'number of species : ', numsp

* -------------------------------
* FINAL READING AND CLOSE FILES
* -------------------------------

* read the line that may exist at end of the file and copy them to
* the output file
8000  CONTINUE
      line=' '
      READ(lin,'(a)',END=8888) inline
      IF (loreply) WRITE(lout,'(a)') inline
c      CALL cleanline(inline,outline(1),lenread)
c      IF(lenstr(line,lenline).eq.0)GOTO 8000
      GOTO 8000
8888  CONTINUE

      CLOSE(lin)
      CLOSE(lout)

* write the CPU time used to achieve the simulation
      tdat(1)=0.
      tdat(2)=0.
!      duree=etime(tdat)
!      WRITE (6,'(a14, 1pe11.4)') 'duree=', duree
!      WRITE (6,'(a14, 1pe11.4)') 'elapsed in checked =', telstd


      STOP '-- no pb found in chemical scheme --'

* -----------
* STOP ERROR
* -----------

9100  CONTINUE
      WRITE(lout,*)'   --error--   end of input file before species',
     &                           ' declaration'
      GOTO 9999
9210  CONTINUE
      WRITE(lout,*)'   --error--   the word species not found'
      GOTO 9999
9220  CONTINUE
      WRITE(lout,*)'   --error--   end of input file while READing',
     &                           ' species'
      GOTO 9999
9400  CONTINUE
      WRITE(lout,*)'   --error--   end of input file before reactions',
     &                           ' declaration'
      GOTO 9999
9410  CONTINUE
      WRITE(lout,*)'   --error--   the word reactions not found'
      GOTO 9999
9420  CONTINUE
      WRITE(lout,*)'   --error--   end of input file while READing',
     &                           ' reactions'
9999  CONTINUE
c
      CLOSE(lin)
      CLOSE(lout)
      STOP '-- stop error--'

      END

