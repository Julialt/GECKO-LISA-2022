*******************************************************************
! NetCDF-specific - returns extra values required for NetCDF o/p
!
* cette routine a pour objet de lire le fichier des               *
* constantes de photolyse tabulees pour dIFferentes               *
* valeurs d'angle zenithal.                                       *
* ATTENTION : la programmation suppose que les frequences         *
* de photolyse sont toutes DOnnees pour le meme jeu               *
* d'angles zenithaux => aucun test ne le controle (a              *
* ajouter dans les versions futures)                              *
*                                                                 *
* NETCDF VERSION: reads ascii *.phot input & assigns chromophores *
*                 to tables AND ALSO saves all photfile input for *
*                 reproduction in NetCDF omnibus output.          *
* NB: ASSUMES THAT EACH CHROMOPHORE TABLE HAS THE SAME SZA VALUES *
*                                                                 *
* EDIT AUTHOR: Julia Lee-Taylor, NCAR, Feb 2018                   *
*                                                                 *
* INPUT :                                                         *
*   lout   : numero de fichier pour ecriture des erreurs          *
* OUTPUT                                                          *
*   numtet      : nombre d'angle zenithaux tabule                 *
*   xang(j)     : valeur des angles zenithaux tabules             *
*   ratpho(i,j) : frequence de photolyse de la reaction i         *
*                 a l'angle xang(j)                               *
*   npos_cf     : numero de la reaction ayant le drapeau /HV/ 1   *
*                 (correspond interativement a O3 -> 2 OH         *
*   coefpho(i,*): valeur des coef. de l'interpolation des         *
*                 frequences de photolyse pour la reaction i      *
*                 (l'interpolation est un polynome de degre 2     *
*                 sur chaque intervalle de frequence de photolyse *
*                 soit 3*numtet coefficient                       *
*******************************************************************
      SUBROUTINE readj4_ncdf(lout,!ncid,
     1                  numsp,chrsp,nt1chromo,chromo1cf,
     1                  chromomedcf,ntmedchromo,chromotopcf,
     2                  numtet,xang,
     3                  rat1pho,coef1pho,ratmedpho,coefmedpho,
     4                  rattoppho,coeftoppho,
     &                  mxjtab,mxsza,llin,
     5                  njtab,nsza,idjtab,jsza,jvref,jnam,jreac)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lout, numsp
      INTEGER  ncid
!      INTEGER   ntchromo
!      INTEGER  chromocf(mchromo)
      CHARACTER(maxlsp) chrsp(maxsp)

      INTEGER  nt1chromo,chromo1cf(mchromo)
      INTEGER  chromomedcf(mmedchromo),ntmedchromo
      INTEGER  chromotopcf(mtopchromo)


* OUTPUT
      INTEGER  numtet
      REAL     xang(maxang)

      REAL     rat1pho(mchromo,maxang),coef1pho(mchromo,nlo)
      REAL     ratmedpho(mmedchromo,maxang),coefmedpho(mmedchromo,nlo)
      REAL     rattoppho(mtopchromo,maxang),coeftoppho(mtopchromo,nlo)

* photolysis input data, for re-output 
      INTEGER  mxjtab,mxsza,llin
!      INTEGER,PARAMETER:: mxjtab=150,mxsza=15,llin=76
      INTEGER  njtab,nsza,idjtab(mxjtab)
      REAL     jsza(mxsza),jvref(mxsza,mxjtab)
      CHARACTER(maxcoe) jnam(mxjtab)
      CHARACTER(llin) jreac(mxjtab)

* LOCAL
      CHARACTER(llin) line
      CHARACTER(4)  keyword
      INTEGER n1test(mchromo),nmedtest(mmedchromo),ntoptest(mtopchromo)
      REAL    val1(maxang),val2(maxang)
      REAL    y(maxang),cpho(nlo)
      LOGICAL lokerr, lofind
      INTEGER nchar
      INTEGER  ii, iunit,isp, nreap, i_val,ierr, i,j,k,l,ij
      REAL     bidon1, bidon2

* for NetCDF
      INTEGER nhdr,ipos
      CHARACTER(8)  attname
*========================================================

      nchar=76
      lokerr=.false.
      n1test(1:mchromo)=0 
      nmedtest(1:mmedchromo)=0
      ntoptest(1:mtopchromo)=0
      
      njtab = 0
      nsza = 0
      nhdr = 0

*-------------------------
* open file of input data
      iunit=12 
      OPEN(iunit,file='jfile.phot',form='formatted',status='old')

*--------------------------------------
* return label to read new set of data
90    CONTINUE
      READ(iunit,'(a4,(a))',end=99) keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'.'.OR.keyword(1:1).EQ.'/'.OR.
     &   keyword(1:1).EQ.'!')GOTO 90

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 100

* read the photolytic data
      ELSE IF (keyword.EQ.'PHOT') THEN
        
        njtab = njtab+1

* check if the species is known in the mechanism
* if species is unknown, data are read and go to the next set
c        CALL akspnum(line,chrsp,numsp,isp)
c        IF (isp.eq.0) THEN
c          lokerr=.true.
c          WRITE(lout,*)' --error--  while reading in jdata, line',line
c          WRITE(lout,*)'            species is unkown'
c           DO k=1,numtet
c             read(iunit,*) bidon1,bidon2
c          ENDDO
c          GOTO 90
c        ENDIF

* read the flag number for the reaction
* assign to idjtab output table

        nreap=i_val(line,1,nchar,2,ierr)
        idjtab(njtab) = nreap 
        IF(ierr.ne.0)then
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading in jdata, line',line
          WRITE(lout,*)'            value nreap not read correctly'
        ENDIF

* read the number of data in the set, 
* assign nsza as size of largest table

        numtet=i_val(line,1,nchar,3,ierr)
        IF(numtet.GT.nsza)nsza=numtet
        IF(ierr.ne.0)then
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading in jdata, line',line
          WRITE(lout,*)'            value numtet not read correctly'     
        ENDIF

* read the data for the set
* assign jsza and jvref output tables

        DO i=1,numtet
           READ(iunit,*) val1(i),val2(i)
           jsza(i) = val1(i)
           jvref(i,njtab) = val2(i)
!          val2(i)=val2(i)*1.54 !!
        ENDDO

* put the dataset to the corresponding table in the mechanism
* a test is made to check that the data are only given one time
* for each reaction

* check if label belong to the "most used" table
        lofind=.false.
        DO ij=1,mtopchromo
          IF (nreap.EQ.chromotopcf(ij)) THEN
            IF (ntoptest(ij).EQ.0) THEN
              DO k=1,numtet
                 xang(k)=val1(k)
                 rattoppho(ij,k)=val2(k)
              ENDDO
              ntoptest(ij)=1
              lofind=.true.
            ELSE
              WRITE(lout,'(a19)') '-- error-- in readj'
              WRITE(lout,'(i5,a25)') nreap,' ident. present + 1 time'
              lokerr=.true.
              lofind=.true.
            ENDIF
          ENDIF
        ENDDO
        IF (lofind) GOTO 90

* check if label belong to the "used only once" tables
        DO ij=1,nt1chromo
          IF (nreap.EQ.chromo1cf(ij)) THEN
            IF (n1test(ij).EQ.0) THEN
              DO k=1,numtet
                 xang(k)=val1(k)
                 rat1pho(ij,k)=val2(k)
              ENDDO
              n1test(ij)=1
              lofind=.true.
            ELSE
              WRITE(lout,'(a19)') '-- error-- in readj'
              WRITE(lout,'(i5,a25)') nreap,' ident. present + 1 time'
              lokerr=.true.
              lofind=.true.
            ENDIF
          ENDIF
        ENDDO
        IF (lofind) GOTO 90

* check if label belong to the "regularly used" table
        DO ij=1,ntmedchromo
          IF (nreap.EQ.chromomedcf(ij)) THEN
            IF (nmedtest(ij).EQ.0) THEN
              DO k=1,numtet
                 xang(k)=val1(k)
                 ratmedpho(ij,k)=val2(k)
              ENDDO
              nmedtest(ij)=1
              lofind=.true.
            ELSE
              WRITE(lout,'(a19)') '-- error-- in readj'
              WRITE(lout,'(i5,a25)') nreap,' ident. present + 1 time'
              lokerr=.true.
              lofind=.true.
            ENDIF
          ENDIF
        ENDDO
        GOTO 90

c        DO 50 ij=1,ntchromo
c          IF (nreap.EQ.chromocf(ij)) THEN
c            IF (ntest(ij).EQ.0) THEN
c              DO k=1,numtet
c                  xang(k)=val1(k)
c                  ratpho(ij,k)=val2(k)
c               ENDDO
c               ntest(ij)=1
c             ELSE
c               WRITE(lout,'(a19)') '-- error-- in readj'
c               WRITE(lout,'(i5,a25)') nreap,' ident. present + 1 time'
c               lokerr=.true.
c             ENDIF 
c          ENDIF
c50      CONTINUE
c        GOTO 90

* keyword unknown : stop
      ELSE
        WRITE(lout,*)' --error--  while reading in readj',keyword
        WRITE(lout,*)'            keyword unkown'
        STOP
      ENDIF

* label to check that the keyword 'END' was found
99    WRITE(lout,*)' --error--  keyword END not found in readj'
      lokerr=.true.

100   CONTINUE
      REWIND(iunit)

*=====================================================
* NetCDF file space definitions 
!! NOW MOVED TO setup_ncdf_op.f90 !!
*=====================================================
      isp = 0

* define dimensions njtab & nsza in NetCDF output file
!      CALL eznc_def_dim(ncid,"njtab",njtab)
!      CALL eznc_def_dim(ncid,"nsza",nsza)
!      CALL eznc_def_dim(ncid,"llin",llin)
! maxcoe is defined in akparameter

* define space for tables in NetCDF output file
* define chromophore name
* (conveniently, "maxcoe" is the correct size for the names)
!      CALL eznc_def_1Dchar(ncid,"jnam","maxcoe","njtab")
!      CALL eznc_def_localatt(ncid,"jnam","title",
!     &                      "j-value species name")
* define chromophore reaction text
!      CALL eznc_def_1Dchar(ncid,"jreac","llin","njtab")
!      CALL eznc_def_localatt(ncid,"jreac","title",
!     &                      "j-value reaction string")
* define chromophore reaction id number
!      CALL eznc_def_1Dint(ncid,"idjtab","njtab")
!      CALL eznc_def_localatt(ncid,"idjtab","title",
!     &               "chromophore id number")
!      CALL eznc_def_localatt(ncid,"idjtab","cross-reference",
!     &               "hvcf in mechanism")
* define values for j-tables
!      CALL eznc_def_1Dreal(ncid,"jsza","nsza")
!      CALL eznc_def_localatt(ncid,"jsza","title",
!     &                      "solar angles for reference j-values")
!      CALL eznc_def_localatt(ncid,"jsza","units",
!     &                      "degrees from zenith")

!      CALL eznc_def_2Dreal(ncid,"jvref","nsza","njtab")
!      CALL eznc_def_localatt(ncid,"jvref","title",
!     &                            "reference j-values")
!      CALL eznc_def_localatt(ncid,"jvref","units","s-1")

* reading input file lines (again)
* return label to read new set of data
190    CONTINUE
      READ(iunit,'(a1,(a))',end=199) keyword,line

* if the line is a header comment, write as attribute
      IF(keyword(1:1).EQ.'/'.AND.line.NE.''.AND.line(1:1).NE.' ')THEN
        nhdr = nhdr + 1
!        IF(nhdr.LT.10)THEN
!          WRITE (attname,"(a7,i1)") "header_",nhdr
!        ELSE
!          WRITE (attname,"(a6,i2)") "header",nhdr
!        ENDIF
!        CALL eznc_def_localatt(ncid,"jvref",attname,line)
        GOTO 190

* if line is a jtable comment (preceded by "PHOT"), save as value "jnam"
* and save preceding line as "jreac"
      ELSE IF(keyword(1:1).EQ.'P'.AND.line(1:3).EQ.'HOT')THEN
        isp = isp+1
        BACKSPACE(iunit)
        BACKSPACE(iunit)
        READ(iunit,'(a2,(a))',end=199) keyword,line
        jreac(isp) = line

        READ(iunit,'(a6,(a))',end=199) keyword,line
        ipos = INDEX(line," ")
        jnam(isp) = line(1:ipos-1)
        GOTO 190

* check for end of file
      ELSE IF(keyword(1:1).EQ.'E'.AND.line(1:2).EQ.'ND') THEN
        GOTO 200

* if not useful text, go back to read next data line
      ELSE
        GOTO 190

      ENDIF !(keyword(1:1).EQ.'/')THEN

* label to check that the keyword 'END' was found
199    WRITE(lout,*)' --error--  keyword END not found in readj'
      lokerr=.true.
200   CONTINUE

*=====================================================
* Chromophore table definitions:
*=====================================================

* check that all photolytic reactions are fed with data
      DO ii=1,mtopchromo
        IF (ntoptest(ii).EQ.0) THEN
          WRITE(lout,'(a16,i6)')'--error--a nreap (top)',chromotopcf(ii)
          WRITE(lout,'(a23)') ' not found in jdata.dat'
          lokerr=.true.
        ENDIF
      ENDDO

      DO ii=1,ntmedchromo
        IF (nmedtest(ii).EQ.0) THEN
          WRITE(lout,'(a16,i6)')'--error--b nreap (med)',chromomedcf(ii)
          WRITE(lout,'(a23)') ' not found in jdata.dat'
          lokerr=.true.
        ENDIF
      ENDDO

      DO ii=1,nt1chromo
        IF (n1test(ii).EQ.0) THEN
          WRITE(lout,'(a16,i6)')'--error--c nreap (1)',chromo1cf(ii)
          WRITE(lout,'(a23)') ' not found in jdata.dat'
          lokerr=.true.
        ENDIF
      ENDDO

* stop if error
      IF (lokerr) STOP 'error in readj_ncdf, check outdat.ou'


* look for the positon of the special reaction : O3 -> 2 OH
* this reaction must be set with the flag reaction 1
c      DO i=1,ntchromo
c        IF (chromocf(i).eq.1) npos_cf1=i
c      ENDDO

* preparation of the interpolating polynomial coefficients      
* for the photolytic rate constants depending on solar position 

* top tables
      DO i=1,mtopchromo
        DO l=1,nlo
          cpho(l)=0.
        ENDDO
        DO j=1,numtet
          y(j)=rattoppho(i,j)
        ENDDO
        CALL splna(numtet,xang,y,cpho,1)                                 
        DO j=1,nlo
          coeftoppho(i,j)=cpho(j)
        ENDDO      
      ENDDO    

* med tables
      DO i=1,ntmedchromo
        DO l=1,nlo
          cpho(l)=0.
        ENDDO
        DO j=1,numtet
          y(j)=ratmedpho(i,j)
        ENDDO
        CALL splna(numtet,xang,y,cpho,1)                                 
        DO j=1,nlo
          coefmedpho(i,j)=cpho(j)
        ENDDO      
      ENDDO    

* 1 table
      DO i=1,nt1chromo
        DO l=1,nlo
          cpho(l)=0.
        ENDDO
        DO j=1,numtet
          y(j)=rat1pho(i,j)
        ENDDO
        CALL splna(numtet,xang,y,cpho,1)                                 
        DO j=1,nlo
          coef1pho(i,j)=cpho(j)
        ENDDO      
      ENDDO    

c      DO 59 i=1,ntchromo
c        DO l=1,nlo
c          cpho(l)=0.
c        ENDDO
c        DO j=1,numtet
c          y(j)=ratpho(i,j)
c        ENDDO
c        call splna (numtet,xang,y,cpho,1)                                 
c        DO j=1,nlo
c          coefpho(i,j)=cpho(j)
c        ENDDO      
c59    CONTINUE      

*=====================================================
* fin de la routine
      CLOSE(iunit)


      END SUBROUTINE readj4_ncdf
