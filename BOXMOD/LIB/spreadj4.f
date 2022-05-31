*******************************************************************
* cette routine a pour objet de lire le fichier des               *
* constantes de photolyse tabulees pour dIFferentes               *
* valeurs d'angle zenithal.                                       *
* ATTENTION : la programmation suppose que les frequences         *
* de photolyse sont toutes DOnnees pour le meme jeu               *
* d'angles zenithaux => aucun test ne le controle (a              *
* ajouter dans les versions futures)                              *
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
      SUBROUTINE readj4(lout,numsp,chrsp,nt1chromo,chromo1cf,
     1                  chromomedcf,ntmedchromo,chromotopcf,
     2                  numtet,xang,
     3                  rat1pho,coef1pho,ratmedpho,coefmedpho,
     4                  rattoppho,coeftoppho)
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lout, numsp
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

* LOCAL
      CHARACTER(76) line
      CHARACTER(4)  keyword
      INTEGER n1test(mchromo),nmedtest(mmedchromo),ntoptest(mtopchromo)
      REAL    val1(maxang),val2(maxang)
      REAL    y(maxang),cpho(nlo)
      LOGICAL lokerr, lofind
      INTEGER nchar
      INTEGER  ii, njdat,isp, nreap, i_val,ierr, i,j,k,l,ij
      REAL     bidon1, bidon2

      nchar=76

      lokerr=.false.
      DO ii=1,mchromo
        n1test(ii)=0
      ENDDO
      DO ii=1,mmedchromo
        nmedtest(ii)=0
      ENDDO
      DO ii=1,mtopchromo
        ntoptest(ii)=0
      ENDDO

* open file of input data
      njdat=12 
      OPEN(njdat,file='jfile.phot',form='formatted',status='old')

* return label to read new set of data
90    CONTINUE
      READ(njdat,'(a4,(a))',end=99) keyword,line

* check if the line is a comment
      IF(keyword(1:1).EQ.'.'.OR.keyword(1:1).EQ.'/'.OR.
     &   keyword(1:1).EQ.'!')GOTO 90

* check end of file
      IF(keyword(1:3).EQ.'END') THEN
        GOTO 100

* read the photolytic data
      ELSE IF (keyword.EQ.'PHOT') THEN

* check is the species is known in the mechanism
* if species is unknown, data are read and go to the next set
c        CALL akspnum(line,chrsp,numsp,isp)
c        IF (isp.eq.0) THEN
c          lokerr=.true.
c          WRITE(lout,*)' --error--  while reading in jdata, line',line
c          WRITE(lout,*)'            species is unkown'
c           DO k=1,numtet
c             read(njdat,*) bidon1,bidon2
c          ENDDO
c          GOTO 90
c        ENDIF

* read the flag number for the reaction
        nreap=i_val(line,1,nchar,2,ierr)
        IF(ierr.ne.0)then
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading in jdata, line',line
          WRITE(lout,*)'            value nreap not read correctly'
        ENDIF

* read the number of data in the set
        numtet=i_val(line,1,nchar,3,ierr)
        IF(ierr.ne.0)then
          lokerr=.true.
          WRITE(lout,*)' --error--  while reading in jdata, line',line
          WRITE(lout,*)'            value numtet not read correctly'     
        ENDIF

* read the data for the set
        DO i=1,numtet
           READ(njdat,*) val1(i),val2(i)
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

* check that all photolytic reactions are feed with data
100   CONTINUE
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
      IF (lokerr) STOP 'error in readj, check outdat.out'


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

* fin de la routine
      CLOSE(njdat)

      END
