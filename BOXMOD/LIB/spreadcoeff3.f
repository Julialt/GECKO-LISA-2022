***********************************************************************
* Cette routine a pour objet de lire le fichier donnant la dependance *
* des coeff. stoechiometriques avec la temperature la temperature     *
* INPUT :                                                             *
*   lout   : numero de fichier pour ecriture des erreurs              *
*                                                                     *
* OUTPUT : (in the common block)                                      *
*   valcoe(i,j,k) : valeur des coef. stoech. a la ieme temperature,   *
*                   pour la kieme reaction cvar et pour la jieme      *
*                   espece                                            *
*   ntype(k)      : type de la kieme reaction (ne sert plus dans      *
*                   cette version - mais laisse au cas ou)            *
*   numcoe(k)     : nbre de coef. stoe. tabule pour la kieme reaction *
*   nopc(k)       : nbre d'espece compteur de la kieme reaction qui   *
*                   doivent etre reevalue a partir des coef. stoe.    *
*                   associes aux operateurs                           *
*   ndatopc(k,j)  : nbre de donnees pour reevaluer le jieme compteur  *
*                   de la kieme reaction                              *
*   nposopc(k,j,l): position (l) des operateurs pour la reevaluation  *
*                   du jieme compteur de la kieme reaction            *
***********************************************************************
      SUBROUTINE readcoeff3(lout, numsp, chrsp, numcvar, cvarcf,
     2                      valcoe,ntype,numcoe,nopc,ndatopc,nposopc)
      USE akparameter_module
      IMPLICIT NONE
c      INCLUDE 'akcommon.h'
c      COMMON /coefvar3/ valcoe(nset,maxcoe,maxcvar),ntype(maxcvar),
c     &                  numcoe(maxcvar),nopc(maxcvar),
c     &                  ndatopc(maxcvar,mopc),nposopc(maxcvar,mopc,mpos)

* INPUT
      INTEGER lout, numcvar
      REAL    cvarcf(maxcvar)

      CHARACTER(maxlsp) chrsp(maxsp)
      INTEGER  numsp

* OUTPUT
      INTEGER ntype(maxcvar), numcoe(maxcvar), nopc(maxcvar)
      INTEGER ndatopc(maxcvar,mopc), nposopc(maxcvar,mopc,mpos)
      REAL    valcoe(nset,maxcoe,maxcvar)

* LOCAL
      CHARACTER(80) line
      INTEGER ntest(maxcvar)
      INTEGER lenli, nfic, memo, isp, idreac, i_val, ierr, ltest
      INTEGER i, ii, j, k, lenstr, ipro, itemp

      lenli=80

* open the file
      nfic=12 
      OPEN(nfic,file='cfile.coe',form='formatted',status='old')

* initialise  
      memo=0
      DO ii=1,numcvar
        ntest(ii)=0
      ENDDO

* read the file
90    CONTINUE 
      READ(nfic,'(a)',end=99) line
      IF (line(1:1).EQ.'/'.or.line(1:1).EQ.'!') GOTO 90       
      IF (line(1:3).EQ.'END') GOTO 100

* check if the species are correctly sorted
      IF (line(1:4).EQ.'SECP') THEN
         CALL akspnum(line(5:lenli),chrsp,numsp,isp)
         IF (isp.EQ.0) THEN
           WRITE(lout,*)'--error-- in file coefficient'
           WRITE(lout,*)'          while reading ',
     &                    line(5:lenstr(line,lenli))
           WRITE(lout,*)'          species is unknown'
           STOP 
         ENDIF
         IF (isp.lt.memo) THEN
           WRITE(lout,*)'--error-- in file coefficient'
           WRITE(lout,*)'because species order not correctly given'
           STOP
         ENDIF
         memo=isp
         GOTO 90
      ENDIF

* fill the common block
      IF (line(1:4).EQ.'IDEN') THEN
         idreac=i_val(line,5,lenli,1,ierr)
         IF (ierr.ne.0) THEN
           WRITE(lout,*)'--error-- in readcoeff '
           WRITE(lout,*)'          while reading coeff 1 of line:'
           WRITE(lout,*) line(1:lenstr(line,lenli))
           STOP
         ENDIF
         ltest=0

* check IF the cvar coefficient is allowed, IF yes => set the value to
* the corresponding tables (STOP IF error is uncountered) :
         DO 20 ii=1,numcvar
           IF (idreac.EQ.nint(cvarcf(ii))) THEN

* ntype : 0=sans type particulier, 1=alcane
             ntype(ii)=i_val(line,5,lenli,2,ierr)
             IF (ierr.ne.0) THEN
                WRITE(lout,*)'--error-- in readcoeff '
                WRITE(lout,*)'          while reading coeff 2 of line:'
                WRITE(lout,*) line(1:lenstr(line,lenli))
                STOP
             ENDIF

* numcoe : nombre de coefficient stoech. cote droit de la react.
*          pour lesquels une parametrisation est DOnne
             numcoe(ii)=i_val(line,5,lenli,3,ierr)
             IF (ierr.ne.0) THEN
                WRITE(lout,*)'--error-- in readcoeff '
                WRITE(lout,*)'          while reading coeff 3 of line:'
                WRITE(lout,*) line(1:lenstr(line,lenli))
                STOP
             ENDIF

* nopc : nombre d'especes compteurs dont les coefficients stoech.
* doivent etre evaluees a partir des coef. stoe. des operateurs 
             nopc(ii)=i_val(line,5,lenli,4,ierr)
             IF (ierr.ne.0) THEN
                WRITE(lout,*)'--error-- in readcoeff '
                WRITE(lout,*)'          while reading coeff 4 of line:'
                WRITE(lout,*) line(1:lenstr(line,lenli))
                STOP
             ENDIF

* ndatopc : nombre de donnees necessaires a la reevaluation des especes
*           compteurs. 
            IF (nopc(ii).ne.0) THEN
              DO j=1,nopc(ii)
                ndatopc(ii,j)=i_val(line,5,lenli,4+j,ierr)
                IF (ierr.ne.0) THEN
                WRITE(lout,*)'--error-- in readcoeff '
                WRITE(lout,*)'          while reading coeff :',4+j
                WRITE(lout,*) line(1:lenstr(line,lenli))
                STOP
                ENDIF
              ENDDO
            ENDIF

* remplissage des position des differentes especes dans la reaction
* necessaire a la reevaluation des compteurs
            DO j=1,nopc(ii)
              READ(nfic,*,err=97) 
     &        (nposopc(ii,j,k),k=1,ndatopc(ii,j)+1)
            ENDDO
            GOTO  98

* stop if error in reading nposopc
97          WRITE(lout,*) '--error-- in readcoeff'
            WRITE(lout,*) '          while reading nposopc coef.'
            WRITE(lout,*) '          for reaction ID number', idreac
            STOP

98          CONTINUE
* remplissage de valcoe  
             DO 40 itemp=1,nset
               READ(nfic,*) (valcoe(itemp,ipro,ii),ipro=1,numcoe(ii))
40           CONTINUE
             IF (ntest(ii).EQ.1) THEN 
                WRITE(lout,'(a21)')'--error-- in readcoeff'
                WRITE(lout,'(i7,a27)')
     &            idreac,' was found more than 1 time'
                STOP
             ELSE
                ntest(ii)=1
                ltest=1
             ENDIF
             GOTO 90
           ENDIF
20       CONTINUE
         IF (ltest.EQ.0) THEN
            DO 17 i=1,5
               READ (nfic,'(a)') line
17           CONTINUE
             GOTO 90
         ENDIF   
      ELSE
         WRITE(lout,*)'--error-- while reading in readcoeff',
     &                  line(1:lenstr(line,lenli))
         STOP
      ENDIF

* keyword end not found 
99    WRITE(lout,*)'--error-- keyword END not found'
      STOP

100   CONTINUE
* check that all cvar reaction were found
      DO 88 ii=1,numcvar
        IF (ntest(ii).EQ.0) THEN
          WRITE(lout,'(a16,i9)')'--error--,ident',nint(cvarcf(ii))
          WRITE(lout,'(a23)')       ' not found in file.coe' 
          STOP
        ENDIF
88    CONTINUE

      CLOSE (nfic)

      END
