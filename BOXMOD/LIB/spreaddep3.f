*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* qui se depose, et de lire les parametres qui y sont associes    *
* INPUT :                                                         *
*   lout   : numero de fichier pour ecriture des erreurs          *
*   iscape : drapeau pour calcul des equilibres thermos           *
* OUTPUT :                                                        *
*   ndepspe  : nombre d'espece qui se depose                      *
*   depnamspe(i) : nom de la ieme espece qui se depose            *
*   depdatspe(i,j) : parametre j de la ieme espece                *
*                    j=1 => D_H2O/D_x                             *
*                    j=2 => constante de Henry                    *
*                    j=3 => facteur de reactivite                 *
*   iddepspe(i)    : numero d'identification de l'espece dans le  *
*                    mecanisme chimique                           *
* JMLT, Feb 2018: deposition arrays now have dimension "mxdep"    *
*******************************************************************
      SUBROUTINE readdep3(lout,chrsp,numsp,iscape,
     1                    ndepspe,depnamspe,depdatspe,iddepspe)
      !$ use OMP_LIB
      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER  lout,iscape, numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) depnamspe(mxdep)
      INTEGER  ndepspe, iddepspe(mxdep)
      REAL     depdatspe(mxdep,3)

* LOCAL
      CHARACTER(80) line
      INTEGER i,j,k,isp

* initialise
      ndepspe=0


      OPEN (12,file='vfile.dep',status='OLD')


* read the head of the file
      DO j=1,10
         READ(12,'(a)') line
      ENDDO

* read the file, end of read occurs when '****' is found
      DO j=1,mxdep
         READ(12,'(a)') line
         IF (line(1:3).eq.'END') GOTO 100
         ndepspe=ndepspe+1
         READ(line,'(1pe7.1,2x,1pe7.1,2x,1pe7.1,2x,a16)',err=10)
     &   (depdatspe(j,k),k=1,3),depnamspe(j)
      ENDDO
      CLOSE (12)

10    WRITE(lout,*) '-error--, while reading depspe.dat'
      WRITE(lout,*) '        , at line : ', j
      STOP

100   CONTINUE

* check if the species are defined in the mechanism
* stop if the species is unknown
!$OMP PARALLEL DO private(i, isp)
      DO i=1,ndepspe
        CALL akspnum(depnamspe(i),chrsp,numsp,isp)
        iddepspe(i)=isp
        IF (isp.eq.0) THEN
          WRITE(lout,*)' --error--  while reading depspe.dat'
          WRITE(lout,*)' species unkown=',depnamspe(i)
          STOP
        ENDIF
      ENDDO
!$OMP END PARALLEL DO

*******************************************************************
* recommence la procedure precedente pour les especes gazeuses    *
* en equilibre thermo                                             *
*******************************************************************
      IF (iscape.eq.1) THEN
        OPEN (12,file='vfileaero.dep',status='OLD')

* read the head of the file
        DO j=1,10
           READ(12,'(a)') line
        ENDDO

* read the file, end of read occurs when '****' is found
        DO j=1,mxdep
           READ(12,'(a)') line
           IF (line(1:4).eq.'****') GOTO 200
           ndepspe=ndepspe+1
           READ(line,'(1pe7.1,2x,1pe7.1,2x,1pe7.1,2x,a16)',err=20)
     &     (depdatspe(j,k),k=1,3),depnamspe(j)
        ENDDO
        CLOSE (12)

20      WRITE(lout,*) '-error--, while reading depspe.dat'
        WRITE(lout,*) '        , at line : ', j
        WRITE(lout,*) line
        STOP

200     CONTINUE

* check if the species are defined in the mechanism
* stop if the species is unknown
        DO i=1,ndepspe
          CALL akspnum(depnamspe(j),chrsp,numsp,isp)
          iddepspe(j)=isp
          IF (isp.eq.0)then
            WRITE(lout,*)' --error--  while reading depspe-aero.dat'
            WRITE(lout,*)' species unkown=',depnamspe(j)
            WRITE(lout,*) line
            STOP
          ENDIF
        ENDDO

      ENDIF

* end of the routine
      END
