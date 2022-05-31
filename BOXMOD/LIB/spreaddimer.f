*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* constituant le (les) compteurs RO2 - schema chimique generateur *
* routine a reecrire par la suite                                 *
*******************************************************************
      SUBROUTINE readdimer(chrsp,numsp,
     1                   ncldimer, numchemdimer,idchemdimer,cdimer)
      USE akparameter_module
      IMPLICIT NONE
      
* INPUT 
      INTEGER  :: numsp
      CHARACTER(maxlsp) :: chrsp(maxsp)

* OUTPUT
      INTEGER  :: ncldimer, numchemdimer(maxdimer)
      INTEGER  :: idchemdimer(mxrdimer,maxdimer)
      REAL     :: cdimer(maxdimer)

* LOCAL
c      CHARACTER(maxlsp) ro2sp(mxrpero,maxro2)
      CHARACTER(maxlsp)   :: dimersp(mxrdimer)
      CHARACTER(maxlsp+1) :: liner
      CHARACTER(20)         :: filnam
      INTEGER              :: i,j,isp,ll,nn, ipos

* set the number of ro2 class (must be .le. maxro2)
      ncldimer=4
      IF (ncldimer.gt.maxdimer) THEN
        WRITE(6,*) '--error--, in spreaddimer.f (read comment)'
        STOP
      ENDIF

* initialize
      DO i=1,mxrdimer
        DO j=1,maxdimer
          idchemdimer(i,j)=0
        ENDDO
      ENDDO      

      DO i=1,maxdimer
        numchemdimer(i)=0
        cdimer(i)=0.
      ENDDO

* -----------------------
* LOOP OVER THE RO2 CLASS
* -----------------------

      DO nn=1,ncldimer

* initialize
        DO i=1,mxrdimer
          dimersp(i)=' '
        ENDDO

* open file
        filnam(1:5)='indat'
        WRITE(filnam(6:6),'(i1)') nn
        filnam(7:)='.dim '
        WRITE(6,*) '    open ',filnam

        OPEN (12, file=filnam,STATUS='OLD')

* read DIM in the class
        i=0
56      READ(12,'(a)',END=59) liner
        IF (liner(1:3).EQ.'***') GOTO 59
        i=i+1
        IF (i.gt.mxrdimer) THEN
          WRITE(6,*) '--error-- number of dimer exceed mxrdimer in file'
          WRITE(6,*)   filnam 
          STOP
        ENDIF 
        ll=index(liner,' ')
        IF (ll.LE.1) THEN
          WRITE(6,*) '--error--, lecture DIMER in file', filnam
          STOP
        ENDIF
        dimersp(i)=liner
        GOTO 56
59      CONTINUE 
        CLOSE(12)
        numchemdimer(nn)=i
      
* get the id for the DIM species and store into idchemdimer table
* the program was changed to increase the search
        ipos=1
        DO 60 i=1,numchemdimer(nn)
          DO j=ipos,maxsp
            IF (dimersp(i).eq.chrsp(j)) THEN
              ipos=j
              idchemdimer(i,nn)=j
              GOTO 60
            ENDIF
          ENDDO

          WRITE(6,*) '--error--,DIMER species unidentified in readdimer'
          WRITE(6,'(a)') dimersp(i)
          WRITE(6,*) 'in file', filnam
          WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE' 
          WRITE(6,*) 'SPECIES MAY NOT BE SORTED' 
          STOP
60      CONTINUE

      ENDDO

* -----------------------
* RETURN
* -----------------------
      
      END 
