*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* constituant le (les) compteurs RO2 - schema chimique generateur *
* routine a reecrire par la suite                                 *
*******************************************************************
      SUBROUTINE readro2sof2_bin(chrsp,numsp,
     1                   nclro2, numchemro2,idchemro2,cro2)
      !$ use OMP_LIB
      USE akparameter_module
      IMPLICIT NONE
      
* INPUT 
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      INTEGER  nclro2, numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
      REAL     cro2(maxro2)

* LOCAL
c      CHARACTER(maxlsp) ro2sp(mxrpero,maxro2)
      CHARACTER(maxlsp) ro2sp(mxrpero)
      CHARACTER(maxlsp+1) liner
      CHARACTER(20)       filnam
      INTEGER i,j,isp,ll,nn, ipos
      LOGICAL    found

* set the number of ro2 class (must be .le. maxro2)
      nclro2=9
      IF (nclro2.gt.maxro2) THEN
        WRITE(6,*) '--error--, in spreadro2.f (read comment)'
        STOP
      ENDIF

* initialize
!$OMP PARALLEL DO private(i,j)
      DO i=1,mxro2cl
        DO j=1,maxro2
          idchemro2(i,j)=0
        ENDDO
      ENDDO   
!$OMP END PARALLEL DO

      DO i=1,maxro2
        numchemro2(i)=0
        cro2(i)=0.
      ENDDO

* -----------------------
* LOOP OVER THE RO2 CLASS
* -----------------------

      DO nn=1,nclro2

* initialize
!$OMP PARALLEL DO private(i)
        DO i=1,mxrpero
          ro2sp(i)=' '
        ENDDO
!$OMP END PARALLEL DO

* open file
        filnam(1:5)='indat'
        WRITE(filnam(6:6),'(i1)') nn
        filnam(7:)='.ro2 '
        WRITE(6,*) '    open ',filnam

        OPEN (12, file=filnam,STATUS='OLD')

* read RO2 in the class
        i=0
56      READ(12,'(a)',END=59) liner
        IF (liner(1:3).EQ.'***') GOTO 59
        i=i+1
        IF (i.gt.mxrpero) THEN
          WRITE(6,*) '--error-- number of RO2 exceed mxrpero in file'
          WRITE(6,*)   filnam 
          STOP
        ENDIF 
        ll=index(liner,' ')
        IF (ll.LE.1) THEN
          WRITE(6,*) '--error--, lecture RO2 in file', filnam
          STOP
        ENDIF
        ro2sp(i)=liner
        GOTO 56
59      CONTINUE 
        CLOSE(12)
        numchemro2(nn)=i
      
* get the id for the RO2 species and store into idchemro2 table
c        DO i=1,numchemro2(nn)
c         ll=INDEX(ro2sp(i),' ')
c         CALL akspnum(ro2sp(i)(1:ll),chrsp,numsp,isp)
c         IF (isp.eq.0) THEN
c           WRITE(6,*) '--error--, RO2 species unidentified in readro2'
c           WRITE(6,'(a)') ro2sp(i)
c           WRITE(6,*) 'in file', filnam
c           STOP
c         ENDIF
c         idchemro2(i,nn)=isp
c        ENDDO

* get the id for the RO2 species and store into idchemro2 table
* the program was changed to increase the search
        ipos=-1

!$OMP PARALLEL DO private(i, found, j, ipos)
        DO i=1,numchemro2(nn)
        
          call akspnum(ro2sp(i),chrsp,numsp,ipos)
!          found = .FALSE.
!          DO j=1,maxsp
!            IF (ro2sp(i).eq.chrsp(j)) THEN
!              ipos=j
!              idchemro2(i,nn)=j
!              found = .TRUE.
!              EXIT
!            ENDIF
!          ENDDO
          IF (ipos == 0) then
            WRITE(6,*) '--error--, RO2 species unidentified in readro2'
            WRITE(6,'(a)') ro2sp(i)
            WRITE(6,*) 'in file', filnam
            WRITE(6,*) 'PLEASE RETURN TO THE PROGRAM AND CHANGE' 
            WRITE(6,*) 'SPECIES MAY NOT BE SORTED' 
            STOP
          endif
          idchemro2(i,nn)=ipos
        ENDDO
!$OMP END PARALLEL DO
        
      ENDDO

* -----------------------
* RETURN
* -----------------------
      
      END 
