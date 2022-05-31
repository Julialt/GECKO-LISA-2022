*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* constituant le (les) compteurs RO2 - schema chimique jenkin     *
* routine a reecrire par la suite                                 *
*******************************************************************
      SUBROUTINE readro2jenk(chrsp,numsp,
     1                   ro2sp,nclro2, numchemro2,idchemro2,cro2)
      USE akparameter_module
      IMPLICIT NONE
      
* INPUT 
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) ro2sp(maxsp,maxro2)
      INTEGER  nclro2, numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
      REAL     cro2(maxro2)

* LOCAL
      CHARACTER(maxlsp) liner
      INTEGER i,j,isp,ll

* set the number of ro2 class
      nclro2=1
      IF (nclro2.gt.maxro2) THEN
        WRITE(6,*) '--error--, in spreadro2.f (read comment)'
        STOP
      ENDIF

* initialize
      DO i=1,mxro2cl
        DO j=1,maxro2
          ro2sp(i,j)=' '
          idchemro2(i,j)=0
        ENDDO
      ENDDO      

      DO i=1,maxro2
        numchemro2(i)=0
        cro2(i)=0.
      ENDDO

* open all file
      OPEN (12, file='indat.ro2',STATUS='OLD')

      i=0
56    READ(12,'(a)',END=59) liner
      IF (liner(1:2).EQ.'**') GOTO 59
      i=i+1
      ll=index(liner,' ')
      IF (ll.LE.1) THEN
        WRITE(6,*) '--error--, lecture RO2'
        STOP
      ENDIF
      ro2sp(i,1)=liner
      GOTO 56
59    CONTINUE 
      CLOSE(12)
      numchemro2(1)=i
      
* get the id for the RO2 species and store into ro2num table
      DO i=1,numchemro2(1)
       ll=INDEX(ro2sp(i,1),' ')
       CALL akspnum(ro2sp(i,1)(1:ll),chrsp,numsp,isp)
       IF (isp.eq.0) THEN
         WRITE(6,*) '--error--, RO2 species unidentified'
         WRITE(6,'(a)') ro2sp(i,1)
         STOP
       ENDIF
       idchemro2(i,1)=isp
      ENDDO

      END 
