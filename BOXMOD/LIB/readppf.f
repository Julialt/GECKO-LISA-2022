***********************************************************************
* subroutine to read existing model output and extract all results    *
* Interpolates to input time & returns values for single time point   *
***********************************************************************
      SUBROUTINE readppa(iotime,iochrsp,nsat,idsat,oconcaer)
      USE akparameter_module
      IMPLICIT NONE

! input
      REAL    iotime
      CHARACTER(maxlsp) iochrsp(maxsp)
      INTEGER nsat, idsat(mxsat)

! output
      !REAL    otemp
      !REAL    oconc(maxsp)
      REAL    oconcaer(mxsat)

* maximum number of data in the initialisation file
      INTEGER,PARAMETER :: mdat=80000

! arrays to read in from initialisation file
      CHARACTER(maxlsp) namextra1,namextra2,achrsp(mxsat)
      REAL    atime(mdat),atemp(mdat)
      REAL    aconc(mdat,mxsat),aconcaer(mdat,mxsat)

! other internals
      INTEGER numsp,numextra,lensp,dummy1,dummy2
      INTEGER i,j,ndat
      REAL    incr

      CHARACTER(15) filename

* -----------------------------
* initialize
* -----------------------------
      namextra1=' '
      namextra2=' '
      !oconc = 0.
      oconcaer = 0.

* read data in the result file
* ----------------------------

      OPEN(21,file='indat.ppa',STATUS='OLD',FORM='UNFORMATTED')

* read parameters of the file and check the size
      READ(21) numsp,numextra,lensp,dummy1,dummy2

      IF (numsp.gt.mxsat) THEN
        WRITE(6,*) '--error--, number of species is gt mxsat'
        STOP
      ENDIF

      IF (lensp.ne.maxlsp) THEN
        WRITE(6,*) '--error--, length of the species not correct'
        STOP
      ENDIF

      IF (numextra.ne.2) THEN
        WRITE(6,*) '--error--, numextra should be 2 (time and temp)'
        STOP
      ENDIF

* read header : 'time', 'temperature' and species names
      READ(21) namextra1,namextra2,(achrsp(i),i=1,numsp)
 
      IF (namextra1(1:5).ne.'TIME ') THEN
        WRITE(6,*) '--error--, first header expected to be TIME'
        STOP
      ENDIF

      IF (namextra2(1:6).ne.'TEMPER') THEN
        WRITE(6,*) '--error--, first header expected to be TEMPERATURE'
        STOP
      ENDIF

* check species identities match those from mechanism

      DO i=1,numsp
        IF(achrsp(i).NE.iochrsp(idsat(i))) THEN
          PRINT*,achrsp(i),iochrsp(idsat(i))
          STOP
        ENDIF
      ENDDO

* read the data in the file
      ndat=0

30    CONTINUE
      ndat=ndat+1
      READ(21,END=50) atime(ndat),atemp(ndat),(aconc(ndat,i),i=1,numsp),
     &                (aconcaer(ndat,i),i=1,numsp)
      IF (ndat.gt.mdat) THEN
        WRITE(6,*) '--error--, number of data exceed mdat'
        STOP
      ENDIF
      IF (atime(ndat).gt.iotime) THEN
        WRITE(7,*) atime(ndat),atemp(ndat),(aconc(ndat,i),i=1,numsp),
     &                (aconcaer(ndat,i),i=1,numsp)
       GOTO 50
      ENDIF

      GOTO 30

50    CONTINUE
      !ndat=ndat-1
      WRITE(6,*) 'number of data in the result file=',ndat
      CLOSE(21)

      IF (atime(ndat).LT.iotime) THEN
        WRITE(6,*)  
     &  '--error--, input data does not encompass desired time'
        WRITE(6,*) atime(ndat),iotime 
        STOP
      ENDIF

* interpolate values for desired time and
* re-index gas-phase values according to idsat(i)

      incr = (iotime-atime(ndat-1))/(atime(ndat)-atime(ndat-1))
      !otemp = atemp(ndat-1) + incr*(atemp(ndat)-atemp(ndat-1))
      DO i=1,numsp 
!        oconc(idsat(i)) = aconc(ndat-1,i) + 
!     &               incr*(aconc(ndat,i)-aconc(ndat-1,i))
        oconcaer(i) = aconcaer(ndat-1,i) + 
     &               incr*(aconcaer(ndat,i)-aconcaer(ndat-1,i))
      ENDDO

      !WRITE(7,*) oconc

      RETURN
      END
