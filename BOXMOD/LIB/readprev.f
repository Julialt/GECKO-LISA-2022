***********************************************************************
* subroutine to read existing model output and extract all results    *
* Interpolates to input time & returns values for single time point   *
* may need updating for 2-box situation!
***********************************************************************
      SUBROUTINE readprev(iotime,iochrsp,idsat,ndim,
     &                    oconc,oconcaer, ibox)
      !$ use OMP_LIB
      USE flags_module,ONLY: soa_fg
      USE akparameter_module
      IMPLICIT NONE

! input
      REAL    iotime
      CHARACTER(maxlsp) iochrsp(maxsp)
      INTEGER  idsat(mxsat),ndim, ibox

! output
      !REAL    otemp
      REAL    oconc(maxsp)
      REAL    oconcaer(mxsat)

* maximum number of data in the initialisation file
!      INTEGER,PARAMETER :: mdat=3000

! arrays to read in from initialisation file
      CHARACTER(maxlsp) namextra1,namextra2,namextra3
      CHARACTER(maxlsp) achrsp(mxsat),gchrsp(maxsp)
      REAL    atime(2),atemp(2), arh(2)
      REAL    gconc(2,maxsp)
      REAL    aconc(2,mxsat),aconcaer(2,mxsat)
      REAL    psat(2,mxsat) ! NB: not passed out of subroutine
      REAL    maer

! other internals
      INTEGER numgas,numsat,numextra,lensp,dummy1,dummy2
      REAL    dummy_wmol(mxsat)
      INTEGER i,j,ndat
      REAL    incr

* -----------------------------
* initialize
* -----------------------------
      namextra1=' '
      namextra2=' '
      namextra3=' '
      oconc = 0.
      oconcaer = 0.
      gconc = 0.
      aconc = 0.
      aconcaer = 0.
      atime = 0.
      atemp = 0.
      arh = 0.
      
! no need to read previous ppa file if soa_fg = 2
! because all data is contained in ppf file
      if (soa_fg .ne. 2) then

* read data in the aerosol result file
* ----------------------------

        if(ibox == 1) then
           OPEN(21,file='indat.ppa',STATUS='OLD',FORM='UNFORMATTED')
        else if (ibox == 2) then
           OPEN(21,file='indat.paa',status='old',form='unformatted')
        endif
           print*,"file opened"

* read parameters of the files and check the size
        READ(21) numsat,ndim,numextra,lensp,dummy1,dummy2
        print*, numsat,ndim,numextra,lensp,dummy1,dummy2

        IF (numsat.gt.mxsat) THEN
          WRITE(6,*) '--error--, number of aer species is gt mxsat'
          STOP
        ENDIF

        WRITE(6,*) 'numsat =',numsat

        IF (lensp.ne.maxlsp) THEN
          WRITE(6,*) '--error--, length of the species not correct'
          STOP
        ENDIF

        IF (numextra.ne.2) THEN
          WRITE(6,*) '--error--, numextra should be 2 (time and temp)'
          STOP
        ENDIF

* read header : 'time', 'temperature' and species names
        READ(21) namextra1,namextra2, (achrsp(i),i=1,numsat),
     &    (dummy_wmol(i),i=1,numsat)

        IF (namextra1(1:5).ne.'TIME ') THEN
          WRITE(6,*) '--error--, first header expected to be TIME'
          STOP
        ENDIF

        IF (namextra2(1:6).ne.'TEMPER') THEN
        WRITE(6,*)'--error--, first header expected to be TEMPERATURE'
          STOP
        ENDIF

* check species identities match those from mechanism

        DO i=1,numsat
          IF(achrsp(i).NE.iochrsp(idsat(i))) THEN
            WRITE(6,*) "--error--, species identity doesn't match mech"          
            WRITE(6,*) achrsp(i),iochrsp(idsat(i))
            STOP
          ENDIF
        ENDDO

* read the data in the file
        ndat=0

30      CONTINUE
! index 2 is current time
! index 1 is previous time
        atime(1) = atime(2)
        atemp(1) = atemp(2)
        aconc(1,:) = aconc(2,:)
        aconcaer(1,:) = aconcaer(2,:)
        psat(1,:) = psat(2,:)
        ndat=ndat+1
        READ(21,END=50) atime(2),atemp(2),
     &                (aconc(2,i),i=1,numsat),
     &                (aconcaer(2,i),i=1,numsat),
     &                (psat(2,i),i=1,numsat),maer
!      IF (ndat.gt.mdat) THEN
!        WRITE(6,*) '--error--, number of data exceed mdat'
!        STOP
!      ENDIF
        IF (atime(2).ge.iotime) THEN
         GOTO 50
        ENDIF

        GOTO 30

50      CONTINUE
        CLOSE(21)
        WRITE(6,*) 'number of times read from aerosol ini file =',ndat

        IF (atime(2).LT.iotime) THEN
          WRITE(6,*)
     &  '--error--, input data does not encompass desired time'
          WRITE(6,*) atime(2),iotime
          STOP
        ENDIF

* interpolate values for desired time and
* re-index gas-phase values according to idsat(i)

        incr = (iotime-atime(1))/(atime(2)-atime(1))
        DO i=1,numsat
          oconcaer(i) = aconcaer(1,i) +
     &               incr*(aconcaer(2,i)-aconcaer(1,i))
        ENDDO
      ENDIF

* ----------------------------
* read data in the gas result file
* ----------------------------
      if(ibox == 1) then
         OPEN(21,file='indat.ppf',STATUS='OLD',FORM='UNFORMATTED')
      else if (ibox == 2) then
         open(21,file='indat.pff',status='old',form='unformatted')
      endif
      READ(21) numgas,numextra,lensp,dummy1,dummy2

      IF (numgas.gt.maxsp) THEN
        WRITE(6,*) '--error--, number of gas species is gt maxsp'
        STOP
      ENDIF

      IF (lensp.ne.maxlsp) THEN
        WRITE(6,*) '--error--, length of the species not correct'
        STOP
      ENDIF

      IF (numextra.ne.3) THEN
        WRITE(6,*) '--error--, numextra should be 2 (time and temp)'
        STOP
      ENDIF

* read header : 'time', 'temperature', 'humidity' and species names
      READ(21) namextra1,namextra2, namextra3,(gchrsp(i),i=1,numgas)

      IF (namextra1(1:5).ne.'TIME ') THEN
        WRITE(6,*) '--error--, first header expected to be TIME'
        STOP
      ENDIF

      IF (namextra2(1:6).ne.'TEMPER') THEN
        WRITE(6,*) '--error--, first header expected to be TEMPERATURE'
        STOP
      ENDIF
      
      IF (namextra3(1:6).ne.'HUMIDI') THEN
        WRITE(6,*) '--error--, first header expected to be HUMIDITY'
        STOP
      ENDIF

* read the data in the file
      ndat=0

60    CONTINUE
      ndat=ndat+1
! index 2 is current time
! index 1 is previous time
      atime(1) = atime(2)
      atemp(1) = atemp(2)
      arh(1) = arh(2)
      gconc(1,:) = gconc(2,:)

      READ(21,END=70) atime(2),atemp(2),arh(2),
     &      (gconc(2,i),i=1,numgas)
!      IF (ndat.gt.mdat) THEN
!        WRITE(6,*) '--error--, number of data exceed mdat'
!        STOP
!      ENDIF
      IF (atime(2).ge.iotime) THEN
       GOTO 70
      ENDIF

      GOTO 60

70    CONTINUE
      CLOSE(21)
      WRITE(6,*) 'number of times read from the gas ini file =',ndat

      IF (atime(2).LT.iotime) THEN
        WRITE(6,*)
     &  '--error--, input data does not encompass desired time'
        WRITE(6,*) atime(2),iotime
        STOP
      ENDIF

* interpolate values for desired time

      incr = (iotime-atime(1))/(atime(2)-atime(1))
!$OMP PARALLEL DO private(i)
      DO i=1,numgas
        oconc(i) = gconc(1,i) +
     &               incr*(gconc(2,i)-gconc(1,i))
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
