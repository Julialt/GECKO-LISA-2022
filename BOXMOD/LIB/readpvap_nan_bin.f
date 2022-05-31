**********************************************************************
* This subroutine read the file in which data for Pvap estimate       *
* are stored (datamyr.sat).                                           *
*                                                                     *
* INPUT :                                                             *
*    -lout                                                            *
*    - chrsp                                                          *
*    - numsp                                                          *
*                                                                     *
* OUTPUT :                                                            *
*    - nsat : number of species in the datamyr.sat file               *
*    - namsat(i) : names of the species for which Pvap can be         *
*                  computed                                           *
*    - Tb(i)  : The boiling point of chem, computed by the Joback     *
*               group contribution method                             *
*    - HBN(i) : Hydrogen bond number                                  *
*    - tau(i) : effective number of torsional bonds                   *
*    - idsat(i) : index of species i in the full list od species      *
*                 in the scheme                                       *
***********************************************************************
      SUBROUTINE readpvap_nan_bin(chrsp,
     &                            numsp,nsat,namsat,
     &                            Tb,dB,idsat)

      USE flags_module
      USE akparameter_module
      USE io_units_module,ONLY: lread
      IMPLICIT NONE
      
* INPUT 
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      REAL               Tb(mxsat)
      REAL               dB(mxsat)
      INTEGER            nsat, idsat(mxsat)

!      CHARACTER(maxlsp) depnamspe(maxsp)
!      INTEGER  ndepspe, iddepspe(maxsp)
!      REAL     depdatspe(maxsp,3)

* LOCAL
      CHARACTER(40)     line
      INTEGER           i, isp
      CHARACTER(6)      temp
      CHARACTER(20)     filnam
      CHARACTER(maxlsp) namtmp

* open the file
* -------------
      nsat=0
      filnam = 'pnan.sat'
      OPEN (lread, file=filnam)

* read the file
* -------------
      READ (lread, '(a6)') temp
      IF (temp.NE.'NANNOO') THEN
        WRITE(6,*) 'wrong choice of pvap SAR'
        STOP
      ENDIF

      DO i=1,mxsat
        dB(i)=0
        READ(lread,'(A40)', END=890) line
        IF (line(1:4).eq.'END ') GOTO 896
        nsat=nsat+1
        READ(line,'(A7,2x,f6.1,2x,f8.4)',err=890) 
     &             namsat(i),Tb(i),dB(i)

! USE AKSPNUM TO GET IDSAT
        IF(soa_fg.EQ.1)THEN
          namtmp(1:1) = "G"
        ELSE
          namtmp(1:1) = "A"
        ENDIF
        namtmp(2:maxlsp) = namsat(i)(2:maxlsp)
        CALL akspnum(namtmp,chrsp,numsp,isp)
        idsat(i) = isp

      ENDDO

890     WRITE (6,*) '--ERROR--, END NOT FOUND IN nannoolal'
        WRITE (6,*) 'error may be that mxsat is underevaluated'
        STOP        'RATE, SHOOT AGAIN'
896     CONTINUE
        CLOSE(lread)

* end of the routine
* -----------------
      END 
