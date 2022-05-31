***********************************************************************
* This subroutine reads the files in which data for Pvap estimates    *
* are stored : SIMPOL method                                          *
*                                                                     *
***********************************************************************
      SUBROUTINE readpvap_sim_bin(chrsp,
     &                        numsp,nsat,namsat,
     &                        bk,simpgroup,idsat)

      USE flags_module
      USE akparameter_module
      USE io_units_module,ONLY: lread
      IMPLICIT NONE
      
      
* INPUT 
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      REAL              simpgroup(mxsat,31),bk(31,4)
      INTEGER           nsat, idsat(mxsat)

!      CHARACTER(maxlsp) depnamspe(maxsp)
!      INTEGER  ndepspe, iddepspe(maxsp)
!      REAL     depdatspe(maxsp,3)

* LOCAL
      CHARACTER(105)    line
      INTEGER           i, j, isp
      CHARACTER(6)      temp
      CHARACTER(20)     filnam
      CHARACTER(maxlsp) namtmp

*----------------------------------------           
*  reading data
*----------------------------------------

      OPEN(lread,file='simpol.dat',form='FORMATTED',status='OLD')
      DO i=1,31
        READ(lread,'(f10.4,2x,f9.6,2x,f11.8,2x,f11.8)') 
     &       bk(i,1),bk(i,2),bk(i,3),bk(i,4)
      ENDDO
      CLOSE(lread)
      
* open the file
* -------------
      nsat=0
      filnam = 'psim.sat'
      OPEN (lread, file=filnam)

* read the file
* -------------
      READ (lread, '(a6)') temp
      IF (temp.NE.'SIMPOL') THEN
        WRITE(6,*) 'wrong choice of pvap SAR'
        STOP
      ENDIF

      DO i=1,mxsat
        READ(lread,'(A105)', END=890) line
        IF (line(1:4).eq.'END ') GOTO 896
        nsat=nsat+1
        READ(line,'(A7,2x,31(f3.0))',err=890) 
     &             namsat(i),(simpgroup(i,j),j=1,31)

! USE AKSPNUM TO GET IDSAT)
        IF(soa_fg.EQ.1)THEN
          namtmp(1:1) = "G"
        ELSE
          namtmp(1:1) = "A"
        ENDIF
        namtmp(2:maxlsp) = namsat(i)(2:maxlsp)
        CALL akspnum(namtmp,chrsp,numsp,isp)
        idsat(i) = isp

      ENDDO
890   WRITE (6,*) '--ERROR--, END NOT FOUND IN simpol'
      WRITE (6,*) 'error may be that mxsat is underevaluated'
      STOP        'RATE, SHOOT AGAIN'
896   CONTINUE
      CLOSE(lread)

* end of the routine
* -----------------
      END SUBROUTINE readpvap_sim_bin
