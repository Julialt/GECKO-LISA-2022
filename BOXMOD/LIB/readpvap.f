***********************************************************************
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
      SUBROUTINE readpvap(chrsp,numsp,nsat,namsat,Tb,HBN,tau,idsat)
      USE akparameter_module
      USE io_units_module,ONLY: lread
      IMPLICIT NONE
      
* INPUT 
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      REAL               Tb(mxsat), HBN(mxsat), tau(mxsat)
      INTEGER            nsat, idsat(mxsat)

!      CHARACTER(maxlsp) depnamspe(maxsp)
!      INTEGER  ndepspe, iddepspe(maxsp)
!      REAL     depdatspe(maxsp,3)

* LOCAL
      CHARACTER(40)     line
      INTEGER            i, j, k
      CHARACTER(6)      temp

* open the file
* -------------
      nsat=0
      OPEN (lread, file='pvap.sat')

* read the file
* -------------
      READ (lread, '(a6)') temp
      IF (temp.NE.'MYRDAL') THEN
        WRITE(6,*) 'wrong choice of pvap SAR'
        STOP
      ENDIF

      DO i=1,mxsat
        READ(lread,'(A40)', END=890) line
          IF (line(1:4).eq.'END ') GOTO 896
          nsat=nsat+1
          READ(line,'(A7,2x,f6.1,2x,f7.4,2x,f6.1)',err=890) 
     &               namsat(i),Tb(i),HBN(i),tau(i)
        ENDDO
890     WRITE (6,*) '--ERROR--, END NOT FOUND IN datamyrdal'
        WRITE (6,*) 'error may be that mxsat is underevaluated'
        STOP        'RATE, SHOOT AGAIN'
896     CONTINUE
        CLOSE(lread)

* set index for corresponding namsat names with chrsp names
* ----------------------------------------------------------
        k=1
        DO 696 i=1,nsat
          DO j=k,numsp
            IF (namsat(i).eq.chrsp(j)) THEN
              idsat(i)=j
              k=j
              GOTO 696
            ENDIF
          ENDDO
          WRITE(6,*) '--error--, the following species in datamyr file'
          WRITE(6,*) '           not in list of species in the scheme'
          STOP 'STOP in readpvap subroutine'
696     CONTINUE


* end of the routine
* -----------------
      END 
