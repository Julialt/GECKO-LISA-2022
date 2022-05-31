* INCA : subroutine getaux
*
* Purpose : 
* read the line that give the auxiliary information of a given
* special reaction (HV, EXTRA, ...). The line that is read is formatted
* as follow :
*    keyword   /numerical_value_1   numerical_value_2 .... /
* Data are separated by a least 1 space (but can be more). The
* line may contain more than one keyword.
* 
*
* INPUT :
*  maxaux    : maximum number of data that can be given between slashes
*              (and maximum number of keyword available)
*  line      : the line to be read
*  lenline   : maximum length of the line
*  lout      : unit number for the output file
*  nauxpar(i): number of of auxiliary info that must be given for the
*              keyword i
*
* OUTPUT
*  lostop    : logical to stop program if error
*  loaux     : logical
*  xauxcf(i,j): jth data for the ith keyword 
*
*
* Called subroutine : getnpar, errline 
* Called Function   : r_val
**********************************************************************
      SUBROUTINE getaux (
     1  maxaux,
     1  line, lenline, lout,
     1  lostop, loaux, 
     1  nauxpar,
     1  xauxcf)
      IMPLICIT NONE

* input
      INTEGER lenline, lout
      INTEGER maxaux
      INTEGER nauxpar(maxaux)
      CHARACTER*(*) line

* ouput
      REAL xauxcf(0:maxaux,maxaux)
      LOGICAL lostop
      LOGICAL loaux

* local
      INTEGER i, j, ii
      INTEGER iterm, ib, ierr,ie
      INTEGER lenmin, iidaux
      INTEGER npar
      REAL    xcoeff
! ASSUMING COMPILATION OPTION real-8
!      DOUBLE PRECISION    xcoeff
      LOGICAL loslash

      REAL r_val
      DOUBLE PRECISION d_val


* initialize
      loslash =.false.
      iterm=lenline
      ib=0
      ierr=0

      IF (.not.loaux) THEN
        loaux=.true.
        DO i=1,maxaux
          DO j=0,maxaux
            xauxcf(j,i)=0.
          ENDDO
        ENDDO
      ENDIF
 
* find the beginning of the line (non blank character). 
* LABEL 4050 IS THE REENTRY POINT (see the goto statement at the end of
* the file)

      ii=0
4050  CONTINUE
      ii=ii+1
      IF (line(ii:ii).eq.' ') THEN
        IF (ii.lt.lenline) THEN
          GOTO 4050
        ELSE
          RETURN
        ENDIF
      ELSE IF (line(ii:ii).ge.'A'.and.line(ii:ii).le.'Z') THEN
        GOTO 4100
      ELSE
        CALL errline(lout,line,1,iterm,ii)
        WRITE(lout,*)'   --error--   unexpected character',
     &               ' in auxiliary information line'
        WRITE(lout,*)
        lostop=.true.
        RETURN
      ENDIF

* find the keyword involved - Set the ID number (iidaux) corresponding
* to the keyword. If keyword not found (lenmin=0) then error.
*
* REMEMBER : link between keywords and number in xauxcf
*     1 = LOW        |   2 = TROE 
*     3 = not used   |   4 = not used
*     5 = HV         |   6 = EXTRA
*     7 = CVAR       |   8 = AOU
*     9 = WOU        |  10 = WIN
*    11 = ISOM
*    (previously and briefly) 11 = AIN        !  12 = ISOM
* ===================================================================

4100  CONTINUE
      lenmin=0
      IF (line(ii:ii+4).eq.'EXTRA') THEN
        lenmin=5
        iidaux=6
      ELSE IF (line(ii:ii+1).eq.'HV') THEN
        lenmin=2
        iidaux=5
      ELSE IF (line(ii:ii+3).eq.'CVAR') THEN
        lenmin=4
        iidaux=7
      ELSE IF (line(ii:ii+3).eq.'TROE') THEN
        lenmin=4
        iidaux=2
      ELSE IF (line(ii:ii+2).eq.'LOW') THEN
        lenmin=3
        iidaux=1
      ELSE IF (line(ii:ii+2).eq.'AOU') THEN
        lenmin=3
        iidaux=8
      ELSE IF (line(ii:ii+2).eq.'WOU') THEN
        lenmin=3
        iidaux=9
      ELSE IF (line(ii:ii+2).eq.'WIN') THEN
        lenmin=3
        iidaux=10
!      ELSE IF (line(ii:ii+2).eq.'AIN') THEN
!        lenmin=3
!        iidaux=12
      ELSE IF (line(ii:ii+3).eq.'ISOM') THEN
        lenmin=4
        iidaux=11
      ENDIF

      IF (lenmin.eq.0) THEN
        CALL errline(lout,line,1,iterm,ii)
        WRITE (lout,*)'   --error--   auxiliary information'
        WRITE (lout,*)'               could not be identified'
        WRITE (lout,*)' '
        lostop=.true.
        RETURN
      ENDIF

* set a flag ("1.") at index 0 of xauxcf.

      IF (xauxcf(0,iidaux).ne.0.0)THEN
        CALL errline(lout,line,1,iterm,ii)
        WRITE(lout,*)'   --error--   the key-words EXTRA, CVAR,',
     &               '               TROE, LOW, or HV can only be'
        WRITE(lout,*)'               given once for each reaction'
        WRITE(lout,*)
        lostop=.true.
      ELSE
        xauxcf(0,iidaux)=1.
      ENDIF
      ii=ii+lenmin

* find the 2 slahes that must be given after any keyword. "ib" is the
* position after the first "/" and ie is the position before the 
* second "/".

      ib=INDEX(line(ii:),'/') 
      IF (ib.eq.0) loslash=.true.
      ib=ii+ib
      ii=ib

      ie=INDEX(line(ii:),'/') 
      IF (ie.eq.0) loslash=.true.
      ie=ii+ie-2 

      IF (loslash) THEN
        WRITE(lout,*)'   --error--   numerical value  must be',
     &               ' given between two slashes. Slashes not found !'
        WRITE(lout,*)
        lostop=.true.
        STOP
        RETURN
      ENDIF

* find the number of data between the slashes (npar). If that number
* is not equal to the one expected then error. For EXTRA reaction
* (iidaux=6), the number of data is free (between 1 and 10). 
      CALL getnpar(line,ib,ie,npar)
      IF (nauxpar(iidaux).ne.npar) THEN
        !IF (iidaux.eq.6.and.npar.ge.1.and.npar.le.10) GOTO 4350
        IF (iidaux.eq.6.and.npar.ge.1.and.npar.le.11) GOTO 4350
        CALL errline(lout,line,1,iterm,ib+1)
        WRITE(lout,*)'   --error--   the number of parameters given',
     &               ' is incorrect'
        WRITE(lout,*)
        lostop=.true.
        STOP
      ENDIF  

* read the numerical value between the slashes
4350  CONTINUE
      DO i=1,npar
        xcoeff=r_val(line,ib,ie,i,ierr)
        !xcoeff=d_val(line,ib,ie,i,ierr)
        !IF(xcoeff.LT.1.e-30.AND.xcoeff.GT.1.e-60)THEN
        !  print*,xcoeff
        !ENDIF
        IF (ierr.ne.0) THEN
          CALL errline(lout,line,1,iterm,ib+1)
          WRITE(lout,*)'   --error--   while reading the',i,
     &                 'th numerical value'
          WRITE(lout,*)'               of',nauxpar(iidaux),
     &                 ' expected values'
          WRITE(lout,*)
          lostop=.true.
          STOP
        ENDIF
        xauxcf(i,iidaux)=xcoeff
      ENDDO

* read next keyword in the line (if any).
      ii=ie+2
      GOTO 4050

      END
