      FUNCTION i_val(line,ibegin,iend,npar,ierr)
      IMPLICIT NONE

      CHARACTER(*) line
      CHARACTER(8)   form
      CHARACTER(2)   itot
      LOGICAL lochar,loapost,loempty
      INTEGER i,j,i_val,ierr,ipar,istart,ibegin,iend,npar,istop

      lochar =.false.
      loapost=.false.
      loempty=.true.

      i_val=0
      ierr=0
      ipar=0

      istart=ibegin
      DO 100 i=ibegin,iend
        IF (line(i:i).EQ.char(39).or.lochar) THEN 
          IF (lochar) THEN 
            IF (line(i:i).EQ.char(39)) THEN 
              IF(loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF
            ELSE 
              IF (loapost) THEN 
                loapost=.false.
                lochar=.false.
                istop=i-2
                ipar=ipar+1
                IF(ipar.EQ.npar)GOTO 900
                istart=i
                GOTO 50
              ELSE
                GOTO 100
              ENDIF
            ENDIF
          ELSE
            istart=i+1
            lochar=.true.
            GOTO 100
          ENDIF
        ENDIF
c
50      IF(line(i:i).EQ.',') THEN 
          IF(i.EQ.1)GOTO 70
          IF(line(i-1:i-1).ne.' ')GOTO 70
          DO 60 j=i-2,1,-1
            IF(line(j:j).EQ.' ')GOTO 60
            IF(line(j:j).EQ.',')GOTO 70
            GOTO 100 
60        CONTINUE
70        loempty=.true.
          istop=i-1
          ipar=ipar+1
          IF(ipar.EQ.npar) GOTO 200
          istart=i+1
          GOTO 100
        ELSE IF(line(i:i).EQ.' ') THEN 
          IF(loempty) THEN 
            istart=i+1
            GOTO 100
          ELSE 
            loempty=.true.
            istop=i-1
            ipar=ipar+1
            IF(ipar.EQ.npar) GOTO 200
            istart=i+1
            GOTO 100
          ENDIF
        ELSE
          loempty=.false.
          GOTO 100
        ENDIF
100   CONTINUE
      IF(lochar)GOTO 900
      IF(loempty)GOTO 800
      istop=iend
      ipar=ipar+1
      IF(ipar.ne.npar)GOTO 800
c
200   CONTINUE
      WRITE(itot,'(i2.2)',err=930) istop-istart+1
      IF(itot.EQ.'00')GOTO 940
      form='(i' // itot // ')'
      read(line(istart:istop),form,err=950) i_val
      return
800   ierr=800
      GOTO 999
900   ierr=900
      GOTO 999
930   ierr=930
      GOTO 999
940   ierr=940
      GOTO 999
950   ierr=950
999   END

* --------------------------------------------------------------

      FUNCTION r_val(line,ibegin,iend,npar,ierr)
      IMPLICIT NONE

      CHARACTER(*) line
      CHARACTER(8)   form
      CHARACTER(2)   itot,idec
      LOGICAL lochar,loapost,loempty,lonumber,lopoint,loeee
      INTEGER i,j,ierr,ipar,istart,ibegin,iend,npar,istop,ipoint,ieee
      REAL    r_val

      lochar  =.false.
      loapost =.false.
      loempty =.true.
      lonumber=.false.
      lopoint =.false.
      loeee   =.false.

      r_val=0.0
      ierr=0
      ipar=0
      istart=ibegin

      DO 100 i=ibegin,iend
        IF (line(i:i).EQ.char(39).or.lochar) THEN 
          IF (lochar) THEN 
            IF (line(i:i).EQ.char(39)) THEN 
              IF (loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF
            ELSE 
              IF (loapost) THEN 
                loapost=.false.
                lochar=.false.
                istop=i-2
                ipar=ipar+1
                IF (ipar.EQ.npar) GOTO 900
                istart=i
                GOTO 50
              ELSE
                GOTO 100
              ENDIF
            ENDIF
          ELSE
            istart=i+1
            lochar=.true.
            GOTO 100
          ENDIF
        ENDIF
c
50      IF (line(i:i).EQ.',') THEN 
          IF (i.EQ.1)GOTO 70
          IF (line(i-1:i-1).ne.' ') GOTO 70
          DO 60 j=i-2,1,-1
            IF (line(j:j).EQ.' ') GOTO 60
            IF (line(j:j).EQ.',') GOTO 70
            GOTO 100 
60        CONTINUE
70        loempty=.true.
          istop=i-1
          ipar=ipar+1
          IF (ipar.EQ.npar) GOTO 200
          istart=i+1
          GOTO 100
        ELSE IF (line(i:i).EQ.' ') THEN 
          IF (loempty) THEN
            istart=i+1
            GOTO 100
          ELSE 
            loempty=.true.
            istop=i-1
            ipar=ipar+1
            IF(ipar.EQ.npar) GOTO 200
            istart=i+1
            GOTO 100
          ENDIF
        ELSE
          loempty=.false.
          GOTO 100
        ENDIF
100   CONTINUE
      IF (lochar) GOTO 900
      IF (loempty) GOTO 800
      istop=iend
      ipar=ipar+1
      IF (ipar.ne.npar) GOTO 800
c
200   DO 300 i=istart,istop
        IF (line(i:i).ge.'0'.and.line(i:i).le.'9') THEN
          lonumber=.true.
          GOTO 300
        ELSE IF (line(i:i).EQ.'+'.or.line(i:i).EQ.'-') THEN
          GOTO 300
        ELSE IF (line(i:i).EQ.'.') THEN
          IF(lopoint) GOTO 930
          lopoint=.true.
          ipoint=i
        ELSE IF (line(i:i).EQ.'e'.or.line(i:i).EQ.'E') THEN
          IF (loeee) GOTO 940
          loeee=.true.
          ieee=i
        ELSE
          GOTO 950
        ENDIF
300   CONTINUE
c
      IF (.not.lonumber) GOTO 960
      WRITE (itot,'(i2.2)',err=970) istop-istart+1
      IF (loeee) THEN
        IF (lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) ieee-ipoint-1
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(e' // itot // '.' // idec // ')'
      ELSE
        IF (lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) istop-ipoint
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(f' // itot // '.' // idec // ')'
      ENDIF
      READ(line(istart:istop),form,err=990) r_val
      RETURN
800   ierr=800
      GOTO 999
900   ierr=900
      GOTO 999
930   ierr=930
      GOTO 999
940   ierr=940
      GOTO 999
950   ierr=950
      GOTO 999
960   ierr=960
      GOTO 999
970   ierr=970
      GOTO 999
980   ierr=980
      GOTO 999
990   ierr=990
999   END

* --------------------------------------------------------------

      SUBROUTINE getnpar(line,ibegin,iend,npar)
      IMPLICIT NONE

      CHARACTER(*) line
      LOGICAL lochar,loapost,loempty
      INTEGER i, j, ibegin, iend, npar, ierr

      lochar =.false.
      loapost=.false.
      loempty=.true.
      ierr=0
      npar=0

      DO 100 i=ibegin,iend
        IF (line(i:i).EQ.char(39).OR.lochar) THEN 
          IF (lochar) THEN 
            IF (line(i:i).EQ.char(39)) THEN 
              IF (loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF
            ELSE 
              IF (loapost) THEN 
                loapost=.false.
                lochar=.false.
                npar=npar+1
                GOTO 50
              ELSE
                GOTO 100
              ENDIF
            ENDIF
          ELSE
            lochar=.true.
            GOTO 100
          ENDIF
        ENDIF
c
50      IF (line(i:i).EQ.',') THEN 
          IF (i.EQ.1)GOTO 70
          IF (line(i-1:i-1).ne.' ')GOTO 70
          DO 60 j=i-2,1,-1
            IF (line(j:j).EQ.' ')GOTO 60
            IF (line(j:j).EQ.',')GOTO 70
            GOTO 100 
60        CONTINUE
70        loempty=.true.
          npar=npar+1
          GOTO 100
        ELSE IF (line(i:i).EQ.' ') THEN 
          IF (loempty) THEN 
            GOTO 100
          ELSE 
            loempty=.true.
            npar=npar+1
            GOTO 100
          ENDIF
        ELSE
          loempty=.false.
          GOTO 100
        ENDIF
100   CONTINUE
      IF (lochar) GOTO 999
      IF (loempty) GOTO 999
      npar=npar+1
999   RETURN
      END

* --------------------------------------------------------------

      SUBROUTINE cleanline(inline,outline,lenline)
      IMPLICIT NONE

      CHARACTER(*) inline,outline
      INTEGER lenline, i
      
      DO i=1,lenline
        outline(i:i)=' '
      ENDDO

      DO 200 i=1,lenline
        IF (inline(i:i).EQ.'!') THEN
          RETURN
        ELSE IF (inline(i:i).GE.'a'.AND.inline(i:i).LE.'z') THEN
          outline(i:i)=char(ichar(inline(i:i))-32)
        ELSE
          outline(i:i)=inline(i:i)
        ENDIF
200   CONTINUE
      END

* --------------------------------------------------------------

      function lenstr(line,maxchr)
      IMPLICIT NONE

      CHARACTER(*) line
      INTEGER  i, maxchr, lenstr
      DO 10 i=maxchr,1,-1
        IF(line(i:i).ne.' ')THEN
          lenstr=i
          RETURN
        ENDIF
10    CONTINUE
      lenstr=0
      END
