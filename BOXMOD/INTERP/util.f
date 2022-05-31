**********************************************************************
* THIS FILES GIVES THE FUNCTIONS AND SUBROUTINES THAT ARE USED
* IN INCA, ESPECIALLY :
*
*      FUNCTION time_diff(time1,time2)
*      INTEGER FUNCTION search(aseek,alist,nlist)
*      SUBROUTINE errline(lout,line,ib,ie,ipos)
*      SUBROUTINE sort(ns,s)
*      SUBROUTINE cleanline (inline, outline, lenline)   
*      SUBROUTINE getnpar(line, ibegin, iend, npar)
*      FUNCTION lenstr(line,maxchr)
*      FUNCTION r_val(line,ibegin,iend,npar,ierr)
**********************************************************************


*
* ======================== FUNCTION time_diff =========================
*

* this fuction return the time difference between time1 and time2
      FUNCTION time_diff(time1,time2)
      IMPLICIT NONE

* input 
      CHARACTER*10 time1,time2

* output      
      REAL time_diff

* internal 
      REAL hr, min, sec, elapsed

      READ(time1(1:2),*)  hr
      READ(time1(3:4),*)  min
      READ(time1(5:10),*) sec

      elapsed = 60.*(min + 60.*hr) + sec

      READ(time2(1:2),*)  hr
      READ(time2(3:4),*)  min
      READ(time2(5:10),*) sec

      time_diff = 60.*(min + 60.*hr) + sec - elapsed
      RETURN
      END

*
* ======================== FUNCTION search  ==========================
*

      INTEGER FUNCTION search(aseek,alist,nlist)
      USE akparameter_module
      IMPLICIT NONE

* input:
      INTEGER            nlist
      CHARACTER*(maxlsp) aseek, alist(nlist)

* internal:
      INTEGER jhi, jlo, jold, j

* initialize:
      search = 0
      jold = 0
      jlo  = 1
      jhi  = nlist + 1

10    j    = (jhi+jlo)/2
       
      IF(j.EQ.jold) GO TO 40
      jold = j
      IF(aseek.GT.alist(j)) GO TO 20
      IF(aseek.EQ.alist(j)) GO TO 30
      jhi  = j
      GO TO 10
      
20    jlo  = j

      go to 10

30    search = j
      RETURN

40    search = -j

* end of SRH5
      RETURN
      END

*
* ======================== FUNCTION sort  =============================
*

* simple bubble sort
      SUBROUTINE sort(ns,s)
      USE akparameter_module
      IMPLICIT NONE

* input/output      
      INTEGER    ns
      CHARACTER*(maxlsp) s(ns)
* internal
      CHARACTER*(maxlsp) store
      INTEGER i,j

20    i=1
30    j=i+1
      IF (s(i).LE.s(j)) GO TO 10
      store = s(j)  
      s(j)  = s(i)
      s(i)  = store
      i     = i-1
      IF (i.EQ.0) GO TO 20
      GO TO 30
10    IF (j.EQ.ns) GO TO 40
      i = i + 1
      GO TO 30
40    RETURN
      END
*
* ===================== SUBROUTINE errline ============================
*

      SUBROUTINE errline(lout,line,ib,ie,ipos)
      IMPLICIT NONE

* input/output
      INTEGER lout
      INTEGER ipos, ib, ie
      CHARACTER*(*) line

* local
      INTEGER i
      CHARACTER*8 form1,form2

* write lines
      WRITE(lout,*)
      WRITE(form1,'(a3,i3.3,a2)') '( a',ie-ib+1,' )'
      WRITE(lout,form1) line(ib:ie)

      IF (ipos.eq.0) THEN
        WRITE(lout,*)
      ELSE
        WRITE(form2,'(a2,i3.3,a3)') '( ',ipos   ,'a )'
        WRITE(lout,form2) (' ',i=1,ipos-1),'^'
      ENDIF

      END

*
* ========================= SUBROUTINE cleanline ==============
*

      SUBROUTINE cleanline (inline, outline, lenline)   
      IMPLICIT NONE

* input/output
      INTEGER lenline
      CHARACTER*(*) inline, outline

* local
      INTEGER i

* initialize
      DO i=1,lenline
        outline(i:i)=' '
      ENDDO

*
      i=INDEX(inline,'!')
      IF (i.eq.0) THEN
        outline(1:lenline)=inline(1:lenline)
        RETURN
      ELSE
        outline(1:i-1)=inline(1:i-1)
        RETURN
      ENDIF
c* change lowercase in uppercase
c      DO 200 i=1,lenline
c        IF(inline(i:i).eq.'!')THEN
c          RETURN
c        ELSE IF(inline(i:i).ge.'a'.and.inline(i:i).le.'z')THEN
c          outline(i:i)=char(ichar(inline(i:i))-32)
c        ELSE
c          outline(i:i)=inline(i:i)
c        ENDIF
c200   CONTINUE

      END


*
* ========================== SUBROUTINE getnpar =================
*

      SUBROUTINE getnpar(line, ibegin, iend, npar)
      IMPLICIT NONE

* input/output
      INTEGER ibegin, iend, npar
      CHARACTER*(*) line

* local
      INTEGER i, j
      LOGICAL lochar, loapost, loempty

* initialize
      lochar =.false.
      loapost=.false.
      loempty=.true.
      npar=0

* 
      DO 100 i=ibegin,iend

        IF(line(i:i).eq.char(39).or.lochar) THEN 
          IF(lochar) THEN 
            IF(line(i:i).eq.char(39)) THEN 
              IF(loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF
            ELSE 
              IF(loapost) THEN 
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

50      CONTINUE
        IF(line(i:i).eq.',') THEN 
          IF(i.eq.1)GOTO 70
          IF(line(i-1:i-1).ne.' ')GOTO 70
          DO 60 j=i-2,1,-1
            IF(line(j:j).eq.' ')GOTO 60
            IF(line(j:j).eq.',')GOTO 70
            GOTO 100 
60        CONTINUE
70        loempty=.true.
          npar=npar+1
          GOTO 100
        ELSE IF(line(i:i).eq.' ') THEN 
          IF(loempty) THEN 
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

      IF(lochar)GOTO 999
      IF(loempty)GOTO 999
      npar=npar+1

999   RETURN

      END

*
* ========================= FUNCTION lenstr ==================
*

      FUNCTION lenstr(line,maxchr)
      IMPLICIT NONE

* input/output
      INTEGER maxchr, lenstr
      CHARACTER*(*) line

* local
      INTEGER i

*
      DO i=maxchr,1,-1
        IF(line(i:i).ne.' ')THEN
          lenstr=i
          RETURN
        ENDIF
      ENDDO

      lenstr=0

      END

*
* ========================  FUNCTION r_val  ======================
*
      FUNCTION r_val(line,ibegin,iend,npar,ierr)
      IMPLICIT NONE
      INCLUDE 'general.h'

* input/output
       INTEGER ibegin, iend, npar, ierr
       REAL r_val
       CHARACTER*(*) line

* local
      INTEGER ieee, ipoint, j
      INTEGER ipar, i, istop
      INTEGER istart
      CHARACTER(LEN=lco+2)   form
      CHARACTER*2   itot,idec
      LOGICAL lochar,loapost,loempty,lonumber,lopoint,loeee

* initialize
      lochar  =.false.
      loapost =.false.
      loempty =.true.
      lonumber=.false.
      lopoint =.false.
      loeee   =.false.

      r_val=0.
      ierr=0
      ipar=0
      istart=ibegin

*
      DO 100 i=ibegin,iend

        IF(line(i:i).eq.char(39).or.lochar) THEN 

          IF(lochar) THEN

            IF(line(i:i).eq.char(39)) THEN 
              IF(loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF

            ELSE 
              IF(loapost) THEN 
                loapost=.false.
                lochar=.false.
                istop=i-2
                ipar=ipar+1
                IF(ipar.eq.npar)GOTO 900
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

50      CONTINUE

        IF(line(i:i).eq.',') THEN 

          IF(i.eq.1) GOTO 70
          IF(line(i-1:i-1).ne.' ') GOTO 70
          DO 60 j=i-2,1,-1
            IF(line(j:j).eq.' ') GOTO 60
            IF(line(j:j).eq.',') GOTO 70
            GOTO 100 
60        CONTINUE
70        CONTINUE
          loempty=.true.
          istop=i-1
          ipar=ipar+1
          IF(ipar.eq.npar) GOTO 200
          istart=i+1
          GOTO 100

        ELSE IF(line(i:i).eq.' ') THEN 
          IF(loempty) THEN
            istart=i+1
            GOTO 100
          ELSE 
            loempty=.true.
            istop=i-1
            ipar=ipar+1
            IF(ipar.eq.npar) GOTO 200
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

200   CONTINUE

      DO 300 i=istart,istop
        IF(line(i:i).ge.'0'.and.line(i:i).le.'9') THEN
          lonumber=.true.
          GOTO 300
        ELSE IF(line(i:i).eq.'+'.or.line(i:i).eq.'-') THEN
          GOTO 300
        ELSE IF(line(i:i).eq.'.') THEN
          IF(lopoint) GOTO 930
          lopoint=.true.
          ipoint=i
        ELSE IF(line(i:i).eq.'E'.or.line(i:i).eq.'e') THEN
          IF(loeee) GOTO 940
          loeee=.true.
          ieee=i
        ELSE
          GOTO 950
        ENDIF
300   CONTINUE

      IF(.not.lonumber) GOTO 960
      WRITE(itot,'(i2.2)',err=970) istop-istart+1

      IF(loeee) THEN
        IF(lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) ieee-ipoint-1
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(E' // itot // '.' // idec // ')'
      ELSE
        IF(lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) istop-ipoint
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(F' // itot // '.' // idec // ')'
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

*
* ========================  FUNCTION d_val  ======================
*

      FUNCTION d_val(line,ibegin,iend,npar,ierr)
      IMPLICIT NONE

* input/output
       INTEGER ibegin, iend, npar, ierr
       DOUBLE PRECISION d_val
       CHARACTER*(*) line

* local
      INTEGER ieee, ipoint, j
      INTEGER ipar, i, istop
      INTEGER istart
      CHARACTER*8   form
      CHARACTER*2   itot,idec
      LOGICAL lochar,loapost,loempty,lonumber,lopoint,loeee

* initialize
      lochar  =.false.
      loapost =.false.
      loempty =.true.
      lonumber=.false.
      lopoint =.false.
      loeee   =.false.

      d_val=0.
      ierr=0
      ipar=0
      istart=ibegin

*
      DO 100 i=ibegin,iend

        IF(line(i:i).eq.char(39).or.lochar) THEN 

          IF(lochar) THEN

            IF(line(i:i).eq.char(39)) THEN 
              IF(loapost) THEN 
                loapost=.false.
                GOTO 100
              ELSE
                loapost=.true.
                GOTO 100
              ENDIF

            ELSE 
              IF(loapost) THEN 
                loapost=.false.
                lochar=.false.
                istop=i-2
                ipar=ipar+1
                IF(ipar.eq.npar)GOTO 900
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

50      CONTINUE

        IF(line(i:i).eq.',') THEN 

          IF(i.eq.1) GOTO 70
          IF(line(i-1:i-1).ne.' ') GOTO 70
          DO 60 j=i-2,1,-1
            IF(line(j:j).eq.' ') GOTO 60
            IF(line(j:j).eq.',') GOTO 70
            GOTO 100 
60        CONTINUE
70        CONTINUE
          loempty=.true.
          istop=i-1
          ipar=ipar+1
          IF(ipar.eq.npar) GOTO 200
          istart=i+1
          GOTO 100

        ELSE IF(line(i:i).eq.' ') THEN 
          IF(loempty) THEN
            istart=i+1
            GOTO 100
          ELSE 
            loempty=.true.
            istop=i-1
            ipar=ipar+1
            IF(ipar.eq.npar) GOTO 200
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

200   CONTINUE

      DO 300 i=istart,istop
        IF(line(i:i).ge.'0'.and.line(i:i).le.'9') THEN
          lonumber=.true.
          GOTO 300
        ELSE IF(line(i:i).eq.'+'.or.line(i:i).eq.'-') THEN
!        ELSE IF(line(i:i).eq.'+'.or.
!     &          (line(i:i).eq.'-'.AND.lco.NE.8)) THEN
          GOTO 300
        ELSE IF(line(i:i).eq.'.') THEN
          IF(lopoint) GOTO 930
          lopoint=.true.
          ipoint=i
        ELSE IF(line(i:i).eq.'E'.or.line(i:i).eq.'e') THEN
          IF(loeee) GOTO 940
          loeee=.true.
          ieee=i
        ELSE
          GOTO 950
        ENDIF
300   CONTINUE

      IF(.not.lonumber) GOTO 960
      WRITE(itot,'(i2.2)',err=970) istop-istart+1

      IF(loeee) THEN
        IF(lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) ieee-ipoint-1
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(E' // itot // '.' // idec // ')'
      ELSE
        IF(lopoint) THEN
          WRITE(idec,'(i2.2)',err=980) istop-ipoint
        ELSE
          WRITE(idec,'(i2.2)',err=980) 0
        ENDIF
        form='(F' // itot // '.' // idec // ')'
      ENDIF

      READ(line(istart:istop),form,err=990) d_val

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
