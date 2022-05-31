* INCA : subroutine getspec
*
* Purpose : 
* Given an input line of the form : 
*     species_name  /weight/  species_name  /weight/  ...
* add the species in the chrsp table and the corresponding molecular
* weight in the wmol table. The line may contains more than one species.
* Molecular weight is optional (set to 0 if not given). Note that if
* no molecular weight is given, the program fail when the species are
* checked (see subroutine chkspec)
* 
*
* INPUT :
*  line      : the line to be read
*  lenline   : maximum length of the line
*  lout      : unit number for the output file
*  maxsp     : maximum number of species
*  maxlsp    : maximum length of the species
*
* OUTPUT
*  numsp     : number of species 
*  lostop    : logical to stop program if error
*  wmol(i)   : molecular weight of the species i
*  chrsp(i)  : name of the species i
*
*
* Called subroutine : errline
* Called Function   : r_val
***********************************************************************
      SUBROUTINE getspec (
     1   line, lenline, lout,
     2   maxsp, maxlsp, numsp,
     3   lostop,
     4   wmol, chrsp)

      IMPLICIT NONE

* input
      INTEGER lenline, lout
      INTEGER maxsp, maxlsp

* output
      INTEGER       numsp
      REAL          wmol(maxsp)
      CHARACTER*(*) chrsp(maxsp)
      CHARACTER*(*) line
      LOGICAL       lostop

* local
      INTEGER i
      INTEGER iterm, ib, ilen, ierr
      LOGICAL loslash,lomemo
      REAL    r_val

* initialize
      loslash=.false.
      lomemo =.false.
      iterm=lenline
      ib=0
      ilen=0
      ierr=0

* --------------
* READ THE LINE 
* --------------
* "i" is the position in the line that is tested. If lomemo is true
* then a name of the species is read. If loslash is true, then 
* molecular weight is read. 

      DO 200 i=1,lenline

* molecular weight between slashes
        IF(line(i:i).eq.'/'.or.loslash)THEN
          IF(lomemo)THEN
            lomemo=.false.
            IF(numsp.ge.maxsp)THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   number of species can',
     &                     'not be greater than maxsp=',maxsp
              WRITE(lout,*)
              lostop=.true.
            ELSE IF (ilen.gt.maxlsp) THEN
              CALL errline(lout,line,1,iterm,i-1)
              WRITE(lout,*)'   --error--   length of species name',
     &                     ' cannot be greater than maxlsp=',maxlsp
              WRITE(lout,*)' => CHECK value of lco in general.h '
              WRITE(lout,*)
              lostop=.true.
            ELSE
              numsp=numsp+1
              chrsp(numsp)=line(ib:i-1)
            ENDIF
          ENDIF
          IF (loslash) THEN
            IF (line(i:i).eq.'/') THEN
              loslash=.false.
              IF (numsp.eq.0) THEN
                CALL errline(lout,line,1,iterm,ib)
                WRITE(lout,*)'   --error--  species name must',
     &                       ' precede molecular weight'
                WRITE(lout,*)
                lostop=.true.
                GOTO 200
              ENDIF
              wmol(numsp)=r_val(line,ib+1,i-1,1,ierr)
              IF (ierr.ne.0) THEN
                CALL errline(lout,line,1,iterm,ib+1)
                WRITE(lout,*)'   --error--   while reading the',
     &                       ' molecular weight (from input file)'
                WRITE(lout,*)
                lostop=.true.
                GOTO 200
              ENDIF
              GOTO 200  
            ELSE
              GOTO 200
            ENDIF
          ELSE
            loslash=.true.
            ib=i
            GOTO 200
          ENDIF
        ENDIF

* species name
        IF (line(i:i).eq.' ') THEN
          IF(lomemo)THEN
            lomemo=.false.
            IF(numsp.ge.maxsp)THEN
              CALL errline(lout,line,1,iterm,ib)
              WRITE(lout,*)'   --error--   number of species can',
     &                     ' not be greater than maxsp=',maxsp
              WRITE(lout,*)
              lostop=.true.
              GOTO 200
            ELSE IF(ilen.gt.maxlsp)THEN
              CALL errline(lout,line,1,iterm,i-1)
              WRITE(lout,*)'   --error--   length of species name',
     &                     ' can not be greater than maxlsp=',maxlsp
              WRITE(lout,*)
              lostop=.true.
              GOTO 200
            ELSE 
              numsp=numsp+1
              chrsp(numsp)=line(ib:i-1)
              GOTO 200
            ENDIF
          ELSE
            GOTO 200
          ENDIF
        ELSE
          IF (lomemo) THEN
            ilen=ilen+1
            GOTO 200
          ELSE
            lomemo=.true.
            ilen=1
            ib=i
            GOTO 200
          ENDIF
        ENDIF
200   CONTINUE          

* ------------
* END OF LINE
* ------------

* If the line does not end with a "blank" then check if the reading 
* of a species was still active. In that case put the species in the
* table
      IF (lomemo) THEN
        lomemo=.false.
        IF(numsp.ge.maxsp)THEN
          CALL errline(lout,line,1,iterm,ib)
          WRITE(lout,*)'   --error--   number of species can',
     &                 ' not be greater than maxsp=',maxsp
          WRITE(lout,*)
          lostop=.true.
        ELSE IF(ilen.gt.maxlsp)THEN
          CALL errline(lout,line,1,iterm,i-1)
          WRITE(lout,*)'   --error--   length of species name',
     &                 ' can not be greater than maxlsp=',maxlsp
          WRITE(lout,*)
          lostop=.true.
        ELSE
          numsp=numsp+1
          chrsp(numsp)=line(ib:i-1)
        ENDIF
      ENDIF

* if a slash was open but is still not closed => error
      IF (loslash) THEN
        CALL errline(lout,line,1,iterm,iterm)
        WRITE(lout,*)'   --error--   molecular weight can not span',
     &               ' two input lines'
        WRITE(lout,*)
        loslash=.false.
        lostop=.true.
      ENDIF

      END
