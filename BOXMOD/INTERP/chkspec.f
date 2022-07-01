* INCA : subroutine chkspec
*
* Purpose : 
*   check that the species given as input are correctly set (no
*   illegal character, no duplicate name, no keyword and molecular
*   weight).
*
* INPUT :
*  lout      : unit number for the output file
*  maxsp     : maximum number of species
*  maxlsp    : maximum length of the species
*  chrsp(i)  : name of the species i
*  wmol(i)   : molecular weight of the species i
!  sortsp(j) : sorted list of the species
!  isortlk(j): link between the sorted list and the unsorted list
! from module sorting:
!  srtid_chrsp : sorted index of chrsp
*
* OUTPUT
*  numsp     : number of species 
*  lostop    : logical to stop program if error
*  wmol(i)   : molecular weight of the species i
*  chrsp(i)  : name of the species i
*
***********************************************************************
      SUBROUTINE chkspec (
     1   lout,
     2   maxsp, numsp, maxlsp,
     3   wmol, chrsp,
     4   loreply, lostop)
!     &   sortsp,isortlk,

      use sorting, only: srtid_chrsp
      IMPLICIT NONE
      INCLUDE 'general.h'

* input
      INTEGER       lout
      INTEGER       maxsp
      INTEGER       numsp, maxlsp
!      INTEGER       isortlk(maxsp)
      REAL          wmol(maxsp)
      CHARACTER*(*) chrsp(maxsp)
!      CHARACTER*(*) sortsp(maxsp)
      LOGICAL       loreply 

* output
      LOGICAL       lostop

* local
      INTEGER       i, j, k
      LOGICAL       loerr
      LOGICAL       locheck(maxsp)


* initialize
* -----------
      DO i=1,numsp
        locheck(i)=.false.
      ENDDO

* check if the name of the species contain illegal character
* ----------------------------------------------------------
* The name of a species must start with A-Z or '(<[{'. All other
* character are allowed, except '+-='

      WRITE(6,*) '        .... check character in species'
      DO 1320 i=1,numsp
        DO 1310 j=1,maxlsp

          IF (chrsp(i)(j:j).ge.'a'.and.chrsp(i)(j:j).le.'z') GOTO 1310
          IF (chrsp(i)(j:j).ge.'A'.and.chrsp(i)(j:j).le.'Z') GOTO 1310
          IF (chrsp(i)(j:j).eq.'('.or .chrsp(i)(j:j).eq.'<'.or.
     &       chrsp(i)(j:j).eq.'['.or .chrsp(i)(j:j).eq.'{') GOTO 1310

          IF (j.eq.1) GOTO 1305
          IF (chrsp(i)(j:j).eq.'=') GOTO 1305
          IF (chrsp(i)(j:j).eq.'+') GOTO 1305
!JMLT, 2021: dashes are ALLOWED in MECHGEN names
          IF (chrsp(i)(j:j).eq.'-'.AND.lco.NE.8) GOTO 1305

          IF (chrsp(i)(j:j).ge.' '.and.chrsp(i)(j:j).le.'~') GOTO 1310

1305      WRITE(lout,*)
          WRITE(lout,*)'   --error--   species name "',chrsp(i),'"',
     &                 ' contains an illegal character'
          WRITE(lout,*)
          locheck(i)=.true.
          lostop=.true.
          GOTO 1320

1310    CONTINUE
1320  CONTINUE


* check if a species is listed more than once
* -------------------------------------------

      WRITE(6,*) '        .... search duplicate species'
!       DO i=1,numsp-1
!         IF (sortsp(i).eq.sortsp(i+1)) THEN
!           WRITE(lout,*)
!           WRITE(lout,*)'   --error--   species "',sortsp(i),'"',
!      &              ' is listed more than once'
!           WRITE(lout,*)
!           locheck(isortlk(i))=.true.
!           locheck(isortlk(i+1))=.true.
!           lostop=.true.
!         ENDIF
!       ENDDO

      DO i=1, numsp -1
!        PRINT*,chrsp(srtid_chrsp(i))
        IF (chrsp(srtid_chrsp(i)) == chrsp(srtid_chrsp(i + 1))) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   species "',
     &      chrsp(srtid_chrsp(i)),'"',' is listed more than once'
          WRITE(lout,*)
          locheck(srtid_chrsp(i))=.true.
          locheck(srtid_chrsp(i+1))=.true.
          lostop=.true.          
        ENDIF
      ENDDO


* check if any illegal name (keyword) is used
* -------------------------------------------

      WRITE(6,*) '        .... search if species=keyword'
      DO i=1,numsp
        loerr=.false.
        IF (chrsp(i)(1:6).eq.'TBODY ')   loerr=.true.
        IF (chrsp(i)(1:3).eq.'HV ')      loerr=.true.
        IF (chrsp(i)(1:8).eq.'FALLOFF ') loerr=.true.
        IF (chrsp(i)(1:5).eq.'CVAR ')    loerr=.true.
        IF (chrsp(i)(1:6).eq.'EXTRA ')   loerr=.true.

        IF (loerr) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   species name "',chrsp(i),'"',
     &                 ' is also a keyword'
          WRITE(lout,*)
          locheck(i)=.true.
          lostop=.true.
        ENDIF
      ENDDO

* check if the molecular weight was set for every species
* -------------------------------------------------------

      WRITE(6,*) '        .... check molecular weight'
      DO  i=1,numsp
        IF (wmol(i).le.0.0) THEN
          WRITE(lout,*)
          WRITE(lout,*)'   --error--   species "',i,'"',
     &                 ' has no molecular weight'
          WRITE(lout,*)
        ENDIF 
      ENDDO 
                
* write the species to the output file
* -------------------------------------

      IF (loreply) THEN
        WRITE(lout,*)
        WRITE(lout,*)' the used species are:'
        DO i=1,numsp
          IF (.not.locheck(i)) WRITE(lout,'(6x,i6,a20,f6.2,$)')
     &                   i,chrsp(i),wmol(i)
          WRITE(lout,*)
        ENDDO
        WRITE(lout,*)
      ENDIF

      END
