!-----------------------------------------------------------------!
! PURPOSE: ascertain whether a given bond is part of any ring
!-----------------------------------------------------------------!
      SUBROUTINE findring(i1,i2,nca,bond,rngflg,ring)
      IMPLICIT NONE
      include 'general.h'

! input: 
      INTEGER  i1,i2,nca
      INTEGER  bond(mca,mca)

! output:
      INTEGER    rngflg       ! 0 = 'no', 1 = 'yes'
      INTEGER    ring(mca)    ! =1 if node participates in current ring

! internal:
      INTEGER    i,j,k        ! current, next, previous nodes
      INTEGER  track(mco,mca)
      INTEGER  trlen(mco)
      INTEGER  ntr

! -----------------------------------------------------
      !print*,'*findring*'

      rngflg=0
      DO i=1,mca
        ring(i)=0
      ENDDO

* find tracks starting at i1
      CALL gettrack(bond,i1,nca,ntr,track,trlen)

* search in the tracks if node i2 is found (at least in beta)
      DO i=1,ntr
        DO j=3,trlen(i)
	  IF (track(i,j).eq.i2) THEN
	    rngflg=1
	    DO k=1,j
	      ring(track(i,k))=1
	    ENDDO
c	    GOTO 100
	  ENDIF
	ENDDO
      ENDDO

100   CONTINUE
      
      END
