!--------------------------------------------------------------------------!
! MASTER MECHANISM - ROUTINE NAME : rjggadd                                 !
!                                                                          !
! PURPOSE: Add ring-joining numerical characters to group strings          !
!          at nodes given by rjg.                                           !
!                                                                          !
! INPUT:                                                                   !
! - nring       : number of rings in chem                                  !
! - rjg(j,2)     : group numbers of ring-joining pairs for each ring        !
!                                                                          !
! IN/OUT:                                                                  !
! - group(i)    : groups at position (node) i                              !
!                                                                          !
!--------------------------------------------------------------------------!
      SUBROUTINE rjgadd(nring,group,rjg)
      IMPLICIT NONE
      INCLUDE 'general.h'

! input:
      INTEGER    nring
      INTEGER    rjg(mri,2)
! in/out:
      CHARACTER(lgr) group(mca)
! internal:
      INTEGER    n,i,ii,j,k
!--------------------------------------------------------------------------!
      !print*,'*rjgadd*'
! loop in reverse order => chars in numerical order if >1 exist at any node
      DO n=nring,1,-1
        DO ii=1,2
          i = rjg(n,ii)
!          IF(group(i)(1:1).EQ.' ') GO TO 10
          IF(i.NE.0)THEN
            IF(group(i)(1:2).EQ.'-O'.OR.group(i)(2:2).EQ.'d')THEN
              j=3
            ELSE
              j=2
            ENDIF
            DO k=lgr,j+1,-1
              group(i)(k:k)=group(i)(k-1:k-1)
            ENDDO
            group(i)(j:j) = digit(n)
          ENDIF
        ENDDO
10      CONTINUE
      ENDDO

      END
