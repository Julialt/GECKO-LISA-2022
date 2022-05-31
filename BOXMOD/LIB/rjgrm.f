!--------------------------------------------------------------------------!
! MASTER MECHANISM - ROUTINE NAME : rjgrm                                 !
!                                                                          !
! PURPOSE: Strip ring-joining numerical characters from group strings.     !
!          Subroutine rjgadd puts them back.                               !
!                                                                          !
! INPUT:                                                                   !
! - nring       : number of rings in chem                                  !
!                                                                          !
! IN/OUT:                                                                  !
! - group(i)    : groups at position (node) i                              !
!                                                                          !
! OUTPUT:                                                                  !
! - rjg(j,2)     : group numbers of ring-joining pairs for each ring        !
!                                                                          !
!--------------------------------------------------------------------------!
      SUBROUTINE rjgrm(nring,group,rjg)
      IMPLICIT NONE
      INCLUDE 'general.h'

! input:
      INTEGER    nring
! in/out:
      CHARACTER(lgr) group(mca)
! output:
      INTEGER    rjg(mri,2)
! internal:
      INTEGER    n,i,j,k

!--------------------------------------------------------------------------!
      !print*,'*rjgrm*'

      DO n=1,nring
        rjg(n,1)=0
        rjg(n,2)=0
        DO i=1,mca
!          IF(group(i)(1:1).EQ.' ') GO TO 10  
! causes problems for rings that have expelled products from branches
          IF(group(i)(1:2).EQ.'-O'.OR.group(i)(2:2).EQ.'d')THEN
            j=3
          ELSE
            j=2
          ENDIF
          IF(group(i)(j:j).EQ.digit(n)) THEN
            DO k=j,lgr-1
              group(i)(k:k)=group(i)(k+1:k+1)
            ENDDO
            IF(rjg(n,1).EQ.0)THEN
              rjg(n,1) = i
            ELSE
              rjg(n,2) = i
            ENDIF
          ENDIF
        ENDDO
10      CONTINUE
      ENDDO
      !DO i=1,10
      !  print*,group(i)
      !ENDDO

      END
