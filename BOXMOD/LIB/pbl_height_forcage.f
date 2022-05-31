      SUBROUTINE pbl_height_forcage!(lout, time, nhd, htop,boxt,boxh,
!     &    height, dhdt)

      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: timemod
      USE forcing_params_module,ONLY: nhd,htop,boxt,boxh,height,dhdt

      IMPLICIT NONE

      INTEGER :: i
!-------------------------------------------------------
* first part: calculate boxes heights and dhdt

      do i=1, nhd
        if (i == nhd) then
! time interval not found
          WRITE (lout,*) '--error--, in pbl_height_forcage -boxheight-'
          WRITE (lout,*) 'upper limit for time not found'
          stop
        endif
        if(boxt(i+1) >= timemod) then
          dhdt=(boxh(i+1)-boxh(i))/(boxt(i+1)-boxt(i))
          height=boxh(i) + dhdt*(timemod-boxt(i))
          exit
        endif
      enddo

!-------------------------------------------------------
      END SUBROUTINE pbl_height_forcage
!-------------------------------------------------------
