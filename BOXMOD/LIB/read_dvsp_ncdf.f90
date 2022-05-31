      SUBROUTINE read_dvsp_ncdf(ncid)
! read diffusion volumes and
! redistribute from mxsat into maxsp space 

      USE akparameter_module,ONLY: mxsat,maxsp
      USE module_data_gecko_main,ONLY: nsat,numsp,dvsp,difvol,difid

      !$ use OMP_LIB
      IMPLICIT NONE

      INTEGER :: i,ncid
! -------------------------------------
      PRINT*,"read_dvsp",numsp
! retrieve diffusion volume (nsat) and species ids (numsp)
      CALL eznc_get_1Dreal(ncid,"difvol",mxsat, &
                                 difvol(1:nsat), &
                                        1,nsat)
      CALL eznc_get_1Dint(ncid,"difid",maxsp, &
                                difid(1:numsp), &
                                      1,numsp)

! redistribute it
!$OMP PARALLEL DO private(i)
       DO i=1,numsp
         IF(difid(i).NE.0)THEN
           dvsp(i)=difvol(difid(i))
         ENDIF
       ENDDO
!$OMP END PARALLEL DO

! -------------------------------------
      END SUBROUTINE read_dvsp_ncdf
! -------------------------------------

