MODULE RATE
  USE PARAMETERS
  IMPLICIT NONE
!
  integer                  :: maxre,mxleft,mxright
  integer                  :: numre
  integer, dimension(:,:), allocatable        ::  idrestoi, idpdstoi, numstoi
  real, dimension(:,:), allocatable           ::  restoicf, pdstoicf
  real(kind=8), dimension(:,:,:), allocatable ::  reacrate

END MODULE RATE
