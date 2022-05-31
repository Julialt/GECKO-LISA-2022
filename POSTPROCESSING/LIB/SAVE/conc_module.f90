MODULE CONC
  USE PARAMETERS
  IMPLICIT NONE
!
  integer                  :: ndat, numsp
  character*(maxlsp), dimension(:), allocatable     ::  chrsp
  real(kind=8),  dimension(:,:), allocatable        ::  concentrations
  real(kind=8),  dimension(:), allocatable          ::  time, temperature, pressure, rh, sza !sec, K, hPa, %
  real(kind=8),  dimension(:), allocatable          ::  seed ! ug m-3
  real(kind=8),  dimension(:,:), allocatable        ::  precursor_lifetime
  integer                                           ::  efold_indices(5)
  integer                                           ::  precursor_indx(mxprecu), nprecu
END MODULE CONC
