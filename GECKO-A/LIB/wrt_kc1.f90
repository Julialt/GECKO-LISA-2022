! BA: Ugly subroutine ... waiting for a nicer one to write rate constant
! for C1 species with OH. Rate consistent with singlec_mech.dat file. 
SUBROUTINE wrt_kc1()
  USE keyparameter, ONLY: kohu
  IMPLICIT NONE
  
  REAL :: koh, tref
  
  tref=298.
!GCH4 + GHO => GCH3O2                          : 2.15E-12  0.   1735 ;
  koh=2.15E-12*exp(-1735./tref)
  WRITE(kohu,'(a7,1x,1pe10.2)') 'CH4    ', koh

!GCH3OH + GHO => 0.85 GCH2O + 0.85 GHO2 + 0.15 GCH3O : 3.10E-12  0. 360 ;
  koh=3.10E-12*exp(-360./tref)
  WRITE(kohu,'(a7,1x,1pe10.2)') 'CH3OH  ', koh

!GCH3OOH + GHO => GHO + GCH2O                  : 1.00E-12  0.   -190 ;
!GCH3OOH + GHO => GCH3O2                       : 1.90E-12  0.   -190 ;
  koh=2.90E-12*exp(190./tref)
  WRITE(kohu,'(a7,1x,1pe10.2)') 'CH3OOH ', koh

!GN01003 + GHO => GH2O + GCH2O + GNO           : 0.301E-12  0.    0. ;
  koh=3.01E-13
  WRITE(kohu,'(a8,1x,1pe10.2)') '!N01003 ', koh   ! commented - not produced or negligible

!GN01001 + GHO => GCH2O + GNO2                 : 4.00E-13  0.    845 ;
  koh=4.00E-13*exp(-845./tref)
  WRITE(kohu,'(a8,1x,1pe10.2)') '!N01001 ', koh   ! commented - not produced or negligible

!GCH2O + GHO => GCO + GHO2                     : 8.60E-12  0.    -20 ;
  koh=8.60E-12*exp(20./tref)
  WRITE(kohu,'(a7,1x,1pe10.2)') 'CH2O   ', koh

!GNO1001 + GHO => GHCOOH + GNO2                : 3.10E-12  0.      0 ;
  koh=3.10E-12
  WRITE(kohu,'(a8,1x,1pe10.2)') '!NO1001  ', koh  ! commented - not produced or negligible

!GHCOOH + GHO => GHO2 + GCO2                   : 4.50E-13  0.      0 ;
  koh=4.50E-13
  WRITE(kohu,'(a7,1x,1pe10.2)') 'HCOOH  ', koh

!GHCOO2H + GHO => GHCOO2                       : 1.90E-12  0.   -190 ;
  koh=1.90E-12*exp(190./tref)
  WRITE(kohu,'(a7,1x,1pe10.2)') 'HCOO2H  ', koh

END SUBROUTINE wrt_kc1
