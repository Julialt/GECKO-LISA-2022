!---------------------------------------------------!
! Define the list of DEFAULTS FOR user-settable key parameters 
! & flags for GECKO-A
! Directory GECKO-A/LIB: accessible to developers
!---------------------------------------------------!  
MODULE keyuser
  IMPLICIT NONE

! manage mechanism size
  INTEGER :: maxgen     ! maximum # of generations allowed

! select the SAR for vapor pressure cutoff: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
  INTEGER :: pvap_sar
  REAL    :: critvp     ! log(Pvap) below which chemistry is ignored

! PARTITIONING FLAGS
!--------------------

! flag to write for gas/particle partitioning "reactions"
  LOGICAL :: g2pfg
! flag to write for gas/wall partitioning "reactions" (chamber simulations)
  LOGICAL :: g2wfg

! CHEMISTRY FLAGS
!-----------------

! "PAM"/"OFR" simulations: H chemistry could be considered
  LOGICAL :: OFRfg
! only NO chemistry for RO2 considered
  LOGICAL :: highnox
! add the RO2+OH reaction in RO2 chemistry
  LOGICAL :: rx_ro2oh
! allow isomer substitution (O=no substitution, 1=allow)
  LOGICAL :: isomerfg
! flag to activate DHF formation
  LOGICAL :: dhffg
! consider cyclic hemi acetal (CHA) species as non volatile
  LOGICAL :: chafg

! ADDITIONAL OUTPUT/DOCUMENTATION FLAGS
!----------------------------------------

! write "operator info" in the output files
  LOGICAL :: wtopefg
! write info about SAR (groups, ...) in dedicated files
  LOGICAL :: sar_info
! write the references for the reactions in the mechanism
  LOGICAL  :: wrtref
! write additional info on screen during generation
  LOGICAL  :: screenfg
  
END MODULE keyuser
!------------------------------------------
SUBROUTINE define_defaults
  USE keyuser, ONLY:maxgen,pvap_sar,critvp,g2pfg,g2wfg,&
                    OFRfg,ofrfg,highnox,rx_ro2oh,isomerfg,dhffg,&
                    chafg,wtopefg,sar_info,wrtref,screenfg
  IMPLICIT NONE

! manage mechanism size
  maxgen=3  ! maximum # of generations allowed
  critvp=-13. ! log(Pvap) below which chemistry is ignored

! PARTITIONING FLAGS
!--------------------

! flag to write for gas/particle partitioning "reactions"
  g2pfg=.TRUE.
! flag to write for gas/wall partitioning "reactions" (chamber simulations)
  g2wfg=.FALSE.
! select the SAR for vapor pressure cutoff: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
  pvap_sar=2

! CHEMISTRY FLAGS
!-----------------

! "PAM"/OFR simulations: H chemistry could be considered
  OFRfg=.FALSE.
  ofrfg=.FALSE.
! only NO chemistry for RO2 considered
  highnox=.FALSE.
! add the RO2+OH reaction in RO2 chemistry
  rx_ro2oh=.FALSE.
! allow isomer substitution (O=no substitution, 1=allow)
  isomerfg=.TRUE.
! flag to activate DHF formation
  dhffg=.FALSE.
! consider cyclic hemi acetal (CHA) species as non volatile
  chafg=.FALSE.
! ADDITIONAL OUTPUT/DOCUMENTATION FLAGS
!----------------------------------------

! write "operator info" in the output files
  wtopefg=.FALSE.      
! write info about SAR (groups, ...) in dedicated files
  sar_info=.FALSE.       
! write the references for the reactions in the mechanism
  wrtref=.TRUE.    
! write additional info on screen during generation
  screenfg=.FALSE.  

END SUBROUTINE define_defaults
