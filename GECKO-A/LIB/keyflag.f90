MODULE keyflag
  IMPLICIT NONE

  REAL,PARAMETER  :: critvp=-13. ! log(Pvap) below which chemistry is ignored
  INTEGER,PARAMETER :: maxgen=4  ! maximum # of generations allowed
  REAL,PARAMETER    :: dT=298.   ! default temperature
  REAL,PARAMETER    :: dM=2.5E19 ! default 3rd body M (molec cm-3)
  
! CHEMISTRY FLAG - SAR SELECTOR
! -----------------------------

! ALKOXY
! R(O.) selector for isomerisation reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
  INTEGER, PARAMETER :: kisomfg=2      
! R(O.) selector for decomposition reaction rate (1=Atkinson 2007; 2=Vereecken 2009)
  INTEGER, PARAMETER :: kdissfg=2      

! OH RATE CONSTANT
! OH addition to alkene (1=Peeters 1997; 2=Ziemann 2009; 3=Jenkin 2016)
  INTEGER, PARAMETER :: kohaddfg=3


! -- FLAGS 

! flag to write for gas/particle partitioning "reactions"
  LOGICAL, PARAMETER :: g2pfg=.TRUE.

! flag to write for gas/wall partitioning "reactions" (chamber simulations)
  LOGICAL, PARAMETER :: g2wfg=.TRUE.

! allow isomer substitution (O=no substitution, 1=allow)
  INTEGER, PARAMETER :: isomerfg=1

! select the SAR for vapor pressure: 1=JR-MY; 2=Nannoolal; 3=SIMPOL-1
  INTEGER, PARAMETER :: pvapfg=2

! write "operator info" in the output files
  INTEGER, PARAMETER :: wtopeflag=0

! write info about SAR (groups, ...) in dedicated files  
  INTEGER, PARAMETER :: sar_info=0

! allow to switch enols to ketones (1: allowed, 0: not allowed)
  INTEGER, PARAMETER :: enolflg=1

! "PAM" simulations: H chemistry could be considered
  LOGICAL, PARAMETER :: lopam=.FALSE.  

! treat RO2+RO2 with the 9 classes (true) of ro2 or only CH3O2 (false)
  LOGICAL, PARAMETER :: multiclass=.TRUE.

! add the RO2+OH reaction in RO2 chemistry  
  LOGICAL, PARAMETER :: rx_ro2oh=.FALSE.

! only NO chemistry for RO2 considered 
  LOGICAL, PARAMETER :: highnox=.FALSE.

! write kOH, kNO3, kO3 for SAR assessment (rate in database not longer used)
! NOTE: if true, run the reactions for the input species only
  LOGICAL, PARAMETER :: losar=.FALSE.  

! flag to activate DHF formation 
  LOGICAL, PARAMETER :: dhffg=.FALSE. 

! consider cyclic hemi acetal (CHA) species as non volatile
  LOGICAL, PARAMETER :: locha=.FALSE.

! OUTPUT FLAGS
! -------------
  LOGICAL,PARAMETER :: screenfg=.FALSE.  ! write additional info on screen during generation
  LOGICAL,PARAMETER :: wrtdhf=.TRUE.     ! write output file with formation enthalpies information
  LOGICAL,PARAMETER :: wrtpvap=.TRUE.    ! write output file with the vapor pressures
  LOGICAL,PARAMETER :: wrthenry=.TRUE.   ! write output file with the Henry's law coefs
  LOGICAL,PARAMETER :: wrtref=.TRUE.     ! write the references for the reactions in the mechanism
  LOGICAL,PARAMETER :: wrttg=.TRUE.      ! write output file with Tg data


END MODULE keyflag
