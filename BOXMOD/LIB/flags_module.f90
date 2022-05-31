!***************************************************************
!  HEADER INFORMATION: flags_module
!  declaration of flags used in box model
!***************************************************************
      MODULE flags_module

      IMPLICIT NONE

      INTEGER :: print_steadystate_fg
! flag for input/output format : 0=NetCDF; 1=binary, 2=both (NetCDF in).
      INTEGER :: iofmt_fg
! flag to use a previously output data file as input
      INTEGER :: prevflag
! flag 1 = developing/Lagrangian time. flag 0 = repeating/Eulerian
      INTEGER :: lagflag
! flag = 1: Henry''s Law deposition routine is called (0 = no deposition)
      INTEGER :: depos_fg
! flag = 1: we have dilution and mixing from outside box
      INTEGER :: mixing_fg
! flag = 1: we have emissions
      INTEGER :: emis_fg
! flag > 0: we emit NO based on soil temperature
! values of flag is multiplicative factor
      REAL :: noem_fg
! flag to consider the formation of soa
! if soa_fg eq 0 run without soa module
! if soa_fg eq 1 run with aerosol module without unifac
! if soa_fg eq 2 run with aerosol module with mass transfer
      INTEGER :: soa_fg
! flag to select the vapor pressure estimation method : 1 for JR-MY, 2 for nannoolal, 3 for SIMPOL.1
      INTEGER :: pvap_fg
! flag to consider dimerisation : 1 with dimer, 0 without
      INTEGER :: dimer_fg
! flag 1 = RO2+RO2 reactions allowed (0 = no RO2+RO2)
      INTEGER :: ro2_fg
! flag 1 = OFR simulation (REQUIRES MECH WITH OFR-SPECIFIC INORG RXNS)
      INTEGER :: OFR_fg
! flag 1 = output all reference j-values at each output time
      INTEGER :: jall_fg
! flag 1 = mtratdyn gas-> aero description (0; mtrat gas -> aero coeff)
      INTEGER :: dyn_fg
! flag 1 = wall losses considered (0; no wall losses)
      INTEGER :: wall_fg
! flag identifying specific chamber
!      1 = Matsunaga/Ziemann, 2 = Ziemann-CU, 3 = Jimenez-CU
      INTEGER :: icham
! flag identifying seed aerosol type (0 = organic, 1 = inorganic)
      INTEGER :: seedtyp_fg
! flag do we allow isopsoa formation? this requires constraining 
! ph, sulfate, nitrates, kappa and naer for inorganic aerosol
      INTEGER :: isopsoa_fg
! simple multiplicative factor for isopsoa uptake
      REAL :: isopsoa_fac      
! flag to include calculation of vbs (Shrivastava et al., 2019)
      INTEGER :: vbs_fg
! flag to include aging of vbs (Shrivastava et al., 2019)
      INTEGER :: vbs_aging_fg
! chamber module flags
      INTEGER :: lightsonoff_fg
      INTEGER :: inject_fg
! thermodynamic
      INTEGER :: iscape

! flag 1 = calculate loss/prod of selected tracers
! WARNING ! this will increase the running time by approx a factor 2
      INTEGER :: tracer_fg
! flag 1 = output instantaneous chemical reaction rates
      INTEGER :: reacrate_fg

! simple multiplicative factors for gas->wall rate
!      REAL :: winfac, weqfac

! FLAGS HARDWIRED IN BOXMOD_MAIN
! print photolysis rates [s-1] at each timestep
! (only bottom box)
      LOGICAL :: printphoto_fg

      END MODULE flags_module
