      SUBROUTINE init_defaults

      USE time_mgmt_module
      USE constraints_module,ONLY : constrained_thing
      USE solver_params_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE flags_module

      IMPLICIT NONE
      !include 'flags.h'

! initialize defaults for some values & flags
! (for production runs, flags are set in script/keyfile)
!====================================================
! MEANING OF FLAGS
!====================================================
! flag for i/o format (0:NetCDF; 1:binary;
!                      2:in=NetCDF, out=binary+NetCDF)
!      INTEGER :: iofmt_fg = 0
! flag to use a previous output data as input
!      INTEGER :: prevflag = 0
! flag 1 = developing/Lagrangian time. flag 0 = repeating/Eulerian
!      INTEGER :: lagflag = 0
! flag = 1: Henry's Law deposition routine is called (0 = no deposition)
!      INTEGER :: depos_fg = 0
! flag to consider the formation of soa
! if soa_fg eq 0 run without soa module
! if soa_fg eq 1 run with aerosol module without unifac
! if soa_fg eq 2 run with aerosol module with mass transfer
!      INTEGER :: soa_fg=2
! flag to select the vapour pressure estimation method : 1 for JR-MY, 2
! for nannoolal, 3 for SIMPOL.1
!      INTEGER :: pvap_fg = 2
! flag to consider dimerisation : 1 with dimer, 0 without
!      INTEGER :: dimer_fg = 0
! flag 1 = RO2+RO2 reactions allowed (0 = no RO2+RO2)
!      INTEGER :: ro2_fg = 1
! flag 1 = OFR simulation (REQUIRES MECH WITH OFR-SPECIFIC INORG RXNS)
!      INTEGER :: OFR_fg = 0
! flag 1 = mtratdyn gas-> aero description (0; mtrat gas -> aero coeff)
!      INTEGER :: dyn_fg=1      !1 mtratdyn gas-> aero description
!      LOGICAL :: wall_fg=0

! initialize subroutine selection flags
      depos_fg = 0
      dimer_fg = 0
      dyn_fg=1
      emis_fg = 0
      icham = 0
      iofmt_fg = 0
      isopsoa_fg = 0
      jall_fg = 1
      lagflag = 0
      mixing_fg = 0
      OFR_fg = 0
      prevflag = 0
      print_steadystate_fg = 1
      printphoto_fg = .true.
      pvap_fg = 2
      ro2_fg = 1
      seedtyp_fg = 0
      soa_fg = 2
      tracer_fg=0
      wall_fg=0
      mixing_fg = 0
      vbs_fg = 0
      vbs_aging_fg = 0
      lightsonoff_fg = 0
      inject_fg = 0
      reacrate_fg = 0

! initialize fixed/varying flags
      dilfix = 0
      noxfix = 0
      presfix = 0
      szafix = 0
      sumcfix = 0
      waterfix = 0

! initialize some other flags
      iscape = 0
      iseas = 0

! initialize variables
      boxt= 0.
      boxh= 0.
      caer= 0.
      cbg = 0.
      cnv = 2.4E9
      conc = 0.
      cro2 = 0.
      cseed = cnv
      eflux = 0.
      gamm = 1.  
      htop = 0.
      idtr = 0
      ibox = 0
      jfac = 1.0
      jftim = 0.
      jfval = 0.
      maer = 0.
      Mp = 426.  ! DOS: default
      nbox = 0
      nhd = 0
      nsd = 0
      njf = 0
      nmx = 0
      noem_fg=0.
      nrh = 0
      nseed = 0
      nskip = 1
      ntk = 0
      ntprint = 0
      ntstep = 1
      nws = 0
      psat = 0.
      psurf = 0.
      qfor = 0.
      rdep = 0.
      rdil = 0.
      rem = 0.
      rex = 0.
      vmix = 0.
      vs = 0.
      rha = 0.
      rhm = 0.
      rhtm = 0.
      rhtim = 0.
      rhval = 0.
      Rp = 0.
      Rpo = 1.25E-5  ! initial particle radius (cm) = 125 nm
      saero = 0.
      sumnox = 0.
      szaval = 0.
      tempm = 0.
      tempa = 0.
      temptm = 0.
      timemod = 0.
      tktim = 0.
      tkval = 0.
      trprod = 0.0
      trloss = 0.0
      tstart = 0.0
      tstop = -1.0
      Vd = 0.
      windm=0.
      winda=0.
      windtm=0.
      wstim=0.
      wsval=0.
      rtol=1e-2
      atol=1e2
      sla = 0.
      slo = 0.
      tz = 0.
      iy = 0
      im = 0
      id = 0

! variables that perhaps ought to be in log file (but aren't yet)
      dtmin = 0.1 ! (parameter for 2step)
      numit = 10  ! (parameter for 2step)
      winfac=1.0    ! simple multiplicative factor for gas->wall rate
      weqfac=1.0    ! simple multiplicative factor for gas->wall rate
      ndim = 0
      isop_fac = 1.0
      mterp_fac = 1.0
      isopsoa_fac = 1.0
      iskip=0
      itout=1 ! (corresponds to initialization time, t0)

      PRINT*,"initialise constraints array"
      cons_spec(:,:)%activefg = .false.
      cons_spec(:,:)%name = ""
      cons_spec(:,:)%unit = ""
      cons_spec(:,:)%index = 0
      cons_spec(:,:)%npoints = 0
      DO i = 1, maxconst
        DO j = 1, mbox
        cons_spec(j, i)%table = 0.0
        ENDDO
      ENDDO

      PRINT*,"initialise emissions array"
      emi_spec(:)%activefg = .false.
      emi_spec(:)%name = ""
      emi_spec(:)%unit = ""
      emi_spec(:)%index = 0
      emi_spec(:)%npoints = 0
      DO i = 1, maxem
        emi_spec(i)%table = 0.0
      ENDDO

      PRINT*, "initialise surface array"
      DO i = 1, msur
        surface_emi(i)%emission(:)%activefg = .false.
        surface_emi(i)%emission(:)%name = ""
        surface_emi(i)%emission(:)%unit = ""
        surface_emi(i)%emission(:)%index = 0
        surface_emi(i)%emission(:)%npoints = 0
        surface_emi(i)%nemis = 0
        DO j = 1, maxem
          surface_emi(i)%emission(j)%table = 0.0
        ENDDO
      ENDDO

      END SUBROUTINE init_defaults

