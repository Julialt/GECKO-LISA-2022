      SUBROUTINE solve_box

!--------------------------------------------------
! Find environmental parameters, 
! => calculate rates, 
! => integrate chemistry
! => call output writing
! Called for each box, each timestep
!--------------------------------------------------

      USE akparameter_module
      USE flags_module,ONLY: depos_fg,OFR_fg,soa_fg,tracer_fg, &
                             lagflag,prevflag,emis_fg,printphoto_fg
      USE io_units_module,ONLY: lpbl,lppa,lpaa
      USE time_mgmt_module,ONLY: itout,idat,iskip,nskip,tstart, &
                                 time,timemod,tlocal,time_save,tout
      USE printphoto_module,ONLY: n_printphoto,idprintphoto,photorates
      USE NetCDF_vars_module,ONLY: ncid_out 
      USE solver_params_module,ONLY: numit,atol,rtol,dtmin,dtmax
      USE forcing_params_module,ONLY: height
      USE module_data_gecko_main,ONLY: numsp,numre,chrsp,ibox,conc, &
                                 sumc,arrhcf,cdimer,isocf,qfor, &
                                 cmeo2,cro2,idmeo2,ido2,id_m,idextra, &
                                 idno,idno2,idno3,idiso,idhv,idfo, &
                                 nself,idrestoi,idpdstoi,numstoi, &
                                 restoicf,pdstoicf, &
                                 idreacro2,idselfreac,idreacdimer, &
                                 inorg_aer,extracf,focf,lpmap, & 
                                 numiso,numhv,nummeo2,numo2, &
                                 numreacdimer,numreacro2, &
                                 ncldimer,nclro2,num_m,numfo,numextra, &
                                 hvfact,id,im,iy,sla,slo,tz,vd,sza, &
                                 coef1pho,coefmedpho,coeftoppho, &
                                 rat1pho,ratmedpho,rattoppho, &
                                 id1chromo,idmedchromo,idtopchromo, &
                                 nt1chromo,ntmedchromo, &
                                 nummedchromo,numtopchromo, &
                                 numtet,xang,szaval,szafix,saero, &
                                 rem,rdep,rex,rdil,noxfix,sumnox, &
                                 temp,water,wmol,maer,ctotaer, &
                                 idtr,ntr,trprod,trloss,cons_spec
      IMPLICIT NONE

      INTEGER :: lbox       ! current box output unit : lppa or lpaa
      REAL    :: maerbox    ! = maer(ibox)
      REAL    :: ctotaerbox ! = ctotaer(ibox)
      REAL,DIMENSION(maxsp) :: cbox ! = conc(1:numsp,ibox)

! mandatory interface for extract_reacrates
      INTERFACE
        SUBROUTINE extract_reacrates(qfor, idreacs, rate_tab )
          REAL, INTENT(IN)                 :: qfor(:)
          INTEGER, ALLOCATABLE,INTENT(IN)  :: idreacs(:)
          REAL, ALLOCATABLE,INTENT(INOUT)  :: rate_tab(:)
        END SUBROUTINE extract_reacrates
      END INTERFACE

!------------------------------------------------------
! set up variable arrays and io units to use for current box

      cbox = conc(:,ibox)
      maerbox = maer(ibox)
      ctotaerbox = ctotaer(ibox)

      SELECT CASE (ibox)
        CASE (1)
          lbox = lppa
        CASE (2)
          lbox = lpaa
        CASE DEFAULT
          lbox = lppa
      END SELECT

! pbl height

      IF(ibox.EQ.1)THEN
        CALL pbl_height_forcage
        !IF (time.EQ.0.OR.iskip.eq.nskip) THEN
        IF (iskip.EQ.0.OR.iskip.eq.nskip) THEN
          !WRITE(lpbl,*) time,",", height
          WRITE(lpbl,*) tout,",", height
        ENDIF
      ENDIF

!------------------------------------------------------
! forcage is called either with repeated day 1 (timemod, "Eulerian") 
! or actual time (time, "Lagrangian"). 
! This currently applies to EMISSIONS ONLY: everything else uses timemod

      IF (lagflag.EQ.1) THEN
        tlocal = time
      ELSE
        tlocal = timemod
      ENDIF

      CALL forcage6(tlocal,cbox,maerbox)

!----------------------------------------------------------------------
! Calculate some collected concentrations (local to current box)

      CALL find_concs(cbox)

!------------------------------------------------------
! Calculations depending on time-varying forcing and 
! either same for both boxes or specific to box 1

      IF(ibox.EQ.1)THEN

! Photolysis rates, adjusted by time-varying factor jfac
! Photolysis rates for OFR mode calculated INSIDE interp5

        CALL interp5( numhv,idhv,hvfact, &
              nt1chromo,id1chromo, &
              ntmedchromo,nummedchromo,idmedchromo, &
              numtopchromo,idtopchromo, &
              timemod, &
              sla,slo,tz,iy,im,id, &
              szafix,szaval, &
              numtet,xang, &
              rat1pho,coef1pho, &
              ratmedpho,coefmedpho, &
              rattoppho,coeftoppho, &
              arrhcf, sza, cbox)

!---------------------------------------------------
! surface characteristics, surface emission and surface deposition

        IF (depos_fg.GT.0.OR.emis_fg .gt. 0) CALL surface_rates

!------------------------------------------------------!
! preserve solver time from box 1 (just in case)
        time_save = time
      ELSE !(i.e. ibox.NE.1)
! reinstate solver time in other box(es) (just in case)
        time = time_save
      ENDIF! (ibox.EQ.1)

!---------------------------------------------------
! compute rate constants

      CALL akkrat6(chrsp,numsp,numre,num_m,numfo,numextra,numo2, &
                   nummeo2,numiso,idiso, &
                   nclro2,numreacro2,idrestoi,id_m,idfo,idextra, &
                   ncldimer,numreacdimer,idreacdimer,cdimer, &
                   ido2,idmeo2,idreacro2, &
                   arrhcf,focf,extracf,isocf,wmol, &
                   vd,cro2,cmeo2, &
                   temp,sumc,ibox,water,height,saero,qfor, &
                   inorg_aer)

!-----------------------------------------------------------------
! find and output photolysis rates [s-1] of all inorganics (box 1 only)

      IF (iskip.EQ.0.OR.iskip.EQ.nskip) THEN
        IF (ibox.EQ.1 .AND. printphoto_fg .AND. n_printphoto.GT.0) THEN
          CALL extract_reacrates(qfor, idprintphoto, photorates)
          CALL write_printphoto
        ENDIF !(ibox.EQ.1) THEN
      ENDIF !(correct iskip)

!----------------------------------------------------------------------
! Calculate saturation vapor pressures (temp-dependent, local to current box) and rate constant for mass transfer (if soa_fg = 2)

      IF (soa_fg.NE.0) THEN
        CALL calc_psat
        IF (soa_fg.EQ.2) CALL soa_dyn(ctotaerbox,maerbox)
      ENDIF

!-----------------------------------------------------------------
! rate of emission, deposition, dilution and exchange 

      CALL phys_rates

!----------------------------------------------------------------
! t0  , READ PREVIOUS CONCENTRATIONS IF SUPPLIED
! t1:n, SOLVE THE GAS-PHASE EQUATIONS 
!----------------------------------------------------------------

      IF (iskip.EQ.0)THEN
        IF(prevflag.EQ.1 .and. time .eq. tstart) THEN
          CALL get_prev(cbox)
          IF(ibox.EQ.1) itout = itout+ idat -1
        ENDIF
      ELSE

        CALL twostep5(maxsp, maxre, numre, mxleft, mxright, mself, &
                      numsp,numstoi,idrestoi,idpdstoi,nself,idselfreac,&
                      time,tout,dtmin,dtmax,qfor,restoicf,pdstoicf, &
                      numit, atol, rtol, &
                      rem, rdep, rex, rdil, &
                      mtr,ntr,idtr,trprod,trloss, &
                      pack(cons_spec(ibox,:)%index, &
                      mask = cons_spec(ibox,:)%activefg), &
                      count(cons_spec(ibox,:)%activefg), &
                      cbox, noxfix(ibox), sumnox, idno, idno2, idno3, &
                      lpmap)

      ENDIF !(time.EQ.tstart) THEN

!----------------------------------------------------------------
! compute soa for equilibrium representation
! soa_fg=1 ==> equilibrium at each time step (here).
! soa_fg=2 ==> update ctotaer and maer for output
!             (mass transfer rates were calculated before integration).
!----------------------------------------------------------------
      IF (soa_fg.EQ.1) CALL soa_equil(cbox,ctotaerbox,maerbox)
      IF (soa_fg.EQ.2) CALL soa_dyn_update(cbox,ctotaerbox,maerbox)

!----------------------------------------------------------------
! output everything

      IF (iskip.EQ.0.OR.iskip.eq.nskip) THEN
        CALL output_vals(cbox,ctotaerbox,maerbox)

! sync data to netcdf file
        CALL sync_ncfile(ncid_out)

!       USER may edit to output specific rates here

!       output tracer production/destruction rates
!       OPTIONAL AND TIME-CONSUMING!
        IF (ibox.EQ.1 .AND. tracer_fg.EQ.1) THEN
          WRITE(30,'(1x,ES10.3,4(1x,ES12.5))') &
          time, &
          trprod(idtr(1)),trloss(idtr(1)), &
          trprod(idtr(2)),trloss(idtr(2))
        ENDIF

      ENDIF

!----------------------------------------------------------------
! update variable arrays to carry forward for current box 

      conc(:,ibox) = cbox
      maer(ibox) = maerbox
      ctotaer(ibox) = ctotaerbox

!------------------------------------------------------
      END SUBROUTINE solve_box
!------------------------------------------------------
