      SUBROUTINE phys_rates
!!-----------------------------------------------------------------
! rates of :
! dilution and exchange between 2 boxes (box-specific)
! emission & deposition (box-specific) 
! dilution is input as whole-box RATE, units s-1
! mixing & subsidence are input as cross-boundary VELOCITY, units cm.s-1
! NB: this subroutine can currently ONLY handle 2 boxes.
! Modification required for > 2 boxes.
!!-----------------------------------------------------------------

      USE flags_module,ONLY: emis_fg,depos_fg,inject_fg
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE module_chamber_tools, only: run_chamber_injections

      IMPLICIT NONE
!------------------------------------------------------
      !PRINT*,"starting phys_rates"

      rem  = 0.
      rdep = 0.
      rdil = 0.
      rex  = 0.

!---------------
      SELECT CASE (ibox)
!---------------
        CASE(1)
! box 1 dilution (loss) & exchange (addition)
          rdil = 0.
          IF (dilfix .eq. 1) THEN
            IF(dilconst.GT.0.) rdil = dilconst
          ELSE 
            IF (dhdt > 0.) rdil =  dhdt/height
          ENDIF

          IF(nbox.eq.2)THEN
            rex = conc(:,2) * (rdil + vmix(1)/height)
          ELSE
            rex = cbg(:) * (rdil + vmix(1)/height)
          ENDIF

          rdil = rdil + vmix(1)/height

! box 1 emission
          IF (emis_fg .GT. 0) rem = eflux/height
          IF (inject_fg .GT. 0) call run_chamber_injections()

! box 1 deposition for gas, aerosol species
        IF (depos_fg .gt. 0) THEN
          DO i=1,ndepspe
            rdep(iddepspe(i)) = vd(i)/height
          ENDDO
          DO i=1, nsat
            rdep(idasat(i)) = inorg_aer(ibox)%vdepaer / height
          ENDDO
        ENDIF

! save conc (box 1, start of timestep) for later if needed by second box
        IF (nbox.GT.1) cbot_sav(:) = conc(:,1)

!---------------
        CASE(2)
! box 2 dilution (loss) & exchange (addition)
          rdil = 0.
          IF (dilfix .eq. 1) THEN
            IF(dilconst.LT.0.) rdil = -dilconst
          ELSE
            IF (dhdt < 0.) rdil = -dhdt/(htop-height)
          ENDIF

! box 2 needs conc(:,1) from the beginning of the timestep!
!         hence use cbot_sav for dilution
          rex  = cbot_sav * rdil    &
               + cbot_sav * vmix(1)/(htop - height) &
               + cbg * vmix(2)/(htop - height)

          rdil = rdil + (vmix(1) + vmix(2))/(htop - height)

! box 2 dilution, add the loss due to atmospheric subsidence
! i.e  gases are "pushed away" by subsiding air from the free trop.
          rdil = rdil + vs / (htop - height)

! box 2 emission is the source coming from atmospheric subsidence
! i.e  gases are slowly subsiding from the free trop.
          rem = cbg * vs / (htop - height)

!---------------
      END SELECT
!---------------

! -------------------------------------
      END SUBROUTINE phys_rates
! -------------------------------------

