      SUBROUTINE surface_rates
!!-----------------------------------------------------------------
! rates of :
! surface emissions
! surface deposition
!!-----------------------------------------------------------------

      USE flags_module,ONLY: emis_fg,noem_fg,depos_fg,isopsoa_fg
      USE io_units_module,ONLY: lout
      USE time_mgmt_module,ONLY: timemod,time,tlocal
      USE akparameter_module
      USE forcing_params_module
      USE module_data_gecko_main

      IMPLICIT NONE

!------------------------------------------------------
      !PRINT*,"starting surface rates"
      
      IF(iseas.eq.0)THEN
         WRITE (lout,*) '--error--'
         WRITE (lout,*) 'SEAS and SURF are required to use DEPO or SEMI'
         STOP
      ENDIF

! Compute surface percent at given time, needed for emissions and
! depositions

        CALL timesurf3(lout,tlocal,mhd,msur,nsd,psurf,surft,xsurf)

! compute surface emissions

        eflux = 0.

        IF (emis_fg .GT. 0) THEN

          CALL surface_emissions(lout,timemod,xsurf,surface_emi,eflux)

          IF (idisop > 0 .AND. isop_fac .GE. 0.0) THEN
            eflux(idisop) = eflux(idisop)*isop_fac
          ENDIF

          IF (mterp_fac .NE. 1.) THEN
            DO i = 1, maxsp
              IF ( chrsp(i) .eq. "GCARENE" .or. &
                   chrsp(i) .eq. "GAPINEN" .or. &
                   chrsp(i) .eq. "GBMYRCN" .or. &
                   chrsp(i) .eq. "GBPINEN" .or. &
                   chrsp(i) .eq. "GTHUJEN" .or. &
                   chrsp(i) .eq. "GCAMPHN" .or. &
                   chrsp(i) .eq. "GLIMONN" .or. &
                   chrsp(i) .eq. "GTPNENE" .or. &
                   chrsp(i) .eq. "GSABINN" .or. &
                   chrsp(i) .eq. "GTPNOLN" .or. &
                   chrsp(i) .eq. "GBOCIMN") then !identify monoterpenes
                eflux(i) = eflux(i)*mterp_fac
              ENDIF
            ENDDO
          ENDIF
!PRINT*,"!     NO emission from tropical rain forest"
! based on old parameterisation by
!  Yienger, J. J., & Levy, H. (1995).
! Empirical model of global soil-biogenic NOχ emissions.
!  Journal of Geophysical Research: Atmospheres,
! 100(D6), 11447–11464. https://doi.org/10.1029/95JD00370
! old but easy to implement...
! according to them, tropical rain forests are unique because researchers
! report no correlation between temperature and emissions
! even at temperatures well below 30°C
! Therefore we set the flux constant with respect to soil moisture condition
! 2.6 ngN m-2 s-1 for wet soil

          IF (noem_fg .GT.0. .AND. ibox .EQ. 1) THEN
            eflux(idno) = eflux(idno) + noem_fg*2.6e-9*6.02e23*1e-4/14  ! convert to molec cm2 s-1
          ENDIF
        ENDIF

! compute deposition for gas and for particles
        IF(depos_fg.GT.0)THEN

          CALL deposition3(timemod,temp,xsurf,Rik,Rlu0k, &
                           Rack,RgsSk,RgsOk,RclSk,RclOk,z0k, &
                           ndepspe,depnamspe,depdatspe,iddepspe, &
                           windm,winda,windtm,nws,wstim,wsval, &
                           iseas,sla,slo,tz,iy,im,id,Vd, ResAerdep)
                           
          if (isopsoa_fg .gt. 1) then
            CALL update_deposition(inorg_aer(ibox), &
                                 temp, ResAerdep, pres)
          endif

        ENDIF

      !PRINT*,"after surface rates"
! -------------------------------------
      END SUBROUTINE surface_rates
! -------------------------------------

