      SUBROUTINE defopvals_ncdf(ncid)


!==================================================================
! PURPOSE: define space for runtime output variables 
!          in NetCDF-format box model output file
! AUTHOR: Julia Lee-Taylor, NCAR, Jan 2018
!==================================================================

      USE netcdf
      USE akparameter_module
      USE forcing_params_module,ONLY: nbox
      USE module_data_gecko_main,ONLY: numaou
      USE flags_module
      IMPLICIT NONE

      INTEGER ncid

!----------------------------------------------------------
!==define runtime outputs in o/p file
      CALL eznc_def_1Dreal(ncid,"time","ntout")
      CALL eznc_def_localatt(ncid,"time","title",
     &              "time since start of simulation series")
      CALL eznc_def_localatt(ncid,"time","units","seconds")

      CALL eznc_def_2Dreal(ncid,"temp","nbox","ntout")
      CALL eznc_def_localatt(ncid,"temp","title",
     &                            "temperature")
      CALL eznc_def_localatt(ncid,"temp","units","Kelvin")

      CALL eznc_def_2Dreal(ncid,"rh","nbox","ntout")
      CALL eznc_def_localatt(ncid,"rh","title",
     &                            "relative humidity")
      CALL eznc_def_localatt(ncid,"rh","units","percent")

      call eznc_def_2Dreal(ncid,"pres","nbox","ntout")
      call eznc_def_localatt(ncid,"pres","title","pressure")
      call eznc_def_localatt(ncid,"pres","units","bar")

      call eznc_def_1Dreal(ncid,"sza","ntout")
      call eznc_def_localatt(ncid,"sza","title",
     &                            "runtime solar zenith angle")
      call eznc_def_localatt(ncid,"sza","units","degrees from zenith")

      IF(jall_fg.EQ.1)THEN
      CALL eznc_def_2Dreal(ncid,"jchrom","nchrom","ntout")
      CALL eznc_def_localatt(ncid,"jchrom","title",
     &              "jvalues by photolysis model species")
      CALL eznc_def_localatt(ncid,"jchrom","units","s-1")
      ENDIF

!---species concentrations---
      CALL eznc_def_3Dreal(ncid,"conc","maxsp","nbox","ntout")
      IF(soa_fg.EQ.2)THEN
      CALL eznc_def_localatt(ncid,"conc","title",
     &                            "concentrations")
      ELSE
      CALL eznc_def_localatt(ncid,"conc","title",
     &                            "gas phase concentrations")
      ENDIF
      CALL eznc_def_localatt(ncid,"conc","congruence","chrsp")
      CALL eznc_def_localatt(ncid,"conc","units","molec cm-3")

!---instantaneous reaction rates---
      IF(reacrate_fg.GT.0)THEN
      CALL eznc_def_3Dreal(ncid,"reacrate","maxre","nbox","ntout")
      CALL eznc_def_localatt(ncid,"reacrate","title",
     &                  "instantaneous reaction rates: all reactions")
      CALL eznc_def_localatt(ncid,"reacrate","units","molec cm-3 s-1")
      ENDIF

!---CH3RO2 concentrations---
      CALL eznc_def_2Dreal(ncid,"cmeo2","nbox","ntout")
      CALL eznc_def_localatt(ncid,"cmeo2","title",
     &                "gas phase concentrations of CH3O2")
      CALL eznc_def_localatt(ncid,"cmeo2","units","molec cm-3")

!---lumped RO2 concentrations---
      CALL eznc_def_3Dreal(ncid,"cro2","maxro2","nbox","ntout")
      CALL eznc_def_localatt(ncid,"cro2","title",
     &                "gas phase concentrations of lumped RO2")
      CALL eznc_def_localatt(ncid,"cro2","congruence","ro2sp")
      CALL eznc_def_localatt(ncid,"cro2","units","molec cm-3")

!---particle phase concentrations---
      IF (soa_fg.GT.0) THEN
        CALL eznc_def_2Dreal(ncid,"ctotaer","nbox","ntout")
        CALL eznc_def_localatt(ncid,"ctotaer","title",
     &         "total particle phase concentration including seed")
        CALL eznc_def_localatt(ncid,"ctotaer","units","molec.cm-3")

!---particle phase mass--
        CALL eznc_def_2Dreal(ncid,"maer","nbox","ntout")
        CALL eznc_def_localatt(ncid,"maer","title",
     &                         "particle phase mass (non-seed)")
        CALL eznc_def_localatt(ncid,"maer","units","ug.m-3")

        IF (soa_fg.EQ.1) THEN
          CALL eznc_def_3Dreal(ncid,"psat","mxsat","nbox","ntout")
          CALL eznc_def_localatt(ncid,"psat","title",
     &                               "saturation vapor pressures")
          CALL eznc_def_localatt(ncid,"psat","congruence","idsat")
          CALL eznc_def_localatt(ncid,"psat","units","atm")

          CALL eznc_def_2Dreal(ncid,"caer","mxsat","ntout")
          CALL eznc_def_localatt(ncid,"caer","title",
     &                            "particle phase concentrations")
          CALL eznc_def_localatt(ncid,"caer","congruence","idsat")
          CALL eznc_def_localatt(ncid,"caer","units","molec cm-3")
        ELSEIF (soa_fg.EQ.2) THEN
          CALL eznc_def_3Dreal(ncid,"psat","maxsp","nbox","ntout")
          CALL eznc_def_localatt(ncid,"psat","title",
     &                               "saturation vapor pressures")
          CALL eznc_def_localatt(ncid,"psat","units","atm")

          CALL eznc_def_2Dreal(ncid,"Rp","nbox","ntout")
          CALL eznc_def_localatt(ncid,"Rp","title",
     &           "time-varying particle radius")
          CALL eznc_def_localatt(ncid,"Rp","units","cm")

          CALL eznc_def_2Dreal(ncid,"mwaer","nbox","ntout")
          CALL eznc_def_localatt(ncid,"mwaer","title",
     &           "particle phase mean molecular mass (including seed)")
          CALL eznc_def_localatt(ncid,"mwaer","units","a.m.u.")

!          CALL eznc_def_1Dreal(ncid,"Cstar298","maxsp")
!          CALL eznc_def_localatt(ncid,"Cstar298","title",
!     &                                "Cstar at 298K")
!          CALL eznc_def_localatt(ncid,"Cstar298","units","ug.m-3")

        ENDIF ! (soa_fg.EQ.1)
      ENDIF ! (soa_fg.GT.0)

!==end define output variables

      END SUBROUTINE defopvals_ncdf
