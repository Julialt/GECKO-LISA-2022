      SUBROUTINE wrtopaer_ncdf(ncid,ibox,itout,nsat,
     &                      maerbox,psat,caer,ctotaerbox)

!==================================================================
! PURPOSE: write runtime gasp-phase output variables 
!          in NetCDF-format box model output file
! AUTHOR: Julia Lee-Taylor, NCAR, 123Jan 2018
!==================================================================

      USE flags_module,ONLY: soa_fg
      USE akparameter_module
      USE forcing_params_module,ONLY: Mp,Rp,cnv
      USE module_data_gecko_main,ONLY: numsp,small
      USE fundamental_consts_module,ONLY: multimass

      IMPLICIT NONE

      INTEGER ncid,itout
      INTEGER ibox,nsat
      REAL    maerbox,ctotaerbox
      REAL,DIMENSION(mxsat)::  caer
      REAL,DIMENSION(maxsp)::  psat
!      REAL    small

      !small =  TINY(1.0)
!----------------------------------------------------------
* write particle phase values for current timestep to netCDF output file
      CALL eznc_put_0Dreal_into2D(ncid,"ctotaer",ctotaerbox,ibox,itout)
      CALL eznc_put_0Dreal_into2D(ncid,"maer",   maerbox   ,ibox,itout)

      IF (soa_fg.EQ.1) THEN
        CALL eznc_put_1Dreal_into3D(ncid,"psat",psat(1:nsat),
     &                                               1,nsat,
     &                                               ibox,itout)
        CALL eznc_put_1Dreal_into3D(ncid,"caer",caer(1:nsat),
     &                                               1,nsat,
     &                                               ibox,itout)
      ELSEIF (soa_fg.EQ.2) THEN
        CALL eznc_put_1Dreal_into3D(ncid,"psat",psat(1:numsp),
     &                                               1,numsp,
     &                                               ibox,itout)
        CALL eznc_put_0Dreal_into2D(ncid,"Rp",Rp,ibox,itout)
        IF (ctotaerbox.GT.small) THEN
        CALL eznc_put_0Dreal_into2D(ncid,"mwaer",
     &               (maerbox/multimass+cnv*Mp)/ctotaerbox,ibox,itout)
        ENDIF
      ENDIF ! (soa_fg.EQ.1)

!==end write output variables

      END SUBROUTINE wrtopaer_ncdf
