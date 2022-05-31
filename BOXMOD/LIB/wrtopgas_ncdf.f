      SUBROUTINE wrtopgas_ncdf(ncid,ibox,itout,cbox)


!==================================================================
! PURPOSE: write runtime gas-phase output variables 
!          in NetCDF-format box model output file
! AUTHOR: Julia Lee-Taylor, NCAR, 123Jan 2018
!==================================================================
      !USE io_units_module,ONLY: lro2
      USE akparameter_module  
      USE flags_module,ONLY: reacrate_fg
      USE time_mgmt_module,ONLY: time
      USE forcing_params_module,ONLY: nbox,sumc,sumcfix,pres,temp,rh
      USE module_data_gecko_main,ONLY: qfor,idrestoi,numre,numsp,sza,
     &                                 cro2,nclro2,cmeo2


      IMPLICIT NONE

      INTEGER ncid,itout,ibox
      REAL,DIMENSION(maxsp)  :: cbox

! calculate instrantaneous rates separately for first and second
! reactants
      INTEGER ire,ileft,isp
      REAL,DIMENSION(maxre) :: reacrate 

!----------------------------------------------------------
      reacrate=0.0

! write gas phase values for current timestep and box to netCDF output file

! NB: time,temp, rh calculated as scalars but written into 2D space
!     cbot etc are calculated as 1D arrays but written into 2D space

      CALL eznc_put_0Dreal_into1D(ncid,"time",time,itout)

      CALL eznc_put_0Dreal_into2D(ncid,"cmeo2",cmeo2,ibox,itout)

      CALL eznc_put_1Dreal_into3D(ncid,"conc",
     &                                  cbox(1:numsp),
     &                                       1,numsp,ibox,itout)

      CALL eznc_put_1Dreal_into3D(ncid,"cro2",
     &                                  cro2(1:nclro2),
     &                                       1,nclro2,ibox,itout)

! DIAGNOSTICS
!      WRITE(61,*) itout,cbox(92551),cbox(92552)
!      WRITE(66,*) itout,cbox(87187),cbox(87188)

! see calculation of sumc in spforcage6.f90
! this calculation is here in case we want to explicitly output pressure
      !pres = sumc(ibox)*8.31446*temp/6.022E+22

      CALL eznc_put_0Dreal_into2D(ncid,"temp",temp,ibox,itout)
      CALL eznc_put_0Dreal_into2D(ncid,"rh"  ,rh  ,ibox,itout)
      CALL eznc_put_0Dreal_into2D(ncid,"pres",pres,ibox,itout)
      CALL eznc_put_0Dreal_into1D(ncid,"sza",sza,itout)
      IF (sumcfix.NE.1) THEN
        CALL eznc_put_1Dreal_into2D(ncid,"sumc",sumc(ibox),
     &                                          ibox,ibox,itout)
      ENDIF

! calculate and output instantaneous reaction rates 
      IF(reacrate_fg.GT.0)THEN
        DO ire=1,numre
          reacrate(ire)=qfor(ire)
          DO ileft=1,mxleft
            isp=(idrestoi(ire,ileft))
            IF(isp.GT.0)THEN
              reacrate(ire)=reacrate(ire)*cbox(isp)
              !PRINT*,ire,ileft,qfor(ire),cbox(isp),reacrate(ire)
            ENDIF
          ENDDO
        ENDDO
        CALL eznc_put_1Dreal_into3D(ncid,"reacrate",
     &                                    reacrate(1:numre),
     &                                             1,numre,ibox,itout)
      ENDIF

!==end write output variables

      END SUBROUTINE wrtopgas_ncdf
