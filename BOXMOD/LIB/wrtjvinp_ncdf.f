      SUBROUTINE wrtjvinp_ncdf(ncid,njtab,nsza,
     &                         mxjtab,mxsza,llin,
     &                         idjtab,jsza,jvref,jnam,jreac)


!==================================================================
! PURPOSE: write input j-value tables 
!          to NetCDF-format box model output file
! AUTHOR: Julia Lee-Taylor, NCAR, 16 Feb 2018
!==================================================================

      USE akparameter_module
      IMPLICIT NONE

      INTEGER ncid,i
* photolysis input data, for writing
      INTEGER  mxjtab,mxsza,llin
      !INTEGER,PARAMETER:: mxjtab=150,mxsza=15,llin=76
      INTEGER  njtab,nsza,idjtab(mxjtab)
      REAL     jsza(mxsza),jvref(mxsza,mxjtab)
      CHARACTER(maxcoe) jnam(mxjtab)
      CHARACTER(llin) jreac(mxjtab)

!----------------------------------------------------------
* write input jvalue tables to netCDF output file

      CALL eznc_put_1Dchar(ncid,"jnam",jnam,maxcoe,1,njtab)

      CALL eznc_put_1Dchar(ncid,"jreac",jreac,llin,1,njtab)

      CALL eznc_put_1Dint(ncid,"idjtab",idjtab(1:njtab),1,njtab)

      CALL eznc_put_1Dreal(ncid,"jsza",jsza(1:nsza),1,nsza)

      CALL eznc_put_2Dreal(ncid,"jvref",jvref(1:nsza,1:njtab),
     &                                        1,nsza,1,njtab)

!==end write jvalue tables

      END SUBROUTINE wrtjvinp_ncdf
