***********************************************************************
* This subroutine reads data for Pvap estimate from the NetCDF link file
*                                                                     *
* INPUT :                                                             *
*    -lout                                                            *
*    - chrsp                                                          *
*    - numsp                                                          *
*                                                                     *
* OUTPUT :                                                            *
*    - nsat : number of species in the datamyr.sat file               *
*    - namsat(i) : names of the species for which Pvap can be         *
*                  computed                                           *
*    - Tb(i)  : The boiling point of chem, computed by the Joback     *
*               group contribution method                             *
*    - HBN(i) : Hydrogen bond number                                  *
*    - tau(i) : effective number of torsional bonds                   *
*    - idsat(i) : index of species i in the full list od species      *
*                 in the scheme                                       *
***********************************************************************
      SUBROUTINE readpvap_nan_ncdf(ncid,chrsp,numsp,nsat,namsat,
     &                              Tb,dB,idsat,satid)

      USE akparameter_module
      USE netcdf

      IMPLICIT NONE

      
* INPUT 
      INTEGER  ncid
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      REAL              Tb(mxsat)
      REAL              dB(mxsat)
      INTEGER           nsat, idsat(mxsat),satid(maxsp)

* LOCAL
      CHARACTER(40)  line
      INTEGER        i, j, k, statcode, varid
      CHARACTER(6)   temp

* NETCDF 
      REAL,DIMENSION(mxsat,2) :: nandat

      Tb(:)=0.
      dB(:)=0.
      nandat(:,:)=0.
      idsat(:)=0

! retrieve "nsat"
      CALL eznc_get_0Dint(ncid,"nsat",nsat)
! retrieve Tb and dB
      CALL eznc_get_2Dreal(ncid,"nandat",mxsat,2,
     $                           nandat(1:nsat,1:2),
     $                                  1,nsat,1,2)

      Tb=nandat(:,1)
      dB=nandat(:,2)

! read namsat names and index to corresponding chrsp names
! if "sat" form of identifiers not found in input file, use "nan" form instead
! ----------------------------------------------------------
      statcode = NF90_INQ_VARID(ncid,"namsat",varid)
      IF (statcode.EQ.-49) THEN

      CALL eznc_get_1Dchar(ncid,"namnan",maxlsp,mxsat,
     $                           namsat,1,nsat)
      CALL eznc_get_1Dint(ncid,"idnan",mxsat,
     $                          idsat(1:nsat),1,nsat)
      CALL eznc_get_1Dint(ncid,"nanid",maxsp,
     $                          satid(1:numsp),1,numsp)
      ELSE
      CALL eznc_get_1Dchar(ncid,"namsat",maxlsp,mxsat,
     $                           namsat,1,nsat)
      CALL eznc_get_1Dint(ncid,"idsat",mxsat,
     $                          idsat(1:nsat),1,nsat)
      CALL eznc_get_1Dint(ncid,"satid",maxsp,
     $                          satid(1:numsp),1,numsp)
      ENDIF

* -----------------
* end of the routine
* -----------------
      END SUBROUTINE readpvap_nan_ncdf
