***********************************************************************
* This subroutine reads data for Pvap estimate from the NetCDF link file
*                                                                     *
* INPUT :                                                             *
*    -lout                                                            *
*    - chrsp                                                          *
*    - numsp                                                          *
*                                                                     *
* OUTPUT :                                                            *
*    - nsat : number of species in the {mech.psim} file               *
*    - namsat(i) : names of the species for which Pvap can be         *
*                  computed                                           *
*    - satid(i) : index of species i in the full list of species      *
*                 in the scheme                                       *
***********************************************************************
      SUBROUTINE readpvap_sim_ncdf(ncid,numsp,nsat,namsat,
     &                              bk,simpgroup,idsat,satid)

      USE akparameter_module
      USE netcdf

      IMPLICIT NONE

* INPUT 
      INTEGER  ncid
      INTEGER  numsp

* OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      REAL              simpgroup(mxsat,31)
      REAL              bk(31,4)
      INTEGER           nsat, idsat(mxsat),satid(maxsp)

* LOCAL
      CHARACTER(40)  line
      INTEGER        statcode, varid


      simpgroup(:,:)=0.
      bk(:,:)=0.
      idsat(:)=0
      satid(:)=0

! retrieve "nsat"
      CALL eznc_get_0Dint(ncid,"nsat",nsat)
! retrieve simpgroup and bk
      CALL eznc_get_2Dreal(ncid,"simpgroup",mxsat,31,
     $                           simpgroup(1:nsat,1:31),
     $                                     1,nsat,1,31)
      CALL eznc_get_2Dreal(ncid,"bk",31,4,
     $                           bk(1:31,1:4),
     $                              1,31,1,4)


! read namsat names and indexes to corresponding chrsp names
! if "sat" form of identifiers not found in input file, use "sim" form instead
! ----------------------------------------------------------
      statcode = NF90_INQ_VARID(ncid,"namsat",varid)
      IF (statcode.EQ.-49) THEN

      CALL eznc_get_1Dchar(ncid,"namsim",maxlsp,mxsat,
     $                           namsat,1,nsat)
      CALL eznc_get_1Dint(ncid,"idsim",mxsat,
     $                          idsat(1:nsat),1,nsat)
      CALL eznc_get_1Dint(ncid,"simid",maxsp,
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
      END SUBROUTINE readpvap_sim_ncdf
