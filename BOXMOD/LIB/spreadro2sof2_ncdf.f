*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* constituant le (les) compteurs RO2 - schema chimique generateur *
* routine a reecrire par la suite                                 *
*******************************************************************
      SUBROUTINE readro2sof2_ncdf(ncid,chrsp,numsp,
     1                            nclro2,numchemro2,idchemro2,cro2)

      USE akparameter_module
      IMPLICIT NONE

      
* INPUT 
      INTEGER  numsp
      INTEGER ncid
      CHARACTER(maxlsp) chrsp(maxsp)
* OUTPUT
      INTEGER  nclro2, numchemro2(maxro2), idchemro2(mxro2cl,maxro2)
      REAL     cro2(maxro2)
* LOCAL
      INTEGER k

* -----------------------
* initialize
      idchemro2(:,:)=0
      numchemro2(:)=0
      cro2(:)=0.

* open link file (now done in calling routine)
!      CALL open_ncdf_readonly("indat.nc",ncid)

* set the number of ro2 class (must be .le. maxro2)
!(could check this value by interrogating NetCDF file)
      nclro2=9
      IF (nclro2.gt.maxro2) THEN
        WRITE(6,*) '--error--, in spreadro2.f (read comment)'
        STOP
      ENDIF

! read numchemro2(maxro2)
      CALL eznc_get_1Dint(ncid,"numchemro2",maxro2,
     $                          numchemro2(1:maxro2),1,maxro2)

      
* get the id for the RO2 species and store into idchemro2 table
      DO k=1,nclro2
        CALL eznc_get_2Dint(ncid,"idchemro2",mxrpero,maxro2,
     $                            idchemro2(1:numchemro2(k),k),
     $                                      1,numchemro2(k),k,k)
      ENDDO

* -----------------------
* RETURN
* -----------------------
      
      END 
