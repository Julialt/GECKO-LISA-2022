*******************************************************************
* Cette routine a pour objet de lire et de definir les especes    *
* qui se depose, et de lire les parametres qui y sont associes    *
* INPUT :                                                         *
*   lout   : numero de fichier pour ecriture des erreurs          *
*   iscape : drapeau pour calcul des equilibres thermos           *
* OUTPUT :                                                        *
*   ndep  : nombre d'espece qui se depose                      *
*   depnam(i) : nom de la ieme espece qui se depose            *
*   depdatspe(i,j) : parametre j de la ieme espece                *
*                    j=1 => D_H2O/D_x                             *
*                    j=2 => constante de Henry                    *
*                    j=3 => facteur de reactivite                 *
*   iddep(i)    : numero d'identification de l'espece dans le  *
*                    mecanisme chimique                           *
*
* NOTE: binary version reads either vfile.dep or vfileaero.dep.
*       This version (so far) only reads one set of dep data
*******************************************************************
      SUBROUTINE readdep3_ncdf(ncid,lout,chrsp,numsp,iscape,
     1                    ndep,depnam,depdat,iddep)

      USE akparameter_module
      IMPLICIT NONE

* INPUT
      INTEGER   ncid
      INTEGER   lout,iscape, numsp
      CHARACTER(maxlsp) chrsp(maxsp)

* OUTPUT
      CHARACTER(maxlsp) depnam(mxdep)
      INTEGER  ndep, iddep(mxdep)
      REAL     depdat(mxdep,3)

* LOCAL
      CHARACTER(80) line
      INTEGER i,j,k,isp

* ----------------------------------------------------------
      depdat(:,:)=0.
      iddep(:)=0

*******************************************************************
* recommence la procedure precedente pour les especes gazeuses    *
* en equilibre thermo                                             *
*******************************************************************
      !IF (iscape.EQ.1) THEN

! retrieve "ndep"
      CALL eznc_get_0Dint(ncid,"ndep",ndep)
! retrieve deposition data
      CALL eznc_get_2Dreal(ncid,"depdat",mxdep,3,
     $                           depdat(1:ndep,1:3),
     $                                  1,ndep,1,3)

! read depnam names and index to corresponding chrsp names
      CALL eznc_get_1Dchar(ncid,"namdep",maxlsp,mxdep,
     $                           depnam,1,ndep)

      CALL eznc_get_1Dint(ncid,"iddep",mxdep,
     $                          iddep(1:ndep),1,ndep)

* check if the species are defined in the mechanism
* stop if the species is unknown
!! PROBABLY UNNECESSARY FOR NETCDF FILE - WAS CHECKED AT TIME OF FILE CREATION
!! AND WE JUST READ IDDEP FROM THE FILE
      DO i=1,ndep
        CALL akspnum(depnam(i),chrsp,numsp,isp)
        iddep(i)=isp
        IF (isp.eq.0) THEN
          WRITE(lout,*)' --error--  while reading depsp.dat'
          WRITE(lout,*)' species unkown=',depnam(i)
          STOP
        ENDIF
      ENDDO

      !ENDIF ! IF(iscape.EQ.1)

* end of the routine
      END
