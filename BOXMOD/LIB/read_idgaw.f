!**********************************************************************
! This subroutine read the file in which data for Pvap estimate       *
! are stored.                                                         *
!                                                                     *
! INPUT :                                                             *
!    -lout                                                            *
!    -chrsp                                                          *
!    -numsp                                                          *
!                                                                     *
! OUTPUT :                                                            *
!    - idgsat(i) : index of species i in the full list of species     *
!                 in the scheme - gas phase                           *
!    - idasat(i) : index of species i in the full list of species     *
!                 in the scheme - particle phase                      *
!    - idwsat(i) : index of species i in the full list of species     *
!                 in the scheme - wall phase                          *
!**********************************************************************
      SUBROUTINE read_idgaw(chrsp,numsp,ndim,nsat,namsat,
     &                      idgsat,idasat,idwsat,iddsat,satid)

      !$ use OMP_LIB

      USE flags_module,ONLY: wall_fg, dimer_fg
      USE akparameter_module
      USE io_units_module,ONLY: lread

      IMPLICIT NONE

! INPUT
      INTEGER  numsp
      CHARACTER(maxlsp) chrsp(maxsp)

! OUTPUT
      CHARACTER(maxlsp) namsat(mxsat)
      INTEGER           nsat,ndim,idgsat(mxsat),idasat(mxsat)
      INTEGER           idwsat(mxsat),iddsat(mxsat),satid(maxsp)

! LOCAL
!      CHARACTER(40)     line
!      CHARACTER(6)      temp
!      REAL               dum1, dum2
!      INTEGER            keof,kspe
      INTEGER            i, j, k
      CHARACTER(maxlsp) aspe,gspe,wspe
!      LOGICAL            found

! ----------------------------------------------
! Binary: use pre-read nsat and namsat from ascii input file
! Netcdf: use pre-read nsat and namsat from netcdf link file 
! ----------------------------------------------------------
! set index for corresponding namsat names with chrsp names
! ----------------------------------------------------------

! -----------------
! gas phase species
! -----------------
!$OMP PARALLEL DO private(i,gspe,j,k)
      satgloop: DO i=1,nsat
        gspe=namsat(i)
        call akspnum(gspe,chrsp,numsp,j)
! do loop exit should occur before if species found -
        if (j == 0) then
          WRITE(6,*) '--error--, the following species in nannoolal '
          WRITE(6,*) '           not in list of species in the scheme'
          WRITE(6,*) 'gas-phase=',gspe
          WRITE(6,*) 'check read_idgaw subroutine'
          STOP 'STOP in read_idgaw subroutine'
        endif
        idgsat(i)=j
        satid(j)= i
      ENDDO satgloop
!$OMP END PARALLEL DO

! -----------------
! particle phase species
! -----------------
!$OMP PARALLEL DO private(i,aspe,j,k)
      sataloop: DO i=1,nsat
        aspe=namsat(i)
        aspe(1:1)="A"
        call akspnum(aspe,chrsp,numsp,j)
! do loop exit should occur if species found -
        if (j == 0) then
          WRITE(6,*) '--error--, the following species in nannoolal '
          WRITE(6,*) '           not in list of species in the scheme'
          WRITE(6,*) 'aerosol-phase=',aspe
          WRITE(6,*) 'check read_idgaw subroutine'
          STOP 'STOP in read_idgaw subroutine'
        endif
        idasat(i)=j
        satid(j)= i
      ENDDO sataloop
!$OMP END PARALLEL DO


! -----------------
! wall phase species
! -----------------
      IF(wall_fg .eq. 1) THEN

        satwloop: DO i=1,nsat
          wspe=namsat(i)
          wspe(1:1)="W"
          call akspnum(wspe,chrsp,numsp,j)
! do loop exit should occur if species found -
          if (j == 0) then
            WRITE(6,*) '--error--, the following species in nannoolal '
            WRITE(6,*) '           not in list of species in the scheme'
            WRITE(6,*) 'wall-phase=',wspe
            WRITE(6,*) 'check read_idgaw subroutine'
            STOP 'STOP in read_idgaw subroutine'
          endif
          idwsat(i) = j
          satid(j)  = i
        ENDDO satwloop
      ENDIF

! DEBUG
!      DO i=1,nsat
!          print*,i,namsat(i),idgsat(i),idasat(i),idwsat(i)
!      ENDDO
! -----------------
! dimer species
! -----------------
      IF (dimer_fg.EQ.0) RETURN
      k=1
      satdloop: DO i=1,nsat
        aspe=namsat(i)
        aspe(1:1)="D"
        DO j=1,numsp
          IF (aspe.eq.chrsp(j)) THEN
            iddsat(i) = j
            satid(j)  = i
            k=j+1
            ndim=ndim+1
            CYCLE satdloop
          ENDIF
        ENDDO
! do loop exit should occur if species found -
!        WRITE(6,*) '--error--, the following species in nannoolal '
!        WRITE(6,*) '           not in list of species in the scheme'
!        WRITE(6,*) 'dimer=',aspe
!        WRITE(6,*) 'check read_idgaw subroutine'
!        STOP 'STOP in read_idgaw subroutine'

      ENDDO satdloop

! end of the routine
! -----------------
      END
