!**********************************************************************
! This subroutine returns the vapor pressure of a list of species     *
! by the nanoonal estimation method                                   *
!                                                                     *
! INPUT :                                                             *
!    - nsat  : number of species for which vapor                      *
!    - temp  : temperature                                            *
!    - Tb(i) : boiling point                                          *
!    - db(i) :                                                        *
!                                                                     *
! OUTPUT :                                                            *
!    - pvap(i) : vapor pressure of the species (atm)                  *
!**********************************************************************
      SUBROUTINE getpvap_nan(nsat,temp,tb,dB,psat)
      !$ use OMP_LIB
      USE akparameter_module
      IMPLICIT NONE

      INTEGER  nsat
      REAL     temp
      REAL     Tb(mxsat)
      REAL     dB(mxsat)
      REAL     psat(mxsat)

! local
      REAL     Trb
      REAL     logPvap
      INTEGER  i

! --------------------
! COMPUTE Pvap
! --------------------
! compute the vapor pressure for temp - unit is atm
!$OMP PARALLEL DO private(i, Trb, logPvap)
      DO i=1,nsat
         Trb = temp/Tb(i)
         logPvap = (4.1012+dB(i))*((Trb-1.)/(Trb-(1./8.)))
         psat(i) = exp(logPvap*log(10.))
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END
