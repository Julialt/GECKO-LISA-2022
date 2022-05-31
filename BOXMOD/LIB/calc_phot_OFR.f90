      SUBROUTINE calc_phot_OFR(cbox,valj,ratact)

      USE OFR_params_module
      USE module_data_gecko_main

      IMPLICIT NONE

      REAL,DIMENSION(maxsp) :: cbox ! = conc(1:numsp,ibox)
      REAL,DIMENSION(maxhv) :: ratact ! = phot rates by rxn
      REAL,DIMENSION(nchrom) :: valj
                                      ! = phot rates by chromophore

!------------------------------------------
! PAM OFR CONDITIONS (assume 1-box simulation)
! Expects the following values for idhv:
!  1) GO3 + HV => GO1D
!  2) GO3 + HV => GO3P
! 10) GO2 + HV => 2.  GO3P
! 11) GH2O + HV => GHO2 + GHO
! NOTE THAT HV RXNS MUST BE SUPPLIED IN CORRECT ORDER !!
!-------------------------------------------------------
! - zero out all photolysis rates
! - calculate jO2 and jH2O from f185 in KEY file, attenuate
! - ATTENUATION: jO3 is a function of [O3], photons
! OFR185: jvalues for all 3 species (OFR185 uses BOTH wavelengths)
! OFR254: jvalues only for O3 (f185 = 0.)
! => assign jh2o and ho2 to valj(_) ONLY if f185 > 0.
! note:   valj(idhv(2)) = 0.  
!------------------------------------------

! zero out all photolysis rates (IF you don't want to use the lamp spectrum)
!      ratact = 0.
!      valj = 0.

! H2O & O2 concentrations .... are now set in find_concs.f90
      conc(idh2o,1) = water(1)
      conc(ido2dic,1) = sumc(1)*0.2

      IF (f185.GT.0.) THEN
!JMLT - calculate jh2o using F185 input from Keyfile
!     - jh2o corresponds to 11th HV reaction
        i = 11
        k=mtopchromo+ntmedchromo+i

        jh2o = sigmh2o_185*f185* &
               EXP(-7.0*(sigmo2_185*cbox(ido2dic) + &
               sigmh2o_185*cbox(idh2o)))

        ratact(i) = jh2o
        valj(k) = jh2o

!JMLT - calculate jo2 using F185 input from Keyfile
!     - jo2 corresponds to 10th HV reaction
        i = 10
        k=mtopchromo+ntmedchromo+i

        jo2 = sigmo2_185*f185* &
              EXP(-7.0*(sigmo2_185*cbox(ido2dic) + &
              sigmh2o_185*cbox(idh2o)))

        ratact(i) = jo2
        valj(k) = jo2

      ENDIF

!JMLT - calculate jo3 using F254 input from Keyfile
!     - jo3 corresponds to 1st HV reaction
      i = 1
      k=mtopchromo+ntmedchromo+i

      jo3 = sigmo3_254*f254 * EXP(-7.0*sigmo3_254*cbox(ido3))
  
      ratact(i) = jo3
      valj(k) = jo3

! DEBUG
!      PRINT*,"after calc_phot_OFR"
!      PRINT*,"valj(h2o), valj(o2), valj(o3)"
!      PRINT*,valj(mtopchromo+ntmedchromo+11), &
!             valj(mtopchromo+ntmedchromo+10), &
!             valj(mtopchromo+ntmedchromo+1)

!------------------------------------------
      END SUBROUTINE calc_phot_OFR
!------------------------------------------
