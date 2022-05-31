      SUBROUTINE defflags_ncdf(ncid)

!===================================================================
! PURPOSE: write boxmodel flags into new o/p file.
!           (NetCDF version of GECKO-A mechanism)
! AUTHOR: Julia Lee-Taylor, NCAR, 18 Jan 2018
!===================================================================

      USE netcdf
      USE akparameter_module
      USE flags_module

      IMPLICIT NONE

      INTEGER ncid

! INPUT FLAGS: from TOP of indat.key
      CALL eznc_def_0Dint(ncid,"depos_fg")
      CALL eznc_def_localatt(ncid,"depos_fg","keyword","DEPO")
      CALL eznc_def_localatt(ncid,"depos_fg","description",
     $ "0: no deposition. 1: Henry deposition")

      CALL eznc_def_0Dint(ncid,"dimer_fg")
      CALL eznc_def_localatt(ncid,"dimer_fg","keyword","DIMR")
      CALL eznc_def_localatt(ncid,"dimer_fg","description",
     $ "0: no dimerisation.    1: yes dimerisation")

      CALL eznc_def_0Dint(ncid,"dyn_fg")
      CALL eznc_def_localatt(ncid,"dyn_fg","keyword","DYNF")
      CALL eznc_def_localatt(ncid,"dyn_fg","description",
     $ "0: mtrat gas -> aero coeff. 1: mtratdyn gas-> aero description")

!      CALL eznc_def_0Dint(ncid,"emis_fg")
!      CALL eznc_def_localatt(ncid,"emis_fg","keyword","EMIS")
!      CALL eznc_def_localatt(ncid,"emis_fg","description",
!     $ "")

      CALL eznc_def_0Dint(ncid,"icham")
      CALL eznc_def_localatt(ncid,"icham","keyword","CHAM")
      CALL eznc_def_localatt(ncid,"icham","description",
     $ "1: Matsunaga/Ziemann, 2: Ziemann-CU, 3: Jimenez-CU ")

!      CALL eznc_def_0Dint(ncid,"iofmt_fg")
!      CALL eznc_def_localatt(ncid,"iofmt_fg","keyword","IFMT")

      CALL eznc_def_0Dint(ncid,"isopsoa_fg")
      CALL eznc_def_localatt(ncid,"isopsoa_fg","keyword","ISOA")
      CALL eznc_def_localatt(ncid,"isopsoa_fg","description",
     $ "0/1: isoprene soa scheme off/on")

      CALL eznc_def_0Dint(ncid,"lagflag")
      CALL eznc_def_localatt(ncid,"lagflag","keyword","LAGR")
      CALL eznc_def_localatt(ncid,"lagflag","description",
     $ "0: repeating/Eulerian. 1: developing/Lagrangian")

!      CALL eznc_def_0Dint(ncid,"mixing_fg")
!      CALL eznc_def_localatt(ncid,"mixing_fg","keyword","MIX")
!      CALL eznc_def_localatt(ncid,"mixing_fg","description",
!     $ "")

      CALL eznc_def_0Dreal(ncid,"noem_fg")
      CALL eznc_def_localatt(ncid,"noem_fg","keyword","NOEM")
      CALL eznc_def_localatt(ncid,"noem_fg","description",
     $ "")

      CALL eznc_def_0Dint(ncid,"OFR_fg")
      CALL eznc_def_localatt(ncid,"OFR_fg","keyword","OFR")
      CALL eznc_def_localatt(ncid,"OFR_fg","description",
     $ " 0: normal mode.   1: OFR mode "//
     $ "(requires OFR-specific mechanism)")
      CALL eznc_def_localatt(ncid,"OFR_fg","note",
     $ "BOXMOD only runs if OFR_fg agrees with original generator flag")

      CALL eznc_def_0Dint(ncid,"print_steadystate_fg")
      CALL eznc_def_localatt(ncid,"print_steadystate_fg","keyword",
     $                                                      "PSSF")
      CALL eznc_def_localatt(ncid,"print_steadystate_fg",
     $                            "description",
     $ "output end-of-run values separately")

      CALL eznc_def_0Dint(ncid,"reacrate_fg")
      CALL eznc_def_localatt(ncid,"reacrate_fg","keyword","RRAT")
      CALL eznc_def_localatt(ncid,"reacrate_fg","description",
     $ "output all reaction rates at each output time")

      CALL eznc_def_0Dint(ncid,"prevflag")
      CALL eznc_def_localatt(ncid,"prevflag","keyword","PREV")
      CALL eznc_def_localatt(ncid,"prevflag","description",
     $ "output from a previous simulation is used as input")

      CALL eznc_def_0Dint(ncid,"pvap_fg")
      CALL eznc_def_localatt(ncid,"pvap_fg","keyword","PVAP")
      CALL eznc_def_localatt(ncid,"pvap_fg","description",
     $ "vapor pressure scheme for gas-particle transfer: "//
     $ "2: Nannoolal. 3: SIMPOL.")

      CALL eznc_def_0Dint(ncid,"ro2_fg")
      CALL eznc_def_localatt(ncid,"ro2_fg","keyword","RO2F")
      CALL eznc_def_localatt(ncid,"ro2_fg","description",
     $ "0/1: RO2+RO2 deactivated/allowed (default = 1)")

      CALL eznc_def_0Dint(ncid,"seedtyp_fg")
      CALL eznc_def_localatt(ncid,"seedtyp_fg","keyword","NVID")
      CALL eznc_def_localatt(ncid,"seedtyp_fg","description",
     $ "0: organic.  1: inorganic")

      CALL eznc_def_0Dint(ncid,"iscape")
      CALL eznc_def_localatt(ncid,"iscape","keyword","SCAP")
      CALL eznc_def_localatt(ncid,"iscape","description",
     $ "0/1 : no/yes thermodynamic representation.")

      CALL eznc_def_0Dint(ncid,"soa_fg")
      CALL eznc_def_localatt(ncid,"soa_fg","keyword","SOAF")
      CALL eznc_def_localatt(ncid,"soa_fg","description",
     $ "0: no soa module. 1: soa but no unifac. 2: soa with mass "//
     $ "transfer")

      CALL eznc_def_0Dint(ncid,"tracer_fg")
      CALL eznc_def_localatt(ncid,"tracer_fg","keyword","?")
      CALL eznc_def_localatt(ncid,"tracer_fg","description",
     $ "")

      CALL eznc_def_0Dint(ncid,"wall_fg")
      CALL eznc_def_localatt(ncid,"wall_fg","keyword","WALL")
      CALL eznc_def_localatt(ncid,"wall_fg","description",
     $ "0: no wall losses.     1: wall losses considered")

      END SUBROUTINE defflags_ncdf
