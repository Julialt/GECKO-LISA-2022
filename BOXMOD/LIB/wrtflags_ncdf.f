      SUBROUTINE wrtflags_ncdf(ncid)

!===================================================================
! PURPOSE: write boxmodel flags into new o/p file.
!           (NetCDF version of GECKO-A mechanism)
! AUTHOR: Julia Lee-Taylor, NCAR, 18 Jan 2018
!===================================================================

      USE flags_module,ONLY: iofmt_fg, prevflag, lagflag, depos_fg,
     &                       mixing_fg, emis_fg, noem_fg, soa_fg, icham,
     &                        pvap_fg, OFR_fg, ro2_fg, dyn_fg, wall_fg, 
     &                        dimer_fg, tracer_fg, seedtyp_fg, iscape,
     &                        print_steadystate_fg, reacrate_fg

      USE akparameter_module

      IMPLICIT NONE

      INTEGER ncid

! INPUT FLAGS
!      CALL eznc_put_0Dint(ncid,"iofmt_fg",iofmt_fg)
      CALL eznc_put_0Dint(ncid,"prevflag",prevflag)
      CALL eznc_put_0Dint(ncid,"lagflag",lagflag)
      CALL eznc_put_0Dint(ncid,"depos_fg",depos_fg)
!      CALL eznc_put_0Dint(ncid,"mixing_fg",mixing_fg)
!      CALL eznc_put_0Dint(ncid,"emis_fg",emis_fg)
      CALL eznc_put_0Dreal(ncid,"noem_fg",noem_fg)
      CALL eznc_put_0Dint(ncid,"soa_fg",soa_fg)
      CALL eznc_put_0Dint(ncid,"pvap_fg",pvap_fg)
      CALL eznc_put_0Dint(ncid,"dimer_fg",dimer_fg)
      CALL eznc_put_0Dint(ncid,"ro2_fg",ro2_fg)
      CALL eznc_put_0Dint(ncid,"dyn_fg",dyn_fg)
      CALL eznc_put_0Dint(ncid,"seedtyp_fg",seedtyp_fg)
      CALL eznc_put_0Dint(ncid,"wall_fg",wall_fg)
      CALL eznc_put_0Dint(ncid,"icham",icham)
      CALL eznc_put_0Dint(ncid,"iscape",iscape)
      CALL eznc_put_0Dint(ncid,"tracer_fg",tracer_fg)
      CALL eznc_put_0Dint(ncid,"print_steadystate_fg",
     $                            print_steadystate_fg)
      CALL eznc_put_0Dint(ncid,"reacrate_fg",reacrate_fg)

       END SUBROUTINE wrtflags_ncdf
