      SUBROUTINE defenvinp_ncdf(ncid)

!==================================================================
! PURPOSE: define space for environmental input variables 
!          in NetCDF-format box model output file
! NOTES: Refer to spreadkey6 for origin of most parameters.
!        Not all argument parameters are used in this subroutine:
!        they are present merely as reminders of purpose
! AUTHOR: Julia Lee-Taylor, NCAR, 123Jan 2018
!==================================================================

      USE netcdf
      USE akparameter_module
      USE flags_module
      USE time_mgmt_module
      USE constraints_module,ONLY : constrained_thing
      USE solver_params_module,ONLY : atol,rtol,dtmin
      USE OFR_params_module
      USE forcing_params_module
      USE module_data_gecko_main
      USE steadystate_vars_module,ONLY: non_zero_species_number

      IMPLICIT NONE

      INTEGER ncid
      INTEGER nemis       ! # emitted spp

!----------------------------------------------------------
      CALL eznc_def_0Dreal(ncid,"tstart")
      CALL eznc_def_localatt(ncid,"tstart","keyword","TSTR")
      CALL eznc_def_localatt(ncid,"tstart","units","s")
      CALL eznc_def_localatt(ncid,"tstart","title", &
                       "time at start of simulation")

      CALL eznc_def_0Dreal(ncid,"tstop")
      CALL eznc_def_localatt(ncid,"tstop","keyword","TSTP")
      CALL eznc_def_localatt(ncid,"tstop","units","s")
      CALL eznc_def_localatt(ncid,"tstop","title", &
                       "time at end of simulation")

      CALL eznc_def_0Dreal(ncid,"tlen")
      CALL eznc_def_localatt(ncid,"tlen","keyword","TLEN")
      CALL eznc_def_localatt(ncid,"tlen","units","s")
      CALL eznc_def_localatt(ncid,"tlen","title", &
                       "length of simulation timestep")

!      CALL eznc_def_0Dint(ncid,"ntstep")
!      CALL eznc_def_localatt(ncid,"ntstep","title", &
!                       "number of simulation timesteps")

!      CALL eznc_def_0Dint(ncid,"ntprint")
!      CALL eznc_def_localatt(ncid,"ntprint","title", &
!                       "number of output timesteps")

      CALL eznc_def_0Dint(ncid,"nskip")
      CALL eznc_def_localatt(ncid,"nskip","keyword","SKIP")
      CALL eznc_def_localatt(ncid,"nskip","title", &
              "number of simulation timesteps between output times")

! NB: ntout is set here because it is a fixed dimension: it must be
! set AFTER it is calculated (in get_envinp), and written in wrtenvinp_ncdf

      !CALL eznc_def_dim(ncid,"ntout",ntout)
      CALL eznc_def_dim_unlim(ncid,"ntout")

!----------------------------------------------------------------
! cbot,ctop are defined elsewhere, since they are also outputs
      CALL eznc_def_0Dreal(ncid,"rtol")
      CALL eznc_def_localatt(ncid,"rtol","keyword","RTOL")
      CALL eznc_def_localatt(ncid,"rtol","title", &
                       "relative solver tolerance")

      CALL eznc_def_0Dreal(ncid,"atol")
      CALL eznc_def_localatt(ncid,"atol","keyword","ATOL")
      CALL eznc_def_localatt(ncid,"atol","title", &
                       "absolute solver tolerance")

      CALL eznc_def_0Dreal(ncid,"dtmin")
      CALL eznc_def_localatt(ncid,"dtmin","keyword","DTMN")
      CALL eznc_def_localatt(ncid,"dtmin","units","s")
      CALL eznc_def_localatt(ncid,"dtmin","title", &
                       "minimum delta-time for solver")

      CALL eznc_def_0Dint(ncid,"nbox")
      CALL eznc_def_localatt(ncid,"nbox","keyword","NBOX")
      CALL eznc_def_localatt(ncid,"nbox","title", &
                       "number of actively-solved boxes (= 1 or 2)")

      CALL eznc_def_0Dreal(ncid,"cnv")
      CALL eznc_def_localatt(ncid,"cnv","keyword","SEED")
      CALL eznc_def_localatt(ncid,"cnv","units","molec cm-3")
      CALL eznc_def_localatt(ncid,"cnv","title", &
                       "concentration of non-volatile seed particles")

      CALL eznc_def_0Dreal(ncid,"Mp")
      CALL eznc_def_localatt(ncid,"Mp","keyword","NVMW")
      CALL eznc_def_localatt(ncid,"Mp","units","a.m.u.")
      CALL eznc_def_localatt(ncid,"Mp","title", &
                 "mean molecular mass of non-volatile seed particles")

      CALL eznc_def_0Dreal(ncid,"Rpo")
      CALL eznc_def_localatt(ncid,"Rpo","keyword","NVRO")
      CALL eznc_def_localatt(ncid,"Rpo","units","cm")
      CALL eznc_def_localatt(ncid,"Rpo","title", &
                       "initial radius of non-volatile seed particles")

      CALL eznc_def_0Dreal(ncid,"gamm")
      CALL eznc_def_localatt(ncid,"gamm","keyword","GAMM")
      CALL eznc_def_localatt(ncid,"gamm","units","none")
      CALL eznc_def_localatt(ncid,"gamm","title", &
                       "bulk activity coeff for SOA")

      CALL eznc_def_1Dreal(ncid,"cbg","maxsp")
      CALL eznc_def_localatt(ncid,"cbg","units","molec cm-3")
      CALL eznc_def_localatt(ncid,"cbg","title", &
                       "background gas phase concentrations")
      CALL eznc_def_localatt(ncid,"cbg","congruence","chrsp")

!----------------------------------------------------------------
! initialized concentrations
      non_zero_species_number = COUNT(conc(:,1) > 1.)
      CALL eznc_def_dim(ncid,"maxinit",non_zero_species_number)

      CALL eznc_def_0Dint(ncid,"ninit")
      CALL eznc_def_localatt(ncid,"ninit","keyword","REAC")
      CALL eznc_def_localatt(ncid,"ninit","title", &
                                  "# of initialized species")

      CALL eznc_def_1Dint(ncid,"idinit","maxinit")
      CALL eznc_def_localatt(ncid,"idinit","title", &
                                  "chrsp id of initialized species")

      CALL eznc_def_1Dchar(ncid,"initnam","maxlsp","maxinit")
      CALL eznc_def_localatt(ncid,"initnam","keyword","REAC: names")
      CALL eznc_def_localatt(ncid,"initnam","title", &
                                  "initialized species names")

      CALL eznc_def_2Dreal(ncid,"initconc","maxinit","nbox")
      CALL eznc_def_localatt(ncid,"initconc","keyword","REAC: concs")
      CALL eznc_def_localatt(ncid,"initconc","units","molec cm-3")
      CALL eznc_def_localatt(ncid,"initconc","title", &
                                  "initialized species concentrations")

      CALL eznc_def_1Dreal(ncid,"initcbg","maxinit")
      CALL eznc_def_localatt(ncid,"initcbg","keyword", &
                                  "REAC: concs (last column)")
      CALL eznc_def_localatt(ncid,"initcbg","units","molec cm-3")
      CALL eznc_def_localatt(ncid,"initcbg","title", &
                                  "background gas phase concentrations")
      CALL eznc_def_localatt(ncid,"initcbg","congruence","initnam")

!----------------------------------------------------------------
! nemis is defined as a variable since defining a dimension might
! create problems for sequential runs with different parameters
      nemis = COUNT(emi_spec(:)%npoints/=0)
      CALL eznc_def_0Dint(ncid,"nemis")
      CALL eznc_def_localatt(ncid,"nemis","keyword","EMIS")
      CALL eznc_def_localatt(ncid,"nemis","title", &
                    "# of emitted species")

      IF(nemis.GT.0)THEN
        CALL eznc_def_1Dint(ncid,"idemis","maxem")
        CALL eznc_def_localatt(ncid,"idemis","title", &
                              "chrsp id of emitted species")

        CALL eznc_def_1Dchar(ncid,"eminam","maxlsp","maxem")
        CALL eznc_def_localatt(ncid,"eminam","keyword","EMIS: names")
        CALL eznc_def_localatt(ncid,"eminam","title", &
                                   "name of emitted species")

        CALL eznc_def_1Dint(ncid,"ntem","maxem")
        CALL eznc_def_localatt(ncid,"ntem","keyword","EMIS: #times")
        CALL eznc_def_localatt(ncid,"ntem","title", &
                                   "# of times in emission input")

        CALL eznc_def_2Dreal(ncid,"emtim","maxinput","maxem")
        CALL eznc_def_localatt(ncid,"emtim","keyword","EMIS: times")
        CALL eznc_def_localatt(ncid,"emtim","units","s")
        CALL eznc_def_localatt(ncid,"emtim","title", &
                                    "times of emission input data")

        CALL eznc_def_2Dreal(ncid,"emval","maxinput","maxem")
        CALL eznc_def_localatt(ncid,"emval","keyword","EMIS: values")
        CALL eznc_def_localatt(ncid,"emval","units","molec cm-2 s-1")
        CALL eznc_def_localatt(ncid,"emval","title", &
                                    "emission input rates")
      ENDIF ! (nemis.GT.0)
      
! surface emissions: one species_data for each surface
      CALL eznc_def_1Dint(ncid, "surf_nemis", "msur")
      CALL eznc_def_localatt(ncid,"surf_nemis","keyword","SEMI")
      CALL eznc_def_localatt(ncid,"surf_nemis","title", &
        "actual # of emitted species for each surface type")
     
! create surface emissions arrays only if any one of them is used
      IF (ANY(surface_emi%nemis > 0)) THEN
        CALL eznc_def_2Dint(ncid,"surf_idemis","maxem", "msur")
        CALL eznc_def_localatt(ncid,"surf_idemis","title", &
          "chrsp id of emitted species for surface i")
        
        CALL eznc_def_2Dchar(ncid,"surf_eminam", &
                            "maxlsp","maxem","msur")
        CALL eznc_def_localatt(ncid,"surf_eminam","keyword", &
                                    "SEMI: names")
        CALL eznc_def_localatt(ncid,"surf_eminam","title", &
           "name of emitted species for surface i")
     
        CALL eznc_def_2Dint(ncid,"surf_ntem","maxem", "msur")
        CALL eznc_def_localatt(ncid,"surf_ntem","keyword", &
                                    "SEMI: #times")
        CALL eznc_def_localatt(ncid,"surf_ntem","title", &
           "# of times in emission input for surface i")   
     
        CALL eznc_def_3Dreal(ncid,"surf_emtim", &
                                  "maxinput","maxem","msur")
        CALL eznc_def_localatt(ncid,"surf_emtim","keyword", &
                                    "SEMI: times")
        CALL eznc_def_localatt(ncid,"surf_emtim","units","s") 
        CALL eznc_def_localatt(ncid,"surf_emtim","title", &
           "times of emission for surface i")

        CALL eznc_def_3Dreal(ncid,"surf_emval", &
                                  "maxinput","maxem", "msur")
        CALL eznc_def_localatt(ncid,"surf_emval","keyword", &
                                    "SEMI: values")
        CALL eznc_def_localatt(ncid,"surf_emval","units", &
                              "molec cm-2 s-1")        
        CALL eznc_def_localatt(ncid,"surf_emval","title", &
            "emission rates for surface i")
      ENDIF

!----------------------------------------------------------------
! ncons is defined as a variable (see above note for nemis)
      
      CALL eznc_def_1Dint(ncid, "ncons","nbox")
      CALL eznc_def_localatt(ncid,"ncons","keyword","CONS")
      CALL eznc_def_localatt(ncid,"ncons","title", &
                    "actual # of constrained species in each box")
     
      IF (ANY(cons_spec%activefg)) THEN
        CALL eznc_def_2Dint(ncid,"idcons","maxconst", "nbox")
        CALL eznc_def_localatt(ncid,"idcons","title", &
                    "chrsp ids of constrained species in each box")
     
        CALL eznc_def_2Dchar(ncid,"consnam",  &
                            "maxlsp","maxconst","nbox")
        CALL eznc_def_localatt(ncid,"consnam","keyword","CONS: names")
        CALL eznc_def_localatt(ncid,"consnam","title", &
                    "names of constrained species in each box")
     
        CALL eznc_def_2Dint(ncid,"ntcons","maxconst","nbox")
        CALL eznc_def_localatt(ncid,"ntcons","keyword","CONS: #times")
        CALL eznc_def_localatt(ncid,"ntcons","title", &
                    "# of times in constraint input in each box")
     
        CALL eznc_def_3Dreal(ncid,"constim",&
                            "maxinput","maxconst","nbox")
        CALL eznc_def_localatt(ncid,"constim","keyword","CONS: times")
        CALL eznc_def_localatt(ncid,"constim","units","s")
        CALL eznc_def_localatt(ncid,"constim","title", &
                    "times of constraint in each box")
        
        CALL eznc_def_3Dreal(ncid,"consval", &
                            "maxinput","maxconst","nbox")
        CALL eznc_def_localatt(ncid,"consval","keyword","CONS: values")
        CALL eznc_def_localatt(ncid,"consval","units","molec cm-3")
        CALL eznc_def_localatt(ncid,"consval","title", &
                    "constrained concentrations in each box")
      ENDIF
      
!----------------------------------------------------------------
      CALL eznc_def_0Dint(ncid,"njf")
      CALL eznc_def_localatt(ncid,"njf","keyword","JFAC: #times")
      CALL eznc_def_localatt(ncid,"njf","title", &
                     "number of j-value adjustment factors")
      IF (njf.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"jftim","mtim")
        CALL eznc_def_localatt(ncid,"jftim","keyword","JFAC: times")
        CALL eznc_def_localatt(ncid,"jftim","units","s")

        CALL eznc_def_1Dreal(ncid,"jfval","mtim")
        CALL eznc_def_localatt(ncid,"jfval","keyword","JFAC: values")
        CALL eznc_def_localatt(ncid,"jfval","title", &
                     "multiplication factor for (all) j-values")
      ENDIF
      CALL eznc_def_0Dint(ncid,"szafix")
      CALL eznc_def_localatt(ncid,"szafix","keyword","SZAF")
      CALL eznc_def_localatt(ncid,"szafix","title", &
                     "flag for constrained solar zenith angle")
      IF (szafix.EQ.1) THEN
        CALL eznc_def_0Dreal(ncid,"szaval")
        CALL eznc_def_localatt(ncid,"szaval","keyword","SZAF: value")
        CALL eznc_def_localatt(ncid,"szaval","units", &
                                             "degrees from zenith")
        CALL eznc_def_localatt(ncid,"szaval","title", &
                     "constrained solar zenith angle")
      ENDIF

!----------------------------------------------------------------
      IF (jo2.GT.0)CALL eznc_def_0Dreal(ncid,"jo2")
      IF (jh2o.GT.0)CALL eznc_def_0Dreal(ncid,"jh2o")

      IF (f185.GT.0)THEN
        CALL eznc_def_0Dreal(ncid,"f185")
        CALL eznc_def_localatt(ncid,"f185","keyword","F185")
        CALL eznc_def_localatt(ncid,"f185","units","photons.cm-2.s-1")
        CALL eznc_def_localatt(ncid,"f185","title", &
            "photon flux at 185nm in the OFR") 
      ENDIF
      IF (f254.GT.0)THEN
        CALL eznc_def_0Dreal(ncid,"f254")
        CALL eznc_def_localatt(ncid,"f254","keyword","F254")
        CALL eznc_def_localatt(ncid,"f254","units","photons.cm-2.s-1")
        CALL eznc_def_localatt(ncid,"f254","title", &
            "photon flux at 254nm in the OFR") 
      ENDIF

!----------------------------------------------------------------
      CALL eznc_def_0Dint(ncid,"htop")
      CALL eznc_def_localatt(ncid,"htop","keyword","HTOP")
      CALL eznc_def_localatt(ncid,"htop","units","cm")
      CALL eznc_def_localatt(ncid,"htop","title", &
                    "height of top of top modeled box")

      CALL eznc_def_0Dint(ncid,"nhd")
      CALL eznc_def_localatt(ncid,"nhd","keyword","HBOX: #times")
      CALL eznc_def_localatt(ncid,"nhd","title", &
          "# of input data points for height of box 1 (mixed layer)")
      IF (nhd.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"boxt","mhd")
        CALL eznc_def_localatt(ncid,"boxt","keyword","HBOX: times")
        CALL eznc_def_localatt(ncid,"boxt","units","s")
        CALL eznc_def_1Dreal(ncid,"boxh","mhd")
        CALL eznc_def_localatt(ncid,"boxh","keyword","HBOX: values")
        CALL eznc_def_localatt(ncid,"boxh","units","cm")
        CALL eznc_def_localatt(ncid,"boxh","title", &
            "height of box 1 (mixed layer)")
      ENDIF

      CALL eznc_def_0Dreal(ncid,"vs")
      CALL eznc_def_localatt(ncid,"vs","keyword","SUBS")
      CALL eznc_def_localatt(ncid,"vs","units","cm.s-1")
      CALL eznc_def_localatt(ncid,"vs","title", &
          "tropospheric subsidence velocity")

      CALL eznc_def_0Dint(ncid,"ndil")
      CALL eznc_def_localatt(ncid,"ndil","keyword","DLTB: #times")
      CALL eznc_def_localatt(ncid,"ndil","title", &
          "# of input data points for rate of dilution of lowest box"//&
          " by background or by box 2")
      IF (ndil.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"diltim","mtim")
        CALL eznc_def_localatt(ncid,"diltim","keyword","DLTB: times")
        CALL eznc_def_localatt(ncid,"diltim","units","s")
        CALL eznc_def_localatt(ncid,"diltim","title", &
          "time points for dilution, lowest box")
        CALL eznc_def_1Dreal(ncid,"dilval","mtim")
        CALL eznc_def_localatt(ncid,"dilval","keyword","DLTB: values")
        CALL eznc_def_localatt(ncid,"dilval","units","s-1")
        CALL eznc_def_localatt(ncid,"dilval","title", &
          "dilution rate, lowest box")
      ENDIF

      CALL eznc_def_0Dint(ncid,"nmx")
      CALL eznc_def_localatt(ncid,"nmx","keyword","MIX: #times")
      CALL eznc_def_localatt(ncid,"nmx","title", &
          "# of input data points for bidrectional mixing velocity")
      IF (nmx.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"mixt","mhd")
        CALL eznc_def_localatt(ncid,"mixt","keyword","MIX: times")
        CALL eznc_def_localatt(ncid,"mixt","units","s")
        CALL eznc_def_localatt(ncid,"mixt","title", &
          "time points for bidirectional mixing velocity")
        CALL eznc_def_2Dreal(ncid,"mixv","nbox","mhd")
        CALL eznc_def_localatt(ncid,"mixv","keyword","MIX: values")
        CALL eznc_def_localatt(ncid,"mixv","units","cm.s-1")
        CALL eznc_def_localatt(ncid,"mixv","title", &
          "bidirectional mixing velocity across each box top")
      ENDIF

!----------------------------------------------------------------
      IF (tempm(1).GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"tempm","nbox")
        CALL eznc_def_localatt(ncid,"tempm","keyword","TEMP: mean")
        CALL eznc_def_localatt(ncid,"tempm","title", &
            "parameter m for temperature sine function")
        CALL eznc_def_1Dreal(ncid,"tempa","nbox")
        CALL eznc_def_localatt(ncid,"tempa","keyword","TEMP: amplitude")
        CALL eznc_def_localatt(ncid,"tempa","title", &
            "parameter a for temperature sine function")
        CALL eznc_def_1Dreal(ncid,"temptm","nbox")
        CALL eznc_def_localatt(ncid,"temptm","keyword", &
                                                  "TEMP: time(maxval)")
        CALL eznc_def_localatt(ncid,"temptm","title", &
            "parameter tm for temperature sine function")

!        CALL eznc_def_2Dreal(ncid,"TEMP","dim3","nbox")
!        CALL eznc_def_localatt(ncid,"TEMP","keyword","TEMP")
!        CALL eznc_def_localatt(ncid,"TEMP","title", &
!            "parameters m,a,tm for temperature sine function")
!        CALL eznc_def_localatt(ncid,"TEMP","units", &
!            "K, K, s")
      ENDIF

      CALL eznc_def_0Dint(ncid,"ntk")
      CALL eznc_def_localatt(ncid,"ntk","keyword","TKTB")
      CALL eznc_def_localatt(ncid,"ntk","title", &
          "# of input data points for temperature")
      IF (ntk.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"tktim","mtim")
        CALL eznc_def_localatt(ncid,"tktim","keyword","TKTB: times")
        CALL eznc_def_localatt(ncid,"tktim","units","s")
        CALL eznc_def_2Dreal(ncid,"tkval","nbox","mtim")
        CALL eznc_def_localatt(ncid,"tkval","keyword","TKTB: T values")
        CALL eznc_def_localatt(ncid,"tkval","units","K")
      ENDIF

!----------------------------------------------------------------
      IF (rhm(1).GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"rhm","nbox")
        CALL eznc_def_localatt(ncid,"rhm","keyword","RHUM: mean")
        CALL eznc_def_localatt(ncid,"rhm","units","%")
        CALL eznc_def_localatt(ncid,"rhm","title", &
            "parameter m for relative humidity sine function")
        CALL eznc_def_1Dreal(ncid,"rha","nbox")
        CALL eznc_def_localatt(ncid,"rha","keyword","RHUM: amplitude")
        CALL eznc_def_localatt(ncid,"rha","units","%")
        CALL eznc_def_localatt(ncid,"rha","title", &
            "parameter a for relative humidity sine function")
        CALL eznc_def_1Dreal(ncid,"rhtm","nbox")
        CALL eznc_def_localatt(ncid,"rhtm","keyword", &
                                                   "RHUM: time(RHmax)")
        CALL eznc_def_localatt(ncid,"rhtm","units","s")
        CALL eznc_def_localatt(ncid,"rhtm","title", &
            "parameter tm for relative humidity sine function")

!        CALL eznc_def_2Dreal(ncid,"RH","dim3","nbox")
!        CALL eznc_def_localatt(ncid,"RH","keyword","RH")
!        CALL eznc_def_localatt(ncid,"RH","units", &
!            "%, %, s")
      ENDIF

      CALL eznc_def_0Dint(ncid,"nrh")
      CALL eznc_def_localatt(ncid,"nrh","keyword","RHTB")
      CALL eznc_def_localatt(ncid,"nrh","title", &
          "# of input data points for relative humidity")
      IF (nrh.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"rhtim","mtim")
        CALL eznc_def_localatt(ncid,"rhtim","keyword","RHTB: times")
        CALL eznc_def_localatt(ncid,"rhtim","units","s")
        CALL eznc_def_2Dreal(ncid,"rhval","nbox","mtim")
        CALL eznc_def_localatt(ncid,"rhval","keyword","RHTB: values")
        CALL eznc_def_localatt(ncid,"rhval","units","%")
      ENDIF

      CALL eznc_def_0Dint(ncid,"npr")
      CALL eznc_def_localatt(ncid,"npr","keyword","PRTB")
      CALL eznc_def_localatt(ncid,"npr","title", &
            "# of input data points for prescribed pressure")
      IF (npr.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"prtim","mtim")
        CALL eznc_def_localatt(ncid,"prtim","keyword","PRTB: times")
        CALL eznc_def_localatt(ncid,"prtim","units","s")
        CALL eznc_def_2Dreal(ncid,"prval","nbox","mtim")
        CALL eznc_def_localatt(ncid,"prval","keyword","PRTB: values")
        CALL eznc_def_localatt(ncid,"prval","units","bar")
      ENDIF

!----------------------------------------------------------------
      CALL eznc_def_0Dint(ncid,"waterfix")
      CALL eznc_def_localatt(ncid,"waterfix","keyword","WATR")
      CALL eznc_def_localatt(ncid,"waterfix","title", &
                     "flag for constrained water vapor")
      IF (waterfix.EQ.1) THEN
        CALL eznc_def_1Dreal(ncid,"water","nbox")
        CALL eznc_def_localatt(ncid,"water","keyword","WATR: values")
        CALL eznc_def_localatt(ncid,"water","units","molec cm-3")
        CALL eznc_def_localatt(ncid,"water","title", &
          "atmospheric water concentration")
      ENDIF

      CALL eznc_def_0Dint(ncid,"sumcfix")
      CALL eznc_def_localatt(ncid,"sumcfix","keyword","NDEN")
      CALL eznc_def_localatt(ncid,"sumcfix","title", &
                    "flag for constrained number density")
      IF (sumcfix.EQ.1) THEN
        CALL eznc_def_1Dreal(ncid,"sumc","nbox")
        CALL eznc_def_localatt(ncid,"sumc","keyword","NDEN: values")
      ELSE
        CALL eznc_def_2Dreal(ncid,"sumc","nbox","ntout")
      ENDIF
      CALL eznc_def_localatt(ncid,"sumc","units","molec cm-3")
      CALL eznc_def_localatt(ncid,"sumc","title",& 
                             "atmospheric number density")

      CALL eznc_def_0Dint(ncid,"presfix")
      CALL eznc_def_localatt(ncid,"presfix","keyword","PRES")
      CALL eznc_def_localatt(ncid,"presfix","title", &
                                "flag for constrained atmos. pressure")
      IF (presfix.EQ.1) THEN
        CALL eznc_def_1Dreal(ncid,"prconst","nbox")
        CALL eznc_def_localatt(ncid,"prconst","keyword","PRES: values")
        CALL eznc_def_localatt(ncid,"prconst","units","bar")
        CALL eznc_def_localatt(ncid,"prconst","title", &
                              "constrained atmospheric pressure")
      ENDIF

      CALL eznc_def_1Dint(ncid,"noxfix","nbox")
      CALL eznc_def_localatt(ncid,"noxfix","keyword","NOx")
      CALL eznc_def_localatt(ncid,"noxfix","title", &
                    "flag for constrained NOx concentration")

!----------------------------------------------------------------
      IF (windm.GT.0) THEN
        CALL eznc_def_0Dreal(ncid,"windm")
        CALL eznc_def_localatt(ncid,"windm","keyword","WIND: mean")
        CALL eznc_def_localatt(ncid,"windm","title", &
            "parameter m for windspeed sine function")
        CALL eznc_def_0Dreal(ncid,"winda")
        CALL eznc_def_localatt(ncid,"winda","keyword","WIND: amplitude")
        CALL eznc_def_localatt(ncid,"winda","title", &
            "parameter a for windspeed sine function")
        CALL eznc_def_0Dreal(ncid,"windtm")
        CALL eznc_def_localatt(ncid,"windtm","keyword", &
                                              "WIND: time(max_wind)")
        CALL eznc_def_localatt(ncid,"windtm","title", &
            "parameter tm for windspeed sine function")
      ENDIF

      CALL eznc_def_0Dint(ncid,"nws")
      CALL eznc_def_localatt(ncid,"nws","keyword","WSTB")
      CALL eznc_def_localatt(ncid,"nws","title", &
          "# of input data points for wind speed")
      IF (nws.GT.0) THEN
        CALL eznc_def_1Dreal(ncid,"wstim","mtim")
        CALL eznc_def_localatt(ncid,"wstim","keyword","WSTB: times")
        CALL eznc_def_localatt(ncid,"wstim","units","s")
        CALL eznc_def_1Dreal(ncid,"wsval","mtim")
        CALL eznc_def_localatt(ncid,"wsval","keyword","WSTB: values")
        CALL eznc_def_localatt(ncid,"wsval","units","m s-1")
      ENDIF

!----------------------------------------------------------------
      CALL eznc_def_0Dreal(ncid,"saero")
      CALL eznc_def_localatt(ncid,"saero","keyword","AERO")
      CALL eznc_def_localatt(ncid,"saero","units","cm2.cm-3")
      CALL eznc_def_localatt(ncid,"saero","title", &
                                              "aerosol surface area")

      CALL eznc_def_0Dint(ncid,"nseed")
      CALL eznc_def_localatt(ncid,"nseed","keyword","SEET")
      CALL eznc_def_localatt(ncid,"nseed","title", &
          "# of input data points for seed aerosol")
      IF (nseed.GT.0) THEN
         CALL eznc_def_1Dreal(ncid,"tseed","mtim")
         CALL eznc_def_localatt(ncid,"tseed","keyword","SEET: times")
         CALL eznc_def_localatt(ncid,"tseed","units","s")
         CALL eznc_def_1Dreal(ncid,"cseed","mtim")
         CALL eznc_def_localatt(ncid,"cseed","keyword","SEET: values")
         CALL eznc_def_localatt(ncid,"cseed","units","molec cm-3")
      ENDIF

      CALL eznc_def_0Dint(ncid,"nsd")
      CALL eznc_def_localatt(ncid,"nsd","keyword","SURF")
      CALL eznc_def_localatt(ncid,"nsd","title", &
          "# input of data points for surface data")
      IF (nsd.GT.0) THEN
         CALL eznc_def_1Dreal(ncid,"surft","mhd")
         CALL eznc_def_localatt(ncid,"surft","keyword","SURF: times")
         CALL eznc_def_localatt(ncid,"surft","units","s")
         CALL eznc_def_2Dreal(ncid,"psurf","mhd","msur")
         CALL eznc_def_localatt(ncid,"psurf","keyword","SURF: values")
         CALL eznc_def_localatt(ncid,"psurf","title", &
                     "fraction of surface type {URB1,CULT,FLEA,FCON}")
         CALL eznc_def_localatt(ncid,"psurf","actual_size", &
                                     "(nsd,msur)")
      ENDIF

      CALL eznc_def_0Dreal(ncid,"isop_fac")
      CALL eznc_def_localatt(ncid,"isop_fac","keyword","ISOP")
      CALL eznc_def_localatt(ncid,"isop_fac","units","none")
      CALL eznc_def_localatt(ncid,"isop_fac","title", &
           "isoprene emission factor")

      CALL eznc_def_0Dreal(ncid,"mterp_fac")
      CALL eznc_def_localatt(ncid,"mterp_fac","keyword","MTER")
      CALL eznc_def_localatt(ncid,"mterp_fac","units","none")
      CALL eznc_def_localatt(ncid,"mterp_fac","title", &
           "monoterpenes emission factor")

!----------------------------------------------------------------
       IF(sla.GT.0)THEN
         CALL eznc_def_0Dreal(ncid,"sla")
         CALL eznc_def_localatt(ncid,"sla","keyword","PPHO: lat")
         CALL eznc_def_localatt(ncid,"sla","units","degrees north")
         CALL eznc_def_localatt(ncid,"sla","title","latitude")
         CALL eznc_def_0Dreal(ncid,"slo")
         CALL eznc_def_localatt(ncid,"slo","keyword","PPHO: lon")
         CALL eznc_def_localatt(ncid,"slo","units","degrees east")
         CALL eznc_def_localatt(ncid,"slo","title","longitude")
         CALL eznc_def_0Dreal(ncid,"tz")
         CALL eznc_def_localatt(ncid,"tz","keyword","PPHO: timezone")
         CALL eznc_def_localatt(ncid,"tz","units","hours to add to UT")
         CALL eznc_def_localatt(ncid,"tz","format","HH")
         CALL eznc_def_0Dint(ncid,"iy")
         CALL eznc_def_localatt(ncid,"iy","keyword","PPHO: year")
         CALL eznc_def_localatt(ncid,"iy","format","yyyy")
         CALL eznc_def_0Dint(ncid,"im")
         CALL eznc_def_localatt(ncid,"im","keyword","PPHO: month")
         CALL eznc_def_localatt(ncid,"im","format","mm")
         CALL eznc_def_0Dint(ncid,"id")
         CALL eznc_def_localatt(ncid,"id","keyword","PPHO: day")
         CALL eznc_def_localatt(ncid,"id","format","dd")
       ENDIF

!----------------------------------------------------------------
!       CALL eznc_def_0Dint(ncid,"iscape")
!       CALL eznc_def_localatt(ncid,"iscape","keyword","SCAP")
       CALL eznc_def_0Dint(ncid,"iseas")
       CALL eznc_def_localatt(ncid,"iseas","keyword","SEAS")

!==end define keyfile inputs

!----------------------------------------------------------------
!==define mask value for concentration arrays
      CALL eznc_def_0Dreal(ncid,"maskval")
      CALL eznc_def_localatt(ncid,"maskval","title", &
          "threshold value for writing model output")

      END SUBROUTINE defenvinp_ncdf
