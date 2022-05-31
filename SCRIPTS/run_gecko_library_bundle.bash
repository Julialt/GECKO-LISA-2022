##!/bin/bash
# run gecko library using bundled submission
# to be run from the SCRIPT directory, otherwise job submissions will fail
#----------------------------------------------------------------------
##### LIBRARY SET #####

## (don't need "precursors" list UNLESS you're generating.)
precursors=(\
"CH3CH3" \
)
#"#CH3Cd(CH3)=CdHCH2CdH=Cd(CH3)CdH=CdH2" \
#"CH3CH2CH3" \
#"#CH3Cd(=CdH2)CdH=CdH2" \
#"CH3Cd(CH3)=CdHCH2CH2Cd(=CdH2)CdH=CdH2" \

#---------------------------------
names=(\
"ethane"\
)
#"ocimene" \
#"propane" \
#"isoprene" \
#"bmyrcene" \
#---------------------------------

#choices=(spinup)
#choices=(spinup gen_sp run_sp) 
#choices=(spinup run_sp) 
#choices=(gen_sp)
choices=(gen_sp run_sp) 
#choices=(run_sp)
#choices=(post_sp)

#---------------------------------
scenarios=(remote remotecontinental continental pollutedcontinental urban)
#scenarios=(remote remotecontinental continental pollutedcontinental)
#scenarios=(remotecontinental continental pollutedcontinental urban)
#scenarios=(remote remotecontinental pollutedcontinental)
#scenarios=(remote urban)
#scenarios=(continental)
#scenarios=(pollutedcontinental)
#scenarios=(remotecontinental)
#scenarios=(remote)
#scenarios=(urban)
#----------------------------------------------------------------------

source cheyenne_scripting_functions.bash
pathfile="setup.dat"

echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
gecko_inp_dir=`  grep "gecko_inp_dir"  ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `

script_directory=$PWD
echo Working in ${script_directory}
if [ ! -e ./run_boxmod_cheyenne.bash ]; then
  echo ERROR: ./run_boxmod_cheyenne.bash is not available in ${script_directory}
  exit
fi

#----------------------------------------------------------------------
### !!!GENERATE SPINUP MECH FIRST SEPARATELY!!! ###
# generate the chemical mechanism for the spinup (only inorganic and background mechanism)
# using the command:
#./generate_scheme_cheyenne.bash -i "CH3CH3" -s settings_spinup -m spinup | grep last_job_id | cut -d= -f2
### !!!END!!! ###

### !!! THEN DO ./run_gecko_library_bash (this script) !!! ###
### The following assumes that the spinup mechanism has been generated beforehand ###

for c in ${!choices[*]}; do
  choice=${choices[c]}
  echo "selecting "$choice
  case "$choice" in
   
### RUN BOXMOD SPINUP SIMULATION ###
  spinup)
    echo "running "$choice

# the submission scripts should only run once the spinup mechanism is generated
    # flag list:
    # a = pathfile
    # c = concfile
    # f = inpdir
    # k = keyfile string
    # l = library_flag
    # m = mechname
    # n = netcdf_flag
    # p = photfile
    # r = runlength (seconds)
    # t = number of threads to use (default=16)

scenar_spinup_job_ids=()

for s in ${!scenarios[*]}; do
  scenar=${scenarios[$s]}
  mech='spinup'
  jobname=sbl_${mech}_${scenar}
  # run only on one processor: too much parallelization overhead for such a small run

# not dependent
      write_cheyenne_monoproc_script\
       ${jobname}.bash\
       ${jobname}\
       "00:01:00" \
       output_${jobname}\
       error_${jobname}\
       eval \
               "./run_boxmod_cheyenne.bash \
                -r 864000 -a setup.dat \
                -u INPUTS/hydrocarbons_library \
                -m ${mech} -k ${mech}_${scenar} \
                -p library.phot -n yes -t 16 "

# -l library flag makes script search for previous steadystate output
#    -> uses spinup in later runs
#    -> spinup itself should NOT use it!

    scenar_spinup_job_ids[$s]=`qsub ${home_dir}/GENERATED_SCRIPTS/${jobname}.bash`

    echo ../submitting spinup script 
    echo ${scenar_spinup_job_ids[$s]}: run spinup ${scenar}

done
;;
#-c ${mech}_${scenar}\

#-------------------------------------#
### GENERATE HYDROCARBON MECHANISMS ###
#-------------------------------------#
  gen_sp)
    echo "running "$choice

    # flag list:
    # a = pathfile
    # i = precursor
    # m = mechname
    # s = settings file
    # t = number of threads to use (default=16)

gen_job_ids=()
for p in ${!precursors[*]}; do # loop over indices of array (note ! before array's name)
  prec=${precursors[$p]}
  mech=${names[$p]}
  jobname=libgen_${mech}
  cheminput='cheminput.dat.'${mech}

  # run only on one processor: too much parallelization overhead for such a small run
# Do Not call generate_scheme_cheyenne_bundle.bash
# Must call launch script *directly* in order to get correct job id.
# Not dependent

      write_cheyenne_monoproc_script\
       ${jobname}.bash\
       ${jobname}\
       "06:00:00" \
       output_${jobname}\
       error_${jobname}\
       ./launch_gen_bundle.bash \
       ${home_dir}/$gecko_version \
       ${scratch_dir}/$gecko_run_dir \
       ${mech} 


    gen_job_ids[$p]=`qsub ${home_dir}/GENERATED_SCRIPTS/${jobname}.bash`

    echo ../submitting generator script 
    echo ${gen_job_ids[$p]}: generate ${mech} scheme

done
;;

#----------------------------------------------------------#
### RUN BOX MODEL WITH SCENARIO SET FOR EACH HYDROCARBON ###
#----------------------------------------------------------#
  run_sp)
    echo "running "$choice

for s in ${!scenarios[*]}; do
  scenar=${scenarios[$s]}
    echo selecting ${scenar}

  for p in ${!precursors[*]}; do
    mech=${names[$p]}
    jobname=sbl_${mech}_${scenar}
    echo run ${jobname}

    # we want the submission script to run only once the mechanism is generated and the spinup for this particular scenario has run
    # submit a small job that will wait until conditions are met to submit the real jobs    
    # flag list:
    # a = pathfile
    # c = concfile
    # f = inpdir
    # k = keyfile string
    # l = library_flag
    # m = mechname
    # n = netcdf_flag (for mechanism INPUT. Set flag in keyfile for NetCDF o/p)
    # p = photfile
    # r = runlength (seconds)
    # t = number of threads to use (default=16)

#------------------------------------------------------------------
# MANUALLY SELECT DEPENDENCY LIST:
     scenar_spinup_job=${scenar_spinup_job_ids[$s]}
     gen_job=${gen_job_ids[$p]}

#---if runs have NO dependency requirement: (e.g. prerequisites already exist)
    if [ ${#choices[@]} = 1 ]&&[ ${choices[0]} = "run_sp" ] ; then
      dependency_exp=none         

#---if mechs already generated, but spinup required:
    elif [ ${#choices[@]} = 2 ]&&[ ${choices[0]} = "spinup" ]&&[ ${choices[1]} = "run_sp" ] ; then
      dependency_exp=afterok:${scenar_spinup_job} 

#---if spinups already done, but mech generation required:
    elif [ ${#choices[@]} = 2 ]&&[ ${choices[0]} = "gen_sp" ]&&[ ${choices[1]} = "run_sp" ] ; then
      dependency_exp=afterok:${gen_job}  

#---if both mech generation and spinup required:
    elif [ ${#choices[@]} = 3 ] ; then
      dependency_exp=afterok:${gen_job}:${scenar_spinup_job} 
    fi

    echo "dependency "$dependency_exp

# END DEPENDENCY SELECTION
#------------------------------------------------------------------
# flags for boxmod run submission              
#------------------------------------------------------------
#    qsub \
#        -N ${name_job}\
#        -A P19010000\
#        -l walltime=00:01:00\
#        -q share\
#        -l select=1:ncpus=1\
#        -j oe\
#        -m abe\
#        -W depend=${dependency_exp}\
#        -o /glade/scratch/$USER/${name_job}\
#        -- ${script_directory}/run_boxmod_cheyenne.bash -l yes -r 864000  -u INPUTS/hydrocarbons_library \
#              -a setup.dat -m ${mech} -p sophie0404_julia_cmv -k ${scenar} -c ${scenar}
#------------------------------------------------------------
# OK to call run_*.bash since we don't need to evaulate job_id
# DO NOT add -d flag to "eval" statement - results in a double-dependency, i
# 2nd "wait" starts AFTER prerequisite job already finished

    if [ "$dependency_exp" = "none" ]; then
      write_cheyenne_monoproc_script\
       ${jobname}.bash\
       ${jobname}\
       "00:01:00" \
       output_${jobname}\
       error_${jobname}\
       eval \
               "./run_boxmod_bundle.bash \
                -l yes -r 864000 -a setup.dat \
                -u INPUTS/hydrocarbons_library \
                -m ${mech} -k ${scenar} -c ${scenar} \
                -f postproc_flags_library.input \
                -p library.phot -n yes -t 16 "
    else          
      write_cheyenne_monoproc_dependentscript\
       ${jobname}.bash\
       ${jobname}\
       "00:01:00" \
       output_${jobname}\
       error_${jobname}\
       ${dependency_exp} \
       eval \
               "./run_boxmod_bundle.bash \
                -l yes -r 864000 -a setup.dat \
                -u INPUTS/hydrocarbons_library \
                -m ${mech} -k ${scenar} -c ${scenar} \
                -f postproc_flags_library.input \
                -p library.phot -n yes -t 16 "
    fi
    
# write submission command into batch script          
    echo ../submitting box model script 
    qsub ${home_dir}/GENERATED_SCRIPTS/${jobname}.bash

  done
done
;;

### RUN POSTPROC SET FOR EACH HYDROCARBON ###
  post_sp)
    echo "running "$choice

for s in ${!scenarios[*]}; do
  scenar=${scenarios[$s]}
    echo scen ${scenar}

  #for p in ${!precursors[*]}; do
  #  mech=${names[$p]}
  for n in ${!names[*]}; do
    mech=${names[$n]}
    runname=${mech}_${scenar}_${scenar}
    jobname=ppr_${runname}
    walltime_post="4:00:00"
    echo run ${runname}

    dependency_exp="none"
    echo "dependency "$dependency_exp

    if [ "$dependency_exp" = "none" ]; then
      write_cheyenne_monoproc_script\
       ${jobname}.bash\
       ${jobname}\
       ${walltime_post}\
       output_postproc_${runname}\
       error_postproc_${runname}\
       ${home_dir}/SCRIPTS/run_postproc_cheyenne.bash\
       ${mech}\
       ${scenar}_${scenar}\
       postproc_flags_library.input
    else          
      write_cheyenne_monoproc_dependentscript\
       ${jobname}.bash\
       ${jobname}\
       ${walltime_post}\
       output_postproc_${runname}\
       error_postproc_${runname}\
       ${dependency_exp} \
       ${home_dir}/SCRIPTS/run_postproc_cheyenne.bash\
       ${mech}\
       ${scenar}_${scenar}\
       postproc_flags_library.input
    fi
             
# write submission command into batch script          
    echo ../submitting postprocessing script 
    qsub ${home_dir}/GENERATED_SCRIPTS/${jobname}.bash
              
  done
done
;;

esac   # end of case structure
done   # end of for-loop through possible cases

#-----------------SPECIES LISTS------------------------
# standard set
#propane \ "CH3CH2CH3" \
#butane \ "CH3CH2CH2CH3" \
#isobutane \ "CH3CH(CH3)CH3" \
#pentane \ "CH3CH2CH2CH2CH3" \
#hexane \ "CH3CH2CH2CH2CH2CH3" \
#heptane \ "CH3CH2CH2CH2CH2CH2CH3" \
#octane \ "CH3CH2CH2CH2CH2CH2CH2CH3" \
#nonane \ "CH3CH2CH2CH2CH2CH2CH2CH2CH3" \
#decane \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#dodecane \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#isoprene \ "#CH3Cd(=CdH2)CdH=CdH2" \
#benzene \ "#mmc1HcHcHcHcHc1H" \
#toluene \ "#mmc1HcHcHcHcHc1CH3" \
#apinene \ "C12HCH2CH(C1(CH3)CH3)CH2CdH=Cd2CH3" \
#bpinene \ "CdH2=Cd1CH2CH2C2HC(CH3)(CH3)C1HC2H2" \
#bmyrcene \ "CH3Cd(CH3)=CdHCH2CH2Cd(=CdH2)CdH=CdH2" \
#limonene \ "CdH2=Cd(CH3)C1HCH2CH2Cd(CH3)=CdHC1H2" \
#ocimene \ "CH3Cd(CH3)=CdHCH2CdH=Cd(CH3)CdH=CdH2" \
#sabinene \ "C1H2Cd(=CdH2)C2HCH2C2(C1H2)CH(CH3)CH3" \

# La ACPD (2016)
#"2-methyloct-1-ene" \ "CdH2=Cd(CH3)CH2CH2CH2CH2CH2CH3" \
#"2-methylbutadec-1-ene" \ "CdH2=Cd(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methylundecane" \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"24-dimethyldecane" \ "CH3CH(CH3)CH2CH(CH3)CH2CH2CH2CH2CH2CH3" \
#"246-trimethylnonane" \ "CH3CH(CH3)CH2CH(CH3)CH2CH(CH3)CH2CH2CH3" \
#"2346-tetramethyloctane" \ "CH3CH(CH3)CH(CH3)CH(CH3)CH2CH(CH3)CH2CH3" \
#"23456-pentamethylheptane" \ "CH3CH(CH3)CH(CH3)CH(CH3)CH(CH3)CH(CH3)CH3" \
#"butadec-6-ene" \ "CH3CH2CH2CH2CH2CH2CdH=CdHCH2CH2CH2CH2CH2CH2CH3" \
#"cyclodecane" \ "C1H2CH2CH2CH2CH2CH2CH2CH2CH2C1H2" \
#"cyclohexane" \ "C1H2CH2CH2CH2CH2C1H2" \
#"cyclopentadecane" \ "C1H2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2C1H2" \
#"dodec-1-ene" \ "CdH2=CdHCH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"heptadecane" \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"heptadec-1-ene" \ "CdH2=CdHCH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"heptadec-8-ene" \ "CH3CH2CH2CH2CH2CH2CH2CdH=CdHCH2CH2CH2CH2CH2CH2CH2CH3" \
#"oct-1-ene" \ "CdH2=CdHCH2CH2CH2CH2CH2CH3" \

# Aumont Faraday (2013)
#"2-methylheptadecane" \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methylnonane" \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"2-methyltridecane" \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methyluncosane" \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"9-10-11-trimethylnonadecane" \ "CH3CH2CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH3" \
#"10-11-dimethyleicosane" \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyldecane" \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyloctadecane" \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyltetradecane" \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"345-trimethylheptane" \ "CH3CH2CH(CH3)CH(CH3)CH(CH3)CH2CH3" \
#"45-dimethyloctane" \ "CH3CH2CH2CH(CH3)CH(CH3)CH2CH2CH3" \
#"567-trimethylundecane" \ "CH3CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH3" \
#"67-dimethyldodecane" \ "CH3CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH3" \
#"789-trimethylpentadecane" \ "CH3CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH3" \
#"89-dimethylhexadecane" \ "CH3CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"docosane"    \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"octadecane"  \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"tetradecane" \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \

# in order of chain length and complexity
#-----------------------------------------
#"propane"           \ "CH3CH2CH3" \
#"butane"            \ "CH3CH2CH2CH3" \
#"isobutane"         \ "CH3CH(CH3)CH3" \
#"pentane"           \ "CH3CH2CH2CH2CH3" \
#"isoprene"          \ "#CH3Cd(=CdH2)CdH=CdH2" \
#"hexane"            \ "CH3CH2CH2CH2CH2CH3" \
#"cyclohexane"       \ "C1H2CH2CH2CH2CH2C1H2" \
#"benzene"           \ "#mmc1HcHcHcHcHc1H" \
#"heptane"           \ "CH3CH2CH2CH2CH2CH2CH3" \
#"toluene"           \ "#mmc1HcHcHcHcHc1CH3" \
#-----------------------------------------
#"octane"               \ "CH3CH2CH2CH2CH2CH2CH2CH3" \
#"oct-1-ene"            \ "CdH2=CdHCH2CH2CH2CH2CH2CH3" \
#"nonane"               \ "CH3CH2CH2CH2CH2CH2CH2CH2CH3" \
#-10----------------------------------------
#"2-methyloct-1-ene"    \ "CdH2=Cd(CH3)CH2CH2CH2CH2CH2CH3" \
#"decane"               \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methylnonane"       \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"45-dimethyloctane"    \ "CH3CH2CH2CH(CH3)CH(CH3)CH2CH2CH3" \
#"345-trimethylheptane" \ "CH3CH2CH(CH3)CH(CH3)CH(CH3)CH2CH3" \
#"cyclodecane"          \ "C1H2CH2CH2CH2CH2CH2CH2CH2CH2C1H2" \
#-10----------------------------------------
#"apinene"           \ "C12HCH2CH(C1(CH3)CH3)CH2CdH=Cd2CH3" \
#"bpinene"           \ "CdH2=Cd1CH2CH2C2HC(CH3)(CH3)C1HC2H2" \
#"bmyrcene"          \ "CH3Cd(CH3)=CdHCH2CH2Cd(=CdH2)CdH=CdH2" \
#"limonene"          \ "CdH2=Cd(CH3)C1HCH2CH2Cd(CH3)=CdHC1H2" \
#"ocimene"           \ "CH3Cd(CH3)=CdHCH2CdH=Cd(CH3)CdH=CdH2" \
#"sabinene"          \ "C1H2Cd(=CdH2)C2HCH2C2(C1H2)CH(CH3)CH3" \
#-12----------------------------------------
#"dodecane"                 \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"dodec-1-ene"              \ "CdH2=CdHCH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methylundecane"         \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"24-dimethyldecane"        \ "CH3CH(CH3)CH2CH(CH3)CH2CH2CH2CH2CH2CH3" \
#"246-trimethylnonane"      \ "CH3CH(CH3)CH2CH(CH3)CH2CH(CH3)CH2CH2CH3" \
#"2346-tetramethyloctane"   \ "CH3CH(CH3)CH(CH3)CH(CH3)CH2CH(CH3)CH2CH3" \
#"23456-pentamethylheptane" \ "CH3CH(CH3)CH(CH3)CH(CH3)CH(CH3)CH(CH3)CH3" \
#-14----------------------------------------
#"tetradecane"              \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methyltridecane"        \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyldecane"         \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"67-dimethyldodecane"      \ "CH3CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH3" \
#"567-trimethylundecane"    \ "CH3CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH3" \
#"butadec-6-ene"            \ "CH3CH2CH2CH2CH2CdH=CdHCH2CH2CH2CH2CH2CH2CH3" \
#-15----------------------------------------
#"cyclopentadecane"         \ "C1H2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2C1H2" \
#"2-methylbutadec-1-ene"    \ "CdH2=Cd(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#-17----------------------------------------
#"heptadecane"              \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"heptadec-1-ene"           \ "CdH2=CdHCH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"heptadec-8-ene"           \ "CH3CH2CH2CH2CH2CH2CH2CdH=CdHCH2CH2CH2CH2CH2CH2CH2CH3" \
#-18----------------------------------------
#"octadecane"               \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methylheptadecane"      \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"89-dimethylhexadecane"    \ "CH3CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyltetradecane"    \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"789-trimethylpentadecane" \ "CH3CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH3" \
#-22----------------------------------------
#"docosane"                    \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"2-methyluncosane"            \ "CH3CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"10-11-dimethyleicosane"      \ "CH3CH2CH2CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"33-diethyloctadecane"        \ "CH3CH2C(CH2CH3)(CH2CH3)CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH2CH3" \
#"9-10-11-trimethylnonadecane" \ "CH3CH2CH2CH2CH2CH2CH2CH2CH(CH3)CH(CH3)CH(CH3)CH2CH2CH2CH2CH2CH2CH2CH3" \
#-------------------------------------------
#-------------------------------------------
