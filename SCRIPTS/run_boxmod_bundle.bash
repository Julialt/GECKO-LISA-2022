#!/bin/bash
#===========================================
# PURPOSE: run boxmodel on cheyenne
# CALLING SYNTAX: ./run_boxmod_bundle.bash -m (mech) -k (keyfile) -v (previous run) -p (photfile) -f (flags_input (postproc))
# (see below for meaning of other flags)
# TESTED: 2020-12-23: JMLT, NCAR
#===========================================

source cheyenne_scripting_functions.bash

#-- default values for some arguments ==:
pathfile="setup.dat"
photfile="sophie0404_julia_cmv.phot"
inpdir="INPUTS"
flags_input="postproc_flags.input"

walltime_compilation=("00:05:00")
walltime_spin=("0:10:00")
walltime_box=("12:00:00")
walltime_post=("3:00:00")

##### Default values for optional arguments
## are we using NetCDF? (user hardwire "yes" or "")
nthreads=16
netcdf_flag="yes"
runlength=3600  # each run will be "runlength" seconds
#runlength=259200  # 3 days # each run will be "runlength" seconds
#runlength=345600  # 4 days # each run will be "runlength" seconds
#runlength=518400  # 6 days # each run will be "runlength" seconds

#===========================================
#== INPUT ARGUMENTS ==
#===========================================
#== general files
# a = file containing source paths (optional: default is setup.dat)
# d = dependencies (optional)
# m = mechanism name (output from generator, input for box model)
#== generator-specific files
# i = cheminput file
# s = settings file for generator
# x = existing mechanism
#== boxmodel-specific files
# c = constrained concentrations file
# k = keyfile (box model input scenario)
# p = phot file (optional: default is sophie0404_julia_cmv)
# v = previous output (created with same mechanism)
# f = file containing postprocessing flags
#== dir paths (optional: obtainable from setup.dat)
# b = box model source path
# h = home directory
# u = directory path for keyfile ("INPUTS/__") default = "INPUTS"
# w = gecko working dir
#== flags
# l = gecko library run (specific rules and input files)
# n = netcdf flag
#== run parameters
# r = what is the duration of individual sub-runs
# t = number of threads to use, default is 16
#--
while getopts "a:b:c:d:f:h:i:k:l:m:n:p:r:s:t:u:v:w:x:z" opt
do
  case "$opt" in
# general files
    a ) pathfile=$OPTARG ;;
    d ) depend=$OPTARG ;;
    m ) mechname=$OPTARG ;;
# generator-specific files
    i ) echo "ERROR! -i is a generator-specific flag!"; exit ;;
    s ) echo "ERROR! -s is a generator-specific flag!"; exit ;;
    x ) echo "ERROR! -x is a generator-specific flag!"; exit ;;
# boxmodel-specific files
    c ) confile=$OPTARG ;;
    k ) keyfile=$OPTARG ;;
    p ) photfile=$OPTARG ;;
    v ) prevout=$OPTARG ;;
    u ) inpdir=$OPTARG ;;
# flags
    l ) library_flag=$OPTARG ;;
    n ) netcdf_flag=$OPTARG ;; 
    f ) flags_input=$OPTARG ;; 
# parameters
    r ) runlength=$OPTARG ;;
    t ) nthreads=$OPTARG ;;
# dir paths, obtainable from setup.dat
#    b ) boxmod_source=$OPTARG ;;
#    h ) home_dir=$OPTARG ;;
#    w ) =$OPTARG ;;
# any other code produces an error
    * ) echo "ERROR! flag not recognized !"; exit ;;
  esac
done

#===========================================
#== SET UP GECKO VERSION AND PATHS==
#===========================================
echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
gecko_data=$home_dir"/"$gecko_version"/DATA"
gecko_outdir=$scratch_dir/$gecko_run_dir
boxmod_source=$home_dir/$boxmod_version
boxmod_outdir=$scratch_dir/$boxmod_run_dir
boxmod_inpdir=$boxmod_source/$inpdir

echo "boxmod source path = "$boxmod_source
echo "boxmod inputs path = "$boxmod_inpdir
echo "boxmod output path = "$boxmod_outdir

#---------check prefix/extension of input data file names -------

# >>keyfile name is used in the script as "indat_${keyfile}.key"
# remove "indat_" prefix from $keyfile if supplied
prefix=${keyfile: 0: 6}
if [ "$prefix" == "indat_" ] ; then
  keyfile=${keyfile##${prefix}}
fi
# remove ".key" extension from $keyfile if supplied
extn=${keyfile: -4}
if [ "$extn" == ".key" ] ; then
  keyfile=${keyfile%${extn}}
fi

# >>confile: if a filename is supplied:
if [ ${#confile} -ne 0 ] ; then
# confile name is used in the script as "init_conc_${confile}.key"
# remove "init_conc_" prefix from $confile if supplied
  prefix=${confile: 0: 10}
  if [ "$prefix" == "init_conc_" ] ; then
    confile=${confile##${prefix}}
  fi
# remove ".key" extension from $confile if supplied
  extn=${confile: -4}
  if [ "$extn" == ".key" ] ; then
    confile=${confile%${extn}}
  fi
fi

# >>photfile name is used in the script as "${photfile}"
# add ".phot" extension to photfile if not already supplied
extn=${photfile: -5}
if [ "$extn" != ".phot" ] ; then
  photfile=$photfile.phot
fi

#---------file existence checks -------
#-- check for boxmod directories
if [ ! -e ${boxmod_source} ] ; then
    echo error, boxmod source directory could not be found
    echo ${boxmod_source}
    exit 1
fi
if [ ! -e ${boxmod_inpdir} ] ; then
    echo error, boxmod inputs directory could not be found
    echo ${boxmod_inpdir}
    exit 1
fi
if [ ! -e ${boxmod_outdir} ] ; then
    echo error, boxmod run directory could not be found
    echo ${boxmod_outdir}
    exit 1
fi

#-- check for input key files (argument, existence)
if [ ${#keyfile} -eq 0 ] ; then
  echo "no input (key) file name supplied" ; exit 2
fi
if [ ! -e ${boxmod_inpdir}/indat_${keyfile}.key ]; then
    echo Input file \'${boxmod_inpdir}/indat_${keyfile}.key\' doesn\'t exist
    exit
fi

#-- check for init_conc files (argument, existence)
if [ ${#confile} -ne 0 ] ; then
  echo "constrained concs file   = "init_conc_$confile.key
  if [ ! -e ${boxmod_inpdir}/init_conc_${confile}.key ]; then
    echo Initial concentrations file \'${boxmod_inpdir}/init_conc_${confile}.key\' doesn\'t exist
    exit
  fi
fi

#-- check for photolysis lookup file existence
if [ ! -e $boxmod_inpdir/$photfile ] ; then
  echo $boxmod_inpdir/$photfile
  echo "photolysis file not found" ; exit 4
fi
echo "photfile = "$photfile

#-- check for mechanism file (argument, existence)
if [ ${#mechname} -eq 0 ] ; then
  echo "no mechanism name supplied" ; exit 3
fi
echo "mechname = "$mechname.mech
gecko_mechdir=$gecko_outdir/$mechname
if [ ! -e ${gecko_mechdir} ] ; then
    echo error, mechanism directory could not be found
    echo ${gecko_mechdir}
    exit 1
fi
echo "mechanism location = "$gecko_mechdir

#-- check for interpreted mech file

if [ -n "${netcdf_flag}" ]; then
  if [ ! -e ${gecko_mechdir}/outdat.nc ]; then
      echo Interpreted mechanism file \'${gecko_mechdir}/outdat.nc\' doesn\'t exist
      exit
  fi
else
  if [ ! -e ${gecko_mechdir}/outdat.akli ]; then
      echo Interpreted mechanism file \'${gecko_mechdir}/outdat.akli\' doesn\'t exist
      exit
  fi
fi

  #---------------- make some links -------

#if [ ${#confile} -eq 0 ] ; then
if [ -z ${confile} ] ; then
  run_inp=${keyfile}
else
  run_inp=${keyfile}_${confile}
fi
run_name=${mechname}_${run_inp}
boxmod_wkdir=${boxmod_outdir}/${run_name}

#==============================================================#
## invoke sequence of cheyenne scripts to run box model:
## write_cheyenne_script
## launch_box_bundle
# can supply up to 4 arguments to the launch script

# look at indat.key to check if we need to decompose into several runs (1 per day)
global_tstart=`grep ^TSTR ${boxmod_inpdir}/indat_${keyfile}.key | cut -d" " -f2`
global_tstop=`grep ^TSTP ${boxmod_inpdir}/indat_${keyfile}.key | cut -d" " -f2`

numruns="$(((${global_tstop}-${global_tstart})/${runlength}+1))"
if [[ "$(((${numruns}-1)*${runlength}))" -eq "$((${global_tstop}-${global_tstart}))" ]] ; then
  let numruns=${numruns}-1
fi

echo ... writing cheyenne script

# loop over multiple runs if needed
counter=0
until [[ ${counter} -eq ${numruns} ]]; do
  let prevcounter=${counter}
  let counter=${counter}+1
# prevout is only used by the first submitted run
# reinitialise to nothing for the following runs
  if [ ${counter} -gt 1 ] ; then
      prevout=""
  fi


  if [ -e ${boxmod_wkdir} ] ; then
    if [[ -z $prevout  && ${counter} -eq 1 ]] ; then   # i.e. $prevout has zero length
      echo Creating ${boxmod_wkdir} ...
      echo WARNING! removing existing directory
      rm -rf ${boxmod_wkdir}
      mkdir ${boxmod_wkdir}
    else
      echo Appending to existing output in dir ${boxmod_wkdir}
    fi
  else
    echo Creating ${boxmod_wkdir} ...
    mkdir ${boxmod_wkdir}
  fi

##-- making links for correct input.output format ---

  if [ -n "${netcdf_flag}" ]; then
    echo linking netcdf link-file ${gecko_mechdir}/outdat.nc
    echo ... to ${boxmod_wkdir}/indat.nc
    ln -s ${gecko_mechdir}/outdat.nc ${boxmod_wkdir}/indat.nc
    echo
  else
    echo linking binary link-file ${gecko_mechdir}/outdat.akli
    echo ... to ${boxmod_wkdir}/indat.li
    ln -s ${gecko_mechdir}/outdat.akli ${boxmod_wkdir}/indat.li
    echo
  fi
 
  #--link some input files 
  if [ -e ${boxmod_inpdir}/*_${keyfile}.input ] ; then
    echo linking user external constraints input files
    for f in `ls ${boxmod_inpdir}/*_${keyfile}.input`
    do
      echo linking $f...
      ln -s $f ${boxmod_wkdir}/
    done
    echo
  fi

  #--link input files
  echo linking photolysis file ${boxmod_inpdir}/${photfile}
  echo to ... ${boxmod_wkdir}/jfile.phot
  ln -s ${boxmod_inpdir}/${photfile} ${boxmod_wkdir}/jfile.phot
  echo

  ##--copy surface properties files 
  #echo copying surface properties files ${boxmod_inpdir}/*.sur
  #echo to ... ${boxmod_wkdir}/*.sur
  #cp ${boxmod_inpdir}/*.sur ${boxmod_wkdir}
  #echo

  echo copying key file ${boxmod_inptdir}/indat_${keyfile}.key...
  echo to ... ${boxmod_wkdir}/indat.key
  cp ${boxmod_inpdir}/indat_${keyfile}.key ${boxmod_wkdir}/indat.key
  echo

#==============================================================#
# -------------- logfile output ---------------
#==============================================================#

# we don't want to change directory in this main script
# what's happening inside the parentheses has no impact outside of them

(cd ${home_dir}
 if [ ${counter} -eq 1 ]; then
 echo '================================'  > $boxmod_wkdir/outdat.README
 echo '= GECKO-A box model simulation =' >> $boxmod_wkdir/outdat.README
 echo '================================' >> $boxmod_wkdir/outdat.README
 echo 'run name      = '${run_name}      >> $boxmod_wkdir/outdat.README
 echo 'user          = '$USER            >> $boxmod_wkdir/outdat.README
 echo 'date          = '`date`           >> $boxmod_wkdir/outdat.README
 echo 'machine       = '$HOSTNAME        >> $boxmod_wkdir/outdat.README
 echo 'model version :'                  >> $boxmod_wkdir/outdat.README
 echo 'GitHub branch = '`git rev-parse --abbrev-ref HEAD` >> $boxmod_wkdir/outdat.README
 echo 'GitHub commit = '`git describe`   >> $boxmod_wkdir/outdat.README
 echo 'model inputs  = indat_'$keyfile'.key' >> $boxmod_wkdir/outdat.README
 if [ ${#confile} -ne 0 ] ; then
 echo 'constraints   = '$confile         >> $boxmod_wkdir/outdat.README
 else
 echo 'no additional constraints file'   >> $boxmod_wkdir/outdat.README
 fi
 echo 'dictionary    = '$mechname        >> $boxmod_wkdir/outdat.README
 echo 'phot file     = '$photfile        >> $boxmod_wkdir/outdat.README
 if [ -n "${prevout}" ]; then
 echo 'continuing from existing output : '$prevout >> $boxmod_wkdir/outdat.README
 else
 echo 'new run, no previous output '     >> $boxmod_wkdir/outdat.README
 fi
 echo '================================' >> $boxmod_wkdir/outdat.README
 echo '(add user comments manually below)'>> $boxmod_wkdir/outdat.README

 fi)
#==============================================================#
# Create reference file git_info.txt in boxmod working directory, to be read by codes.
#==============================================================#
  echo creating gitinfo file ${boxmod_outdir}/indat.gitinfo

(cd ${home_dir}
 echo `git rev-parse --abbrev-ref HEAD` > $boxmod_wkdir/indat.gitinfo 
 echo `git describe` >> $boxmod_wkdir/indat.gitinfo )

#==============================================================#
# check if END statement already exists in indat.key and remove it
# then remove (all) blank lines from file

  sed -i s/^END//g ${boxmod_wkdir}/indat.key
  sed -i '/^$/d' ${boxmod_wkdir}/indat.key

  let starttime="$((${global_tstart}+(${counter}-1)*${runlength}))"
  echo starttime=$starttime

  if [[ ${counter} -gt 1 ]] ; then

## launch script sets previous run directory if we're in a subrun >1
## for subrun = 1, it will be created if needed later
    #boxmod_prevdir=$boxmod_outdir/${run_name}_${prevcounter}_${numruns}
    #boxmod_prevdir=$boxmod_outdir/${run_name}
    #echo "previous output in = "$boxmod_prevdir
    
# if PREV present, make sure is = 1
    if grep -q "PREV" ${boxmod_wkdir}/indat.key ; then
      sed -i 's/PREV 0/PREV 1/g' ${boxmod_wkdir}/indat.key
# if PREV missing, insert flag at beginning of indat file
    else
      sed -i '1 i\PREV 1' ${boxmod_wkdir}/indat.key
    fi

# should consider test depending on timestep length
    sed -i "s/^TSTR [0-9]\+/TSTR ${starttime}/" ${boxmod_wkdir}/indat.key
  fi

  let stoptime="$((${global_tstart}+(${counter})*${runlength}))"
  if [ $stoptime -gt $global_tstop ]; then
    echo stoptime=$global_tstop
  else
    echo stoptime=$stoptime
  fi 

  sed -i "s/^TSTP [0-9]\+/TSTP ${stoptime}/" ${boxmod_wkdir}/indat.key

  if [[ ${counter} -eq ${numruns} ]]; then
    sed -i "s/^TSTP [0-9]\+/TSTP ${global_tstop}/" ${boxmod_wkdir}/indat.key
  fi

  if [[ ${counter} -eq 1 ]]; then
    if [ -n "${prevout}" ]; then
      boxmod_prevdir=$boxmod_outdir/$prevout
      echo "previous output in = "$boxmod_prevdir

  # first run: if ${prevout} provided, set PREV = 1
      if grep -q "PREV" ${boxmod_wkdir}/indat.key ; then
        sed -i 's/PREV 0/PREV 1/g' ${boxmod_wkdir}/indat.key
      else
        sed -i '1 i\PREV 1' ${boxmod_wkdir}/indat.key
      fi

    else # if no ${prevout}, default in model is PREV 0
# replace any existing "PREV 1" line (for documentation purposes)
      sed -i 's/PREV 1/PREV 0/g' ${boxmod_wkdir}/indat.key

    fi
  fi

  if [[ -n "${prevout}" || ${counter} -gt 1 ]]; then
  #if [[ -n "${prevout}" || ${counter} -eq 1 ]]; then
    if [ ! -e ${boxmod_prevdir} ] ; then
      echo error, previous output directory could not be found
      echo ${boxmod_prevdir}
      exit 1
    fi
    #-- check for previous output file outdat.nc
    if [ -n "${prevout}" ]; then
      if [ ! -e ${boxmod_prevdir}/outdat.nc ]; then
        echo previous output file \'${boxmod_prevdir}/outdat.nc\' doesn\'t exist
        exit
      fi
    fi
    
    echo linking existing output ${boxmod_prevdir}/outdat.nc
    echo ... to ${boxmod_wkdir}/prevdat.nc
    ln -sf ${boxmod_prevdir}/outdat.nc ${boxmod_wkdir}/prevdat.nc
    echo

    # ln BOXMOD exec from previous dir
    # if the compilation has not finished yet
    # the link will point to nothing until it is done
    # it is fine since the boxmod runs start only after
    # compilation is done
    # and there is a check at that stage to check that the executable exists
# do not do here since we check this in launch script, and compile first time around
    #if [ ${boxmod_prevdir} -ne ${boxmod_wkdir} ]; then
    #  if [ ${counter} -gt 1 ]; then
    #if [ ${counter} -eq 1 ]; then
    #    echo copying executable ${boxmod_prevdir}/BOXMOD...
    #    echo ... to ${boxmod_wkdir}/BOXMOD
    #    ln -s ${boxmod_prevdir}/BOXMOD ${boxmod_wkdir}/BOXMOD
    #    echo
    #  fi
    #fi
    
  fi

# these initial concentrations will be overwritten by previous run if needs,
# but we need this file for emissions


#==============================================================#
## section for auto-adding precursor concs to keyfile (LIBRARY)
#==============================================================#
  if [ -n "${library_flag}" ] ; then
    echo Finding primary species codes...
    if [ -e ${gecko_mechdir}/${mechname}.prec ]; then
      ln -sf ${gecko_mechdir}/${mechname}.prec ${gecko_mechdir}/findname_output
    fi
    primary=`cat ${gecko_mechdir}/findname_output | cut -f 1 -d":"`
    primary=`echo ${primary} | xargs`  # trick to trim variable
    nc=`cat ${gecko_mechdir}/findname_output | cut -f 3 -d":"`
    if [ $nc -eq 0 ] ; then
      initconc=0
    else
      initconc=`echo 25000000/$nc | bc`
    fi

    echo "REAC G${primary}  ${initconc}  0.00E+00  0.00E+00" >> ${boxmod_wkdir}/indat.key
    echo ${primary}, ${initconc}
    echo added to ${boxmod_workingdir}/indat.key
  fi

  if [ -e ${boxmod_inpdir}/init_conc_${confile}.key ]; then
    echo adding initial concentrations ${boxmod_inpdir}/init_conc_${confile}.key...
    echo to ... ${boxmod_wkdir}/indat.key
    cp ${boxmod_inpdir}/init_conc_${confile}.key ${boxmod_wkdir}/init_conc.key
    cat ${boxmod_wkdir}/init_conc.key >> ${boxmod_wkdir}/indat.key
    echo
  fi

#==============================================================#
## section for auto-adding spinup concs to keyfile (LIBRARY)
#==============================================================#
  if [ -n "${library_flag}" ]; then
    if [ ! -e ${boxmod_outdir}/spinup_spinup_${keyfile}/USEROUT/steadystate.key ]; then
      echo no steadystate.key file found to initialize library run in 
      echo ${boxmod_outdir}/spinup_spinup_${keyfile}
      exit
    fi
    echo adding steady-state spinup...
    echo ${boxmod_outdir}/spinup_spinupi${keyfile}/USEROUT/steadystate.key...
    echo to ... ${boxmod_wkdir}/indat.key
      cp ${boxmod_outdir}/spinup_spinup_${keyfile}/USEROUT/steadystate.key ${boxmod_wkdir}/init_steadystate.key
      cat ${boxmod_wkdir}/init_steadystate.key >> ${boxmod_wkdir}/indat.key
    echo   
  fi

#==============================================================#
## make sure that indat.key has an "END" line
#==============================================================#
  echo END >> ${boxmod_wkdir}/indat.key

  cp ${boxmod_wkdir}/indat.key ${boxmod_wkdir}/indat.key${counter}

#==============================================================#
# generate scripts
#==============================================================#
# no need for himem scripts, boxmod runs rarely take more that 2000MB memory (tested for 800000 species)
  if [ $stoptime -gt $global_tstop ]; then
    echo ${counter} tstart=$starttime, tstop=$global_tstop
  else
    echo ${counter} tstart=$starttime, tstop=$stoptime
  fi
  if [[ ${counter} -eq 1 ]]; then
    write_cheyenne_script\
        bun_${run_name}_${counter}_${numruns}.bash\
        bun_${run_name}\
        ${walltime_box}\
        ${nthreads}\
        output_bun_${run_name}_${counter}_${numruns}\
        error_bun_${run_name}_${counter}_${numruns}\
        ${home_dir}/SCRIPTS/launch_box_bundle.bash\
        ${home_dir}/SCRIPTS/${pathfile}\
        ${mechname}\
        ${run_inp}\
        ${counter} \
        ${flags_input}
  else
    write_cheyenne_dependentscript\
        bun_${run_name}_${counter}_${numruns}.bash\
        bun_${run_name}\
        ${walltime_box}\
        ${nthreads}\
        output_bun_${run_name}_${counter}_${numruns}\
        error_bun_${run_name}_${counter}_${numruns}\
        afterany:${previous_job_id} \
        ${home_dir}/SCRIPTS/launch_box_bundle.bash\
        ${home_dir}/SCRIPTS/${pathfile}\
        ${mechname}\
        ${run_inp}\
        ${counter} \
        ${flags_input}
  fi

  echo ... submitting cheyenne script
  current_wkdir=`pwd`
  cd ${boxmod_wkdir} 
  previous_job_id=`qsub ${home_dir}/GENERATED_SCRIPTS/bun_${run_name}_${counter}_${numruns}.bash`
  cd ${current_wkdir}

done

# print last boxmod job if needed by other scripts
echo last_job_id=$previous_job_id

