#!/bin/bash

#===========================================
#==SCRIPT TO GENERATE MECHANISM FROM MULTIPLE PRECURSORS==
# PURPOSE: run generator on cheyenne
# TESTED: CAN BUILD ONTO A PREVIOUSLY-EXISTING MECHANISM ==
# CALLING SYNTAX: 
# ./generate_scheme_cheyenne.bash -m (mechname) -i (cheminput) -s (settingsfile) -x (existing_mech)
# (see below for meaning of other flags)
#===========================================
#--intended to run from scratch in 3 stages:
#--1) 1st precursor
#--2) more precursors, sequentially, with output each time
#--3) post-process: create pvap & Henry files; add W, A reactions
#=========================================================

source cheyenne_scripting_functions.bash

module unload netcdf
module unload intel
module load gnu
module load netcdf

#-- default values for some arguments ==:
pathfile="setup.dat"
settings="settings_default"
cheminput="cheminput.dat"

walltime_compilation="00:05:00"
walltime_gecko="12:00:00"


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
# p = phot file (optional: default is sophie_0404_julia_cmv)
# v = previous output (created with same mechanism)
#== dir paths (optional: obtainable from setup.dat)
# b = box model source path
# h = home directory
# w = gecko working dir
#--
while getopts "a:b:c:d:h:i:k:m:p:s:v:w:x:" opt
do
  case "$opt" in
# general files
    a ) pathfile=$OPTARG ;;
    d ) depend=$OPTARG ;;
    m ) mechname=$OPTARG ;;
# generator-specific files
    i ) cheminput=$OPTARG ;;
    s ) settings=$OPTARG ;;
    x ) existing=$OPTARG ;;
# boxmodel-specific files
    c ) echo "ERROR! -c is a boxmod-specific flag!"; exit ;;
    k ) echo "ERROR! -k is a boxmod-specific flag!"; exit ;;
    p ) echo "ERROR! -p is a boxmod-specific flag!"; exit ;;
    v ) echo "ERROR! -v is a boxmod-specific flag!"; exit ;;
# dir paths, obtainable from setup.dat
#    b ) boxmod_source=$OPTARG ;;
#    h ) home_dir=$OPTARG ;;
#    w ) =$OPTARG ;;
# any other code produces an error
    * ) echo "ERROR! flag not recognized !"; exit ;;
  esac
done

#====================================
#== SET UP GECKO VERSION AND PATHS ==
#====================================

echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
gecko_inp_dir=`  grep "gecko_inp_dir"  ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
gecko_source=$home_dir/$gecko_version
gecko_inputs=$home_dir/$gecko_inp_dir
gecko_outdir=$scratch_dir/$gecko_run_dir
boxmod_source=$home_dir/$boxmod_version

echo "GECKO source = "$gecko_source
echo "input path   = "$gecko_inputs
echo "output path  = "$gecko_outdir

#===========================================
#== SET UP PRECURSORS & MECH NAME ==

echo "existing=" ${existing}
if [ -z ${existing} ] ; then
    echo existing is zero length
    prevflag="n"
else
    echo existing is non zero length
    prevflag="y"
    existing_mech=${existing}
fi
echo "prevflag=" $prevflag

if [ ! -e ${gecko_inputs}/$cheminput ] ; then
    echo "no file of that name: using "$cheminput" as precursor"
    # create file cheminput."precu" in INPUTS directory
    cheminput_file=cheminput.dat.$cheminput
    echo $cheminput > tmp1.txt
    echo END > tmp2.txt
    cat tmp1.txt tmp2.txt > $gecko_inputs/$cheminput_file
    rm tmp*.txt
else
  cheminput_file=$cheminput
  echo $cheminput_file
fi

if [ -z ${mechname}  ] ; then
  echo "no mech name supplied" ; exit 2 ; fi

#===========================================
# DIRECTORY CHECKS AND LINKS
#===========================================

if [ ! -e ${gecko_source} ] ; then
	echo error, gecko source directory could not be found
	echo ${gecko_source}
	echo check gecko_source variable in $0
	exit ; fi

if [ "$prevflag" == "y" ] && [ ! -e ${gecko_outdir}/${existing_mech} ] ; then
        echo existing mechanism \'${gecko_outdir}/${existing_mech}\' doesn\'t exist
        echo ; fi

if [ ! -e ${gecko_outdir} ] ; then
	echo error, gecko run directory could not be found
	echo ${gecko_outdir}
	echo check gecko_outdir variable in $0
	exit ; fi

if [ ! -e ${gecko_inputs}/${settings} ]; then
	echo Settings file \'${gecko_inputs}/${settings}\' doesn\'t exist
	exit ; fi

if [ ! -e ${gecko_inputs}/${cheminput_file} ]; then
	echo Cheminput file \'${gecko_inputs}/${cheminput_file}\' doesn\'t exist
	exit ; fi

if [ ! -e ${gecko_outdir}/DATA ]; then
	echo "creating (updating) link to DATA in "${gecko_outdir}
	ln -s ${gecko_source}/DATA ${gecko_outdir}/DATA
        fi

#===========================================
# RUN THE GENERATOR
#===========================================

#== directory for generated output ==

gecko_wkdir=${gecko_outdir}/${mechname}

if [ "$prevflag" == "y" ] ; then
  gecko_prevdir=${gecko_outdir}/${existing}
  if [ ! -e ${gecko_prevdir} ] ; then
        echo WARNING! previous output dir ${gecko_prevdir}
        echo          does NOT exist. Stopping !!
        exit 99
  elif [ ${gecko_prevdir} == ${gecko_wkdir} ] ;then
        echo adding to previous output in same directory
        echo ${gecko_wkdir}
  fi
elif [ -e ${gecko_wkdir} ] ; then
  echo CAUTION: OVERWRITING DIRECTORY ${gecko_wkdir}
  rm -R ${gecko_wkdir}
fi
mkdir ${gecko_wkdir}

cd ${gecko_wkdir}

echo ------------------------------------------
echo We are now in directory ...... ; pwd
echo ------------------------------------------

#!!! TODO:  compile as part of job launch_generator_cheyenne.bash !!!

echo "submitting job to compile GECKO on cheyenne"
      write_cheyenne_monoproc_script\
        compile_gecko.bash\
        compile_gecko\
        ${walltime_compilation}\
        output_gecko_compilation\
        error_gecko_compilation\
        $home_dir/SCRIPTS/compile_gecko.bash ${gecko_source} ${gecko_wkdir}
    echo
   # submit on hold until all dependent jobs are submitted (-h option)
    compilation_job_id=`qsub -h ${home_dir}/GENERATED_SCRIPTS/compile_gecko.bash`
    echo $compilation_job_id

echo linking executable files...
echo ${gecko_source}/RUN/cm 
#if [ ! -e ${gecko_source}/RUN/cm ] ; then
  cd ${gecko_source}
  make cm
  cd ${gecko_wkdir}
#fi
ln -s ${gecko_source}/RUN/cm ${gecko_wkdir}/cm

echo Preparing input files...
cp ${gecko_inputs}/${settings}  ${gecko_wkdir}/userparams.input
echo  ${gecko_wkdir}/userparams.input
cp ${gecko_inputs}/${cheminput_file}  ${gecko_wkdir}/cheminput.dat
echo  ${gecko_wkdir}/cheminput.dat

#--Preparing cheminput.dat :
#  delete lines after first 'END' from working copy
#  and find # of precursors in file (non-blank lines not starting with *"

nend=$(grep -nm1 ^END cheminput.dat |awk '{ VAR =+ $1} END {print VAR}')
head -n $nend cheminput.dat > tmp.txt ; mv tmp.txt cheminput.dat
numpre=$(grep -c "^[^*E]" cheminput.dat)
echo "File cheminput.dat contains "$numpre" precursors"

echo "prevflag = "$prevflag

if [ "$prevflag" == "y" ] ; then

# dictionary: test for: 
#             1) renamed existing dict
#             2) existing but un-renamed dict (= fort.7)
#             3) crashed dict in working directory (= existing.dict)
# result: file "existing.dict" to be ingested by generator

  echo "looking for ... ${gecko_outdir}/${existing_mech}/userdat.cheminput"
  if [ -e  ${gecko_outdir}/${existing_mech}/userdat.cheminput ] ; then 
    echo " Copying existing cheminput > existing.cheminput"
    cp ${gecko_outdir}/${existing_mech}/userdat.cheminput existing.cheminput
  fi

  echo "looking for ... ${gecko_outdir}/${existing_mech}/${existing_mech}.dict"
  if [ -e  ${gecko_outdir}/${existing_mech}/${existing_mech}.dict ] ; then 
    echo " Copying "${existing_mech}".dict > existing.dict"
    cp ${gecko_outdir}/${existing_mech}/${existing_mech}.dict ${gecko_wkdir}/existing.dict
  else
    echo "looking for ... ${gecko_outdir}/${existing_mech}/fort.7"
    if [ -e ${gecko_outdir}/${existing_mech}/fort.7 ] ; then
      echo " Copying fort.7 > existing.dict"
      cp ${gecko_outdir}/${existing_mech}/fort.7 ${gecko_wkdir}/existing.dict
    else 
      echo "looking for ... ${gecko_outdir}/${existing_mech}/existing.dict"
      if [ -e ${gecko_outdir}/${existing_mech}/existing.dict ] ; then
        echo " Copying  existing.dict > existing.dict"
        cp ${gecko_outdir}/${existing_mech}/existing.dict ${gecko_wkdir}/existing.dict
      else
        echo "cannot find an existing dictionary"
        exit
      fi
    fi
  fi

# mechanism: use existing.mech if available, otherwise construct a new one.
# result: file "existing.mech" to be ingested by generator
  echo "looking for ... "${gecko_outdir}/${existing_mech}"/existing.rxns"
  if [ -e  ${gecko_outdir}/${existing_mech}/existing.rxns ] ; then

    echo " Copying existing.rxns > existing.rxns"
    cp ${gecko_outdir}/${existing_mech}/existing.rxns ${gecko_wkdir}/existing.rxns
  else

    if [ -e  ${gecko_outdir}/${existing_mech}/${existing_mech}.mech ] ; then
      echo " constructing file existing.rxns from existing mechanism"

#   "nend1" = location of first 'END' label
      nend1=$(grep -nm1 END ${gecko_outdir}/${existing_mech}/${existing_mech}.mech |awk '{ VAR =+ $1} END {print VAR}')

#   "nend2" = # of lines in mech file (location of final 'END' label)
      nend2=$(wc -l ${gecko_outdir}/${existing_mech}/${existing_mech}.mech |awk '{ VAR =+ $1} END {print VAR}')

#   "nlin" = # of lines to keep (deleting 3 x header lines, 1 x final 'END') line
      let nlin=$nend2-$nend1-3-1
      echo $nend2" - "$nend1" = "$nlin 

#   write previous.mech to tmp.rxns, excluding final line (= 'END')
      sed '$d' ${gecko_outdir}/${existing_mech}/${existing_mech}.mech > ${gecko_wkdir}/tmp.rxns

#   write tmp.rxns to existing.rxns, excluding top list to first 'END'
      tail -n $nlin ${gecko_wkdir}/tmp.rxns > ${gecko_wkdir}/existing.rxns

#   delete tmp.rxns
      rm ${gecko_wkdir}/tmp.rxns

#   "nend3" = # of unwanted extra 'END' labels (stop execution if > 0)
      nend3=$(grep -c 'END' ${gecko_wkdir}/existing.rxns |awk '{ VAR =+ $1} END {print VAR}')
      echo $nend3" extra END labels"
      if [ $nend3 != 0 ] ; then
        echo "what to do about END statements in file ${gecko_wkdir}/existing.rxns?"
        exit
      fi

    else
      echo "cannot find an existing reaction list"
      exit
    fi
  fi

# log file: test for:
#           1) renamed existing log file
#           2) existing but un-renamed log file
#           3) crashed log file in working directory (= existing.log, concatenate to )
# result: file "existing.log" to be concatenated to scheme.log by post-processor

  echo " Copying scheme.log > existing.log"
  if [ -e  ${gecko_outdir}/${existing_mech}/${existing_mech}.log ] ; then
    cp ${gecko_outdir}/${existing_mech}/${existing_mech}.log ${gecko_wkdir}/existing.log
  else
    if [ ! -e  ${gecko_wkdir}/existing.log ] ; then
      cp ${gecko_outdir}/${existing_mech}/scheme.log ${gecko_wkdir}/existing.log
    else
      cat ${gecko_outdir}/${existing_mech}/scheme.log ${gecko_wkdir}/existing.log > ${gecko_wkdir}/tmp.log
      mv ${gecko_wkdir}/tmp.log ${gecko_wkdir}/existing.log
    fi
  fi

fi # prevflag = "y"

#==LOOP RUNS OVER PRECURSORS

counter=0
until [ $counter -eq $numpre ]; do
  let counter=$counter+1
  #TODO: here we can tune walltime as a function of the precursor
  #walltime="01:00"
#== define wallclock time
#  let cnum=$( sed -n $counter','$counter'p' ${gecko_wkdir}/cheminput.dat | awk -F\C '{print NF-1}')
#  if [ $cnum -lt 10 ] ; then
#    walltime="00:30"
#  else
#    walltime="03:00"
#  fi
  if [ $counter -eq 1 ] ; then
    if [ -z "$depend" ] ; then
      write_cheyenne_monoproc_dependentscript gen_${mechname}_${counter}_${numpre}.bash \
        gen_${mechname}_${counter}_${numpre}\
        ${walltime_gecko}\
        output_gen_${mechname}_${counter}_${numpre}\
        error_gen_${mechname}_${counter}_${numpre}\
        afterok:$compilation_job_id \
        $home_dir/SCRIPTS/launch_generator_cheyenne.bash $counter $numpre ${mechname} ${prevflag} ${gecko_source} ${gecko_wkdir} 
    else
# forced NO COMPILE
      write_cheyenne_monoproc_script gen_${mechname}_${counter}_${numpre}.bash \
        gen_${mechname}_${counter}_${numpre}\
        ${walltime_gecko}\
        output_gen_${mechname}_${counter}_${numpre}\
        error_gen_${mechname}_${counter}_${numpre}\
        $home_dir/SCRIPTS/launch_generator_cheyenne.bash $counter $numpre ${mechname} ${prevflag}
      #write_cheyenne_monoproc_dependentscript gen_${mechname}_${counter}_${numpre}.bash \
      #  gen_${mechname}_${counter}_${numpre}\
      #  ${walltime_gecko}\
      #  output_gen_${mechname}_${counter}_${numpre}\
      #  error_gen_${mechname}_${counter}_${numpre}\
      #  afterok:${depend}\
      #  $home_dir/SCRIPTS/launch_generator_cheyenne.bash $counter $numpre ${mechname} ${prevflag}
    fi
  else
    let prevcounter=$counter-1
    write_cheyenne_monoproc_dependentscript gen_${mechname}_${counter}_${numpre}.bash \
      gen_${mechname}_${counter}_${numpre}\
      ${walltime_gecko}\
      output_gen_${mechname}_${counter}_${numpre}\
      error_gen_${mechname}_${counter}_${numpre}\
      afterok:$previous_job_id \
      $home_dir/SCRIPTS/launch_generator_cheyenne.bash $counter $numpre ${mechname} ${prevflag}
  fi

  echo submitting $counter:$numpre
  previous_job_id=`qsub $home_dir/GENERATED_SCRIPTS/gen_${mechname}_${counter}_${numpre}.bash`
  echo --------------------------------------

#-- end run loop
done

#===========================================
# CREATE & COLLECT X-FILES & INTERP INPUT
#===========================================

walltime="01:00:00"

write_cheyenne_monoproc_dependentscript pgn_${mechname}.bash \
  pgn_${mechname}\
  $walltime\
  output_postgen_${mechname}\
  error_postgen_${mechname}\
  afterok:$previous_job_id \
  $home_dir/SCRIPTS/gen_package_cheyenne.bash $gecko_source $gecko_wkdir $mechname $boxmod_source

echo submitting postprocess
last_job_id=`qsub $home_dir/GENERATED_SCRIPTS/pgn_${mechname}.bash`
echo -----------------------------------

qrls $compilation_job_id

echo last_job_id=$last_job_id

exit
