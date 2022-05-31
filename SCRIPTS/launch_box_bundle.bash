#!/bin/bash
## UNDER DEVELOPMENT for some arguments ==:

# this script must be launched from the $boxmod_wkdir folder
pathfile="setup.dat"

pathfile=$1
mechname=$2
run_inp=$3
counter=$4
flags_input=$5

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" ${pathfile} | awk '{print $3}' `

gecko_outdir=$scratch_dir/$gecko_run_dir
gecko_mechdir=$gecko_outdir/$mechname
boxmod_source=$home_dir/$boxmod_version
boxmod_outdir=$scratch_dir/$boxmod_run_dir

run_name=${mechname}_${run_inp}
boxmod_wkdir=${boxmod_outdir}/${run_name}

## set previous run directory to working dir for subruns >1
## so that we append results to existing file.
if [[ ${counter} -gt 1 ]] ; then
  boxmod_prevdir=$boxmod_wkdir
  echo linking existing output ${boxmod_prevdir}/outdat.nc
  echo ... to ${boxmod_wkdir}/prevdat.nc
  ln -sf ${boxmod_prevdir}/outdat.nc ${boxmod_wkdir}/prevdat.nc
  echo
fi

cd $boxmod_wkdir

cp indat.key${counter} indat.key
#ln -s ${gecko_mechdir}/outdat.nc ${boxmod_wkdir}/indat.nc

#----------------------------------#
# CHECK IF compiled codes exist
# YES -> no action required
# NO -> compile the box model
# re: compile_boxmod.bash
#----------------------------------#

if [ ! -e ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ] ; then

  echo "boxmodel executable -> compiling"
  rm -rf ${boxmod_wkdir}/WORK_BOXMOD
  mkdir ${boxmod_wkdir}/WORK_BOXMOD
  mkdir ${boxmod_wkdir}/WORK_BOXMOD/LIB
  echo

  echo linking...
  for f in `ls ${boxmod_source}/LIB/*.[fh]*`
  do
    cp $f ${boxmod_wkdir}/WORK_BOXMOD/LIB/
  done

  mkdir ${boxmod_wkdir}/WORK_BOXMOD/OBJ
  cp ${boxmod_source}/OBJ/Makefile ${boxmod_source}/OBJ/template.mk ${boxmod_wkdir}/WORK_BOXMOD/OBJ/
  mkdir ${boxmod_wkdir}/WORK_BOXMOD/PROG
  cp ${boxmod_source}/PROG/boxmod_main.f90 ${boxmod_wkdir}/WORK_BOXMOD/PROG/boxmod_main.f90

# link to the relevant akparameter.h file
  echo akparameter...
  cp ${gecko_mechdir}/akparameter_module.f90 ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter_module.f90

# we don't want to change directory in this main script
# what's happening inside the parentheses has no impact outside of them
  echo Making...
  ( cd ${boxmod_wkdir}/WORK_BOXMOD/OBJ && exec make ../PROG/BOXMOD )

  # final check
  if [ ! -e ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ] ; then

  # if no executable, generate error notice and quit
    echo error, could not find boxmodel executable
    echo ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD
    echo Check compilation!
    exit 1

  else
    cp ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ${boxmod_wkdir}/BOXMOD
  #dereference link
  #  cp ./BOXMOD ./BOXMOD_tmp
  #  mv ./BOXMOD_tmp ./BOXMOD

  # run boxmod
    export OMP_NUM_THREADS=$NCPUS
    ./BOXMOD > $TMPDIR/stdout.$PBS_JOBID
  fi

else # (previous version of BOXMOD does exist)
  cp ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ${boxmod_wkdir}/BOXMOD
#dereference link
  #cp ./BOXMOD ./BOXMOD_tmp
  #mv ./BOXMOD_tmp ./BOXMOD

# run boxmod
  export OMP_NUM_THREADS=$NCPUS
  ./BOXMOD > $TMPDIR/stdout.$PBS_JOBID

fi

#----------------------------------#
# check for results files

if [ -s ${boxmod_wkdir}/outdat.nc ] ; then
  echo box model has run successfully
  echo check output ${boxmod_wkdir}/outdat.nc
fi

#* tidy up
cat outdat.README outdat.out fort.11 > tmp.txt
mv tmp.txt outdat.README

rm fort.*
rm dummy.*

#rm outdat.out

## delete local copy of executable
rm BOXMOD

## delete work directory 
## comment this out if you want to avoid re-compiling every iteration
#rm -R WORK_BOXMOD

## delete symbolic links
#find -lname '*' -delete

## delete zero-length files
find ./ -type f -size 0 -delete

## delete temporary files
#rm *.gitinfo

#===========================================
# RUN POST_PROCESSOR SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================
echo ../running postprocessing script
  #qsub  ${home_dir}/SCRIPTS/run_postproc_cheyenne.bash ${mechname} ${run_inp} ${flags_input}
cd ${home_dir}/SCRIPTS
./run_postproc_cheyenne.bash ${mechname} ${run_inp} ${flags_input}



