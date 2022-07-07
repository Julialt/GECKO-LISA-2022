#!/bin/bash
#========================================================
#==SCRIPT TO PACKAGE GENERATED MECHANISM FOF BOX MODEL
#========================================================

echo ''
echo '------------------------'
echo In script gen_package_local.bash
echo '------------------------'

#--default value
pathfile=$1
mech=$2

#====================================
#== SET UP GECKO VERSION AND PATHS ==
#====================================

#--possible input arguments : incompatible with above input args
#while getopts "a:m" opt
#do
#  case "$opt" in
#    a ) pathfile=$OPTARG ;;
#    m ) mech=$OPTARG ;;
#  esac
#done

echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`     grep "home_dir"      ${pathfile} | awk '{print $3}' `
scratch_dir=`  grep "scratch_dir"   ${pathfile} | awk '{print $3}' `
gecko_version=`grep "gecko_version" ${pathfile} | awk '{print $3}' `
gecko_inp_dir=`grep "gecko_inp_dir" ${pathfile} | awk '{print $3}' `
gecko_run_dir=`grep "gecko_run_dir" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
gecko_source=${home_dir}/${gecko_version}
gecko_outdir=$scratch_dir/$gecko_run_dir
gecko_wkdir=${gecko_outdir}/${mech}

#===========================================
# RUN AKPARAMETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

echo ''
echo '------------------------'
echo create akparameter file...
echo '------------------------'

cd ${home_dir}/SCRIPTS

./write_akparameter.bash ${mech}
rm ${gecko_wkdir}/akparameter.h

#==========================================
# CREATE GITINFO FILE > $gecko_wkdir
# test for git repository or not
#==========================================
if [ ! -e ${gecko_wkdir}/${mech}.gitinfo ] ; then
 if [ ! -d .git ]; then
  echo 'creating gitinfo file'
## GitHub branch & commit:
  echo `git rev-parse --abbrev-ref HEAD` > ${gecko_wkdir}/${mech}.gitinfo
  echo `git describe` >> ${gecko_wkdir}/${mech}.gitinfo
 else
## IF NOT WITHIN A GIT REPOSITORY READ INFO FROM VERSION_NOTES DIRECTORY
  echo copying gitinfo file
  cp ${gecko_source}/VERSION_NOTES/README_gitinfo ${gecko_wkdir}/${mech}.gitinfo
 fi
fi
ln -sf ${gecko_wkdir}/${mech}.gitinfo ${gecko_wkdir}/indat.gitinfo

#===========================================
# RUN INTERPRETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

echo ''
echo '------------------------'
echo run interpreter...
echo '------------------------'

cd ${home_dir}/SCRIPTS
./run_interp_local.bash -m ${mech}
exit
#===========================================
# tidy up
#===========================================

cd ${gecko_wkdir}

# We should find a way to remove these only
# if we are certain the interpreter was successful
# existence of outdat.nc is not good enough
#rm X*
#rm fort.*
#rm indat.*

echo ''
echo ------------------------------------------
echo delete outstanding temporary files ...
echo ------------------------------------------
echo final directory listing ......
ls 
echo ------------------------------------------

exit


