#!/bin/bash
#=================================================================
#==SCRIPT TO PACKAGE TRANSLATED MECHGEN INPUT FOR GECKO BOX MODEL
#=================================================================

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
# CREATE & COLLECT X-FILES & INTERP INPUT
#===========================================

echo ''
echo '------------------------'
echo Linking input files
echo And copying dummy files
echo '------------------------'

if [ ! -e ${gecko_wkdir}/fort.7 ] ; then
  ln -s ${gecko_wkdir}/${mech}'.dict' ${gecko_wkdir}/fort.7 ; fi

#-----------------------------
if [ ! -d .git ]; then
# check if indat.gitinfo already exists and is non-zero
  if [ ! -s indat.gitinfo ]; then
# Create reference file git_info.txt in output directory, to be read by codes.
# GitHub branch, commit:
    echo "** creating gitinfo file **"
    echo `git rev-parse --abbrev-ref HEAD` > ${gecko_wkdir}/indat.gitinfo
    echo `git describe`   >> ${gecko_wkdir}/indat.gitinfo
  else 
    echo "** using existing gitinfo file **"
  fi
else
# IF NOT WITHIN A GIT REPOSITORY READ INFO FROM VERSION_NOTES DIRECTORY
  echo "** copying gitinfo file **"
  cp ${home_dir}/VERSION_NOTES/README_gitinfo ${gecko_wkdir}/indat.gitinfo
fi

cp ${gecko_wkdir}/indat.gitinfo ${gecko_wkdir}/${mech}.gitinfo
#-----------------------------

cp ${home_dir}/GECKO-A/INPUTS/settings_default ${gecko_wkdir}/userdat.settings

echo ''
echo '------------------------'
echo Finding Precursor Names
echo '------------------------'

cd ${gecko_source}
#make findname_cmv

ln -s ${gecko_source}/RUN/findname_cmv ${gecko_wkdir}/findname_cmv

let nlin=`wc -l ${gecko_wkdir}/'userdat.cheminput' | awk '{print $1}'`-1
head -n $nlin ${gecko_wkdir}/'userdat.cheminput' > ${gecko_wkdir}/lstprim.out

cd ${gecko_wkdir}
./findname_cmv > findname.out
rm findname_cmv

rm lstprim.out
mv findname.out ${mech}'.prec'
echo 'output => '${mech}'.prec'

echo ''
echo '------------------------'
echo Counting RO2s...
echo '------------------------'

if [ -e ${gecko_source}/RUN/COMPTEUR/compteur ] ; then
  rm ${gecko_source}/RUN/COMPTEUR/compteur ; fi

cd ${gecko_source}
#make compteur
ln -s ${gecko_source}/RUN/COMPTEUR/compteur ${gecko_wkdir}/compteur

cd ${gecko_wkdir}
./compteur

rm compteur
#rm fort*

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
# TODO: add test for git repository or not
#==========================================
#if [ ! -d .git ]; then
## GitHub branch:
#echo `git rev-parse --abbrev-ref HEAD` > ${gecko_wkdir}/${mechname}.gitinfo
## GitHub commit:
#echo `git describe`   >> ${gecko_wkdir}/${mechname}.gitinfo
#else
## IF NOT WITHIN A GIT REPOSITORY READ INFO FROM VERSION_NOTES DIRECTORY
#  cp ${gecko_source}/VERSION_NOTES/README_gitinfo ${mech}.gitinfo
#fi

#===========================================
# RUN INTERPRETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

echo ''
echo '------------------------'
echo run interpreter...
echo '------------------------'

cd ${home_dir}/SCRIPTS
./run_interp_local.bash -m ${mech}

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


