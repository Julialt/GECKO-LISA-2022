#!/bin/bash

## are we using NetCDF? (user hardwire "yes" or "")
netcdf_flag="yes"

## get paths from setup.dat ##
pathfile="setup.dat"

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `

echo 'home dir = ' $home_dir
echo 'gecko_version = ' $gecko_version
echo 'scratch dir = ' $scratch_dir
echo 'gecko_run_ dir = ' $gecko_run_dir
echo 'boxmod_version = ' $boxmod_version

while getopts "a:b:c:d:h:i:k:m:p:s:v:w:x" opt
do
  case "$opt" in
    m ) mech=$OPTARG ;;
# any other code produces an error
    * ) echo "ERROR! flag not recognized !"; exit ;;
  esac
done

# construct paths and report to screen.
mech_name=$mech
echo mech_name = ${mech}

gecko_source=$home_dir"/"$gecko_version"/LIB"
echo 'gecko_source = ' $gecko_source

gecko_data=$home_dir"/"$gecko_version"/DATA"
echo 'gecko_data = ' $gecko_data

boxmod_dir=$home_dir"/"$boxmod_version"/LIB"
echo 'boxmod_dir = ' $boxmod_dir

work_dir=$scratch_dir"/"$gecko_run_dir"/"$mech_name
echo 'work dir = ' $work_dir

interp_dir=$home_dir"/"$boxmod_version"/INTERP"
echo 'interp dir = ' $interp_dir

# check for existence of mechanism, shared files

if [ ! -e $work_dir'/'$mech_name'.mech' ] ; then
	echo 'error, could not find file '$work_dir/$mech_name
	exit ; fi

if [ ! -e $work_dir'/akparameter_module.f90' ] ; then
	echo 'error, could not find file '$work_dir/$mech_name
	exit ; fi

if [ ! -e $boxmod_dir'/general.h' ] ; then
	echo 'error, could not find file general.h'
	exit ; fi

if [ ! -e $gecko_source'/common.h' ] ; then
        echo 'error, could not find file common.h'
	exit ; fi

# Create reference file git_info.txt in mechanism directory, to be read by codes.
# GitHub branch:
#echo `git rev-parse --abbrev-ref HEAD` > $work_dir/${mech}.gitinfo
# GitHub commit:
#echo `git describe`   >> $work_dir/${mech}.gitinfo

# Here we need to compile the interpreter because it depends on akparameter_module.f90
if [ -e ${work_dir}/WORK_INTERP ] ; then 
 rm -rf ${work_dir}/WORK_INTERP ; fi
mkdir ${work_dir}/WORK_INTERP

ln -s ${work_dir}/akparameter_module.f90 ${work_dir}/WORK_INTERP/akparameter_module.f90
ln -s ${boxmod_dir}/sorting_module.f90 ${work_dir}/WORK_INTERP/sorting_module.f90
ln -s ${boxmod_dir}/general.h ${work_dir}/WORK_INTERP/general.h
ln -s ${gecko_source}/common.h ${work_dir}/WORK_INTERP/common.h
ln -s ${interp_dir}/*.f* ${work_dir}/WORK_INTERP/
ln -s ${interp_dir}/template.mk ${work_dir}/WORK_INTERP/template.mk
ln -s ${interp_dir}/Makefile ${work_dir}/WORK_INTERP/Makefile

# which makefile depends on NetCDF_flag
#if [ -n ${netcdf_flag} ]; then
#  ln -s ${interp_dir}/makefile_ncdf.intp.local ${work_dir}/WORK_INTERP/makefile
#else
#  ln -s ${interp_dir}/makefile.intp.local ${work_dir}/WORK_INTERP/makefile
#fi

# compile the interpreter (after linking to akparameter.h)
cd ${work_dir}/WORK_INTERP
make all

interp_exe="inca"

if [ ! -e ${work_dir}/WORK_INTERP/${interp_exe} ] ; then
	echo error, could not find interpreter
	echo ${work_dir}/WORK_INTERP/${interp_exe}
	echo maybe check compilation?
	exit
fi
cp ${work_dir}/WORK_INTERP/${interp_exe} ${work_dir}/intp.exe

echo ------------------------------------------
echo We are now in directory ......
pwd
echo ....interpreting $mech_name
echo ------------------------------------------

ln -sf ${gecko_data}/simpol.dat ${work_dir}/simpol.dat
ln -sf ${mech_name}.mech ${work_dir}/indat.mech
ln -sf ${mech_name}.dict ${work_dir}/indat.dict
ln -sf ${mech_name}.prec ${work_dir}/indat.prec
ln -sf ${mech_name}.difv ${work_dir}/indat.difv
ln -sf ${mech_name}.pnan ${work_dir}/indat.pnan
ln -sf ${mech_name}.psim ${work_dir}/indat.psim
ln -sf ${mech_name}.Henry ${work_dir}/indat.Henry
ln -sf ${mech_name}.kNO3 ${work_dir}/indat.kNO3
ln -sf ${mech_name}.kO3 ${work_dir}/indat.kO3
ln -sf ${mech_name}.kOH ${work_dir}/indat.kOH
ln -sf ${mech_name}.gitinfo ${work_dir}/indat.gitinfo
ln -sf userdat.settings ${work_dir}/userparams.input

echo linking RO2 files...
ln -sf ${work_dir}/XP1O2 ${work_dir}/indat1.ro2
ln -sf ${work_dir}/XP2O2 ${work_dir}/indat2.ro2
ln -sf ${work_dir}/XP3O2 ${work_dir}/indat3.ro2
ln -sf ${work_dir}/XS1O2 ${work_dir}/indat4.ro2
ln -sf ${work_dir}/XS2O2 ${work_dir}/indat5.ro2
ln -sf ${work_dir}/XS3O2 ${work_dir}/indat6.ro2
ln -sf ${work_dir}/XT1O2 ${work_dir}/indat7.ro2
ln -sf ${work_dir}/XT2O2 ${work_dir}/indat8.ro2
ln -sf ${work_dir}/XACO3 ${work_dir}/indat9.ro2
echo

cd ${work_dir}

echo ------------------------------------------
#echo ! ! We are now in ${workdir} ! !
echo We are now in directory ......
pwd
echo ------------------------------------------

## free up memory to avoid initial segfaults:
ulimit -s unlimited

## run interpreter
./intp.exe
#nohup ./intp.exe

echo ------------------------------------------
echo tidy up...
echo ------------------------------------------

mv outdat.akoi outdat.report

rm -R WORK_INTERP
rm intp.exe
rm dummy*
rm userparams.input

if [ ${netcdf_flag} ]; then
  rm outdat.akli ; fi
if [ -e bidon ] ; then
     rm bidon ; fi

# remove symbolic links
find -lname '*' -delete

echo ------------------------------------------
echo Completed mechanism is in directory ......
pwd
echo ------------------------------------------
echo directory listing ......
ls
echo ------------------------------------------
