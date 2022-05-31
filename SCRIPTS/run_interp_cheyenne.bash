#!/bin/bash
# script gets called by gen_package_cheyenne.bash
# requires presence of X-files
# so CANNOT be run interactively as a standalone:
# requires computeur to be run FIRST.

## are we using NetCDF? (user hardwire "yes" or "")
netcdf_flag="yes"

#-- default values for paths ==:
pathfile="setup.dat"
echo pathfile = $pathfile

if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

## get paths from setup.dat ##

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

if [ $# -ne 1 ]; then
  echo "enter mechanism name"
  read mechname
    if [ ${#mechname} -eq 0 ] ; then echo "no value supplied" ; exit 2 ; fi
else
  mechname=$1
fi


work_dir=$scratch_dir"/"$gecko_run_dir"/"$mechname
gecko_source=$home_dir"/"$gecko_version"/LIB"
gecko_data=$home_dir"/"$gecko_version"/DATA"
boxmod_dir=$home_dir"/"$boxmod_version"/LIB"
interp_dir=$home_dir"/"$boxmod_version"/INTERP"

interp_exe="inca"

echo "working dir (i/o) = "$work_dir
echo "interp_dir = "$interp_dir

#------------------------------#
## check for file existence   ##
#------------------------------#

if [ ! -e $work_dir'/'$mechname'.mech' ] ; then
	echo error, could not find mech file
	echo $work_dir/$mechname
	exit
fi

if [ ! -e $work_dir'/akparameter_module.f90' ] ; then
	echo error, could not find akparameter file
	echo $work_dir/$mechname
	exit
fi

if [ ! -e $boxmod_dir'/general.h' ] ; then
	echo error, could not find file general.h
	exit
fi

if [ ! -e $gecko_source'/common.h' ] ; then
	echo error, could not find file common.h
	exit
fi

module load git

# Create reference file git_info.txt in mechanism directory, to be read by codes.
# GitHub branch:
echo `git rev-parse --abbrev-ref HEAD` > $work_dir/$mech.gitinfo
# GitHub commit:
echo `git describe`   >> $work_dir/$mech.gitinfo
ln -s ${work_dir}/$mech.gitinfo ${work_dir}/indat.gitinfo

# Here we need to compile the interpreter because it depends on akparameter.h
rm -rf ${work_dir}/WORK_INTERP
mkdir ${work_dir}/WORK_INTERP
ln -s ${interp_dir}/*.f ${work_dir}/WORK_INTERP/
ln -s ${interp_dir}/*.f90 ${work_dir}/WORK_INTERP/
ln -s ${work_dir}/akparameter_module.f90 ${work_dir}/WORK_INTERP/akparameter_module.f90
ln -s ${boxmod_dir}/sorting_module.f90 ${work_dir}/WORK_INTERP/sorting_module.f90
ln -s ${boxmod_dir}/general.h ${work_dir}/WORK_INTERP/general.h
ln -s ${gecko_source}/common.h ${work_dir}/WORK_INTERP/common.h
ln -s ${interp_dir}/template.mk ${work_dir}/WORK_INTERP/template.mk
ln -s ${interp_dir}/Makefile ${work_dir}/WORK_INTERP/Makefile

# which makefile depends on NetCDF_flag
#if [ -n "${netcdf_flag}" ]; then
#  ln -s ${interp_dir}/makefile_ncdf.intp ${work_dir}/WORK_INTERP/makefile
#else
#  ln -s ${interp_dir}/makefile.intp ${work_dir}/WORK_INTERP/makefile
#fi

cd ${work_dir}/WORK_INTERP/
make

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
echo ....interpreting $mechname
echo ------------------------------------------


ln -s ${gecko_data}/simpol.dat ${work_dir}/simpol.dat
ln -s ${mechname}.mech ${work_dir}/indat.mech
ln -s ${mechname}.dict ${work_dir}/indat.dict
ln -s ${mechname}.prec ${work_dir}/indat.prec
ln -s ${mechname}.difv ${work_dir}/indat.difv
ln -s ${mechname}.pnan ${work_dir}/indat.pnan
ln -s ${mechname}.psim ${work_dir}/indat.psim
ln -s ${mechname}.Henry ${work_dir}/indat.Henry
ln -s ${mechname}.kNO3 ${work_dir}/indat.kNO3
ln -s ${mechname}.kOH  ${work_dir}/indat.kOH
ln -s ${mechname}.kO3  ${work_dir}/indat.kO3
ln -s userdat.settings ${work_dir}/userparams.input

echo linking RO2 files...
ln -s ${work_dir}/XP1O2 ${work_dir}/indat1.ro2
ln -s ${work_dir}/XP2O2 ${work_dir}/indat2.ro2
ln -s ${work_dir}/XP3O2 ${work_dir}/indat3.ro2
ln -s ${work_dir}/XS1O2 ${work_dir}/indat4.ro2
ln -s ${work_dir}/XS2O2 ${work_dir}/indat5.ro2
ln -s ${work_dir}/XS3O2 ${work_dir}/indat6.ro2
ln -s ${work_dir}/XT1O2 ${work_dir}/indat7.ro2
ln -s ${work_dir}/XT2O2 ${work_dir}/indat8.ro2
ln -s ${work_dir}/XACO3 ${work_dir}/indat9.ro2
echo

cd ${work_dir}

echo ------------------------------------------
#echo ! ! We are now in ${workdir} ! !
echo We are now in directory ......
pwd
echo ------------------------------------------
./intp.exe

echo ------------------------------------------
echo tidy up...
echo ------------------------------------------

mv outdat.akoi outdat.report

rm bidon
rm intp.exe
rm dummy*
rm userparams.input
if [ ${netcdf_flag} ]; then
  rm outdat.akli ; fi
rm -rf ${work_dir}/WORK_INTERP

# remove symbolic links
find -lname '*' -delete

echo ------------------------------------------
echo Completed mechanism is in directory ......
pwd
echo ------------------------------------------
echo directory listing ......
ls
echo ------------------------------------------

