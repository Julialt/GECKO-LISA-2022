#!/bin/bash
#===========================================
# script: run_postproc_cheyenne.bash
# intent: to be called by other Cheyenne scripts
# syntax: ./run_postproc_cheyenne.bash {mech} {case} {flagfile}
# (DO NOT use "-x" flags in the invoking command)
# tested: JMLT, NCAR, 2020-12-23
#===========================================

#===========================================
#== SET UP GECKO VERSION AND PATHS==
#===========================================
# SET HERE ... (by uncommenting following 3 lines)
#home_dir="/home/julial/GECKO/GIT_COPY"
#gecko_version="GECKO-A"
#run_dir="GECKO_SCRATCH"

# ... OR READ IN FROM file "setup.dat"
if [ ! -e "setup.dat" ] ; then
    echo error, gecko setup file could not be found
    echo file paths must be set in script
    exit
fi
home_dir=`      grep "home_dir"      setup.dat | awk '{print $3}' `
gecko_version=` grep "gecko_version" setup.dat | awk '{print $3}' `
gecko_inp_dir=` grep "gecko_inp_dir" setup.dat | awk '{print $3}' `
scratch_dir=`   grep "scratch_dir"   setup.dat | awk '{print $3}' `
gecko_run_dir=` grep "gecko_run_dir" setup.dat | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" setup.dat | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" setup.dat | awk '{print $3}' `
postproc_version=` grep "postproc_version" setup.dat | awk '{print $3}' `
postproc_run_dir=` grep "postproc_run_dir" setup.dat | awk '{print $3}' `

#--construct paths and report to screen.
gecko_outdir=$scratch_dir/$gecko_run_dir
boxmod_outdir=$scratch_dir/$boxmod_run_dir
postproc_source=$home_dir/$postproc_version
postproc_outdir=$scratch_dir/$postproc_run_dir
gecko_source=$home_dir/$gecko_version

# 2 arguments
if [ $# -lt 2 ]; then
    echo We need at least 2 input parameters
    exit
fi

mechname=$1
casename=$2
flags_input=$3

# if invoked from command line will need to reload modules for netcdf to work
    module unload intel
    module unload netcdf
    module load gnu
    module load netcdf

gecko_workingdir=${gecko_outdir}/${mechname}
if [ ! -e ${gecko_workingdir} ] ; then
    echo error, gecko output directory could not be found
    echo ${gecko_workingdir}
    echo check gecko_outdir, mechname variables
    echo in $0
    exit
fi

boxmod_workingdir=${boxmod_outdir}/${mechname}_${casename}
if [ ! -e ${boxmod_workingdir} ] ; then
    echo error, boxmod run directory could not be found
    echo ${boxmod_workingdir}
    echo check boxmod_outdir, mechname and casename variables
    echo in $0
    exit
fi

postproc_workingdir=$postproc_outdir/${mechname}/${casename}
echo Creating ${postproc_workingdir}...
if [ -e ${postproc_workingdir} ] ; then
  echo WARNING! directory already exists !
else
  mkdir -p ${postproc_workingdir}
fi

if [ -n "${flags_input}" ] ; then
  echo linking custom postprocessing options ${postproc_source}/INPUTS/${flags_input}
  ln -sf ${postproc_source}/INPUTS/${flags_input} ${postproc_workingdir}/postproc_flags.input
  echo ... to ${postproc_workingdir}/postproc_flags.input
else 
  echo linking default postprocessing options ${postproc_source}/INPUTS/postproc_flags.input
  ln -sf ${postproc_source}/INPUTS/postproc_flags.input ${postproc_workingdir}/postproc_flags.input
  echo ... to ${postproc_workingdir}/postproc_flags.input
fi


echo linking boxmod results file ${boxmod_workingdir}/outdat.nc

if [ -e ${boxmod_workingdir}/outdat.nc ] ; then
  ln -s ${boxmod_workingdir}/outdat.nc ${postproc_workingdir}/outdat.nc
else
  echo FILE NOT FOUND
  exit 99
fi

#if [ -e ${boxmod_workingdir}/fort.13 ] ; then
#  mv ${boxmod_workingdir}/fort.13 ${boxmod_workingdir}/outdat.ppf
#fi
#ln -s ${boxmod_workingdir}/outdat.ppf ${postproc_workingdir}/outdat.ppf
#echo ... to ${postproc_workingdir}/outdat.ppf

#if [ -e ${boxmod_workingdir}/fort.17 ] ; then
#  mv ${boxmod_workingdir}/fort.17 ${boxmod_workingdir}/outdat.ppa
#fi
#ln -s ${boxmod_workingdir}/outdat.ppa ${postproc_workingdir}/outdat.ppa
#echo ... to ${postproc_workingdir}/outdat.ppa

#if [ -e ${boxmod_workingdir}/outdat.pff ] ; then
#   echo linking binary results file ${boxmod_workingdir}/outdat.pff
#   ln -s ${boxmod_workingdir}/outdat.pff ${postproc_workingdir}/outdat.pff
#   echo ... to ${postproc_workingdir}/outdat.pff
#fi

if [ -e ${boxmod_workingdir}/outdat.pbl ] ; then
   echo linking pbl data file ${boxmod_workingdir}/outdat.pbl
   cp ${boxmod_workingdir}/outdat.pbl ${postproc_workingdir}/outdat.pbl
   echo ... to ${postproc_workingdir}/outdat.pbl
elif [ -e ${boxmod_workingdir}/USEROUT/outdat.pbl ] ; then
   echo linking pbl data file ${boxmod_workingdir}/USEROUT/outdat.pbl
   ln -s ${boxmod_workingdir}/USEROUT/outdat.pbl ${postproc_workingdir}/outdat.pbl
   echo ... to ${postproc_workingdir}/outdat.pbl
fi

echo linking dictionary from ${gecko_workingdir}/${mechname}.dict 
ln -s ${gecko_workingdir}/${mechname}.dict ${postproc_workingdir}/dictionary
echo ... to ${postproc_workingdir}/dictionary
echo linking pvap data from ${gecko_workingdir}/${mechname}.sat 
ln -s ${gecko_workingdir}/${mechname}.pnan  ${postproc_workingdir}/pvap.dat
echo ... to ${postproc_workingdir}/pvap.dat

echo linking henry data from ${gecko_workingdir}/${mechname}.Henry 
ln -s ${gecko_workingdir}/${mechname}.Henry ${postproc_workingdir}/henry.dat
echo ... to ${postproc_workingdir}/henry.dat

cd ${gecko_workingdir}
echo counting primary species codes...

ln -s ${mechname}.prec findname_output
nprecu=`wc -l findname_output | cut -f1 -d" "`

echo linking data from ${gecko_workingdir}/${mechname}.prec
echo -n " &userinput
  precursor_codes(1:$nprecu)=" > ${postproc_workingdir}/userinput.nml
for i in `cat findname_output | cut -f 1 -d":"` ; do echo -n \"G$i\"" " >> ${postproc_workingdir}/userinput.nml ; done
echo "," >> ${postproc_workingdir}/userinput.nml
cat ${postproc_workingdir}/postproc_flags.input >> ${postproc_workingdir}/userinput.nml

echo ${primary}
echo and creating ${postproc_workingdir}/user_input.nml

echo compiling postprocessor...
#mkdir -p  ${postproc_workingdir}/POSTPROC_CODE/OBJ ${postproc_workingdir}/POSTPROC_CODE/RUN
mkdir -p  ${postproc_workingdir}/POSTPROC_CODE/OBJ
cp -r $postproc_source/LIB $postproc_source/RUN ${postproc_workingdir}/POSTPROC_CODE/
cp $postproc_source/OBJ/Makefile $postproc_source/OBJ/template.mk ${postproc_workingdir}/POSTPROC_CODE/OBJ

# compile the gecko code just in case geckolib.a does not exist
if [ ! -e $gecko_source/OBJ/geckolib.a ] ; then
  echo cannot find $gecko_source/OBJ/geckolib.a needed to compile the postprocessor
  ( cd $gecko_source && make cm)
  #exit
else # check date on existing geckolib.a: replace if old
  today=$(date --iso-8601=date )
  #echo $today
  compiledate=$(stat -c %y $gecko_source'/OBJ/geckolib.a' | cut -c1-10)
  #echo $compiledate
  if [ $compiledate == $today ]; then
    echo "geckolib.a file is up-to-date"
  else
    echo "making new geckolib.a file"
    ( cd $gecko_source && rm 'OBJ/geckolib.a' && make cm )
  fi
fi
cp  $gecko_source/OBJ/geckolib.a ${postproc_workingdir}/POSTPROC_CODE/OBJ

cd ${postproc_workingdir}/POSTPROC_CODE/OBJ
make ../RUN/postproc 

echo linking postproc executable ${postproc_workingdir}/POSTPROC_CODE/RUN/postproc
ln -s ${postproc_workingdir}/POSTPROC_CODE/RUN/postproc ${postproc_workingdir}/postproc
echo ... to ${postproc_workingdir}/postproc

echo running postproc
cd ${postproc_workingdir}
./postproc

##--tidy up
rm ${gecko_workingdir}/findname_output
rm -R POSTPROC_CODE
## remove symbolic links
find -lname '*' -delete

