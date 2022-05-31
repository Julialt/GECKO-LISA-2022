#!/bin/bash

#defaults 
flagfile="postproc_flags.input"

#===========================================
#== SET UP INPUT FILES==
#===========================================
while getopts "f:k:m:r:" opt
do
  case "$opt" in
    m ) mech=$OPTARG ;;
    k ) keyfile=$OPTARG ;;
    f ) flagfile=$OPTARG ;;
    r ) runname=$OPTARG ;;
# dir paths, obtainable from setup.library
#    b ) boxmod_source=$OPTARG ;;
#    h ) home_dir=$OPTARG ;;
#    w ) =$OPTARG ;;
  esac
done

# verify input arguments
if [ $# -lt 2 ]; then
    echo We need at least 2 input params: -m mech, -k keyfile
    echo .... argument -f flagfile is optional
    exit
fi

#-------------------
# >>keyfile name is originally in the form "indat_${keyfile}.key"
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

mech_name=${mech}
initfile=${keyfile}
flags_input="INPUTS/"${flagfile}
runname=$mech"_"$keyfile

#-------------------
# verify inputs
echo " "
echo "Inputs: -k keyfile -m mech (-f flagfile)"
echo ".... flagfile = "$flagfile
echo ".... keyfile = "$initfile 
echo ".... mech = "$mech 
echo ".... runname = "$runname
#-------------------


#===========================================
#== SET UP GECKO VERSION AND PATHS==
#===========================================
# ... IF SET HERE
#home_dir="/home/julial/GECKO/GIT_COPY"
#gecko_version="GECKO-A"
#run_dir="GECKO_SCRATCH"

# ... IF READ IN FROM file "setup.library"
if [ ! -e "setup.library" ] ; then
    echo error, gecko setup file could not be found
    echo file paths must be set in script
    exit
fi
home_dir=`      grep "home_dir"      setup.library | awk '{print $3}' `
gecko_version=` grep "gecko_version" setup.library | awk '{print $3}' `
gecko_inp_dir=` grep "gecko_inp_dir" setup.library | awk '{print $3}' `
scratch_dir=`   grep "scratch_dir"   setup.library | awk '{print $3}' `
gecko_run_dir=` grep "gecko_run_dir" setup.library | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" setup.library | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" setup.library | awk '{print $3}' `
postproc_version=` grep "postproc_version" setup.library | awk '{print $3}' `
postproc_run_dir=` grep "postproc_run_dir" setup.library | awk '{print $3}' `

#--construct paths and report to screen.
gecko_outdir=$scratch_dir/$gecko_run_dir
boxmod_outdir=$scratch_dir/$boxmod_run_dir
postproc_source=$home_dir/$postproc_version
postproc_outdir=$scratch_dir/$postproc_run_dir
gecko_source=$home_dir/$gecko_version

gecko_workingdir=${gecko_outdir}/${mech_name}
echo $gecko_workingdir
if [ ! -e ${gecko_workingdir} ] ; then
    echo error, gecko output directory could not be found
    echo ${gecko_workingdir}
    echo check gecko_outdir, mech_name variables
    echo in $0
    exit
fi

boxmod_workingdir=${boxmod_outdir}/${runname}
echo $boxmod_workingdir
if [ ! -e ${boxmod_workingdir} ] ; then
    echo error, boxmod run directory could not be found
    echo ${boxmod_workingdir}
    echo check boxmod_outdir, mech_name and initfile variables
    echo in $0
    exit
fi

postproc_workingdir=$postproc_outdir/${mech_name}/${initfile}
echo Creating ${postproc_workingdir}...
if [ -e ${postproc_workingdir} ] ; then
    #echo WARNING! removing existing directory
    echo WARNING! directory already exists !
    #rm -rf ${postproc_workingdir}
fi
mkdir -p ${postproc_workingdir}

if [ -e "${postproc_source}/${flags_input}" ] ; then
  echo linking custom postprocessing options ${postproc_source}/${flags_input}
  ln -s ${postproc_source}/${flags_input} ${postproc_workingdir}/postproc_flags.input
  echo ... to ${postproc_workingdir}/postproc_flags.input
else
    echo error, postprocessing inputs file could not be found
    echo check filename and try again
    exit
fi
#else 
#  echo linking default postprocessing options ${postproc_source}INPUTS/postproc_flags.input
#  ln -s ${postproc_source}/INPUTS/postproc_flags.input ${postproc_workingdir}/postproc_flags.input
#  echo ... to ${postproc_workingdir}/postproc_flags.input
#fi

echo linking boxmod results file ${boxmod_workingdir}/outdat.nc

if [ -e ${boxmod_workingdir}/outdat.nc ] ; then
  ln -s ${boxmod_workingdir}/outdat.nc ${postproc_workingdir}/outdat.nc
  echo ... to ${postproc_workingdir}/outdat.nc
else
  echo FILE NOT FOUND
  exit 99
fi

if [ -e ${boxmod_workingdir}/outdat.pbl ] ; then
   echo linking pbl data file ${boxmod_workingdir}/outdat.pbl
   ln -s ${boxmod_workingdir}/outdat.pbl ${postproc_workingdir}/outdat.pbl
   echo ... to ${postproc_workingdir}/outdat.pbl
elif [ -e ${boxmod_workingdir}/USEROUT/outdat.pbl ] ; then
   echo linking pbl data file ${boxmod_workingdir}/USEROUT/outdat.pbl
   ln -s ${boxmod_workingdir}/USEROUT/outdat.pbl ${postproc_workingdir}/outdat.pbl
   echo ... to ${postproc_workingdir}/outdat.pbl
fi

echo linking dictionary from ${gecko_workingdir}/${mech_name}.dict 
ln -s ${gecko_workingdir}/${mech_name}.dict ${postproc_workingdir}/dictionary
echo ... to ${postproc_workingdir}/dictionary

echo linking pvap data from ${gecko_workingdir}/${mech_name}.pnan 
ln -s ${gecko_workingdir}/${mech_name}.pnan  ${postproc_workingdir}/pvap.dat
echo ... to ${postproc_workingdir}/pvap.dat

echo linking henry data from ${gecko_workingdir}/${mech_name}.Henry 
ln -s ${gecko_workingdir}/${mech_name}.Henry ${postproc_workingdir}/henry.dat
echo ... to ${postproc_workingdir}/henry.dat

cd ${gecko_workingdir}
echo counting primary species codes...

#ln -s ${mech_name}.prec findname_output
nprecu=`wc -l findname_output | cut -f1 -d" "`

echo linking data from ${gecko_workingdir}/${mech_name}.prec
echo -n " &userinput
  precursor_codes(1:$nprecu)=" > ${postproc_workingdir}/userinput.nml
for i in `cat findname_output | cut -f 1 -d":"` ; do echo -n \"G$i\"" " >> ${postproc_workingdir}/userinput.nml ; done
echo "," >> ${postproc_workingdir}/userinput.nml
cat ${postproc_workingdir}/postproc_flags.input >> ${postproc_workingdir}/userinput.nml

echo ${primary}
echo and creating ${postproc_workingdir}/user_input.nml

echo compiling postprocessor...
mkdir -p  ${postproc_workingdir}/POSTPROC_CODE/OBJ
cp -r $postproc_source/LIB $postproc_source/RUN ${postproc_workingdir}/POSTPROC_CODE/
cp $postproc_source/OBJ/Makefile $postproc_source/OBJ/template.mk ${postproc_workingdir}/POSTPROC_CODE/OBJ

# compile the gecko code just in case geckolib.a does not exist
( cd $gecko_source && make cm)
if [ ! -e $gecko_source/OBJ/geckolib.a ] ; then
  echo cannot find $gecko_source/OBJ/geckolib.a needed to compile the postprocessor
  exit
fi
cp  $gecko_source/OBJ/geckolib.a ${postproc_workingdir}/POSTPROC_CODE/OBJ

cd ${postproc_workingdir}/POSTPROC_CODE/OBJ
make ../RUN/postproc 

echo linking postproc executable ${postproc_workingdir}/POSTPROC_CODE/RUN/postproc
ln -s ${postproc_workingdir}/POSTPROC_CODE/RUN/postproc ${postproc_workingdir}/postproc
echo ... to ${postproc_workingdir}/postproc

echo running postproc
cd ${postproc_workingdir}
pwd

./postproc

# copy over runtime-produced ascii files
if [ -e ${boxmod_workingdir}/USEROUT ] ; then
  cp ${boxmod_workingdir}/USEROUT/outdat.jvals ${postproc_workingdir}
  cp ${boxmod_workingdir}/USEROUT/outdat.ro2 ${postproc_workingdir}
  cp ${boxmod_workingdir}/USEROUT/outdat.pbl ${postproc_workingdir}
else 
  cp ${boxmod_workingdir}/outdat.jvals ${postproc_workingdir}
  cp ${boxmod_workingdir}/outdat.ro2 ${postproc_workingdir}
  cp ${boxmod_workingdir}/outdat.pbl ${postproc_workingdir}
fi

#--tidy up
#rm ${gecko_workingdir}/findname_output
rm -R ${postproc_workingdir}/POSTPROC_CODE
rm -R ${postproc_workingdir}/postproc
# remove symbolic links
find -lname '*' -delete


#--copy to ftp (optional)
cp * /ftp/user/julial
