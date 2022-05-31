#!/bin/bash

netcdf_flag="yes"

module unload netcdf
module unload intel
module load gnu
module load netcdf

# launch to cheyenne with arguments
# may be called (currently by generate_scheme_cheyenne.bash)
# or run interactively from SCRIPTS dir, with all appropriate arguments

gecko_source=$1
gecko_wkdir=$2
mechname=$3
boxmod_source=$4

#--------------------
# supercomputer-specific stuff here
#--------------------

cd ${gecko_wkdir}

echo ''
echo '------------------------'
echo Linking input files
echo '------------------------'

if [ ! -e ${gecko_wkdir}/fort.7 ] ; then
  ln -s ${gecko_wkdir}/${mechname}'.dict' ${gecko_wkdir}/fort.7 ; fi

echo ''
echo '------------------------'
echo Finding Precursor Names
echo '------------------------'

if [ ! -f ${gecko_source}'/RUN/findname_cmv' ] ; then
  echo 'did not find findname_cmv in '${gecko_source}'/RUN'
  ls ${gecko_source}'/RUN'
  cd ${gecko_source}
  make findname_cmv
else
  echo 'found findname_cmv executable'
fi

cd ${gecko_wkdir}
ln -s ${gecko_source}/RUN/findname_cmv ${gecko_wkdir}/findname_cmv

let nlin=`wc -l ${gecko_wkdir}/'userdat.cheminput' | awk '{print $1}'`-1
head -n $nlin ${gecko_wkdir}/'userdat.cheminput' > ${gecko_wkdir}/lstprim.out

#./findname_cmv > findname.out
#mv findname.out ${mechname}'.prec'
./findname_cmv > ${mechname}'.prec'
# remove lines resulting from duplicate precursors
sed -i '/^ : /d' ${mechname}'.prec'

echo 'output => '${mechname}'.prec'
cp ${mechname}'.prec' findname_output

rm findname_cmv
rm lstprim.out 

echo''
echo '------------------------'
echo Counting RO2s...
echo '------------------------'

if [ ! -f ${gecko_source}'/RUN/COMPTEUR/compteur' ] ; then
  echo 'did not find compteur at '${gecko_source}'/RUN/COMPTEUR'
  ls ${gecko_source}'/RUN/COMPTEUR'
  cd ${gecko_source}
  make compteur
else
  echo 'found compteur executable'
fi
cd ${gecko_wkdir}
ln -s ${gecko_source}/RUN/COMPTEUR/compteur ${gecko_wkdir}/compteur

./compteur
rm compteur

#===========================================
# RUN AKPARAMETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

echo ''
echo '------------------------'
echo create akparameter file...
echo '------------------------'

cd ${gecko_source}/../SCRIPTS
./write_akparameter.bash ${mechname}

rm ${gecko_wkdir}/akparameter.h


#===========================================
# RUN INTERPRETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

echo ------------------------------------------
echo We are now in directory ...... ; pwd
echo ------------------------------------------

echo run interpreter...
./run_interp_cheyenne.bash ${mechname}

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
echo delete outstanding temporary files ...
echo ------------------------------------------
echo final directory listing ......
echo ------------------------------------------
ls

exit

