#!/bin/bash
# purpose: script to be submitted as a job to cheyenne to compile the boxmodel for
# the given mechanism
# also takes care of archiving sources and mechanism
# inputs: $1 = mechanism directory
#         $2 = boxmod source code directory

gecko_mechdir=$1
boxmod_source=$2
boxmod_wkdir=$3
gecko_outdir=$4

module list
echo "boxmodel executable -> compiling"
echo
rm -rf ${boxmod_wkdir}/WORK_BOXMOD
mkdir ${boxmod_wkdir}/WORK_BOXMOD
mkdir ${boxmod_wkdir}/WORK_BOXMOD/LIB

echo linking...
for f in `ls ${boxmod_source}/LIB/*.[fh]*`
do
  cp $f ${boxmod_wkdir}/WORK_BOXMOD/LIB/
done
#cp ${boxmod_source}/Makefile ${boxmod_wkdir}/WORK_BOXMOD/Makefile
mkdir ${boxmod_wkdir}/WORK_BOXMOD/OBJ
cp ${boxmod_source}/OBJ/Makefile ${boxmod_source}/OBJ/template.mk ${boxmod_wkdir}/WORK_BOXMOD/OBJ/
mkdir ${boxmod_wkdir}/WORK_BOXMOD/PROG
cp ${boxmod_source}/PROG/boxmod_main.f90 ${boxmod_wkdir}/WORK_BOXMOD/PROG/boxmod_main.f90

# link to the relevant akparameter.h file
echo akparameter...
#rm ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter.h
#cp ${gecko_mechdir}/akparameter.h ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter.h
cp ${gecko_mechdir}/akparameter_module.f90 ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter_module.f90

  echo Making...
# we don't want to change directory in this main script
# what's happening inside the parentheses has no impact outside of them
( cd ${boxmod_wkdir}/WORK_BOXMOD/OBJ && exec make ../PROG/BOXMOD )

if [ ! -e ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ] ; then
    echo error, could not find boxmodel executable
    echo ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD
    echo maybe check compilation?
    exit_status=1
else
	mv ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ${boxmod_wkdir}/BOXMOD

# package boxmod source files for future reference
#  	echo 'Creating source code tar file'
#  	mv ${boxmod_wkdir}/WORK_BOXMOD ${boxmod_wkdir}/BOXMOD_CODE
# remove OBJ files from BOXMOD_CODE; create tar file with INTERP, PROG and LIB dirs
#  	( cd ${boxmod_wkdir}/BOXMOD_CODE/OBJ && rm * && cd ../ && cp -R ${boxmod_source}/INTERP . && rm INTERP/inca && cd ${boxmod_wkdir} && tar -cvf boxmod_source_backup.tar BOXMOD_CODE )

# remove working directory BOXMOD_CODE
#  	rm -R ${boxmod_wkdir}/BOXMOD_CODE

# package mechanism files for future reference
#  	echo 'Creating mechanism tar file'
#  make temporary copy of mech dir, remove executable, binary, and link files; create tar file

#  	( cd ${boxmod_wkdir} && cp -R ${gecko_mechdir} GECKO_MECH_BACKUP && cd GECKO_MECH_BACKUP && rm BOXMOD && rm cm && rm compteur && rm outdat.ak* && find . -maxdepth 1 -type l -exec rm -f {} \; && cd  ${boxmod_wkdir} && tar -cvf gecko_mech_backup.tar GECKO_MECH_BACKUP )
#  	mv ${boxmod_wkdir}/gecko_mech_backup.tar ${boxmod_wkdir}/.
#  	rm -R ${boxmod_wkdir}/GECKO_MECH_BACKUP

	exit_status=0
fi

exit ${exit_status}

