#!/bin/bash
# purpose: script to be submitted as a job to compile the generator for
# the given mechanism
# inputs: $1 = gecko directory

gecko_source=$1
gecko_wkdir=$2

mkdir ${gecko_wkdir}/GECKO_CODE
cp -r ${gecko_source}/* ${gecko_wkdir}/GECKO_CODE
cd ${gecko_wkdir}/GECKO_CODE
make clean
make cm
make compteur
make findname_cmv

