#!/bin/bash
## UNDER DEVELOPMENT
# this script must be launched from the $boxmod_wkdir folder

boxmod_wkdir=$1

cd $boxmod_wkdir

if [ -e ./BOXMOD ]; then
  #dereference link
  cp ./BOXMOD ./BOXMOD_tmp
  mv ./BOXMOD_tmp ./BOXMOD

  # run boxmod
  export OMP_NUM_THREADS=$NCPUS
  ./BOXMOD > $TMPDIR/stdout.$PBS_JOBID
else
  echo couldnt find BOXMOD executable
  echo check compilation
  exit 1
fi

# check for results files
#if [ -s ${boxmod_wkdir}/outdat.ppf ] ; then
#  echo box model has run successfully
#  echo check output ${boxmod_wkdir}/outdat.ppf
#fi
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

## delete executable
rm BOXMOD

## delete work directory 
## comment this out for testing, if you want to avoid re-compiling every time
rm -R WORK_BOXMOD

## delete symbolic links
#find -lname '*' -delete

## delete zero-length files
#find ./ -type f -size 0 -delete

## delete temporary files
#rm *.gitinfo

#-----------
## move output into subdirectory
mkdir USEROUT
mv *.* USEROUT
cd USEROUT

# ... and move standard output back to main run directory
#     (leaving non-standard output in subdirectory)
mv outdat.nc ../
mv outdat.README ../
mv indat.key ../
#-----------


