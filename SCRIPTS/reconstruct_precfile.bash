#!/bin/bash
#---------------------------
# script to extract precursor formulae
# from outdat.nc IF IT IS PRESENT
# (can be either a mech file or a run file)
# Recreates userdat.cheminput and <mechname>.prec
#---------------------------
# Calling syntax:
# ./reconstruct_precfile.bash -m {mechname} 
# => finds mech dir and rebuilds *.prec and userdat.cheminput
# ./reconstruct_precfile.bash -r {runname} 
# => finds run dir and rebuilds *.prec and userdat.cheminput IN THE RUN DIR
#---------------------------

#-- default values for some arguments ==:
pathfile="setup.dat"

spinup_fg=0
while getopts "a:m:r:" opt
do
  case "$opt" in
# pathfile =usually setup.dat
    a ) pathfile=$OPTARG ;;
    m ) mech=$OPTARG ;;
    r ) runname=$OPTARG ;;
  esac
done

echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" ${pathfile} | awk '{print $3}' `

boxmod_source=$home_dir/$boxmod_version

if [ "$mech" != "" ] ; then
  outname=${mech}
  outdir=$scratch_dir/$gecko_run_dir/${outname}
else
  # use boxmod run name
  outname=${runname}
  outdir=$scratch_dir/$boxmod_run_dir/${outname}
fi

echo "outdir = "$outdir
cd $outdir

# find if precursor chem list exists in outdat.nc
ncdump -h  outdat.nc > tmp.hdr
nlin=`grep precnam tmp.hdr > tmp.txt | cut -b -1`

# if precursor list present in netcdf output, do the following:
if [ "$nlin" != "0" ]; then
  # send the precursor data lines to a new file (w/o header)
  # strip out initial "
  ncdump -v precnam  outdat.nc |sed -e '1,/data:/d' -e '$d' > tmp.prec
  grep '"' tmp.prec | awk '{print $1}' > $outname'.naminput'
  sed -i 's/\"//g' $outname'.naminput'

  ncdump -v precchem  outdat.nc |sed -e '1,/data:/d' -e '$d' > tmp.prec
  grep '"' tmp.prec | awk '{print $1}' > $outname'.cheminput'
  sed -i 's/\"//g' $outname'.cheminput'

  # different approach needed for precnc: it's an integer, not a character.
  ncdump -v precnc  outdat.nc |sed -e '1,/data:/d' -e '$d' > tmp.prec
  tail -n 1 tmp.prec > $outname'.ncinput'
  sed -i 's/precnc =//g' $outname'.ncinput'
  sed -i 's/;//g' $outname'.ncinput'
  sed -i 's/ //g' $outname'.ncinput'
  sed -i 's/,/\n/g' $outname'.ncinput'

  # find numer of precursors (for info only)
  nprec=`wc -l $outname'.naminput' | cut -f1 -d" "` ; echo "nprec = "$nprec

# write precfile (overwriting previous versions)
echo -n "" > $outname'.prec'
counter=0
for i in `cat  $outname'.naminput' ` 
do  
  let counter=${counter}+1
  nc=`head -n ${counter} $outname'.ncinput'   | tail -n 1`
  chem=`head -n ${counter} $outname'.cheminput' | tail -n 1`
  echo  " "$i": "${chem}":     "${nc} >> $outname'.prec'
done

# if precursor list NOT present in netcdf output, do the following:
else
  echo "precursor list not present in NetCDF output : Need ascii version!"
fi

# tidy up
echo END >> $outname'.cheminput'
mv $outname'.cheminput' userdat.cheminput

rm tmp.hdr
rm tmp.txt
rm tmp.prec
rm *.naminput
rm *.ncinput
