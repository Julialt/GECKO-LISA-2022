#!/bin/bash
#===========================================
# run boxmodel from command line
# CALLING SYNTAX: ./run_boxmod_local.bash -m (mech) -k (keyfile) -v (previous run) -p (photfile)
# STILL TO DO: 
#            - implement sequential runs
# SCRIPT DOES NOT:
#                - auto-invoke postprocessing. Must be done interactively
#===========================================

#-- default values for some arguments ==:
pathfile="setup.dat"
photfile="O3_300DU.phot"
inpdir="INPUTS"
nthreads=1

#===========================================
#== INPUT ARGUMENTS ==
#===========================================
#== general files
# a = file containing source paths (optional: default is setup.dat)
# d = dependencies (optional)
# m = mechanism name (output from generator, input for box model)
#== generator-specific files
# i = cheminput file
# s = settings file for generator
# x = existing mechanism
#== boxmodel-specific files
# c = constrained concentrations file
# k = keyfile (box model input scenario)
# p = phot file (optional: default is sophie0404_julia_cmv)
# u = directory path for keyfile ("INPUTS/__") default = "INPUTS"
# v = previous output (created with same mechanism)
#== dir paths (optional: obtainable from setup.dat)
# b = box model source path
# h = home directory
# w = gecko working dir
#--

## HARDWIRE FOR DEBUG ##
#mechname="pentane"
#keyfile="indat_benchtest_pentane_NOx_2box"
#prevout="pentane_benchtest_pentane_NOx_2box"
#photfile="BEACHON_280DU_0tauaer"
## END DEBUG SCETION ##

spinup_fg=0
while getopts "a:b:c:d:f:h:i:k:m:p:s:t:u:v:w:x:" opt
do 
  case "$opt" in
# general files
    a ) pathfile=$OPTARG ;;   
    d ) depend=$OPTARG ;;   
    m ) mechname=$OPTARG ;;   
# generator-specific files
#    i ) cheminput=$OPTARG ;;   
#    s ) settings=$OPTARG ;;   
#    x ) existing=$OPTARG ;;   
# boxmodel-specific files
    c ) confile=$OPTARG ;;   
    k ) keyfile=$OPTARG ;;   
    p ) photfile=$OPTARG ;; 
    v ) prevout=$OPTARG ;;   
    u ) inpdir=$OPTARG ;;
# dir paths, obtainable from setup.dat
#    b ) boxmod_source=$OPTARG ;;   
#    h ) home_dir=$OPTARG ;;   
#    w ) =$OPTARG ;;   
#run option
    t ) nthreads=$OPTARG ;;
  esac
done

#===========================================
#== SET UP GECKO VERSION AND PATHS==
#===========================================
echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
gecko_version=`  grep "gecko_version"  ${pathfile} | awk '{print $3}' `
boxmod_version=` grep "boxmod_version" ${pathfile} | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
gecko_data=$home_dir"/"$gecko_version"/DATA"
gecko_outdir=$scratch_dir/$gecko_run_dir
boxmod_source=$home_dir/$boxmod_version

if [  ${#inpdir} -eq 0 ] ; then
    inpdir="INPUT"
fi
boxmod_inpdir=$boxmod_source/$inpdir
boxmod_outdir=$scratch_dir/$boxmod_run_dir

echo "boxmod source path = "$boxmod_source
echo "boxmod inputs path = "$boxmod_inpdir
echo "boxmod output path = "$boxmod_outdir

#--previous output path, if relevant
if [ -n "${prevout}" ]; then
  boxmod_prevdir=$boxmod_outdir/$prevout
  echo "previous output in = "$boxmod_prevdir
fi

#---------check prefix/extension of input data file names -------

# >>keyfile name is used in the script as "indat_${keyfile}.key"
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

# >>confile: if a filename is supplied:
if [ ${#confile} -ne 0 ] ; then
# confile name is used in the script as "init_conc_${confile}.key"
# remove "init_conc_" prefix from $confile if supplied
  prefix=${confile: 0: 10}
  if [ "$prefix" == "init_conc_" ] ; then
    confile=${confile##${prefix}}
  fi
# remove ".key" extension from $confile if supplied
  extn=${confile: -4}
  if [ "$extn" == ".key" ] ; then
    confile=${confile%${extn}}
  fi
fi

# >>photfile name is used in the script as "${photfile}"
# add ".phot" extension to photfile if not already supplied
extn=${photfile: -5}
if [ "$extn" != ".phot" ] ; then
  photfile=$photfile.phot
fi

# set value of netcdf_flag based on flag in keyfile
ifmt=`grep "IFMT" ${boxmod_inpdir}/indat_${keyfile}.key | awk '{print $2}' `
if [ "$ifmt" = 1 ] ; then
  netcdf_flag=""
else
  netcdf_flag="yes"
fi

#---------file existence checks -------
#-- check for boxmod directories
if [ ! -e ${boxmod_source} ] ; then
    echo error, boxmod source directory could not be found
    echo ${boxmod_source}
    exit 1
fi
if [ ! -e ${boxmod_inpdir} ] ; then
    echo error, boxmod inputs directory could not be found
    echo ${boxmod_inpdir}
    exit 1
fi
if [ ! -e ${boxmod_outdir} ] ; then
    echo error, boxmod run directory could not be found
    echo ${boxmod_outdir}
    exit 1
fi
if [ -n "${prevout}" ]; then 
  if [ ! -e ${boxmod_prevdir} ] ; then
    echo error, previous output directory could not be found
    echo ${boxmod_prevdir}
    exit 1
  fi
fi

#-- check for input key files (argument, existence)
if [ ${#keyfile} -eq 0 ] ; then
  echo "no input (key) file name supplied" ; exit 2 
fi
if [ ! -e ${boxmod_inpdir}/indat_${keyfile}.key ]; then
    echo Input file \'${boxmod_inpdir}/indat_${keyfile}.key\' doesn\'t exist
    exit
fi

#-- check for init_conc files (argument, existence)
if [ ${#confile} -ne 0 ] ; then
  echo "constrained concs file   = "init_conc_$confile.key
  if [ ! -e ${boxmod_inpdir}/init_conc_${confile}.key ]; then
    echo Initial concentrations file \'${boxmod_inpdir}/init_conc_${confile}.key\' doesn\'t exist
    exit
  fi
fi

#-- check for photolysis lookup file existence
if [ ! -e $boxmod_inpdir/$photfile ] ; then
  echo $boxmod_inpdir/$photfile
  echo "photolysis file not found" ; exit 4 
fi
echo "photfile = "$photfile

#-- check for mechanism file (argument, existence)
if [ ${#mechname} -eq 0 ] ; then
  echo "no mechanism name supplied" ; exit 3 
fi
echo "mechname = "$mechname.mech
gecko_mechdir=$gecko_outdir/$mechname
if [ ! -e ${gecko_mechdir} ] ; then
    echo error, mechanism directory could not be found
    echo ${gecko_mechdir}
    exit 1
fi
echo "mechanism location = "$gecko_mechdir

#-- check for interpreted mech file (binary = *.akli; NetCDF = *.nc)
if [ ! ${netcdf_flag} ]; then 
  if [ ! -e ${gecko_mechdir}/outdat.akli ]; then
      echo Interpreted mechanism file \'${gecko_mechdir}/outdat.akli\' doesn\'t exist
      exit
  fi
fi

if [ ${netcdf_flag} ]; then 
  if [ ! -e ${gecko_mechdir}/outdat.nc ]; then
    echo Interpreted mechanism file \'${gecko_mechdir}/outdat.nc\' doesn\'t exist
    exit
  fi
fi

#-- check for previous output file, outdat.nc
if [ -n "${prevout}" ]; then 
  if [ -e ${boxmod_prevdir}/outdat.nc ]; then
    continue
  else
    echo previous output file \'${boxmod_prevdir}/outdat.nc\' doesn\'t exist
    exit
  fi
fi

echo

#---------------- make some links -------

run_name=${mechname}_${keyfile}
if [ ${#confile} -ne 0 ] ; then
  run_name=${run_name}_${confile}
fi

boxmod_wkdir=${boxmod_outdir}/${run_name}

#if [ ! -e ${boxmod_wkdir} ] ; then
#    echo WARNING! removing existing directory
#    rm -rf ${boxmod_wkdir}
#    mkdir ${boxmod_wkdir}
#fi

if [ ! -e ${boxmod_wkdir} ] ; then
    echo Creating ${boxmod_wkdir}...
    mkdir ${boxmod_wkdir}
fi


if [ ${netcdf_flag} ]; then
  echo linking NetCDF link-file ${gecko_mechdir}/outdat.nc
  echo ... to ${boxmod_wkdir}/indat.nc
  ln -s ${gecko_mechdir}/outdat.nc ${boxmod_wkdir}/indat.nc
  echo
else
  echo linking binary link-file ${gecko_mechdir}/outdat.akli
  echo ... to ${boxmod_wkdir}/indat.li
  ln -s ${gecko_mechdir}/outdat.akli ${boxmod_wkdir}/indat.li
  echo
fi

if [ -n "${prevout}" ]; then 

#  # NetCDF
#  if [ ${netcdf_flag} ]; then
    echo linking existing output ${boxmod_prevdir}/outdat.nc
    echo ... to ${boxmod_wkdir}/prevdat.nc
    ln -sf ${boxmod_prevdir}/outdat.nc ${boxmod_wkdir}/prevdat.nc
    echo

#  else
#    # Binary particle
#    if [ -e ${boxmod_prevdir}/outdat.ppa ]; then
#      echo linking existing output ${boxmod_prevdir}/outdat.ppa
#      echo ... to ${boxmod_wkdir}/indat.ppa
#      ln -s ${boxmod_prevdir}/outdat.ppa ${boxmod_wkdir}/indat.ppa
#      echo
#    fi
#    # Binary gas
#    if [ -e ${boxmod_prevdir}/outdat.ppf ]; then
#      echo linking existing output ${boxmod_prevdir}/outdat.ppf
#      echo ... to ${boxmod_wkdir}/indat.ppf
#      ln -s ${boxmod_prevdir}/outdat.ppf ${boxmod_wkdir}/indat.ppf
#      echo
#    fi
#  fi
fi

#---------------- executable -----------------
# if no existing executable, go ahead and compile one
# using the current mechanism's akparameter.h file

#if [ ! -f ${boxmod_wkdir}'/BOXMOD' ] ; then  #(BOXMOD NOT FOUND)
#    echo "could not find boxmodel executable -> compiling"
    echo "compiling boxmodel executable"
    echo
    if [ ! -e ${boxmod_wkdir}/WORK_BOXMOD ] ; then
      mkdir ${boxmod_wkdir}/WORK_BOXMOD
      mkdir ${boxmod_wkdir}/WORK_BOXMOD/LIB
      mkdir ${boxmod_wkdir}/WORK_BOXMOD/OBJ
      mkdir ${boxmod_wkdir}/WORK_BOXMOD/PROG
    fi

    echo linking program files...
    for f in `ls ${boxmod_source}/LIB/*.[fh]`
    do
      ln -sf $f ${boxmod_wkdir}/WORK_BOXMOD/LIB/
    done
    ln -s ${boxmod_source}/LIB/*.f90 ${boxmod_wkdir}/WORK_BOXMOD/LIB/
    ln -s ${boxmod_source}/PROG/boxmod_main.f ${boxmod_wkdir}/WORK_BOXMOD/PROG/boxmod_main.f
    ln -s ${boxmod_source}/PROG/boxmod_main.f90 ${boxmod_wkdir}/WORK_BOXMOD/PROG/boxmod_main.f90

# link to the relevant akparameter.h file
    echo linking akparameter...
#    ln -sf ${gecko_mechdir}/akparameter.h ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter.h
    ln -sf ${gecko_mechdir}/akparameter_module.f90 ${boxmod_wkdir}/WORK_BOXMOD/LIB/akparameter_module.f90

    echo linking Makefile...
      ln -sf ${boxmod_source}/OBJ/Makefile ${boxmod_wkdir}/WORK_BOXMOD/OBJ/Makefile
      ln -sf ${boxmod_source}/OBJ/template.mk ${boxmod_wkdir}/WORK_BOXMOD/OBJ/template.mk

    echo Making...
    cd ${boxmod_wkdir}/WORK_BOXMOD/OBJ
    make all #BOXMOD

    if [ ! -e ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ] ; then
       echo error, could not find boxmodel executable
       echo ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD
       echo maybe check compilation?
       exit
    fi

    cp ${boxmod_wkdir}/WORK_BOXMOD/PROG/BOXMOD ${boxmod_wkdir}/BOXMOD
    echo DONE MAKE...
    ls -l ${boxmod_wkdir}/BOXMOD

## end "if executable exists"
#fi

# stop script after make (for testing)
# exit

#-------------------------------------------------------------------------
if [ ! -d .git ]; then
# Create reference file git_info.txt in output directory, to be read by codes.
# GitHub branch, commit:
  echo "** creating gitinfo file **"
  echo `git rev-parse --abbrev-ref HEAD` > ${boxmod_wkdir}/indat.gitinfo
  echo `git rev-parse HEAD`   >> ${boxmod_wkdir}/indat.gitinfo
else
# IF NOT WITHIN A GIT REPOSITORY READ INFO FROM VERSION_NOTES DIRECTORY
  echo "** copying gitinfo file **"
  cp ${home_dir}/VERSION_NOTES/README_gitinfo ${boxmod_wkdir}/indat.gitinfo
fi
#-------------------------------------------------------------------------

if [ ! ${netcdf_flag} ]; then
  echo building RO2 files...
  cd ${gecko_source}
  make compteur
  ln -s ${gecko_source}/RUN/COMPTEUR/compteur ${gecko_mechdir}/compteur
  cd ${gecko_mechdir}
  ./compteur
  rm compteur
  echo linking RO2 files...
  ls ${gecko_mechdir}

  ln -s ${gecko_mechdir}/XP1O2 ${boxmod_wkdir}/indat1.ro2
  ln -s ${gecko_mechdir}/XP2O2 ${boxmod_wkdir}/indat2.ro2
  ln -s ${gecko_mechdir}/XP3O2 ${boxmod_wkdir}/indat3.ro2
  ln -s ${gecko_mechdir}/XS1O2 ${boxmod_wkdir}/indat4.ro2
  ln -s ${gecko_mechdir}/XS2O2 ${boxmod_wkdir}/indat5.ro2
  ln -s ${gecko_mechdir}/XS3O2 ${boxmod_wkdir}/indat6.ro2
  ln -s ${gecko_mechdir}/XT1O2 ${boxmod_wkdir}/indat7.ro2
  ln -s ${gecko_mechdir}/XT2O2 ${boxmod_wkdir}/indat8.ro2
  ln -s ${gecko_mechdir}/XACO3 ${boxmod_wkdir}/indat9.ro2
  echo
fi

echo copying file ${boxmod_inpdir}/indat_${keyfile}.key...
echo to ... ${boxmod_wkdir}/indat.key
cp ${boxmod_inpdir}/indat_${keyfile}.key ${boxmod_wkdir}/indat.key
echo

#==============================================================#
# check if END statement already exists in indat.key and remove it
# then remove (all) blank lines from file

sed -i s/^END//g ${boxmod_wkdir}/indat.key
sed -i '/^$/d' ${boxmod_wkdir}/indat.key

#  for previous o/p: make sure PREV = 1 is present
if [ -n "${prevout}" ]; then 
    if grep -q "PREV" ${boxmod_wkdir}/indat.key ; then
      sed -i 's/PREV 0/PREV 1/g' ${boxmod_wkdir}/indat.key
    else
      sed -i '1 i\PREV 1' ${boxmod_wkdir}/indat.key
    fi
fi

if [ -e ${boxmod_inpdir}/init_conc_${confile}.key ]; then
  echo adding initial concentrations ${boxmod_inpdir}/init_conc_${confile}.key...
  echo to ... ${boxmod_wkdir}/indat.key
  cp ${boxmod_inpdir}/init_conc_${confile}.key ${boxmod_wkdir}/init_conc.key
  cat ${boxmod_wkdir}/init_conc.key >> ${boxmod_wkdir}/indat.key
  echo
fi

if [ -e ${boxmod_inpdir}/init_steadystate_spinup_${keyfile}.key ] ; then
  echo adding stead-state spinup ${boxmod_inpdir}/init_steadystate_spinup_${keyfile}.key...
  echo to ... ${boxmod_wkdir}/indat.key
    cp ${boxmod_inpdir}/init_steadystate_spinup_${keyfile}.key ${boxmod_wkdir}/init_steadystate.key
    cat ${boxmod_wkdir}/init_steadystate.key >> ${boxmod_wkdir}/indat.key
  echo
fi

if [ -e ${boxmod_inpdir}/*_${keyfile}.input ] ; then
  echo linking user external constraints input files
  for f in `ls ${boxmod_inpdir}/*_${keyfile}.input`
  do
    echo linking $f...
    ln -s $f ${boxmod_wkdir}/
  done
  echo
fi

echo END >> ${boxmod_wkdir}/indat.key
 
#--link input files
echo linking photolysis file ${boxmod_inpdir}/${photfile}
echo to ... ${boxmod_wkdir}/jfile.phot
ln -s ${boxmod_inpdir}/${photfile} ${boxmod_wkdir}/jfile.phot
echo
echo linking diffusion volume file ${gecko_mechdir}/${mechname}.difv ...
echo to ... ${boxmod_wkdir}/difv.dat
ln -s ${gecko_mechdir}/${mechname}.difv ${boxmod_wkdir}/difv.dat
echo
echo linking Nannoolal pvap file ${gecko_mechdir}/${mechname}.pnan ...
echo to ... ${boxmod_wkdir}/pnan.sat
ln -s ${gecko_mechdir}/${mechname}.pnan ${boxmod_wkdir}/pnan.sat
echo
echo linking SIMPOL pvap file ${gecko_mechdir}/${mechname}.psim ...
echo to ... ${boxmod_wkdir}/psim.sat
ln -s ${gecko_mechdir}/${mechname}.psim ${boxmod_wkdir}/psim.sat
echo
echo linking SIMPOL data file ${gecko_data}/simpol.dat ...
echo to ... ${boxmod_wkdir}/simpol.dat
ln -s ${gecko_data}/simpol.dat ${boxmod_wkdir}/simpol.dat
echo
echo linking Henry file ${gecko_mechdir}/${mechname}.Henry ...
echo to ... ${boxmod_wkdir}/Henry.dat
ln -s ${gecko_mechdir}/${mechname}.Henry ${boxmod_wkdir}/Henry.dat
echo
#echo linking gitinfo file ${gecko_mechdir}/${mechname}.gitinfo ...
#echo to ... ${boxmod_wkdir}/indat.gitinfo
#ln -s ${gecko_mechdir}/${mechname}.gitinfo ${boxmod_wkdir}/indat.gitinfo
#echo
 
#==============================================================#
# -------------- logfile output ---------------
#==============================================================#

cd ${home_dir}

 echo '================================'  > $boxmod_wkdir/outdat.README
 echo '= GECKO-A box model simulation =' >> $boxmod_wkdir/outdat.README
 echo '================================' >> $boxmod_wkdir/outdat.README
 echo 'run name      = '${run_name}      >> $boxmod_wkdir/outdat.README
 echo 'user          = '$USER            >> $boxmod_wkdir/outdat.README
 echo 'date          = '`date`           >> $boxmod_wkdir/outdat.README
 echo 'machine       = '$HOSTNAME        >> $boxmod_wkdir/outdat.README
 echo 'model version :'                  >> $boxmod_wkdir/outdat.README
 echo 'GitHub branch = '`git rev-parse --abbrev-ref HEAD` >> $boxmod_wkdir/outdat.README
 echo 'GitHub commit = '`git rev-parse HEAD`   >> $boxmod_wkdir/outdat.README
 echo 'model inputs  = indat_'$keyfile'.key' >> $boxmod_wkdir/outdat.README
 if [ ${#confile} -ne 0 ] ; then
 echo 'constraints   = '$confile         >> $boxmod_wkdir/outdat.README
 else
 echo 'constraints   = (no file)'        >> $boxmod_wkdir/outdat.README
 fi
 echo 'dictionary    = '$mechname        >> $boxmod_wkdir/outdat.README
 echo 'phot file     = '$photfile        >> $boxmod_wkdir/outdat.README
 if [ -n "${prevout}" ]; then 
 echo 'continuing from existing output : '$prevout >> $boxmod_wkdir/outdat.README
 else
 echo 'new run, no previous output '     >> $boxmod_wkdir/outdat.README
 fi
 echo '================================' >> $boxmod_wkdir/outdat.README
 echo '(add user comments manually below)'>> $boxmod_wkdir/outdat.README

#==============================================================#
# remove (all) blank lines from indat.key
sed -i '/^$/d' ${boxmod_wkdir}/indat.key

cd ${boxmod_wkdir}
echo ------------------------------------------
echo ! ! We are now in ${boxmod_wkdir} ! !
echo       ... running box model ...
echo ------------------------------------------

## free up memory to avoid initial segfaults:
ulimit -s unlimited

## run box model:
## option 1) standard run, output to terminal
OMP_NUM_THREADS=${nthreads} ./BOXMOD

## option 2) "nohup" = run in background, safer if connection times out
#OMP_NUM_THREADS=${nthreads} nohup ./BOXMOD

# output results filenames to screen
# binary output file
if [ -s ${boxmod_wkdir}/outdat.ppf ] ; then
  echo check output ${boxmod_wkdir}/outdat.ppf
fi
# netCDF output file
if [ -s ${boxmod_wkdir}/outdat.nc ] ; then
  echo check output ${boxmod_wkdir}/outdat.nc
fi
# no output file!!
if [ ! -e ${boxmod_wkdir}/outdat.ppf ] && [ ! -e ${boxmod_wkdir}/outdat.nc ]; then
  echo no output: run failed
fi
# logging file with run-time output
echo check logfile ${boxmod_wkdir}/outdat.README

if [ ${spinup_fg} -eq 1 ] && [ -e ${boxmod_wkdir}/steadystate.key ] ; then
  echo creating link to steadystate file
  cp ${boxmod_wkdir}/steadystate.key ${boxmod_inpdir}/init_steadystate_spinup_${keyfile}.key
fi

#-----------
#exit 99
#-----------
## TIDY UP ##

cat outdat.README outdat.out fort.11 > tmp.txt
mv tmp.txt outdat.README

#rm fort.*
rm dummy.*
rm outdat.out

# delete executable
rm BOXMOD

## delete work directory 
## comment this out for testing, if you want to avoid re-compiling every time
rm -R WORK_BOXMOD

# delete symbolic links
# delete zero-length files
# delete temporary files
find -lname '*' -delete
find ./ -type f -size 0 -delete
rm indat.gitinfo

# move output into subdirectory
mkdir USEROUT
mv *.* USEROUT
cd USEROUT

#-----------
# ... and move standard output back to main run directory
#     (leaving non-standard output in subdirectory)
mv outdat.nc ../
mv outdat.README ../
mv indat.key ../
#-----------






