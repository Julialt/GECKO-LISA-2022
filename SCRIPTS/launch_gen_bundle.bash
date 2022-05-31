#!/bin/bash
#----------------------------
# launch_gen_bundle.bash
#----------------------------
# this script must be launched from the $gecko_wkdir folder
# launch to cheyenne with arguments:
# inputs: $1 = gecko directory
# $2 = working directory inside mech dir
# $4 = mechname
# $5 = $numpre
# $5 = prevflag

gecko_source=$1
gecko_wkdir=$2
mechname=$3
numpre=$4
prevflag=$5

ref_dir=${gecko_wkdir}/spinup 

#----------------------------------#
# CHECK IF compiled codes exist
# YES -> link to them
# NO -> compile the generator
# re: compile_gecko.bash
#----------------------------------#

if [ -e ${ref_dir}/cm ] ; then
  ln -s ${ref_dir}/cm ${gecko_wkdir}/cm
else
  mkdir ${gecko_wkdir}/GECKO_CODE
  cp -r ${gecko_source}/* ${gecko_wkdir}/GECKO_CODE
  cd ${gecko_wkdir}/GECKO_CODE
  make clean
  make cm
fi

# the following assumes pre-compiled compteur and findname_cmv exist

ln -s ${gecko_source}/RUN/COMPTEUR/compteur ${gecko_wkdir}/compteur
ln -s ${gecko_source}/RUN/findname_cmv ${gecko_wkdir}/findname_cmv

#----------------------------------#
# link files etc 
# re: generate_scheme_cheyenne.bash
#----------------------------------#

ln -s ${gecko_wkdir}/GECKO_CODE/RUN/cm ${gecko_wkdir}/cm

#----------------------------------#
# generate the scheme
# re: generate_scheme_cheyenne.bash
#----------------------------------#
cd ${gecko_wkdir}

## single-submission way ...
# (edit 3 lines out out for multi-submission)
counter=0
until [ $counter -eq $numpre ]; do
  let counter=$counter+1

#----- either single- or multi-submission
echo starting $counter:$numpre

echo constructing management file
echo " &manage" > manage.input
# flag to read previous output : all except first iteration
if [ $counter -eq 1 ] && [ "$prevflag" != "y" ] ; then
  echo "     prevflag = 0" >> manage.input
else
  echo "     prevflag = 1" >> manage.input
  # make sure existing dict and mech are available
  if [ ! -e existing.dict ] ; then
    if [ ! -e fort.7 ] ; then
      cp fort.7 existing.dict
    else
      echo "no previous dictionary available => exiting"
      exit
    fi
  fi
  if [ ! -e existing.rxns ] ; then
    echo "no previous reaction list available => exiting"
    exit
  fi
# remove 'END' from pre-existing net reaction rate files
  if [ -e ./${mechname}.kOH ] ; then 
    nend=$(grep -c END ./${mechname}.kOH | awk '{ VAR += $1} END {print VAR}')
    if [ $nend > 0 ] ; then 
      head -n -1 ./${mechname}.kOH > tmp.txt 
      mv tmp.txt ./${mechname}.kOH 
    fi 
  fi

  if [ -e ./${mechname}.kO3 ] ; then 
    nend=$(grep -c END ./${mechname}.kO3 | awk '{ VAR += $1} END {print VAR}')
    if [ $nend > 0 ] ; then 
      head -n -1 ./${mechname}.kO3 > tmp.txt 
      mv tmp.txt ./${mechname}.kO3 
    fi 
  fi

  if [ -e ./${mechname}.kNO3 ] ; then 
     nend=$(grep -c END ./${mechname}.kNO3 | awk '{ VAR += $1} END {print VAR}')
    if [ $nend > 0 ] ; then 
       head -n -1 ./${mechname}.kNO3 > tmp.txt 
       mv tmp.txt ./${mechname}.kNO3 
    fi 
  fi

fi

# flag for post-processing : final iteration only
if [ $counter -ne $numpre ] ; then
  echo "     postflag = 0" >> manage.input
else
  echo "     postflag = 1" >> manage.input
fi
echo " /" >> manage.input

# run generator
./cm

if [ $? -eq 0 ] ; then
  echo generator has run successfully
else
  echo Mechanism generation failed
  exit
fi

# create continuation files for next iteration

if [ ! -e existing.rxns ] ; then
  cat fort.21 fort.17 > ${mechname}.mech
  mv fort.17 existing.rxns
else
  cat existing.rxns fort.17 > tmp.txt ; mv tmp.txt existing.rxns
  cat fort.21 existing.rxns > ${mechname}.mech
fi
cp fort.7 existing.dict

#== if existing.rxns ends with 'END', strip the last line
#-- (except for the last time around)

if [ $counter -ne $numpre ] ; then
  nend=$(grep -c END ./existing.rxns | awk '{ VAR += $1} END {print VAR}')
  if [ $nend > 0 ] ; then
    head -n -1  ./existing.rxns > tmp.txt
    mv tmp.txt ./existing.rxns
  fi
fi

# accumulate list of lumped species
  if [ -e fort.48 ] ; then
    if [ -e $mechname.lump ] ; then
      cat $mechname.lump fort.48 > tmp.txt ; mv tmp.txt $mechname.lump
    else
      mv fort.48 $mechname.lump ; fi
  fi

# accumulate net reaction rate lists
  if [ -e fort.70 ] ; then
    if [ -e $mechname.kO3 ] ; then
      cat $mechname.kO3 fort.70 > tmp.txt ; mv tmp.txt $mechname.kO3
    else
      mv fort.70 $mechname.kO3 ; fi
  fi

  if [ -e fort.71 ] ; then
    if [ -e $mechname.kNO3 ] ; then
      cat $mechname.kNO3 fort.71 > tmp.txt ; mv tmp.txt $mechname.kNO3
    else
      mv fort.71 $mechname.kNO3 ; fi
  fi

  if [ -e fort.72 ] ; then
    if [ -e $mechname.kOH ] ; then
      cat $mechname.kOH fort.72 > tmp.txt ; mv tmp.txt $mechname.kOH
    else
      mv fort.72 $mechname.kOH ; fi
  fi

#-- end run loop (single submission: else comment out 3 lines)
  echo done $counter:$numpre
  echo ---------------------------------------
done

# last time around: add 'END' to net reaction rate files
if [ $counter -eq $numpre ] ; then
    echo 'END' >> $mechname.kOH
    echo 'END' >> $mechname.kO3
    echo 'END' >> $mechname.kNO3
fi

#----------------------------------#
# post-process the scheme
# re: postprocess_generator.bash
#----------------------------------#
home_dir=${gecko_source}/..
netcdf_flag="yes"
#----------------------------------#

# supercomputer-specific stuff here
source $home_dir/SCRIPTS/cheyenne_scripting_functions.bash
#--------------------
echo 'rename files for mechanism ...'

mv difvol.dat $mechname.difv
mv pnan.dat   $mechname.pnan
mv psim.dat   $mechname.psim
mv dHeatf.dat $mechname.dHeatf
mv henry.dat  $mechname.Henry
if [ -e existing.log ] ; then
  cat existing.log scheme.log > $mechname.log
  rm scheme.log
  rm existing.log
else
  mv scheme.log $mechname.log
fi
mv fort.7 $mechname.dict

if [ ! -e $mechname.mech ] ; then
  mv fort.17 $mechname.mech
fi

if [ -e fort.48 ] ; then
  if [ ! -e $mechname.lump ] ; then
    sort -u fort.48 > $mechname.lump
  else
    sort -u fort.48 >> $mechname.lump
  fi
fi


#cd ${gecko_wkdir}
#!echo Counting RO2s...
#./compteur

#echo Finding primary species codes...
#./findname_cmv > findname_output

mv userparams.input userdat.settings

mv lstprim.out      userdat.cheminput
sed 's/^ *//g' -i userdat.cheminput
sed 's/ $*//g' -i userdat.cheminput
echo 'END' >> userdat.cheminput

echo tidy up ...

rm fort.*
rm dummy.*
rm sortlist.*
rm manage.input
rm existing.*

ln -s $mechname.dict fort.7
ln -s $mechname.mech fort.17

#==========================================
# RUN AKPARAMETER.H GENERATION
#==========================================

#cd ${home_dir}/SCRIPTS
#echo -----------------------------------------
#echo run akparameter generation
#./write_akparameter.bash ${mechname}

#===========================================
# RUN INTERPRETER SCRIPT in SCRIPTS > $gecko_wkdir
#===========================================

#echo ------------------------------------------
#echo We are now in directory ...... ; pwd
#echo ------------------------------------------
#
#echo run interpreter...
#./run_interp_automatic.bash ${mechname}
#
#cd $home_dir/SCRIPTS

#==========================================
# RUN PACKAGING ROUTINE
#==========================================

cd ${home_dir}/SCRIPTS
echo -----------------------------------------
echo package mechanism for box model....
# local version of command
#./gen_package_local.bash ${pathfile} ${mech}
# cheyenne version of command
#./gen_package_cheyenne.bash ${gecko_source} ${gecko_wkdir} &{mechname} &{boxmod_source}
./gen_package_cheyenne.bash ${gecko_source} ${gecko_wkdir} ${mechname} ${boxmod_source}
