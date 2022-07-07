#==SCRIPT TO GENERATE MECHANISM FROM MULTIPLE PRECURSORS==
#===========================================
#--intended to run from scratch in 3 stages:
#--1) 1st precursor
#--2) more precursors, sequentially, with output each time
#--3) post-process: create pvap & Henry files; add W, A reactions
#=========================================================
#== READ IN GECKO VERSION AND PATHS FROM FILE setup.dat ==
#=========================================================

#-- default values for some arguments ==:
pathfile="setup.dat"
settings="settings_default"
cheminput="cheminput.dat"

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
# p = phot file (optional: default is sophie_0404_julia_cmv)
# v = previous output (created with same mechanism)
#== dir paths (optional: obtainable from setup.dat)
# b = box model source path
# h = home directory
# w = gecko working dir
#--
# z = dummy argument to make sure all desired args are accepted
#--
while getopts "a:b:c:d:h:i:k:m:p:s:v:w:x:z" opt
do
  case "$opt" in
# general files
    a ) pathfile=$OPTARG ;;
    d ) depend=$OPTARG ;;
    m ) mech=$OPTARG ;;
# generator-specific files
    i ) cheminput=$OPTARG ;;
    s ) settingsfile=$OPTARG ;;
    x ) existing=$OPTARG ;;
# boxmodel-specific files
    c ) echo "ERROR! -c is a boxmod-specific flag!"; exit ;;
    k ) echo "ERROR! -k is a boxmod-specific flag!"; exit ;;
    p ) echo "ERROR! -p is a boxmod-specific flag!"; exit ;;
    v ) echo "ERROR! -v is a boxmod-specific flag!"; exit ;;
# any other code produces an error
    * ) echo "ERROR! flag not recognized !"; exit ;;
  esac
done

#====================================
#== SET UP GECKO VERSION AND PATHS ==
#====================================

echo pathfile = $pathfile
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

# z = dummy argument to make sure all desired args are accepted
home_dir=`     grep "home_dir"      ${pathfile} | awk '{print $3}' `
gecko_version=`grep "gecko_version" ${pathfile} | awk '{print $3}' `
gecko_inp_dir=`grep "gecko_inp_dir" ${pathfile} | awk '{print $3}' `
gecko_run_dir=`grep "gecko_run_dir" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
gecko_source=$home_dir/$gecko_version
gecko_inputs=$home_dir/$gecko_inp_dir
gecko_outdir=$home_dir/$gecko_run_dir

echo "GECKO source = "$gecko_source
echo "input path   = "$gecko_inputs
echo "output path  = "$gecko_outdir
echo "cheminput = "$cheminput
echo "mech = "$mech

if [ -z ${mech} ] ; then
    echo no mech supplied. Aborting.
    exit
fi
#===========================================
#== SET UP PRECURSORS & MECH NAME ==

echo "existing = "$existing
if [ -z ${existing} ] ; then
    echo existing is zero length
    prevflag="n"
else
    echo existing is non zero length
    prevflag="y"
    existing_mech=${existing}
fi
echo "prevflag=" $prevflag

#===========================================
# DIRECTORY CHECKS AND LINKS
#===========================================

if [ ! -e ${gecko_source} ] ; then
	echo error, gecko source directory could not be found
	echo ${gecko_source}
	echo check gecko_source variable in $0
	exit ; fi

if [ ! -e ${gecko_outdir} ] ; then
	echo error, gecko run directory could not be found
	echo ${gecko_outdir}
	echo check gecko_outdir variable in $0
	exit ; fi

if [ $prevflag=="y" ] && [ ! -e ${gecko_outdir}/${existing} ]; then
	echo existing mechanism \'${gecko_outdir}/${existing}\' doesn\'t exist
	exit ; fi

if [ ! -e ${gecko_inputs}/${settings} ]; then
	echo Settings file \'${gecko_inputs}/${settings}\' doesn\'t exist
	exit ; fi

if [ ! -e ${gecko_inputs}/${cheminput} ]; then
	echo Cheminput file \'${gecko_inputs}/${cheminput}\' doesn\'t exist
	exit ; fi

#===========================================
# RUN THE GENERATOR
#===========================================
# make link to settings file BEFORE compiling
echo linking to settings file...
cd ${gecko_dir}/LIB
ln -s ${gecko_dir}/INPUTS/$settingsfile keyflag.f90

# PARIS version compiles in OBJ and doesn't require executable name
cd ${gecko_source}/OBJ
#echo NOT Compiling GECKO-A...
echo Compiling GECKO-A...
echo .... in directory ...... ; pwd
#make clean
make 

#== directory for generated output ==
#== make links to executable and PARIS structure ==
gecko_wkdir=${gecko_outdir}/${mech}

if [ "$prevflag" == "y" ] ; then
  gecko_prevdir=${gecko_outdir}/${existing}
  if [ ! -e ${gecko_prevdir} ] ; then
	echo WARNING! previous output dir ${gecko_prevdir} 
        echo          does NOT exist. Stopping !!
        exit 99
  elif [ ${gecko_prevdir} == ${gecko_wkdir} ] ;then
	echo adding to previous output in same directory
        echo ${gecko_wkdir}
  fi
elif [ -e ${gecko_wkdir} ] ; then
  echo CAUTION: OVERWRITING DIRECTORY ${gecko_wkdir}
  rm -R ${gecko_wkdir} 
#  rm -R ${gecko_outdir}/tmp
#  mv ${gecko_wkdir} ${gecko_outdir}/tmp
fi
mkdir ${gecko_wkdir} 
# Build temporary working dir structure mirroring PARIS RUN dir structure 
mkdir ${gecko_wkdir}/OUT

cd ${gecko_wkdir}
echo ------------------------------------------
echo We are now in directory ...... ; pwd
echo ------------------------------------------
#-------------------------------------------------------------------------
# Create reference file git_info.txt in output directory, to be read by codes.
# GitHub branch, commit:
if [ ! -d .git ]; then 
  echo "** creating gitinfo file **"
  echo `git rev-parse --abbrev-ref HEAD` > ${mech}.gitinfo
  echo `git describe`   >> ${mech}.gitinfo
else
# IF NOT WITHIN A GIT REPOSITORY READ INFO FROM VERSION_NOTES DIRECTORY
  echo "** copying gitinfo file **"
  cp ${gecko_source}/VERSION_NOTES/README_gitinfo ${mech}.gitinfo
fi
#-------------------------------------------------------------------------

echo linking executable files...
ln -s ${gecko_source}/RUN/cm ${gecko_wkdir}/cm
if [ ! -e ${gecko_outdir}/DATA ]; then
  ln -s ${gecko_source}/DATA ${gecko_outdir}/DATA 
fi

echo Preparing input files...
cp ${gecko_inputs}/${settings}  ${gecko_wkdir}/userparams.input
echo  ${gecko_wkdir}/userparams.input
cp ${gecko_inputs}/${cheminput}  ${gecko_wkdir}/cheminput.dat
echo  ${gecko_wkdir}/cheminput.dat

#--Preparing cheminput.dat :
#  delete lines after first 'END' from working copy
#  and find # of precursors in file (non-blank lines not starting with *)

# find FIRST line with "<newline>END"
nend=$(grep -nm1 "^END" cheminput.dat |awk '{ VAR =+ $1} END {print VAR}')
head -n $nend cheminput.dat > tmp.txt ; mv tmp.txt cheminput.dat
numpre=$(grep -c "^[^*E]" cheminput.dat)
echo "File cheminput.dat contains "$numpre" precursors"

echo "prevflag = "$prevflag
if [ "$prevflag" == "y" ] ; then

 echo " Copying existing cheminput > existing.cheminput"
  cp ${gecko_outdir}/${existing}/userdat.cheminput existing.cheminput

 echo " Copying existing dictionary > existing.dict"
  cp ${gecko_outdir}/${existing}/${existing}.dict existing.dict

 echo "Copying existing mechanism > existing.rxns"
  nend1=$(grep -nm1 END ${gecko_outdir}/${existing}/${existing}.mech |awk '{ VAR =+ $1} END {print VAR}')
  nend2=$(wc -l ${gecko_outdir}/${existing}/${existing}.mech |awk '{ VAR =+ $1} END {print VAR}')
  let nlin=$nend2-$nend1-3
  tail -n $nlin ${gecko_outdir}/${existing}/${existing}.mech > existing.rxns

fi

#==LOOP RUNS OVER PRECURSORS

counter=0
until [ $counter -eq $numpre ]; do
  let counter=$counter+1
  echo starting $counter:$numpre

  echo constructing management file
  echo " &manage" > manage.input
# flag to read previous output : all except first iteration
  if [ $counter -eq 1 ] && [ "$prevflag" != "y" ] ; then
    echo "     prevflag = 0" >> manage.input
  else
    echo "     prevflag = 1" >> manage.input
    # make sure existing dict and mech are available
    if [ -e OUT/dictionary.out ] ; then
      mv OUT/dictionary.out existing.dict
    fi
    if [ ! -e existing.dict ] ; then
      echo "no previous dictionary available => exiting"
      exit
    fi
    if [ ! -e $mech.mech ] ; then
      echo "no previous reaction list available => exiting"
      exit
    fi
  fi
# flag for post-processing : final iteration only
  if [ $counter -ne $numpre ] ; then
    echo "     postflag = 0" >> manage.input
  else
    echo "     postflag = 1" >> manage.input
  fi
  echo " /" >> manage.input

#================================================#
# RUN GENERATOR
#================================================#

# DEBUG : use this version of execution statment #:
#  ./cm > out.txt 
# nohup ./cm > out.txt
# END DEBUG #

  ./cm 
#  nohup ./cm 


# check that a temporary dictionary file generated in this iteration
  if [ -e OUT/dictionary.out ] ; then
	echo generator has run successfully
  else
	echo Mechanism generation failed
	exit 17
  fi
#================================================#
# create continuation files for next iteration
#-- existing.rxns

  if [ ! -e $mech.mech ] ; then
    echo "no existing rxns"
    cp OUT/reactions.dum existing.rxns
    #--
    cp OUT/reactionswithcom.dum existing.rxns.notes
  else
    echo "yes existing rxns"

#== if existing.rxns ends with 'END', strip the last line
#-- (except for the last time around)

    if [ $counter -ne $numpre ] ; then
      nend=$(grep -c END ${gecko_wkdir}/existing.rxns | awk '{ VAR += $1} END {print VAR}')
      if [ $nend > 0 ] ; then
        head -n -1  ${gecko_wkdir}/existing.rxns > tmp.txt
        mv tmp.txt ${gecko_wkdir}/existing.rxns
        rm 0
      fi
      #--
      nend=$(grep -c END ${gecko_wkdir}/existing.rxns.notes | awk '{ VAR += $1} END {print VAR}')
      if [ $nend > 0 ] ; then
        head -n -1  ${gecko_wkdir}/existing.rxns.notes > tmp.txt
        mv tmp.txt ${gecko_wkdir}/existing.rxns.notes
        rm 0
      fi
#== and strip header 4 lines 
#== and concatenate reactions.dum onto existing.rxns
#-- (except for the first time around)
      if [ $counter -ne 1 ] ; then
        tail -n -4  ${gecko_wkdir}/OUT/reactions.dum > tmp.txt
        mv tmp.txt ${gecko_wkdir}/OUT/reactions.dum
        cat existing.rxns OUT/reactions.dum > tmp.txt ; mv tmp.txt existing.rxns
        #--
        tail -n -4  ${gecko_wkdir}/OUT/reactionswithcom.dum > tmp.txt
        mv tmp.txt ${gecko_wkdir}/OUT/reactionswithcom.dum
        cat existing.rxns.notes OUT/reactionswithcom.dum > tmp.txt ; mv tmp.txt existing.rxns.notes
      fi
    fi
  fi
  
  mv existing.rxns $mech.mech
  mv existing.rxns.notes $mech.mech.notes

#================================================#
# concatenate and rename some output files:
#-------------------------
# list of lumped species (not in current PARIS output)
  if [ -e fort.48 ] ; then
    if [ -e $mech.lump ] ; then
      cat $mech.lump fort.48 > tmp.txt ; mv tmp.txt $mech.lump
    else
      mv fort.48 $mech.lump 
    fi
  fi

#-------------------------
# net reaction rates
# (new format)
  if [ -e OUT/kicovi.dat ] ; then
    cp OUT/kicovi.dat $mech.kOH
  fi
  if [ -e OUT/kicovj.dat ] ; then
    cp OUT/kicovj.dat $mech.kNO3
  fi

#-------------------------
# warnings log
  if [ -e OUT/warning.out ] ; then
    if [ -e $mech.warnings ] ; then
      cat $mech.warnings OUT/warning.out > tmp.out 
      mv tmp.out $mech.warnings
    else
      cp OUT/warning.out $mech.warnings
    fi
  fi

#-------------------------
#-- end run loop
  echo done $counter:$numpre
  echo ---------------------------------------
done

#================================================#
# step into mech dir to rearrange files
#cd ${gecko_wkdir}
#---------------------------------------

if [ -e existing.cheminput ] ; then
  echo '  --precursors from previously existing mechanism--' 
  echo '  --precursors from previously existing mechanism--' > tmp.txt
  cat OUT/scheme.log tmp.txt existing.cheminput > tmp.log
  mv tmp.log scheme.log
  rm tmp.txt
fi


# NB: may (later) wish to retain OUT directory to keep things neater in $mech
echo ''
echo '------------------------'
echo  rename output files  ...
echo '------------------------'

cp OUT/dictionary.out $mech.dict
cp OUT/dHeatf.dat     $mech.dHeatf
cp OUT/difvol.dat     $mech.difv
# question: do we want henry.dat or henry.dep?
cp OUT/henry.dat      $mech.henry
cp OUT/henry_gromhe.dat $mech.henry_gromhe
cp OUT/maxyield.dat   $mech.maxyield
cp OUT/pvap.myr.dat   $mech.pmyr
cp OUT/pvap.nan.dat   $mech.pnan
cp OUT/pvap.sim.dat   $mech.psim
cp OUT/scheme.log     $mech.log
cp OUT/size.dum       $mech.size
cp OUT/Tg.dat         $mech.Tg
cp OUT/gasspe.dum     $mech.sp_gas
cp OUT/partspe.dum    $mech.sp_part
cp OUT/wallspe.dum    $mech.sp_wall
cp OUT/pero1.dat      $mech.pero1
cp OUT/pero2.dat      $mech.pero2
cp OUT/pero3.dat      $mech.pero3
cp OUT/pero4.dat      $mech.pero4
cp OUT/pero5.dat      $mech.pero5
cp OUT/pero6.dat      $mech.pero6
cp OUT/pero7.dat      $mech.pero7
cp OUT/pero8.dat      $mech.pero8
cp OUT/pero9.dat      $mech.pero9

if [ -e fort.48 ] ; then
  cp fort.48 $mech.lump ; fi

cat OUT/listprimary.dat $mech.prec > tmp.txt
mv tmp.txt  $mech.prec

# move inputs to "userdat" files
mv userparams.input userdat.settings
mv cheminput.dat    userdat.cheminput

# remove END statement  and blank lines from existing.cheminput
sed -i 's/^END//g' existing.cheminput
sed -i '/^$/d' existing.cheminput
# and add to START of userdat.cheminput
cat existing.cheminput userdat.cheminput > tmp.txt
mv tmp.txt userdat.cheminput

# remove leading blanks from cheminput file 
# so it can be reused 
sed 's/^ *//g' -i userdat.cheminput
sed 's/ $*//g' -i userdat.cheminput

echo ''
echo '------------------------'
echo tidy up ...
echo '------------------------'

rm cm
rm manage.input
rm ${gecko_outdir}/DATA

## Optional for testing
#echo 'package disabled: edit script to reinstate'
#exit
## End Optional

#===========================================
# INVOKE GENERATOR PACKAGING ROUTINES
#===========================================

echo ''
echo '------------------------'
echo 'package mechanism for box model....'

## Optional Interactive version
#read -rsp $'   (press any key to continue)' -n1 key
#echo 
## End Interactive version

echo '------------------------'

cd ${home_dir}/SCRIPTS
./gen_package_local.bash ${pathfile} ${mech}

exit
