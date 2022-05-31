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
    s ) settings=$OPTARG ;;
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

echo Compiling GECKO-A...
cd ${gecko_source}
echo .... in directory ...... ; pwd
make cm

#== directory for generated output ==

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

cd ${gecko_wkdir}
echo ------------------------------------------
echo We are now in directory ...... ; pwd
echo ------------------------------------------
#-------------------------------------------------------------------------
if [ ! -d .git ]; then 
# Create reference file git_info.txt in output directory, to be read by codes.
# GitHub branch, commit:
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
    if [ -e fort.7 ] ; then
      mv fort.7 existing.dict
    fi
    if [ ! -e existing.dict ] ; then
      echo "no previous dictionary available => exiting"
      exit
    fi
    if [ ! -e existing.rxns ] ; then
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
  if [ -e fort.7 ] ; then
	echo generator has run successfully
  else
	echo Mechanism generation failed
	exit 17
  fi
#================================================#
# create continuation files for next iteration
#-- existing.rxns

  if [ ! -e existing.rxns ] ; then
    mv fort.17 existing.rxns
  else

#== if existing.rxns ends with 'END', strip the last line
#-- (except for the last time around)
#  if [ $counter -ne $numpre ] ; then
    nend=$(grep -c END ${gecko_wkdir}/existing.rxns | awk '{ VAR += $1} END {print VAR}')
    if [ $nend > 0 ] ; then
      head -n -1  ${gecko_wkdir}/existing.rxns > tmp.txt
      mv tmp.txt ${gecko_wkdir}/existing.rxns
      rm 0
    fi
#  fi

    cat existing.rxns fort.17 > tmp.txt ; mv tmp.txt existing.rxns
  fi

  cat fort.21 existing.rxns > ${mech}.mech

#================================================#
# concatenate and rename some output files:
# list of lumped species
  if [ -e fort.48 ] ; then
    if [ -e $mech.lump ] ; then
      cat $mech.lump fort.48 > tmp.txt ; mv tmp.txt $mech.lump
    else
      mv fort.48 $mech.lump ; fi
  fi

# net reaction rates
  fns=(70 71 72)
  for n in ${!fns[*]}; do
   fn=${fns[n]}
   case "$fn" in 
    70) 
       ksp=kO3 ;;
    71) 
       ksp=kNO3 ;;
    72) 
       ksp=kOH ;;
   esac

   if [ -e fort.$fn ] ; then
    if [ -e $mech'.'$ksp ] ; then
      nlin=$(wc -l $mech'.'$ksp |awk '{ VAR =+ $1} END {print VAR}')
      nend=$(grep -nm1 END $mech'.'$ksp |awk '{ VAR =+ $1} END {print VAR}')
      if [ $nend > 0 ] ; then
      let nlin=$nend-1
      fi
      head -n $nlin $mech'.'$ksp > $mech'.tmp'
      cat $mech'.tmp' 'fort.'$fn > tmp.txt ; mv tmp.txt $mech'.'$ksp
    else
      mv 'fort.'$fn $mech'.'$ksp ; fi
   fi
  done

  if [ -e warnings.out ] ; then
    if [ -e fort.50 ] ; then
      if [ -e $mech.warnings ] ; then
        cat $mech.warnings warnings.out fort.50 > tmp.out 
        mv tmp.out $mech.warnings
      else
        cat warnings.out fort.50 >  $mech.warnings
      fi
    else
      if [ -e $mech.warnings ] ; then
        cat $mech.warnings warnings.out > tmp.out 
        mv tmp.out $mech.warnings
      else
        mv warnings.out $mech.warnings
      fi
    fi
  fi

#-- end run loop
  echo done $counter:$numpre
  echo ---------------------------------------
done

#================================================#

if [ -e existing.cheminput ] ; then
  echo '  --precursors from previously existing mechanism--' 
  echo '  --precursors from previously existing mechanism--' > tmp.txt
  cat scheme.log tmp.txt existing.cheminput > tmp.log
  mv tmp.log scheme.log
  rm tmp.txt
fi

cd ${gecko_wkdir}
echo ''
echo '------------------------'
echo rename output files  ...
echo '------------------------'

mv difvol.dat  $mech.difv
mv pnan.dat    $mech.pnan
mv psim.dat    $mech.psim
mv dHeatf.dat  $mech.dHeatf
mv henry.dat   $mech.Henry
mv scheme.log  $mech.log
mv fort.7 $mech.dict
if [ -e fort.48 ] ; then
  mv fort.48 $mech.lump ; fi

mv lstprim.out      userdat.cheminput
mv userparams.input userdat.settings

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

# add 'END' to some output files
echo 'END' >> $mech.kOH 
echo 'END' >> $mech.kO3 
echo 'END' >> $mech.kNO3 
echo 'END' >> userdat.cheminput

echo ''
echo '------------------------'
echo tidy up ...
echo '------------------------'
rm cm

## Optional for testing
#echo 'package disabled: edit script to reinstate'
#exit
## End Optional

rm fort.*
rm sortlist.*
rm manage.input
#rm existing.*
rm cheminput.dat # has been copied to userdat.cheminput
if [ -e chaname.dat ] ; then 
 rm chaname.dat ; fi
rm ${gecko_outdir}/DATA
if [ -e warnings.out ] ; then 
 rm warnings.out  # has been copied to ${mech}.warnings
fi

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
