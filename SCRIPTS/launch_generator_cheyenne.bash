#!/bin/bash

# this script must be launched from the $gecko_wkdir folder
# launch to cheyenne with arguments:
# $1 = $counter
# $2 = $numpre

counter=$1
numpre=$2
mechname=$3
prevflag=$4

#---------------------------------
#compile if required (does not work yet)
#------------------------------

#./compile_gecko.bash
#echo linking executable files...

#---------------------------------
# loop through required species
#---------------------------------
echo starting $counter:$numpre

echo constructing management file
echo " &manage" > manage.input
# flag to read previous output : all except first iteration
if [ $counter -eq 1 ] && [ "$prevflag" != "y" ] ; then
  echo "     prevflag = 0" >> manage.input
else
  echo "     prevflag = 1" >> manage.input
  # make sure existing dict and mech are available
  if [  -e fort.7 ] ; then
     cp fort.7 existing.dict
  fi
  if [ ! -e existing.dict ] ; then
     echo "no previous dictionary available => exiting"
     exit
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

echo "------------------------------------"
echo "list of all output files immediately after run"
ls 
echo "------------------------------------"

# check that a temporary dictionary file was generated in this iteration
  if [ -e fort.7 ] ; then
        echo generator has run successfully
  else
        echo Mechanism generation failed
        exit 17
  fi

# create continuation files for next iteration

if [ ! -e existing.rxns ] ; then
  mv fort.17 existing.rxns
else
  cat existing.rxns fort.17 > tmp.txt ; mv tmp.txt existing.rxns
fi

cat fort.21 existing.rxns > ${mechname}.mech
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

  if [ -e warnings.out ] ; then
    if [ -e fort.50 ] ; then
      if [ -e $mechname.warnings ] ; then
        cat $mechname.warnings warnings.out fort.50 > tmp.out
        mv tmp.out $mechname.warnings
      else
        cat warnings.out fort.50 >  $mechname.warnings
      fi
    else
      if [ -e $mechname.warnings ] ; then
        cat $mechname.warnings warnings.out > tmp.out
        mv tmp.out $mechname.warnings
      else
        mv warnings.out $mechname.warnings
      fi
    fi
  fi

if [ -e existing.cheminput ] ; then
  echo '  --precursors from previously existing mechanism--' > tmp.txt
  cat scheme.log tmp.txt existing.cheminput > tmp.log
  mv tmp.log scheme.log
  rm tmp.txt
fi

# last time around: add 'END' to net reaction rate files
# and rename output files  ...
if [ $counter -eq $numpre ] ; then

mv difvol.dat  $mechname.difv
mv pnan.dat    $mechname.pnan
mv psim.dat    $mechname.psim
mv dHeatf.dat  $mechname.dHeatf
mv henry.dat   $mechname.Henry
mv scheme.log  $mechname.log
mv fort.7      $mechname.dict
mv lstprim.out      userdat.cheminput
mv userparams.input userdat.settings

# add existing precursors to END of logfile
#if [ -e existing.cheminput ] ; then
#  echo '  --precursors from previously existing mechanism--' > tmp.txt
#  cat scheme.log tmp.txt existing.cheminput > tmp.log
#  mv tmp.log scheme.log
#  rm tmp.txt
#fi

# add existing precursors to START of userdat.cheminput
# -> remove END statement line from existing.cheminput
sed -i '/^END/d' existing.cheminput
cat existing.cheminput userdat.cheminput > tmp.txt
mv tmp.txt userdat.cheminput
# remove leading blanks from cheminput file
# so it can be reused
sed 's/^ *//g' -i userdat.cheminput
sed 's/ $*//g' -i userdat.cheminput

echo 'END' >> $mechname.kOH
echo 'END' >> $mechname.kO3
echo 'END' >> $mechname.kNO3
echo 'END' >> userdat.cheminput

# ... and tidy up
rm cm
rm fort.*
rm sortlist.*
rm manage.input
rm existing.*
rm cheminput.dat
rm warnings.out
rm  ../DATA

fi

