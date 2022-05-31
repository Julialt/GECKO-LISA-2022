#!/bin/bash

#PRO diagnostics
#==================================================================
#PURPOSE:
#Process a completed GECKO-A mechanism and return information on:
# - mechanism and dictionary size => for akparameter.h file setup
# - names of precursors => helps interpretation of boxmodel results
#INPUTS:
# - (where {file} is the name of the mechanism)
# - cheminput.dat
# - {file}.mech
# - {file}.dict
# - {file}.prec
# - {file}.pvap
# - {file}.Henry
# - X* files
# - (all files should be located in the same directory)
#OUTPUT:
# - file.info
#RUN LOCATION: dir 'WORK' containing results dir {FILE}
#==================================================================

if [ $# -ne 1 ]; then
  echo "enter mechanism name"
  read mech
    if [ ${#mech} -eq 0 ] ; then echo "no value supplied" ; exit 2 ; fi
else
  mech=$1
fi

pathfile="setup.dat"
if [ ! -e $pathfile ] ; then
    echo path setup file could not be found ; exit 1 ; fi

home_dir=`       grep "home_dir"       ${pathfile} | awk '{print $3}' `
scratch_dir=`    grep "scratch_dir"    ${pathfile} | awk '{print $3}' `
gecko_run_dir=`  grep "gecko_run_dir"  ${pathfile} | awk '{print $3}' `
boxmod_run_dir=` grep "boxmod_run_dir" ${pathfile} | awk '{print $3}' `

#--construct paths and report to screen.
work_dir=$scratch_dir/$gecko_run_dir/$mech
echo 'work dir = ' $work_dir

#--------------------------------------------------
infile=${work_dir}'/cheminput.dat'
dictfile=${work_dir}'/'${mech}'.dict'
mechfile=${work_dir}'/'${mech}'.mech'
kOHfile=${work_dir}'/'${mech}'.kOH'
kO3file=${work_dir}'/'${mech}'.kO3'
kNO3file=${work_dir}'/'${mech}'.kNO3'
precfile=${work_dir}'/'${mech}'.prec'
pvapfile=${work_dir}'/'${mech}'.pnan'
Henryfile=${work_dir}'/'${mech}'.Henry'
generalfile=${work_dir}'/general.h'

echo linking RO2 files...
ln -sf ${work_dir}/XP1O2 ${work_dir}/indat1.ro2
ln -sf ${work_dir}/XP2O2 ${work_dir}/indat2.ro2
ln -sf ${work_dir}/XP3O2 ${work_dir}/indat3.ro2
ln -sf ${work_dir}/XS1O2 ${work_dir}/indat4.ro2
ln -sf ${work_dir}/XS2O2 ${work_dir}/indat5.ro2
ln -sf ${work_dir}/XS3O2 ${work_dir}/indat6.ro2
ln -sf ${work_dir}/XT1O2 ${work_dir}/indat7.ro2
ln -sf ${work_dir}/XT2O2 ${work_dir}/indat8.ro2
ln -sf ${work_dir}/XACO3 ${work_dir}/indat9.ro2

ln -sf  ${home_dir}/GECKO-A/LIB/general.h ${generalfile}

ifill='INTEGER,PARAMETER ::'
lfill='LOGICAL,PARAMETER ::'

echo 'creating akparameter file for mechanism '${mech}
echo '!akparameter.h file for mechanism '${mech} > $work_dir'/akparameter.h'

# create a module file in parallel
echo '!akparameter_module file for mechanism '${mech} > $work_dir'/akparameter_module.f90'
echo '      MODULE akparameter_module' >> $work_dir'/akparameter_module.f90'
echo '                               ' >> $work_dir'/akparameter_module.f90'
echo '      IMPLICIT NONE            ' >> $work_dir'/akparameter_module.f90'
echo '                               ' >> $work_dir'/akparameter_module.f90'
#--------

#maxlsp=8
#maxlsp=10
# read line with lco from general.h
# EXCLUDE lines deactivated with !PARAMETER
maxlsp=`grep -v "\!P" ${generalfile} | grep lco= |cut -d= -f2 | cut -d\) -f1`
let "maxlsp = maxlsp +2"
echo '! max length of species names' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxlsp = '${maxlsp} >> $work_dir'/akparameter.h'

maxreac_char=90
echo '! max length of a printed reaction' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxreac_char = ' ${maxreac_char} >> $work_dir'/akparameter.h'

#maxsp=`wc -l ${dictfile} | awk '{print $1}'`
#get the first occurence of END, this corresponds to the end of the species list
maxsp=`grep -n END ${mechfile} -m1 | cut -d: -f1`
# remove 2 to account for SPECIES and END
let "maxsp = maxsp -2"
echo '! max number of species' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxsp = '${maxsp} >> $work_dir'/akparameter.h'

mxleft=2
echo '! max number of reactants in a reaction' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxleft = '${mxleft} >> $work_dir'/akparameter.h'

mxright=4
echo '! max number of products in a reaction' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxright = '${mxright} >> $work_dir'/akparameter.h'

maxre=`grep -c '=>' ${mechfile} | awk '{print $1}'`
echo '! max number of reactions' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxre = '${maxre} >> $work_dir'/akparameter.h'

max_m=`grep -c '+ M' ${mechfile} | awk '{print $1}'`
echo '! max number of reactions with "M"' >> $work_dir'/akparameter.h'
echo '      '$ifill' max_m = '${max_m} >> $work_dir'/akparameter.h'

maxfo=`grep -c 'TROE' ${mechfile} | awk '{print $1}'`
echo '! max number of fall-off reactions ' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxfo = '${maxfo} >> $work_dir'/akparameter.h'

nrhvorg=$(grep -c '+HV' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
nrhvinorg=$(grep -c '+ HV' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
let maxhv=$((nrhvorg+nrhvinorg))
echo '! max number of reactions with "HV"' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxhv = '${maxhv} >> $work_dir'/akparameter.h'

nrcvar=$(grep -c 'CVAR' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
let maxcvar=$nrcvar+1
echo '! max number of reactions with "CVAR"' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxcvar = '${maxcvar} >> $work_dir'/akparameter.h'

let maxextra=$(grep -c '+ EXTRA' ${mechfile} | awk '{print $1}')+$(grep -c '+EXTRA' ${mechfile} | awk '{print $1}')
echo '! max number of reactions with "EXTRA"' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxextra = '${maxextra} >> $work_dir'/akparameter.h'

maxo2=`grep -c 'OXYGEN' ${mechfile} | awk '{print $1}'`
echo '! max number of reactions with "OXYGEN"' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxo2 = '${maxo2} >> $work_dir'/akparameter.h'

nriso=$(grep -c 'ISOM' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
let maxiso=$nriso+1
echo '! max number of ISOMERIZATIONS' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxiso = '${maxiso} >> $work_dir'/akparameter.h'

# search on 'A******' name in list of species at start of mech
# note: requires FIXED FORMAT species list
#maxt=`grep -c '^A...............\/' ${mechfile} | awk '{print $1}'`
#on the species list (maxsp first lines of indat.mech), count species not starting with G
# grep -v option inverts the selection
#maxt=`head -${maxsp} ${mechfile} | grep -vc ^G`

echo '! max number of species undergoing phase equilibrium' >> $work_dir'/akparameter.h'
maxt=1
nt=$(grep -c '+AIN' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; maxt=$(( nt > maxt ? nt : maxt ))
nt=$(grep -c '+AOU' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; maxt=$(( nt > maxt ? nt : maxt ))
nt=$(grep -c '+WIN' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; maxt=$(( nt > maxt ? nt : maxt ))
nt=$(grep -c '+WOU' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; maxt=$(( nt > maxt ? nt : maxt ))
echo '      '$ifill' maxt = '${maxt} >> $work_dir'/akparameter.h'

#maxaux=12
maxaux=11
echo '! max number of different types of auxiliary information' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxaux = '${maxaux} >> $work_dir'/akparameter.h'

maxro2=9
echo '! max number of different classes of RO2' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxro2 = '${maxro2} >> $work_dir'/akparameter.h'

# defining mxro2cl as maximum n(PEROx) + 1 in 'X' files
echo '! max number of PEROx in a class' >> $work_dir'/akparameter.h'
mxro2cl=1
np=$(wc -l ${work_dir}'/indat1.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat2.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat3.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat4.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat5.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat6.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat7.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat8.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
np=$(wc -l ${work_dir}'/indat9.ro2' | awk '{VAR =+ $1} {print VAR}') ; mxro2cl=$(( np > mxro2cl ? np : mxro2cl ))
echo '      '$ifill' mxro2cl = '${mxro2cl} >> $work_dir'/akparameter.h'

# defining mxrpero as maximum n(PEROx) in mech file +1 (end-of-file marker)
echo '! max number of reactions with PEROx' >> $work_dir'/akparameter.h'
mxrpero=1
np=$(grep -c 'MEPERO' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO1' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO2' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO3' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO4' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO5' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO6' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO7' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO8' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))
np=$(grep -c 'PERO9' ${mechfile} | awk '{VAR =+ $1} {print VAR}') ; mxrpero=$(( $np > $mxrpero ? $np : $mxrpero ))

# if there are more PERO in a class than available reactions, substitute mxro2cl into mxrpero
mxrpero=$(( $mxro2cl > $mxrpero ? $mxro2cl : $mxrpero ))

let mxrpero=$mxrpero   
echo '      '$ifill' mxrpero = '${mxrpero} >> $work_dir'/akparameter.h'

maxdimer=4
echo '! max number of different classes of dimer' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxdimer = '${maxdimer} >> $work_dir'/akparameter.h'

ndim0=1 # dummy value to avoid defining a zero-dimension array
ndim1=$(grep -c 'DIM_1' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
ndim2=$(grep -c 'DIM_2' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
ndim3=$(grep -c 'DIM_3' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
ndim4=$(grep -c 'DIM_4' ${mechfile} | awk '{VAR =+ $1} {print VAR}')
ardim=($ndim0 $ndim1 $ndim2 $ndim3 $ndim4)
IFS=$'\n'
mxrdimer=$(echo "${ardim[*]}" | sort -nr | head -n1)
echo '! max number of dimers in a given class' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxrdimer = '${mxrdimer} >> $work_dir'/akparameter.h'

maxcoe=30
echo '! max number of variable coefficients in CVAR type reaction' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxcoe = '${maxcoe} >> $work_dir'/akparameter.h'

nset=4
echo '! max number of data set in function of temp. in CVAR type reaction' >> $work_dir'/akparameter.h'
echo '      '$ifill' nset = '${nset} >> $work_dir'/akparameter.h'

maxang=12
echo '! max number of angles in "HV" function' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxang = '${maxang} >> $work_dir'/akparameter.h'

echo '! -------------------------- DATA FOR CHROMOPHORE ---------------' >> $work_dir'/akparameter.h'

grep "HV/" ${mechfile} > grep.txt
grep "HV /" ${mechfile} >> grep.txt
cat grep.txt | cut -d'/' -f 2 | sort -n | uniq > uniq.txt
mchromo=`wc -l uniq.txt | awk '{print $1}'`
rm grep.txt
rm uniq.txt
echo '! max # of different types of chromophore' >> $work_dir'/akparameter.h'
echo '      '$ifill' mchromo = '${mchromo} >> $work_dir'/akparameter.h'

mtopchromo=10 #formerly 7
echo '! number of chromophores  to be stored in "most used"  chromophore' >> $work_dir'/akparameter.h'
echo '!     (the "top" tables)' >> $work_dir'/akparameter.h'
echo '      '$ifill' mtopchromo = '${mtopchromo} >> $work_dir'/akparameter.h'

let msptopchromo=$maxhv
# note: this sets an absolute maximum. We really want a max PER CHROMOPHORE
echo '! max # of species that can be stored in "most used"  chromophore' >> $work_dir'/akparameter.h'
echo '!     (the "top" tables)' >> $work_dir'/akparameter.h'
echo '      '$ifill' msptopchromo = '${msptopchromo} >> $work_dir'/akparameter.h'

mmedchromo=500
echo '! number of chromophores  to be stored in "regularly used"  chromophore' >> $work_dir'/akparameter.h'
echo '!     (the "med" tables)' >> $work_dir'/akparameter.h'
echo '      '$ifill' mmedchromo = '${mmedchromo} >> $work_dir'/akparameter.h'

let mspmedchromo=$maxhv
# note: this sets an absolute maximum. We really want a max PER CHROMOPHORE
echo '! max # of species that can be stored in "regularly used"  chromophore' >> $work_dir'/akparameter.h'
echo '!     (the "med" tables)' >> $work_dir'/akparameter.h'
echo '      '$ifill' mspmedchromo = '${mspmedchromo} >> $work_dir'/akparameter.h'

echo '! -------------------------- ' >> $work_dir'/akparameter.h'

nlo='maxang*3'
echo '! max coefficient for interpolation of the photolytic frequencies' >> $work_dir'/akparameter.h'
echo '      '$ifill' nlo = '${nlo} >> $work_dir'/akparameter.h'

mbox=2
echo '! max number of boxes in the model' >> $work_dir'/akparameter.h'
echo '      '$ifill' mbox = '${mbox} >> $work_dir'/akparameter.h'

mhd=60
echo '! max number of data to compute mixing height' >> $work_dir'/akparameter.h'
echo '      '$ifill' mhd = '${mhd} >> $work_dir'/akparameter.h'

msur=4
echo '! max number of isurface types' >> $work_dir'/akparameter.h'
echo '      '$ifill' msur = '${msur} >> $work_dir'/akparameter.h'

mopc=4
echo '! max number of counting species for which stoe. coff. need to be' >> $work_dir'/akparameter.h'
echo '!     evaluated in CVAR application' >> $work_dir'/akparameter.h'
echo '      '$ifill' mopc = '${mopc} >> $work_dir'/akparameter.h'

mpos=6
echo '! max number of positions used to evaluate stoe. coff. for counting' >> $work_dir'/akparameter.h'
echo '!     species from the operator species' >> $work_dir'/akparameter.h'
echo '      '$ifill' mpos = '${mpos} >> $work_dir'/akparameter.h'


# need to add the number of criegee intermediates to mself
# because each of them potentially has a self reaction
ncriegee=`grep -c '^4' ${dictfile}`
let mself=$ncriegee*2+20
echo '! max number of self reactions' >> $work_dir'/akparameter.h'
echo '      '$ifill' mself = '${mself} >> $work_dir'/akparameter.h'

echo '! -------------------------- ' >> $work_dir'/akparameter.h'

if [ -e ${kOHfile} ]; then
  mxkOH=`wc -l ${kOHfile} | awk '{print $1}'`
else
  mxkOH=0
fi
echo '! max number of species for which kOH is evaluated' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxkOH = '${mxkOH} >> $work_dir'/akparameter.h'

if [ -e ${kO3file} ]; then
  mxkO3=`wc -l ${kO3file} | awk '{print $1}'`
else
  mxkO3=0
fi
echo '! max number of species for which kO3 is evaluated' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxkO3 = '${mxkO3} >> $work_dir'/akparameter.h'

if [ -e ${kNO3file} ]; then
  mxkNO3=`wc -l ${kNO3file} | awk '{print $1}'`
else
  mxkNO3=0
fi
echo '! max number of species for which kNO3 is evaluated' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxkNO3 = '${mxkNO3} >> $work_dir'/akparameter.h'

mxsat=`wc -l ${pvapfile} | awk '{print $1}'`
echo '! max number of species for which Psat is evaluated' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxsat = '${mxsat} >> $work_dir'/akparameter.h'

mxdep=`wc -l ${Henryfile} | awk '{print $1}'`
let mxdep=$mxdep-`grep -c '*' ${Henryfile}`
echo '! max number of species for which Vdep is evaluated' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxdep = '${mxdep} >> $work_dir'/akparameter.h'

echo '! -------------------------- ' >> $work_dir'/akparameter.h'

mxprec=`wc -l ${precfile} | awk '{print $1}'`
echo '! max number of precursors' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxprec = '${mxprec} >> $work_dir'/akparameter.h'

mes=60
echo '! max number of emitted species' >> $work_dir'/akparameter.h'
echo '      '$ifill' mes = '${mes} >> $work_dir'/akparameter.h'
echo '      '$ifill' maxem = '${mes} >> $work_dir'/akparameter.h'

mtim=50
echo '! max number of emitted times' >> $work_dir'/akparameter.h'
echo '      '$ifill' mtim = '${mtim} >> $work_dir'/akparameter.h'

mtr=100
echo '! max number of spp assessing prod & loss rates' >> $work_dir'/akparameter.h'
echo '      '$ifill' mtr = '${mtr} >> $work_dir'/akparameter.h'

maxconst=20
echo '! max number of species that can be constrained' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxconst = '${maxconst} >> $work_dir'/akparameter.h'

maxinput=1500
echo '! max number of lines for datapoints for constrained species' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxinput = '${maxinput} >> $work_dir'/akparameter.h'


nhyd=`grep -c 'HYD' ${mechfile} | awk '{print $1}'`
let maxhyd=$nhyd/2
echo '! max number of hydration reactions ' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxhyd = '${maxhyd} >> $work_dir'/akparameter.h'

nacid=`grep -c 'ACID' ${mechfile} | awk '{print $1}'`
let maxacid=$nacid/2
echo '! max number of acid-base reactions ' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxacid = '${maxacid} >> $work_dir'/akparameter.h'

nohaq=`grep -c 'Aq' ${mechfile} | awk '{print $1}'`
let mxohaq=$nohaq/2
echo '! max number of aqueous phase reactions ' >> $work_dir'/akparameter.h'
echo '      '$ifill' mxohaq = '${mxohaq} >> $work_dir'/akparameter.h'

ntr=`grep -c 'HENRY' ${mechfile} | awk '{print $1}'`
let maxtr=$ntr/2
echo '! max number of aqueous mass transfer reactions ' >> $work_dir'/akparameter.h'
echo '      '$ifill' maxtr = '${maxtr} >> $work_dir'/akparameter.h'

echo '! -------------------------- ' >> $work_dir'/akparameter.h'
#wall_fg='.FALSE.'
#echo '! flag to allow wall partitioning ' >> $work_dir'/akparameter.h'
#echo '      '$lfill' wall_fg = '${wall_fg} >> $work_dir'/akparameter.h'

echo '! structure for external data used to constrain
! species emissions, concentrations, deposition...
! with possibility for using an input file' >> $work_dir'/akparameter.h'
echo '      TYPE species_data ' >> $work_dir'/akparameter.h'
echo '        LOGICAL              :: activefg ' >> $work_dir'/akparameter.h'
echo '        CHARACTER(LEN=maxlsp):: name ' >> $work_dir'/akparameter.h'
echo '        CHARACTER(LEN=10)    :: unit ' >> $work_dir'/akparameter.h'
echo '        INTEGER              :: index,npoints ' >> $work_dir'/akparameter.h'
echo '        REAL                 :: table(maxinput,2) ' >> $work_dir'/akparameter.h'
echo '      END TYPE ' >> $work_dir'/akparameter.h'

echo '!structure for storing surface data (emissions for now
! TODO later: add deposition data...' >> $work_dir'/akparameter.h'
echo '      TYPE surface_data ' >> $work_dir'/akparameter.h'
echo '        INTEGER              :: nemis ' >> $work_dir'/akparameter.h'
echo '        TYPE(species_data)   :: emission(maxem) ' >> $work_dir'/akparameter.h'
echo '      END TYPE ' >> $work_dir'/akparameter.h'

# -----
cat $work_dir'/akparameter_module.f90' $work_dir'/akparameter.h' > $work_dir'/akparameter.tmp'
mv $work_dir'/akparameter.tmp' $work_dir'/akparameter_module.f90'
echo '                               ' >> $work_dir'/akparameter_module.f90'
echo '      END MODULE akparameter_module' >> $work_dir'/akparameter_module.f90'
# -----

# tidy up
rm ${work_dir}/*.ro2

exit
