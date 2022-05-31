#!/bin/bash

#PRO diagnostics
#==================================================================
#PURPOSE:
#Process a completed GECKO-A mechanism NETCDF FILE
# and reconstruct the relevant akparameter_module.f90 file
#INPUTS:
# - outdat.nc
#OUTPUT:
# - akparameter_module.f90
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

#--construct paths and report to screen.
work_dir=$scratch_dir/$gecko_run_dir/$mech
echo 'work dir = ' $work_dir

infile=${work_dir}'/outdat.nc'
echo 'infile = ' $infile
#--------------------------------------------------

ifill='INTEGER,PARAMETER ::'
lfill='LOGICAL,PARAMETER ::'

echo 'creating akparameter file for mechanism '${mech}
echo '!akparameter.h file for mechanism '${mech} > $work_dir'/akparameter.h'

# create a module file in parallel
# just write header now, 
# then add the rest at the end, using 'cat'
echo '!akparameter_module file for mechanism '${mech} > $work_dir'/akparameter_module.f90'
echo '!(reconstructed) '${mech} > $work_dir'/akparameter_module.f90'
echo '      MODULE akparameter_module' >> $work_dir'/akparameter_module.f90'
echo '                               ' >> $work_dir'/akparameter_module.f90'
echo '      IMPLICIT NONE            ' >> $work_dir'/akparameter_module.f90'
echo '                               ' >> $work_dir'/akparameter_module.f90'
#--------
#output function
write_output () {
  dimval=0
  tmpval=$(ncdump -h ${infile} | grep -m 1 ${dimnam} | awk '{print $3}') ; dimval=$(( tmpval > dimval ? tmpval : dimval ))
  echo $comment                               >> $work_dir'/akparameter.h'
  echo '      '$ifill' '$dimnam' = '${dimval} 
  echo '      '$ifill' '$dimnam' = '${dimval} >> $work_dir'/akparameter.h'
}
#--------

dimnam='maxlsp'
comment='! max length of species names' 
write_output

dimnam='maxreac_char'
comment='! max length of a printed reaction'
write_output

dimnam='maxsp'
comment='! max number of species'
write_output

dimnam='mxleft'
comment='! max number of reactants in a reaction'
write_output

dimnam='mxright'
comment='! max number of products in a reaction'
write_output

dimnam='maxre'
comment='! max number of reactions'
write_output

dimnam='max_m'
comment='! max number of reactions with "M"'
write_output

dimnam='maxfo'
comment='! max number of fall-off reactions'
write_output

dimnam='maxhv'
comment='! max number of reactions with "HV"'
write_output

dimnam='maxcvar'
comment='! max number of reactions with "CVAR"'
write_output

dimnam='maxextra'
comment='! max number of reactions with "EXTRA"'
write_output

dimnam='maxo2'
comment='! max number of reactions with "OXYGEN"'
write_output

dimnam='maxiso'
comment='! max number of ISOMERIZATIONS'
write_output

dimnam='maxt'
comment='! max number of species undergoing phase equilibrium'
write_output

dimnam='maxaux'
comment='! max number of different types of auxiliary information'
write_output

dimnam='maxro2'
comment='! max number of different classes of RO2'
write_output

dimnam='mxro2cl'
comment='! max number of PEROx in a class'
write_output

dimnam='mxrpero'
comment='! max number of reactions with PEROx'
write_output

dimnam='maxdimer'
comment='! max number of different classes of dimer'
write_output

dimnam='mxrdimer'
comment='! max number of dimers in a given class'
write_output

dimnam='maxcoe'
comment='! max number of variable coefficients in CVAR type reaction'
write_output

dimnam='nset'
comment='! max number of data sets in function of temp. in CVAR type reaction'
write_output

dimnam='maxang'
comment='! max number of angles in "HV" function'
write_output
echo '! -------------------------- DATA FOR CHROMOPHORE ---------------' >> $work_dir'/akparameter.h'

dimnam='mchromo'
comment='! max number of different types of chromophore'
write_output

dimnam='mtopchromo'
comment='! number of chromophores to be stored in "most used" chromophore'
write_output

dimnam='msptopchromo'
comment='! max # of species that can be stored in "most used" chromophore'
write_output

dimnam='mmedchromo'
comment='! number of chromophores to be stored in "regularly used" chromophore'
write_output

dimnam='mspmedchromo'
comment='! max # of species that can be stored in "regularly used" chromophore'
write_output
echo '! -------------------------- ' >> $work_dir'/akparameter.h'

dimnam='mspmedchromo'
comment='! max # of species that can be stored in "regularly used" chromophore'
write_output

dimnam='nlo'
comment='! max coefficient for interpolation of the photolytic frequencies'
write_output

dimnam='mbox'
comment='! max number of boxes in the model'
write_output

dimnam='mhd'
comment='! max number of data to compute mixing height'
write_output

dimnam='msur'
comment='! max number of isurface types'
write_output

dimnam='mopc'
comment='! max number of counting species for which stoe. coff. need to be'
echo '!     evaluated in CVAR application' >> $work_dir'/akparameter.h'
write_output

dimnam='mpos'
comment='! max number of positions used to evaluate stoe. coff. for counting'
echo '!     species from the operator species' >> $work_dir'/akparameter.h'
write_output

dimnam='mself'
comment='! max number of self reactions'
write_output
echo '! -------------------------- ' >> $work_dir'/akparameter.h'

dimnam='mxkOH'
comment='! max number of species for which kOH is evaluated'
write_output

dimnam='mxkO3'
comment='! max number of species for which kO3 is evaluated'
write_output

dimnam='mxkNO3'
comment='! max number of species for which kNO3 is evaluated'
write_output

dimnam='mxsat'
comment='! max number of species for which Psat is evaluated'
write_output

dimnam='mxdep'
comment='! max number of species for which Vdep is evaluated'
write_output
echo '! -------------------------- ' >> $work_dir'/akparameter.h'

dimnam='mes'
comment='! max number of emitted species'
write_output

dimnam='maxem'
comment='! ALSO max number of emitted species'
write_output

dimnam='mtim'
comment='! max number of emitted times'
write_output

dimnam='mtr'
comment='! max number of spp assessing prod & loss rates'
write_output

dimnam='maxconst'
comment='! max number of species that can be constrained'
write_output

dimnam='maxinput'
comment='! max number of lines for datapoints for constrained species'
write_output

dimnam='maxhyd'
comment='! max number of hydration reactions'
write_output

dimnam='maxacid'
comment='! max number of acid-base reactions'
write_output

dimnam='mxohaq'
comment='! max number of aqueous phase reactions'
write_output

dimnam='maxtr'
comment='! max number of aqueous mass transfer reactions'
write_output

echo '! -------------------------- ' >> $work_dir'/akparameter.h'

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

exit
