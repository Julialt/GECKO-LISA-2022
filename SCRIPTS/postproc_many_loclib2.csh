#!/bin/csh
# PURPOSE: postprocess multiple library files on Nitrogen
set scriptdir = /ur/julial/GECKO/GIT_COPY/GECKO-A/SCRIPTS
set opdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS/LIBRARY_RUNS
set postdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS/POST_SCRATCH
set webdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS
set flagfile =  postproc_flags_library.input

# 10-11-dime,,,
foreach mech ( "propane" "butane" "isobutane" "pentane" "isoprene" "hexane" "cyclohexane" "benzene" "heptane" "toluene" "octane" "oct-1-ene" "nonane" "2-methyloct-1-ene" "decane" "2-methylnonane" "45-dimethyloctane" "345-trimethylheptane" "cyclodecane" "apinene" "bpinene" "bmyrcene" "limonene" "ocimene" g   "sabinene"      "dodecane"             "dodec-1-ene"          "2-methylundecane"     "24-dimethyldecane"    "246-trimethylnonane"  "2346-tetramethyloctane" "23456-pentamethylheptane" "tetradecane"          "2-methyltridecane"    "33-diethyldecane"     "67-dimethyldodecane"  "567-trimethylundecane" "butadec-6-ene"        "cyclopentadecane"     "2-methylbutadec-1-ene" "heptadecane"          "heptadec-1-ene"       "heptadec-8-ene"       "octadecane"           "2-methylheptadecane"  "89-dimethylhexadecane" "33-diethyltetradecane" "789-trimethylpentadecane" "docosane"             "2-methyluncosane"     "10-11-dimethyleicosane" "33-diethyloctadecane" "9-10-11-trimethylnonadecane" )

  cd $opdir
  tar -xvf $opdir/$mech.tar
  cd $scriptdir

  foreach scenario ("remote" "remotecontinental" "continental" "pollutedcontinental") 
    echo './run_postproc_library.bash -m '$mech' -k indat_'$scenario'_'$scenario'_1_1.key -f postproc_flags_library_extended.input'
    .run_postproc_library.bash -m $mech -k 'indat_'$scenario'_'$scenario'_1_1.key' -f postproc_flags_.input
    cp $postdir'/'$mech'/'$scenario'_'$scenario'_1_1/*.csv' $webdir'/'$mech'/'$scenario'/.'
  end

  rm -R $mech'*1_1'

end
