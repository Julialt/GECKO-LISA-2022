#!/bin/csh
# PURPOSE: postprocess multiple library files on Nitrogen
#-----------------------------------------------
set scriptdir = /ur/julial/GECKO/GIT_COPY/GECKO-A/SCRIPTS
set opdir =  /net/modeling1/data14a/GECKO-A/WEB/2020_LIBRARY/BOXMOD_RUNS
set tmpdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS/LIBRARY_RUNS
#set opdir =  /ur/julial/GECKO/GIT_COPY/GECKO-A/BOXMOD_RUNS
set postdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS/POST_SCRATCH
set webdir =  /acomstaff/julial/GECKO_WEBSITE_OUTPUTS
set flagfile =  postproc_flags_library.input

  #cd $opdir
  #tar -xvf $mech.tar


#-----------------------------------------------
#-----------------------------------------------
#foreach combo ( "E11" "F11" "G11" "H11" "I11" "E20" "F20" "G20" "H20" "I20" )
#  set col = `echo $combo | awk '{print substr($0,1,1)}'`
#  set row = `echo $combo | awk '{print substr($0,2,2)}'`
#  if ( $col == "D" ) set scenario = "remote"
#  if ( $col == "F" ) set scenario = "continental"
#  if ( $col == "G" ) set scenario = "pollutedcontinental"
#  if ( $col == "H" ) set scenario = "urban"
#
#  if ( $row == "12" ) set mech = "10-11-dimethyleicosane"


foreach mech ( "10-11-dimethyleicosane" ) #"23456-pentamethylheptane" "2346-tetramethyloctane" "246-trimethylnonane" "24-dimethyldecane" "2-methylnonane" "2-methyluncosane" "2-methylundecane" "2-methyltridecane" "2-methylheptadecane" "2-methylbutadec-1-ene" "2-methyloct-1-ene" "2-methylundecane" "33-diethyldecane" "33-diethyloctadecane" "33-diethyltetradecane" "45-dimethyloctane" "345-trimethylheptane" "45-dimethyloctane" "567-trimethylundecane" "67-dimethyldodecane" "789-trimethylpentadecane" "89-dimethylhexadecane"  "9-10-11-trimethylnonadecane" "butane" "butadec-6-ene" "isobutane" "cyclodecane" "cyclohexane" "cyclopentadecane" "decane" "dodecane" "dodec-1-ene" "docosane" "heptane" "heptadecane" "heptadec-1-ene" "heptadec-8-ene" "hexane" "isobutane" "nonane" "octane" "oct-1-ene"  "octadecane" "pentane" "propane" "tetradecane"  "benzene" "toluene" "apinene" "bpinene" "bmyrcene" "limonene" "ocimene" "sabinene" )
  foreach scenario ( "remote" "remotecontinental" "continental" "pollutedcontinental" "urban" ) 

  #set flagfile = 'postproc_flags.input'
  #set flagfile = 'postproc_flags.input.test'
  #set flagfile = 'postproc_flags_library.'$mech
  #set flagfile = 'postproc_flags_library.massspec'
  #set flagfile = 'postproc_flags_library.pvap'
  #set flagfile = 'postproc_flags_library.bubble'
  #set flagfile = 'postproc_flags_library.ocimene'
  #set flagfile = 'postproc_flags_library.henry'

#-----------------------------------------------
# create local link to box model output
    mkdir $tmpdir'/'$mech'_'$scenario'_'$scenario
    cd $tmpdir'/'$mech'_'$scenario'_'$scenario
# switch to correct version directory
    if ( -d $opdir/$mech'_'$scenario'_'$scenario'_1_1' ) then
      echo $mech'/'$scenario'_1_1'
      ln -sf $opdir/$mech'_'$scenario'_'$scenario'_1_1'/outdat.nc outdat.nc
    else
      echo $mech'/'$scenario' simple'
      ln -sf $opdir/$mech'_'$scenario'_'$scenario/outdat.nc outdat.nc
    endif
    mkdir $postdir'/'$mech'/'$scenario'_'$scenario

# run postprocessor
    cd $scriptdir
      echo './run_postproc_library.bash -m '$mech' -k indat_'$scenario'_'$scenario'.key -f '${flagfile}
           ./run_postproc_library.bash -m $mech -k 'indat_'$scenario'_'$scenario'.key' -f ${flagfile}

#    end
    cd $postdir'/'$mech'/'$scenario'_'$scenario
#    cp *.csv $webdir'/'$mech'/'$scenario

    #cd $opdir
    #rm -R $mech"_"$scenario"_"$scenario

  end

end
exit
#-----------------------------------------------
"10-11-dimethyleicosane" 
"23456-pentamethylheptane" "2346-tetramethyloctane" 
"246-trimethylnonane" "24-dimethyldecane" 
"2-methylnonane" "2-methyluncosane" "2-methylundecane" "2-methyltridecane" "2-methylheptadecane" "2-methylbutadec-1-ene" "2-methyloct-1-ene" "2-methylundecane"
"33-diethyldecane" "33-diethyloctadecane" "33-diethyltetradecane" "45-dimethyloctane" 
"345-trimethylheptane" 
"45-dimethyloctane"
"567-trimethylundecane" 
"67-dimethyldodecane" 
"789-trimethylpentadecane" 
"89-dimethylhexadecane"  
"9-10-11-trimethylnonadecane" 
"butane" "butadec-6-ene" "isobutane" 
"cyclodecane" "cyclohexane" "cyclopentadecane" 
"decane" "dodecane" "dodec-1-ene" "docosane" 
"heptane" "heptadecane" "heptadec-1-ene" "heptadec-8-ene"
"hexane" 
"isobutane" "isoprene"
"nonane" 
"octane" "oct-1-ene"  "octadecane" 
"pentane" "propane" 
"tetradecane"  
"benzene" "toluene" 
"apinene" "bpinene" "bmyrcene" "limonene" "ocimene" "sabinene"
