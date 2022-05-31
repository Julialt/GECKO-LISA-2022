#!/bin/csh
# PURPOSE: postprocess multiple library files on Nitrogen
#-----------------------------------------------
set scriptdir = /ur/julial/GECKO/GIT_COPY/GECKO-A/SCRIPTS
set opdir     = /ur/julial/GECKO/GIT_COPY/GECKO-A/BOXMOD_RUNS
set postdir   = /ur/julial/GECKO/GIT_COPY/GECKO-A/POST_SCRATCH
set flagfile =  postproc_flags.input.MG

#foreach species ( "propane" "butane" "2Mfuran" "isoprene" "chexane" "octane" "apinene")
#foreach species ( "isoprene" )
#  set mech = 'GEK-'$species'_fixi'

#foreach species ( "PROPANE" )
foreach species ( "PROPANE" "N-C4" "2MFURAN" "ISOPRENE" "CYCC6" "N-C8" "APINENE")
  set mech = 'FM-'$species'_fixi'

  if ( $species == "propane"  ) set date = "211105"
  if ( $species == "butane"   ) set date = "211117"
  if ( $species == "2Mfuran"  ) set date = "211115"
  if ( $species == "isoprene" ) set date = "211203"
  if ( $species == "chexane"  ) set date = "211117"
  if ( $species == "octane"   ) set date = "211117"
  if ( $species == "apinene"  ) set date = "211203"

  if ( $species == "PROPANE"  ) set date = "211104"
  if ( $species == "N-C4"     ) set date = "211118"
  if ( $species == "2MFURAN"  ) set date = "211202"
  if ( $species == "ISOPRENE" ) set date = "211207"
  if ( $species == "CYCC6"    ) set date = "211202"
  if ( $species == "N-C8"     ) set date = "211202"
  if ( $species == "APINENE"  ) set date = "211207"

#-----------------------------------------------
  set scenario = "FIXI_"$date
#-----------------------------------------------
# create local link to box model output
    mkdir $postdir'/'$mech'/'$scenario
    cd $postdir'/'$mech'/'$scenario

# run postprocessor
    cd $scriptdir
      echo './run_postproc_local.bash -m '$mech' -k indat_'$scenario'.key -f '${flagfile}
           ./run_postproc_local.bash -m $mech -k 'indat_'$scenario'.key' -f ${flagfile}

    cd $postdir'/'$mech'/'$scenario

  end

end
exit
#-----------------------------------------------
