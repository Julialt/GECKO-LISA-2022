#!/bin/csh
# PURPOSE: postprocess multiple library files on Nitrogen
#-----------------------------------------------
#set scriptdir = /ur/julial/GECKO/GIT_COPY/GECKO-A/SCRIPTS
#set opdir =  /ur/julial/GECKO/GIT_COPY/GECKO-A/BOXMOD_RUNS
set keyfile =  indat_ch2_W_N_211221.key
set photfile = CU_UV_190215.phot
#set flagfile = postproc_flags_library.input

#foreach nc ( 1 2 3 4 5 ) 
#  set mech = 'dodecan'$nc'one_3g_gaw_210521'
foreach nc ( 1 2 3 4 5 ) 
  set mech = 'dec'$nc'nitr_2g_gaw_DHF_211217'

# run box model
  echo './run_boxmod_local.bash -m '${mech}' -k '${keyfile}' -p '${photfile}
  ./run_boxmod_local.bash -m ${mech} -k ${keyfile} -p ${photfile}

end
exit
