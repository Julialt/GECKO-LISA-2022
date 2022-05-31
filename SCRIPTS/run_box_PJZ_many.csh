#!/bin/csh

set mdate = '220406'
set msuff = 'gaw_DHF_'$mdate
set photfile = 'CU_UV_190215.phot'
set kdate = '220406'

foreach ng ( '1g' )

#decanols
set X = "O"
#set keyfile = 'indat_ch2_W_'$X'_'$kdate'.key'
set keyfile = 'indat_ch2_W_'$X'_'$kdate'.key'
foreach carbon ( '1' '2' '3' '4' '5' )
  set species = 'decan-'$carbon'-ol'

  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$ng'_'$msuff
  echo 'mech = '$mech

  echo './run_boxmod_local.bash -m '$mech' -k '$keyfile' -p '$photfile
  bash ./run_boxmod_local.bash -m $mech -k $keyfile -p $photfile
  echo " "

end

#nitrates
set X = "N"
set keyfile = 'indat_ch2_W_'$X'_'$kdate'.key'
foreach carbon ( '1' '2' '3' '4' '5' )
  set species = 'dec'$carbon'nitr'

  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$ng'_'$msuff
  echo 'mech = '$mech

  echo './run_boxmod_local.bash -m '$mech' -k '$keyfile' -p '$photfile
  bash ./run_boxmod_local.bash -m $mech -k $keyfile -p $photfile
  echo " "

end

# decanones
set X = "K"
set keyfile = 'indat_ch2_W_'$X'_'$kdate'.key'
foreach carbon ( '2' '3' '4' '5' '6')
 set species = 'dodecan-'$carbon'-one'
   
  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$ng'_'$msuff
  echo 'mech = '$mech

  echo './run_boxmod_local.bash -m '$mech' -k '$keyfile' -p '$photfile
  bash ./run_boxmod_local.bash -m $mech -k $keyfile -p $photfile
  echo " "
end

end

#--------------------------------------------------
exit
#--------------------------------------------------

