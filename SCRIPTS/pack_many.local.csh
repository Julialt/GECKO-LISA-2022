#!/bin/csh

set date = '220406'
set suff = '1g_gaw_DHF_'$date
set settingsfile = 'settings_PJZcham_1g_DHFon'


foreach carbon ( '1' '2' '3' '4' '5' )
 foreach species ( 'decan-'$carbon'-ol' 'dec'$carbon'nitr')

  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$suff
  echo 'mech = '$mech

  echo './gen_package_local.bash setup.dat '$mech
  echo " "
  bash ./gen_package_local.bash setup.dat ${mech} 

 end
end

foreach carbon ( '2' '3' '4' '5'  '6' )
  set species = 'dodecan-'$carbon'-one'

  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$suff
  echo 'mech = '$mech

  echo './gen_package_local.bash setup.dat '$mech
  echo " "
  bash ./gen_package_local.bash setup.dat ${mech} 

end
exit

foreach carbon ( '1' '2' '3' '4' '5' )
  set species = 'decan-'$carbon'-ol'
   
  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$suff
  echo 'mech = '$mech

  echo './gen_package_local.bash setup.dat '$mech
  echo " "
  bash ./gen_package_local.bash setup.dat ${mech} 

end
exit

foreach carbon ( '1' '2' '3' )#'4' '5' )
  set species = 'dec'$carbon'nitr'

  set input = 'cheminput.'$species
  echo "inpt = "$input

  set mech = $species'_'$suff
  echo 'mech = '$mech

  echo './gen_package_local.bash setup.dat '$mech
  echo " "
  bash ./gen_package_local.bash setup.dat ${mech} 

end
exit
   

