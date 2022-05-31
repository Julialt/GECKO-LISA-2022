#!/bin/csh

set topdir = $argv[1]
set file = $argv[2]
set version1 = GECKO-A_200430
set version2 = GECKO-A_200506
set dir = BOXMOD/$topdir

#echo
#echo "dir = "$dir
echo '< '$version1
echo '> '$version2
#echo
diff -b $version1/$dir/$file $version2/$dir/$file
#echo

exit
# exit optional to short-circuit the copying section
echo 'copy '$file 
echo 'from '$version2'/'$dir' (>)'
echo '  to '$version1'/'$dir' (<)   ?' 
set yesno = $<

if ( $yesno == 'y' ) then
  echo 'copying '$file' to '$version1/$dir/$file
  cp $version2/$dir/$file $version1/$dir
else
  echo 'ok, I wont copy it'
endif

