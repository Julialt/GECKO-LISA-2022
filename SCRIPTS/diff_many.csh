#!/bin/csh

set topdir = $argv[1]

set version1 = GECKO-A_200430
set dir = BOXMOD/$topdir
set listfile = $topdir.list

ls -1 $version1/$dir > $listfile

set v = `cat $listfile`
@ i = 1
while ( $i <= $#v )
    echo " "
    echo "----------------------"
    echo $v[$i]
    ./diff.csh $topdir $v[$i]

    echo " "                      >> diff.out
    echo "----------------------" >> diff.out
    echo $v[$i]                   >> diff.out
    ./diff.csh $topdir $v[$i]     >> diff.out

    echo -n hit any key to continue
    set answer = $<

    @ i = $i + 1
end


