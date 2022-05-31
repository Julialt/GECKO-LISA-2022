#!/bin/csh

set file=$argv[1]

set dir1 = LIB
set dir2 = LIB_SY

diff -b $dir1/$file $dir2/$file
