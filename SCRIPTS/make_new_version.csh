#!/bin/csh
# =================================================== #
# PURPOSE: make new copy of GECKO-A for distribution
# =================================================== #

if ($#argv != 1) goto error1

set version=$1

set srcdir='/ur/julial/GECKO/GIT_COPY/GECKO-A'
set newdir='/ur/julial/GECKO/TUTORIAL/'$version
mkdir $newdir
cd $newdir 

# 1: code & output directories
foreach dir (BOXMOD BOXMOD_RUNS GECKO-A GECKO_SCRATCH POSTPROCESSING POST_SCRATCH)
  mkdir $dir 
end

# 2: copy-everything directories
foreach dir (SCRIPTS TUV_V4.2_GECKO VERSION_NOTES)
  cp -R $srcdir/$dir  $newdir
end

# 3: model section directory
# 4: code & model data directories
# 5: input and compilation directories
set dir=GECKO-A
  foreach subdir (DATA LIB RUN)
    cp -R $srcdir/$dir/$subdir  $newdir/$dir #/$subdir
  end
  foreach subdir (INPUTS OBJ)
    mkdir  $newdir/$dir/$subdir
  end

set dir=BOXMOD
  foreach subdir (INTERP LIB PROG)
    cp -R $srcdir/$dir/$subdir  $newdir/$dir #/$subdir
  end
  foreach subdir (INPUTS OBJ)
    mkdir  $newdir/$dir/$subdir
  end

set dir=POSTPROCESSING
  foreach subdir (LIB RUN)
    cp -R $srcdir/$dir/$subdir  $newdir/$dir #/$subdir
  end
  foreach subdir (INPUTS OBJ)
    mkdir  $newdir/$dir/$subdir
  end

# 6: reference files
foreach dir (BOXMOD GECKO-A POSTPROCESSING)
  cp -R $srcdir/$dir/INPUTS/EXAMPLES $newdir/$dir/INPUTS
  switch ($dir)
    case 'BOXMOD':
      foreach file (README_BOXMOD_INDAT)
        cp $srcdir/$dir/INPUTS/$file $newdir/$dir/INPUTS
      end
    breaksw
    case 'GECKO-A':
      foreach file (cheminput.master cheminput.test settings_default settings_test)
        cp $srcdir/$dir/INPUTS/$file $newdir/$dir/INPUTS
      end
    breaksw
    case 'POSTPROCESSING':
      foreach file (postproc_flags.input userinput.nml.example README_POSTPROCESSING)
        cp $srcdir/$dir/INPUTS/$file $newdir/$dir/INPUTS
      end
    breaksw
  endsw
end
exit

error1:
  echo "ERROR: you must supply a version name as an argument"
  exit 1
