#!/bin/csh

#=====================================================================
# PURPOSE: link to EITHER binary or NetCDF versions of subroutines
#=====================================================================

# input argument = "bin" or "ncdf"
set mode = $argv[1]

foreach file (read_idgaw spakinit9 spreadro2sof2 readpvap_nan spreaddep3)
  # save a copy of previous version in case of accidents
  cp $file.f $file.f.save
  # assign filename for makefile to version with desired i/o format 
  cp $file'_'$mode.f $file.f
end
