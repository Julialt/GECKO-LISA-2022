#!/bin/csh

foreach file( 	   calc_Cstar298.f90 calc_psat.f90 defopvals_ncdf.f get_cheminp.f90 get_prev.f90 gitadd.csh manage_initial.f90 module_data_gecko_main.f90 mtdyn_rp.f mtrat.f open_op_files.f90 read_idgaw.f readpvap_nan_bin.f readpvap_nan_ncdf.f readpvap_sim.f soa_dyn.f90 sorting_module.f90 spakinit9_bin.f spakinit9_ncdf.f spprint_bin.f spreaddep3_ncdf.f write_op_hdrs.f90)

echo "git add "$file
git add $file

end

