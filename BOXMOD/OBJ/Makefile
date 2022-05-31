# Makefile created by mkmf 19.3.1


MK_TEMPLATE = template.mk
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: ../PROG/BOXMOD
akparameter_module.o: ../LIB/akparameter_module.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/akparameter_module.f90
apply_constraints.o: ../LIB/apply_constraints.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/apply_constraints.f
bilan_species.o: ../LIB/bilan_species.f ../LIB/general.h akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/bilan_species.f
boxmod_main.o: ../PROG/boxmod_main.f90 flags_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../PROG/boxmod_main.f90
calc_Cstar298.o: ../LIB/calc_Cstar298.f90 flags_module.o akparameter_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/calc_Cstar298.f90
calc_phot_OFR.o: ../LIB/calc_phot_OFR.f90 forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/calc_phot_OFR.f90
calc_psat.o: ../LIB/calc_psat.f90 flags_module.o akparameter_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/calc_psat.f90
cdcase.o: ../LIB/cdcase.f ../LIB/general.h ../LIB/organic.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/cdcase.f
chemmap.o: ../LIB/chemmap.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/chemmap.f
cleanup.o: ../LIB/cleanup.f90 flags_module.o module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/cleanup.f90
constraints_module.o: ../LIB/constraints_module.f90 akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/constraints_module.f90
defenvinp_ncdf.o: ../LIB/defenvinp_ncdf.f90 akparameter_module.o flags_module.o module_data_gecko_main.o constraints_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/defenvinp_ncdf.f90
defflags_ncdf.o: ../LIB/defflags_ncdf.f akparameter_module.o flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/defflags_ncdf.f
defopvals_ncdf.o: ../LIB/defopvals_ncdf.f akparameter_module.o forcing_params_module.o module_data_gecko_main.o flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/defopvals_ncdf.f
extract_reacrates.o: ../LIB/extract_reacrates.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/extract_reacrates.f
find_concs.o: ../LIB/find_concs.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/find_concs.f90
findring.o: ../LIB/findring.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/findring.f
fit.o: ../LIB/fit.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/fit.f
flags_module.o: ../LIB/flags_module.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/flags_module.f90
forcing_params_module.o: ../LIB/forcing_params_module.f90 constraints_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/forcing_params_module.f90
general_module.o: ../LIB/general_module.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/general_module.f90
get_cheminp.o: ../LIB/get_cheminp.f90 sorting_module.o flags_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/get_cheminp.f90
get_envinp.o: ../LIB/get_envinp.f90 flags_module.o module_data_gecko_main.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/get_envinp.f90
get_prev.o: ../LIB/get_prev.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/get_prev.f90
getpvap_nan.o: ../LIB/getpvap_nan.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/getpvap_nan.f
gettrack.o: ../LIB/gettrack.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/gettrack.f
init_defaults.o: ../LIB/init_defaults.f90 module_data_gecko_main.o constraints_module.o forcing_params_module.o flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/init_defaults.f90
inorganic_aer_module.o: ../LIB/inorganic_aer_module.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/inorganic_aer_module.f90
iter4.o: ../LIB/iter4.f prodloss_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/iter4.f
manage_initial.o: ../LIB/manage_initial.f90 ../LIB/flags.h module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/manage_initial.f90
manage_loop.o: ../LIB/manage_loop.f90 module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/manage_loop.f90
map_indices.o: ../LIB/map_indices.f90 module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/map_indices.f90
map_sp_ids.o: ../LIB/map_sp_ids.f90 flags_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/map_sp_ids.f90
module_chamber_tools.o: ../LIB/module_chamber_tools.f90 akparameter_module.o module_data_gecko_main.o flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/module_chamber_tools.f90
module_data_gecko_main.o: ../LIB/module_data_gecko_main.f90 akparameter_module.o inorganic_aer_module.o prodloss_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/module_data_gecko_main.f90
mtdyn_rp.o: ../LIB/mtdyn_rp.f flags_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/mtdyn_rp.f
mtrat.o: ../LIB/mtrat.f flags_module.o akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/mtrat.f
mtrat_ba.o: ../LIB/mtrat_ba.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/mtrat_ba.f
ncutil.o: ../LIB/ncutil.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/ncutil.f
newdt.o: ../LIB/newdt.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/newdt.f
nodmap.o: ../LIB/nodmap.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/nodmap.f
open_io_files.o: ../LIB/open_io_files.f90 flags_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/open_io_files.f90
open_op_files.o: ../LIB/open_op_files.f90 flags_module.o module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/open_op_files.f90
output_vals.o: ../LIB/output_vals.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/output_vals.f90
pbl_height_forcage.o: ../LIB/pbl_height_forcage.f module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/pbl_height_forcage.f
phys_rates.o: ../LIB/phys_rates.f90 flags_module.o akparameter_module.o forcing_params_module.o module_data_gecko_main.o module_chamber_tools.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/phys_rates.f90
printreaction.o: ../LIB/printreaction.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/printreaction.f
printsteadystate.o: ../LIB/printsteadystate.f90 akparameter_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/printsteadystate.f90
prodloss_module.o: ../LIB/prodloss_module.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/prodloss_module.f90
read_dvsp_bin.o: ../LIB/read_dvsp_bin.f90 flags_module.o module_data_gecko_main.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/read_dvsp_bin.f90
read_dvsp_ncdf.o: ../LIB/read_dvsp_ncdf.f90 akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/read_dvsp_ncdf.f90
read_idgaw.o: ../LIB/read_idgaw.f flags_module.o akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/read_idgaw.f
readkeyflags.o: ../LIB/readkeyflags.f flags_module.o module_data_gecko_main.o forcing_params_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readkeyflags.f
readppf.o: ../LIB/readppf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readppf.f
readprev.o: ../LIB/readprev.f flags_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readprev.f
readprev_ncdf.o: ../LIB/readprev_ncdf.f ../LIB/flags.h akparameter_module.o module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/readprev_ncdf.f
readpvap.o: ../LIB/readpvap.f akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readpvap.f
readpvap_nan_bin.o: ../LIB/readpvap_nan_bin.f ../LIB/flags.h akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/readpvap_nan_bin.f
readpvap_nan_ncdf.o: ../LIB/readpvap_nan_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readpvap_nan_ncdf.f
readpvap_sim_bin.o: ../LIB/readpvap_sim_bin.f ../LIB/flags.h akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/readpvap_sim_bin.f
readpvap_sim_ncdf.o: ../LIB/readpvap_sim_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/readpvap_sim_ncdf.f
rjgadd.o: ../LIB/rjgadd.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/rjgadd.f
rjgrm.o: ../LIB/rjgrm.f ../LIB/general.h
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I../LIB	../LIB/rjgrm.f
setup_GECKO.o: ../LIB/setup_GECKO.f90 flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/setup_GECKO.f90
setup_ncdf_op.o: ../LIB/setup_ncdf_op.f90 flags_module.o module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/setup_ncdf_op.f90
setup_printphoto.o: ../LIB/setup_printphoto.f90 module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/setup_printphoto.f90
soa_dyn.o: ../LIB/soa_dyn.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/soa_dyn.f90
soa_dyn_update.o: ../LIB/soa_dyn_update.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/soa_dyn_update.f90
soa_equil.o: ../LIB/soa_equil.f90 module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/soa_equil.f90
soapartition.o: ../LIB/soapartition.f90 akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/soapartition.f90
solve_box.o: ../LIB/solve_box.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/solve_box.f90
sorting_module.o: ../LIB/sorting_module.f90 akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/sorting_module.f90
spakinit9_bin.o: ../LIB/spakinit9_bin.f akparameter_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakinit9_bin.f
spakinit9_ncdf.o: ../LIB/spakinit9_ncdf.f akparameter_module.o flags_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakinit9_ncdf.f
spakkextra4.o: ../LIB/spakkextra4.f flags_module.o module_data_gecko_main.o akparameter_module.o inorganic_aer_module.o vbs_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakkextra4.f
spakkfo.o: ../LIB/spakkfo.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakkfo.f
spakkrat6.o: ../LIB/spakkrat6.f flags_module.o akparameter_module.o inorganic_aer_module.o module_data_gecko_main.o vbs_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakkrat6.f
spakspnum.o: ../LIB/spakspnum.f sorting_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spakspnum.f
spchge_cvar4.o: ../LIB/spchge_cvar4.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spchge_cvar4.f
spdatsurf.o: ../LIB/spdatsurf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spdatsurf.f
spdeposition3.o: ../LIB/spdeposition3.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spdeposition3.f
spforcage6.o: ../LIB/spforcage6.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o constraints_module.o inorganic_aer_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spforcage6.f90
spgetemi4.o: ../LIB/spgetemi4.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spgetemi4.f
spinterp5.o: ../LIB/spinterp5.f akparameter_module.o flags_module.o module_data_gecko_main.o forcing_params_module.o module_chamber_tools.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spinterp5.f
spprint_bin.o: ../LIB/spprint_bin.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spprint_bin.f
spreadcoeff3.o: ../LIB/spreadcoeff3.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadcoeff3.f
spreaddep3.o: ../LIB/spreaddep3.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreaddep3.f
spreaddep3_ncdf.o: ../LIB/spreaddep3_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreaddep3_ncdf.f
spreaddimer.o: ../LIB/spreaddimer.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreaddimer.f
spreadj4.o: ../LIB/spreadj4.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadj4.f
spreadj4_ncdf.o: ../LIB/spreadj4_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadj4_ncdf.f
spreadkey6.o: ../LIB/spreadkey6.f90 flags_module.o akparameter_module.o module_data_gecko_main.o constraints_module.o forcing_params_module.o module_chamber_tools.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadkey6.f90
spreadro2jenk.o: ../LIB/spreadro2jenk.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadro2jenk.f
spreadro2sof2_bin.o: ../LIB/spreadro2sof2_bin.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadro2sof2_bin.f
spreadro2sof2_ncdf.o: ../LIB/spreadro2sof2_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spreadro2sof2_ncdf.f
spsolar.o: ../LIB/spsolar.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spsolar.f
spsplna.o: ../LIB/spsplna.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spsplna.f
spsplnb.o: ../LIB/spsplnb.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/spsplnb.f
sptimesurf3.o: ../LIB/sptimesurf3.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/sptimesurf3.f
sputil.o: ../LIB/sputil.f
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/sputil.f
surface_emissions.o: ../LIB/surface_emissions.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/surface_emissions.f
surface_rates.o: ../LIB/surface_rates.f90 flags_module.o module_data_gecko_main.o akparameter_module.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/surface_rates.f90
surfdep_wesely89.o: ../LIB/surfdep_wesely89.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/surfdep_wesely89.f
twostep5.o: ../LIB/twostep5.f prodloss_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/twostep5.f
vbs_module.o: ../LIB/vbs_module.f90 module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/vbs_module.f90
wall_Huang2018.o: ../LIB/wall_Huang2018.f90 module_chamber_tools.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wall_Huang2018.f90
wall_test.o: ../LIB/wall_test.f90 flags_module.o constraints_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wall_test.f90
wall_Krech2016_J.o: ../LIB/wall_Krech2016_J.f90 flags_module.o constraints_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wall_Krech2016_J.f90
wall_Krech2016_Z.o: ../LIB/wall_Krech2016_Z.f90 flags_module.o constraints_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wall_Krech2016_Z.f90
wall_MatZie2010.o: ../LIB/wall_MatZie2010.f90 forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wall_MatZie2010.f90
write_op_hdrs.o: ../LIB/write_op_hdrs.f90 flags_module.o module_data_gecko_main.o forcing_params_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/write_op_hdrs.f90
write_printphoto.o: ../LIB/write_printphoto.f90 module_data_gecko_main.o flags_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/write_printphoto.f90
wrtenvinp_ncdf.o: ../LIB/wrtenvinp_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wrtenvinp_ncdf.f
wrtflags_ncdf.o: ../LIB/wrtflags_ncdf.f flags_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wrtflags_ncdf.f
wrtjvinp_ncdf.o: ../LIB/wrtjvinp_ncdf.f akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wrtjvinp_ncdf.f
wrtopaer_ncdf.o: ../LIB/wrtopaer_ncdf.f flags_module.o akparameter_module.o forcing_params_module.o module_data_gecko_main.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wrtopaer_ncdf.f
wrtopgas_ncdf.o: ../LIB/wrtopgas_ncdf.f forcing_params_module.o akparameter_module.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	../LIB/wrtopgas_ncdf.f
./manage_loop.f90: ../LIB/manage_loop.f90
	cp ../LIB/manage_loop.f90 .
./setup_printphoto.f90: ../LIB/setup_printphoto.f90
	cp ../LIB/setup_printphoto.f90 .
./apply_constraints.f: ../LIB/apply_constraints.f
	cp ../LIB/apply_constraints.f .
./defenvinp_ncdf.f90: ../LIB/defenvinp_ncdf.f90
	cp ../LIB/defenvinp_ncdf.f90 .
./spreadro2sof2_ncdf.f: ../LIB/spreadro2sof2_ncdf.f
	cp ../LIB/spreadro2sof2_ncdf.f .
./spforcage6.f90: ../LIB/spforcage6.f90
	cp ../LIB/spforcage6.f90 .
./spinterp5.f: ../LIB/spinterp5.f
	cp ../LIB/spinterp5.f .
./spakkrat6.f: ../LIB/spakkrat6.f
	cp ../LIB/spakkrat6.f .
./readpvap.f: ../LIB/readpvap.f
	cp ../LIB/readpvap.f .
./spreaddep3_ncdf.f: ../LIB/spreaddep3_ncdf.f
	cp ../LIB/spreaddep3_ncdf.f .
./defopvals_ncdf.f: ../LIB/defopvals_ncdf.f
	cp ../LIB/defopvals_ncdf.f .
./map_indices.f90: ../LIB/map_indices.f90
	cp ../LIB/map_indices.f90 .
./mtrat_ba.f: ../LIB/mtrat_ba.f
	cp ../LIB/mtrat_ba.f .
./get_cheminp.f90: ../LIB/get_cheminp.f90
	cp ../LIB/get_cheminp.f90 .
./ncutil.f: ../LIB/ncutil.f
	cp ../LIB/ncutil.f .
./rjgrm.f: ../LIB/rjgrm.f
	cp ../LIB/rjgrm.f .
./read_dvsp_ncdf.f90: ../LIB/read_dvsp_ncdf.f90
	cp ../LIB/read_dvsp_ncdf.f90 .
./write_printphoto.f90: ../LIB/write_printphoto.f90
	cp ../LIB/write_printphoto.f90 .
./surface_rates.f90: ../LIB/surface_rates.f90
	cp ../LIB/surface_rates.f90 .
./setup_GECKO.f90: ../LIB/setup_GECKO.f90
	cp ../LIB/setup_GECKO.f90 .
./soa_equil.f90: ../LIB/soa_equil.f90
	cp ../LIB/soa_equil.f90 .
./spreadkey6.f90: ../LIB/spreadkey6.f90
	cp ../LIB/spreadkey6.f90 .
./findring.f: ../LIB/findring.f
	cp ../LIB/findring.f .
./spsplnb.f: ../LIB/spsplnb.f
	cp ../LIB/spsplnb.f .
./readprev.f: ../LIB/readprev.f
	cp ../LIB/readprev.f .
./wrtflags_ncdf.f: ../LIB/wrtflags_ncdf.f
	cp ../LIB/wrtflags_ncdf.f .
./spsplna.f: ../LIB/spsplna.f
	cp ../LIB/spsplna.f .
./akparameter_module.f90: ../LIB/akparameter_module.f90
	cp ../LIB/akparameter_module.f90 .
./sputil.f: ../LIB/sputil.f
	cp ../LIB/sputil.f .
./cdcase.f: ../LIB/cdcase.f
	cp ../LIB/cdcase.f .
./vbs_module.f90: ../LIB/vbs_module.f90
	cp ../LIB/vbs_module.f90 .
./newdt.f: ../LIB/newdt.f
	cp ../LIB/newdt.f .
./readprev_ncdf.f: ../LIB/readprev_ncdf.f
	cp ../LIB/readprev_ncdf.f .
./spreadj4_ncdf.f: ../LIB/spreadj4_ncdf.f
	cp ../LIB/spreadj4_ncdf.f .
./readkeyflags.f: ../LIB/readkeyflags.f
	cp ../LIB/readkeyflags.f .
./wall_test.f90: ../LIB/wall_test.f90
	cp ../LIB/wall_test.f90 .
./wall_Krech2016_J.f90: ../LIB/wall_Krech2016_J.f90
	cp ../LIB/wall_Krech2016_J.f90 .
./wall_Krech2016_Z.f90: ../LIB/wall_Krech2016_Z.f90
	cp ../LIB/wall_Krech2016_Z.f90 .
./wall_Huang2018.f90: ../LIB/wall_Huang2018.f90
	cp ../LIB/wall_Huang2018.f90 .
./wall_MatZie2010.f90: ../LIB/wall_MatZie2010.f90
	cp ../LIB/wall_MatZie2010.f90 .
./module_data_gecko_main.f90: ../LIB/module_data_gecko_main.f90
	cp ../LIB/module_data_gecko_main.f90 .
./spgetemi4.f: ../LIB/spgetemi4.f
	cp ../LIB/spgetemi4.f .
./open_op_files.f90: ../LIB/open_op_files.f90
	cp ../LIB/open_op_files.f90 .
./spdeposition3.f: ../LIB/spdeposition3.f
	cp ../LIB/spdeposition3.f .
./printreaction.f: ../LIB/printreaction.f
	cp ../LIB/printreaction.f .
./iter4.f: ../LIB/iter4.f
	cp ../LIB/iter4.f .
./map_sp_ids.f90: ../LIB/map_sp_ids.f90
	cp ../LIB/map_sp_ids.f90 .
./boxmod_main.f90: ../PROG/boxmod_main.f90
	cp ../PROG/boxmod_main.f90 .
./spchge_cvar4.f: ../LIB/spchge_cvar4.f
	cp ../LIB/spchge_cvar4.f .
./spakinit9_ncdf.f: ../LIB/spakinit9_ncdf.f
	cp ../LIB/spakinit9_ncdf.f .
./cleanup.f90: ../LIB/cleanup.f90
	cp ../LIB/cleanup.f90 .
./sptimesurf3.f: ../LIB/sptimesurf3.f
	cp ../LIB/sptimesurf3.f .
./rjgadd.f: ../LIB/rjgadd.f
	cp ../LIB/rjgadd.f .
./forcing_params_module.f90: ../LIB/forcing_params_module.f90
	cp ../LIB/forcing_params_module.f90 .
./mtrat.f: ../LIB/mtrat.f
	cp ../LIB/mtrat.f .
./wrtenvinp_ncdf.f: ../LIB/wrtenvinp_ncdf.f
	cp ../LIB/wrtenvinp_ncdf.f .
./surfdep_wesely89.f: ../LIB/surfdep_wesely89.f
	cp ../LIB/surfdep_wesely89.f .
./wrtopgas_ncdf.f: ../LIB/wrtopgas_ncdf.f
	cp ../LIB/wrtopgas_ncdf.f .
./manage_initial.f90: ../LIB/manage_initial.f90
	cp ../LIB/manage_initial.f90 .
./init_defaults.f90: ../LIB/init_defaults.f90
	cp ../LIB/init_defaults.f90 .
./constraints_module.f90: ../LIB/constraints_module.f90
	cp ../LIB/constraints_module.f90 .
./readpvap_nan_ncdf.f: ../LIB/readpvap_nan_ncdf.f
	cp ../LIB/readpvap_nan_ncdf.f .
./spakkextra4.f: ../LIB/spakkextra4.f
	cp ../LIB/spakkextra4.f .
./extract_reacrates.f: ../LIB/extract_reacrates.f
	cp ../LIB/extract_reacrates.f .
./mtdyn_rp.f: ../LIB/mtdyn_rp.f
	cp ../LIB/mtdyn_rp.f .
./general_module.f90: ../LIB/general_module.f90
	cp ../LIB/general_module.f90 .
./wrtjvinp_ncdf.f: ../LIB/wrtjvinp_ncdf.f
	cp ../LIB/wrtjvinp_ncdf.f .
./readpvap_sim_ncdf.f: ../LIB/readpvap_sim_ncdf.f
	cp ../LIB/readpvap_sim_ncdf.f .
./soa_dyn_update.f90: ../LIB/soa_dyn_update.f90
	cp ../LIB/soa_dyn_update.f90 .
./readpvap_nan_bin.f: ../LIB/readpvap_nan_bin.f
	cp ../LIB/readpvap_nan_bin.f .
./open_io_files.f90: ../LIB/open_io_files.f90
	cp ../LIB/open_io_files.f90 .
./wrtopaer_ncdf.f: ../LIB/wrtopaer_ncdf.f
	cp ../LIB/wrtopaer_ncdf.f .
./write_op_hdrs.f90: ../LIB/write_op_hdrs.f90
	cp ../LIB/write_op_hdrs.f90 .
./spreadro2sof2_bin.f: ../LIB/spreadro2sof2_bin.f
	cp ../LIB/spreadro2sof2_bin.f .
./chemmap.f: ../LIB/chemmap.f
	cp ../LIB/chemmap.f .
./getpvap_nan.f: ../LIB/getpvap_nan.f
	cp ../LIB/getpvap_nan.f .
./general.h: ../LIB/general.h
	cp ../LIB/general.h .
./fit.f: ../LIB/fit.f
	cp ../LIB/fit.f .
./spsolar.f: ../LIB/spsolar.f
	cp ../LIB/spsolar.f .
./read_idgaw.f: ../LIB/read_idgaw.f
	cp ../LIB/read_idgaw.f .
./nodmap.f: ../LIB/nodmap.f
	cp ../LIB/nodmap.f .
./module_chamber_tools.f90: ../LIB/module_chamber_tools.f90
	cp ../LIB/module_chamber_tools.f90 .
./pbl_height_forcage.f: ../LIB/pbl_height_forcage.f
	cp ../LIB/pbl_height_forcage.f .
./solve_box.f90: ../LIB/solve_box.f90
	cp ../LIB/solve_box.f90 .
./spakspnum.f: ../LIB/spakspnum.f
	cp ../LIB/spakspnum.f .
./spreadj4.f: ../LIB/spreadj4.f
	cp ../LIB/spreadj4.f .
./setup_ncdf_op.f90: ../LIB/setup_ncdf_op.f90
	cp ../LIB/setup_ncdf_op.f90 .
./inorganic_aer_module.f90: ../LIB/inorganic_aer_module.f90
	cp ../LIB/inorganic_aer_module.f90 .
./gettrack.f: ../LIB/gettrack.f
	cp ../LIB/gettrack.f .
./calc_psat.f90: ../LIB/calc_psat.f90
	cp ../LIB/calc_psat.f90 .
./bilan_species.f: ../LIB/bilan_species.f
	cp ../LIB/bilan_species.f .
./get_envinp.f90: ../LIB/get_envinp.f90
	cp ../LIB/get_envinp.f90 .
./readpvap_sim_bin.f: ../LIB/readpvap_sim_bin.f
	cp ../LIB/readpvap_sim_bin.f .
./readppf.f: ../LIB/readppf.f
	cp ../LIB/readppf.f .
./twostep5.f: ../LIB/twostep5.f
	cp ../LIB/twostep5.f .
./spakkfo.f: ../LIB/spakkfo.f
	cp ../LIB/spakkfo.f .
./get_prev.f90: ../LIB/get_prev.f90
	cp ../LIB/get_prev.f90 .
./spakinit9_bin.f: ../LIB/spakinit9_bin.f
	cp ../LIB/spakinit9_bin.f .
./prodloss_module.f90: ../LIB/prodloss_module.f90
	cp ../LIB/prodloss_module.f90 .
./flags_module.f90: ../LIB/flags_module.f90
	cp ../LIB/flags_module.f90 .
./spprint_bin.f: ../LIB/spprint_bin.f
	cp ../LIB/spprint_bin.f .
./calc_Cstar298.f90: ../LIB/calc_Cstar298.f90
	cp ../LIB/calc_Cstar298.f90 .
./spreaddep3.f: ../LIB/spreaddep3.f
	cp ../LIB/spreaddep3.f .
./calc_phot_OFR.f90: ../LIB/calc_phot_OFR.f90
	cp ../LIB/calc_phot_OFR.f90 .
./sorting_module.f90: ../LIB/sorting_module.f90
	cp ../LIB/sorting_module.f90 .
./output_vals.f90: ../LIB/output_vals.f90
	cp ../LIB/output_vals.f90 .
./defflags_ncdf.f: ../LIB/defflags_ncdf.f
	cp ../LIB/defflags_ncdf.f .
./soapartition.f90: ../LIB/soapartition.f90
	cp ../LIB/soapartition.f90 .
./read_dvsp_bin.f90: ../LIB/read_dvsp_bin.f90
	cp ../LIB/read_dvsp_bin.f90 .
./spreadcoeff3.f: ../LIB/spreadcoeff3.f
	cp ../LIB/spreadcoeff3.f .
./flags.h: ../LIB/flags.h
	cp ../LIB/flags.h .
./surface_emissions.f: ../LIB/surface_emissions.f
	cp ../LIB/surface_emissions.f .
./spreaddimer.f: ../LIB/spreaddimer.f
	cp ../LIB/spreaddimer.f .
./spreadro2jenk.f: ../LIB/spreadro2jenk.f
	cp ../LIB/spreadro2jenk.f .
./printsteadystate.f90: ../LIB/printsteadystate.f90
	cp ../LIB/printsteadystate.f90 .
./find_concs.f90: ../LIB/find_concs.f90
	cp ../LIB/find_concs.f90 .
./soa_dyn.f90: ../LIB/soa_dyn.f90
	cp ../LIB/soa_dyn.f90 .
./spdatsurf.f: ../LIB/spdatsurf.f
	cp ../LIB/spdatsurf.f .
./organic.h: ../LIB/organic.h
	cp ../LIB/organic.h .
./phys_rates.f90: ../LIB/phys_rates.f90
	cp ../LIB/phys_rates.f90 .
SRC = ../LIB/surfdep_wesely89.f ../LIB/get_cheminp.f90 ../LIB/findring.f ../LIB/write_op_hdrs.f90 ../LIB/calc_psat.f90 ../LIB/flags_module.f90 ../LIB/calc_Cstar298.f90 ../LIB/spreadro2sof2_bin.f ../LIB/getpvap_nan.f ../LIB/soapartition.f90 ../LIB/sorting_module.f90 ../LIB/fit.f ../LIB/defflags_ncdf.f ../LIB/readpvap_nan_bin.f ../LIB/init_defaults.f90 ../LIB/module_chamber_tools.f90 ../LIB/get_envinp.f90 ../LIB/wrtjvinp_ncdf.f ../LIB/spsplna.f ../LIB/spreadro2sof2_ncdf.f ../LIB/readpvap_sim_ncdf.f ../LIB/open_io_files.f90 ../LIB/spgetemi4.f ../LIB/wrtflags_ncdf.f ../LIB/spakspnum.f ../LIB/extract_reacrates.f ../LIB/readppf.f ../LIB/spreadcoeff3.f ../LIB/spreaddimer.f ../LIB/vbs_module.f90 ../LIB/readprev.f ../LIB/soa_equil.f90 ../LIB/spakinit9_ncdf.f ../LIB/spdeposition3.f ../LIB/wrtopgas_ncdf.f ../LIB/spchge_cvar4.f ../LIB/soa_dyn.f90 ../LIB/spprint_bin.f ../PROG/boxmod_main.f90 ../LIB/manage_initial.f90 ../LIB/apply_constraints.f ../LIB/iter4.f ../LIB/readkeyflags.f ../LIB/nodmap.f ../LIB/ncutil.f ../LIB/open_op_files.f90 ../LIB/chemmap.f ../LIB/bilan_species.f ../LIB/solve_box.f90 ../LIB/surface_rates.f90 ../LIB/mtdyn_rp.f ../LIB/map_indices.f90 ../LIB/get_prev.f90 ../LIB/phys_rates.f90 ../LIB/twostep5.f ../LIB/wrtopaer_ncdf.f ../LIB/wrtenvinp_ncdf.f ../LIB/read_dvsp_bin.f90 ../LIB/mtrat.f ../LIB/sputil.f ../LIB/forcing_params_module.f90 ../LIB/cleanup.f90 ../LIB/wall_test.f90 ../LIB/wall_Krech2016_Z.f90 ../LIB/defenvinp_ncdf.f90 ../LIB/spreaddep3_ncdf.f ../LIB/read_dvsp_ncdf.f90 ../LIB/wall_Huang2018.f90 ../LIB/surface_emissions.f ../LIB/readprev_ncdf.f ../LIB/readpvap_sim_bin.f ../LIB/mtrat_ba.f ../LIB/setup_ncdf_op.f90 ../LIB/module_data_gecko_main.f90 ../LIB/newdt.f ../LIB/spakkrat6.f ../LIB/spinterp5.f ../LIB/constraints_module.f90 ../LIB/wall_MatZie2010.f90 ../LIB/spreadkey6.f90 ../LIB/spreaddep3.f ../LIB/output_vals.f90 ../LIB/manage_loop.f90 ../LIB/spdatsurf.f ../LIB/akparameter_module.f90 ../LIB/write_printphoto.f90 ../LIB/spakkfo.f ../LIB/defopvals_ncdf.f ../LIB/spforcage6.f90 ../LIB/rjgadd.f ../LIB/printsteadystate.f90 ../LIB/printreaction.f ../LIB/spreadj4_ncdf.f ../LIB/spsplnb.f ../LIB/wall_Krech2016_J.f90 ../LIB/calc_phot_OFR.f90 ../LIB/map_sp_ids.f90 ../LIB/general_module.f90 ../LIB/inorganic_aer_module.f90 ../LIB/spreadro2jenk.f ../LIB/soa_dyn_update.f90 ../LIB/spsolar.f ../LIB/cdcase.f ../LIB/gettrack.f ../LIB/readpvap_nan_ncdf.f ../LIB/setup_GECKO.f90 ../LIB/readpvap.f ../LIB/sptimesurf3.f ../LIB/prodloss_module.f90 ../LIB/pbl_height_forcage.f ../LIB/spreadj4.f ../LIB/find_concs.f90 ../LIB/spakinit9_bin.f ../LIB/setup_printphoto.f90 ../LIB/spakkextra4.f ../LIB/rjgrm.f ../LIB/read_idgaw.f ../LIB/general.h ../LIB/flags.h ../LIB/organic.h
OBJ = surfdep_wesely89.o get_cheminp.o findring.o write_op_hdrs.o calc_psat.o flags_module.o calc_Cstar298.o spreadro2sof2_bin.o getpvap_nan.o soapartition.o sorting_module.o fit.o defflags_ncdf.o readpvap_nan_bin.o init_defaults.o module_chamber_tools.o get_envinp.o wrtjvinp_ncdf.o spsplna.o spreadro2sof2_ncdf.o readpvap_sim_ncdf.o open_io_files.o spgetemi4.o wrtflags_ncdf.o spakspnum.o extract_reacrates.o readppf.o spreadcoeff3.o spreaddimer.o vbs_module.o readprev.o soa_equil.o spakinit9_ncdf.o spdeposition3.o wrtopgas_ncdf.o spchge_cvar4.o soa_dyn.o spprint_bin.o boxmod_main.o manage_initial.o apply_constraints.o iter4.o readkeyflags.o nodmap.o ncutil.o open_op_files.o chemmap.o bilan_species.o solve_box.o surface_rates.o mtdyn_rp.o map_indices.o get_prev.o phys_rates.o twostep5.o wrtopaer_ncdf.o wrtenvinp_ncdf.o read_dvsp_bin.o mtrat.o sputil.o forcing_params_module.o cleanup.o wall_test.o wall_Krech2016_Z.o defenvinp_ncdf.o spreaddep3_ncdf.o read_dvsp_ncdf.o wall_Huang2018.o surface_emissions.o readprev_ncdf.o readpvap_sim_bin.o mtrat_ba.o setup_ncdf_op.o module_data_gecko_main.o newdt.o spakkrat6.o spinterp5.o constraints_module.o wall_MatZie2010.o spreadkey6.o spreaddep3.o output_vals.o manage_loop.o spdatsurf.o akparameter_module.o write_printphoto.o spakkfo.o defopvals_ncdf.o spforcage6.o rjgadd.o printsteadystate.o printreaction.o spreadj4_ncdf.o spsplnb.o wall_Krech2016_J.o calc_phot_OFR.o map_sp_ids.o general_module.o inorganic_aer_module.o spreadro2jenk.o soa_dyn_update.o spsolar.o cdcase.o gettrack.o readpvap_nan_ncdf.o setup_GECKO.o readpvap.o sptimesurf3.o prodloss_module.o pbl_height_forcage.o spreadj4.o find_concs.o spakinit9_bin.o setup_printphoto.o spakkextra4.o rjgrm.o read_idgaw.o
OFF = ../LIB/manage_loop.f90 ../LIB/setup_printphoto.f90 ../LIB/apply_constraints.f ../LIB/defenvinp_ncdf.f90 ../LIB/spreadro2sof2_ncdf.f ../LIB/spforcage6.f90 ../LIB/spinterp5.f ../LIB/spakkrat6.f ../LIB/readpvap.f ../LIB/spreaddep3_ncdf.f ../LIB/defopvals_ncdf.f ../LIB/map_indices.f90 ../LIB/mtrat_ba.f ../LIB/get_cheminp.f90 ../LIB/ncutil.f ../LIB/rjgrm.f ../LIB/read_dvsp_ncdf.f90 ../LIB/write_printphoto.f90 ../LIB/surface_rates.f90 ../LIB/setup_GECKO.f90 ../LIB/soa_equil.f90 ../LIB/spreadkey6.f90 ../LIB/findring.f ../LIB/spsplnb.f ../LIB/readprev.f ../LIB/wrtflags_ncdf.f ../LIB/spsplna.f ../LIB/akparameter_module.f90 ../LIB/sputil.f ../LIB/cdcase.f ../LIB/vbs_module.f90 ../LIB/newdt.f ../LIB/readprev_ncdf.f ../LIB/spreadj4_ncdf.f ../LIB/readkeyflags.f ../LIB/wall_Krech2016_J.f90 ../LIB/module_data_gecko_main.f90 ../LIB/spgetemi4.f ../LIB/open_op_files.f90 ../LIB/spdeposition3.f ../LIB/printreaction.f ../LIB/iter4.f ../LIB/map_sp_ids.f90 ../PROG/boxmod_main.f90 ../LIB/spchge_cvar4.f ../LIB/spakinit9_ncdf.f ../LIB/cleanup.f90 ../LIB/sptimesurf3.f ../LIB/rjgadd.f ../LIB/forcing_params_module.f90 ../LIB/mtrat.f ../LIB/wrtenvinp_ncdf.f ../LIB/surfdep_wesely89.f ../LIB/wrtopgas_ncdf.f ../LIB/manage_initial.f90 ../LIB/init_defaults.f90 ../LIB/constraints_module.f90 ../LIB/readpvap_nan_ncdf.f ../LIB/spakkextra4.f ../LIB/extract_reacrates.f ../LIB/wall_Krech2016_Z.f90 ../LIB/mtdyn_rp.f ../LIB/general_module.f90 ../LIB/wrtjvinp_ncdf.f ../LIB/readpvap_sim_ncdf.f ../LIB/soa_dyn_update.f90 ../LIB/readpvap_nan_bin.f ../LIB/open_io_files.f90 ../LIB/wrtopaer_ncdf.f ../LIB/write_op_hdrs.f90 ../LIB/spreadro2sof2_bin.f ../LIB/chemmap.f ../LIB/getpvap_nan.f ../LIB/general.h ../LIB/fit.f ../LIB/spsolar.f ../LIB/read_idgaw.f ../LIB/nodmap.f ../LIB/module_chamber_tools.f90 ../LIB/pbl_height_forcage.f ../LIB/solve_box.f90 ../LIB/spakspnum.f ../LIB/spreadj4.f ../LIB/setup_ncdf_op.f90 ../LIB/inorganic_aer_module.f90 ../LIB/gettrack.f ../LIB/calc_psat.f90 ../LIB/bilan_species.f ../LIB/get_envinp.f90 ../LIB/readpvap_sim_bin.f ../LIB/readppf.f ../LIB/twostep5.f ../LIB/spakkfo.f ../LIB/get_prev.f90 ../LIB/spakinit9_bin.f ../LIB/prodloss_module.f90 ../LIB/flags_module.f90 ../LIB/spprint_bin.f ../LIB/calc_Cstar298.f90 ../LIB/spreaddep3.f ../LIB/calc_phot_OFR.f90 ../LIB/sorting_module.f90 ../LIB/wall_Huang2018.f90 ../LIB/output_vals.f90 ../LIB/defflags_ncdf.f ../LIB/soapartition.f90 ../LIB/read_dvsp_bin.f90 ../LIB/spreadcoeff3.f ../LIB/flags.h ../LIB/surface_emissions.f ../LIB/spreaddimer.f ../LIB/spreadro2jenk.f ../LIB/printsteadystate.f90 ../LIB/find_concs.f90 ../LIB/soa_dyn.f90 ../LIB/spdatsurf.f ../LIB/organic.h ../LIB/wall_MatZie2010.f90 ../LIB/phys_rates.f90
clean: neat
	-rm -f .../PROG/BOXMOD.cppdefs $(OBJ) ../PROG/BOXMOD

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

../PROG/BOXMOD: $(OBJ)
	$(LD) $(OBJ) -o ../PROG/BOXMOD  $(LDFLAGS)

