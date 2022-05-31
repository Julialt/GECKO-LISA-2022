# Setup options for the postprocessing program

The postprocessor will run calculations based on the flags selected in the input namelist.
The user can create a namelist like the `postproc_flags.input` and provide it as input to `./run_postprocessing.bash`
If no option is turned on, the program will not bother reading input files.

Here is a breakdown of the meaning of each available option.

Options starting with flag_ are turned on and off with `.true.` and `.false.`
e.g.: `flag_carbonchain = .true.`
array inputs need to provide the array indices after the name
e.g., note `(1:2)` following the name of the variable: `selected_species(1:2) = "GAR0084", "GO3"`

output file names are indicated below. if extracting for multiple boxes, filenames are prefixed with `bottom` and `top`

## Namelist Inputs

- `input_type`: 
   - what kind of file is the program reading. No reason to set this to anything other than "netcdf"
- `flag_amsfactors`
   - output file: `ams_factors_aer_ug.csv`
   - outputs mass contribution of identified mass factors in De Sà et al., 2018 (MO-OOA, LO-OOA, IEPOX-SOA, ADOA, BBOA, HOA)
   - _may be broken_
- `flag_atomratios`
   - output files: `atom_ratios_gas.csv`,`atom_ratios_aer.csv`
   - output the time evolution of `O/C`, `H/C` and `N/C` ratios for the gas and aerosol phases
- `flag_bubble`
   - output file: `bubble_mass_lifetime.csv`
   - outputs 2D mass distribution as a function of average carbon oxidation state (OSc) and pvap (atm)
   - if only one precursor, output timesteps for the first 5 lifetimes of the precursors
   - in other situation, the 5 output timesteps are equidistant times between start and end
   - the resolution of the vapor pressure spectrum is 1 atm, hardcoded at the moment but could be a user input one day
   - the resolution of the OSc spectrum is 0.2, hardcoded at the moment but could be a user input one day
- `flag_carbonchain`:
   - output files: `carbon_chain_gas_ppbC_time.csv`, `carbon_chain_aer_ppbC_time.csv`
   - outputs the distribution of length (C1,C2,C3...) for the gas and aerosol phase in ppbC
- `flag_chochonfreq`
   - output files: `chochonfreq_CHO.csv`,`chochonfreq_CHNO.csv`,`chochonfreq_CHOS.csv`,`chochonfreq_CHNOS.csv`
   - outputs frequency of species with CHO, CHNO, CHOS, CHNOS elemental composition as a function of number of carbon atoms
   - _may be broken_
- `flag_chon`:
   - output files: `top_chon_ppbc_gas_time.csv`, `top_chon_ppbc_aer_time.csv`
   - output the mixing ratios of the most important elemental formulas( e.g.: C5H5O2, etc...) in gas and aerosol phase
- `flag_contributingspecs`:
   - output files: `contributing_species_molec_gas_time.csv`,`contributing_species_molec_aer_time.csv`
   - output for each time step, how many species are needed to represent 30% and 90% on a molecular basis in gas and aerosol phase
- `flag_dbeai` 
   - output files: `dbe.csv`,`ai.csv`
   - outputs mass weighted averages of double bond equivalent and aromaticity index in every phase
- `flag_dyn_filter`
   - output files: `gas_filtered_*.csv`, `aer_filtered_*.csv`
   - outputs number concentration, mass concentration and mixing ratios of sum of species defined by filters listed in `gas_dyn_filter` and `aer_dyn_filter`
- `flag_elementscontrib`
   - output files: `elements_contribution_gas_time.csv`, `elements_contribution_aer_time.csv`
   - output the contribution of individual elements (C, H, O, N) on an atom basis
- `flag_entropy`
   - output file: `aerosol_entropy.csv`
   - outputs the first and second order informational entropy in aerosol particles  according to Riemer and West (2013)
- `flag_functions`:
   - output files: `rofc_gas_time.csv`,`rofc_aer_time.csv`,`functions_mass_gas_time.csv`,`functions_mass_aer_time.csv`
   - output the distribution of organic functionalities as number of functions per carbon atom (`rofc_` files) 
     and the total mass of species bearing each function (double counting happens, `functions_mass_` files)
- `flag_henry`
   - output files: `Henry_Matm_distribution_ppbC_lifetime.csv`, `Henry_Matm_distribution_ug_lifetime.csv`
   - output the ppbC and mass distribution of Henry's law constants (in M/atm) in every phase
   - if only one precursor, output timesteps for the first 5 lifetimes of the precursors
   - in other situation, the 5 output timesteps are equidistant times between start and end
   - the resolution of the Henry's law constants is 2 M/atm, hardcoded at the moment but could be a user input one day
- `flag_massspectrum`
   - output file: `mass_spectrum_lifetime.csv`
   - outputs molecular mass number distribution in every phase
   - if only one precursor, output timesteps for the first 5 lifetimes of the precursors
   - in other situation, the 5 output timesteps are equidistant times between start and end
   - the resolution of the mass pressure spectrum is 1 g/mol, hardcoded at the moment but could be a user input one day
- `flag_nitrates`
   - output file: `nitrate_func_gas.csv`, `nitrate_func_aer.csv`
   - output mass distribution and functions number distribution of nitrates in gas and aerosol phases
   - nitrates are classified as primary, secondary and tertiary
- `flag_NO3_exposure`
   - output file: `NO3_exposure.csv`
   - outputs cumulative exposure to NO3 in molec cm-3 s
- `flag_NO3R`
   - output file: `NO3_reactivity.csv`, `NO3_reac_prods.csv`, `NO3_top_spc_reac_time.csv`
   - outputs timeseries of total and non-precursor NO3 reactivity (s-1)
   - outputs timeseries of reactivities of top `n_topspecies` species with respect to NO3 (s-1)
- `flag_O3R`
   - output file: `O3_reactivity.csv`, `O3_reac_prods.csv`, `O3_top_spc_reac_time.csv`
   - outputs timeseries of total and non-precursor O3 reactivity (s-1)
   - outputs timeseries of reactivities of top `n_topspecies` species with respect to O3 (s-1)
- `flag_OHR`
   - output file: `OH_reactivity.csv`, `OH_reac_prods.csv`, `OH_top_spc_reac_time.csv`
   - outputs timeseries of total and non-precursor OH reactivity (s-1)
   - outputs timeseries of reactivities of top `n_topspecies` species with respect to OH (s-1)
- `flag_oh_exposure`
   - output file: `oh_exposure.csv`
   - outputs cumulative exposure to oh radicals in molec cm-3 s
- `flag_phasedist`:
   - output files: `carb_dist_time.csv`, `mass_dist_time.csv`
   - output the total mixing ratio (in ppbC) and the total mass concentration(µg m-3) in each phase (GAS, AEROSOL, PRECURSOR, CO-CO2, INORG, ...)
- `flag_phasestate`
   - output file: `aerosol_glass_transition_temp.csv`
   - estimates dry glass transition temperature in the aeroosol phase and viscosity (Pa s) following Shiraiwa et al (2017).
   - updated following DeRieux et al., (2018)
   - This is a crude estimate but it gives an idea of the phase state of the aerosol particles
- `flag_pvap`
   - output files:`pvap_atm_distribution_ppbC_lifetime.csv`,`pvap_atm_distribution_ug_lifetime.csv`
   - output the ppbC and mass distribution of vapor pressures (in atm) in every phase
   - if only one precursor, output timesteps for the first 5 lifetimes of the precursors
   - in other situation, the 5 output timesteps are equidistant times between start and end
   - the resolution of the vapor pressure spectrum is 2 atm, hardcoded at the moment but could be a user input one day
- `flag_rates`:
   - output file: `top_rates_{rate_species}_time.csv`
   - this will output the `n_toprates` reaction rate timeseries for the production and loss of the specified `rate_species`
   - units are molec cm-3 s-1
- `rate_species`:
   - list the ONE species for which you want to output rates if `flag_rates` is on. Need to give dictionary code name.
   - the output file contains code names of up to two reactants in the header
- `flag_selected`:
   - output file: `selected_species_molec_time.csv`
   - this will output the number concentration of species listed in the `selected_species` array
- `selected_species`:
   - output file: `selected_species_molec_time.csv`
   - list the species you want to output if `flag_selected` is on. Need to give dictionary code name.
   - the output file contains both code names and formulas in header
- `flag_soayield`:
   - output file: `soa_yield_time.csv`
   - output the soa yield for the simulation. Makes sense if there is only on precursor and no emission
- `flag_topspec`:
   - output files:`top_species_ppbc_gas_time.csv`, `top_species_ppbc_aer_time.csv`
   - output the `n_topspecies` most important species on a ppbC basis, for the gas and aerosol phases
- `n_toprates`: (integer)
   - number of reactionss to account for `flag_rates` flag
- `n_topspecies`: (integer)
   - number of species to account for `flag_topspec` flag
- `gas_dyn_filter` and `aer_dyn_filter`: (array of character strings)
   - output files: `gas_filtered_*.csv`, `aer_filtered_*.csv`
   - output number concentration, mass concentration and mixing ratios of sum of species defined by these filters
   - each variable contains up to 100 different filters (one output file for each
   - each filter is a character string defined as follows
     - three sections are defined in the string, separated by `;` and `!`. Each section is optional but these 2 characters need to be present
     - the first section (beginning to `;`) defines atomic contraints: list atoms (CHNO) followed by number:
         - `C5H7` selects species with 5 carbon atoms and 7 hydrogen atoms. The number of nitrogen and oxygen atoms is not constrained in this example
         - `C5H7N0` selects species with 5 carbon atoms, 7 hydrogen atoms and no nitrogen atom. The number of oxygen atoms is not constrained in this example
     - the second section (between `;` and `!`)selects a minimum number of organic functions defined by their code(ABDEFGHKLNOPRTUVX1234S, defined in the [**list of functions**](https://github.com/NCAR/GECKO-A/tree/ML2019/POSTPROCESSING#list-of-functions)
       - `NN` selects species that have 2 nitrates functions or more, other functions are not contrained
       - `NNO` selects species that have 2 nitrates functions or more and 1 hydroxy function or more, other functions are not contrained
     - the last section (from `!` to end) selects maximum (not inclusive) number of organic functions defined by their code(ABDEFGHKLNOPRTUVX1234S, defined in the [**list of functions**](https://github.com/NCAR/GECKO-A/tree/ML2019/POSTPROCESSING#list-of-functions)
       - `NN` selects species with 0 or 1 nitrate function, other functions are not contrained
       - `NNO` selects species with 0 or 1 nitrate function and 0 hydroxy function, other functions are not contrained
   - the three sections are combined to establish the filter that will be applied. Examples: 
     - `C5;NN!`: selects species with exactly 5 carbon atoms and at least 2 nitrate functions. Note there is no constraint on the max number of any function
     - `;N!NN`: selects species with at least 1 nitrate function and less than 2 nitrate function, i.e. species with only one nitrate function (but can have any other function". Note the absence of atom constraints
     - `C5;!OO`: selects species with exactly 5 carbon atoms and less than 2 hydroxy function (0 or 1). Note the absence of a minimum number of any organic function
     - `C6;OH!HH`: selects species with exactly 5 carbon atom, at least one hydroxy and one hydroperoxide functions and less that 2 hydroperoxide functions
## List of functions  
 - A: carboxylic acid -CO(OH)
 - B: carboxylate -CO(Om)
 - D: aldehyde -CHO
 - E: ether -O-
 - F: fluoride -F
 - G: percarboxylic acid -CO(OOH)
 - H: hydroperoxide -OOH
 - K: ketone -CO-
 - L: cloride -Cl
 - N: nitrate -ONO2
 - O: hydroxy -OH
 - P: peroxyacynnitrate -CO(OONO2)
 - R: aromatic cycle
 - T: aliphatic cycle
 - U: double bond
 - V: nitroso -NO2
 - X: ketene -C=C=O
 - 1: alcoxy radical -O.
 - 2: peroxy radical -OO.
 - 3: acylperoxy radical -CO(OO.)
 - 4: criegee radical -C.(OO.)
 - S: sulfate -OSO3

