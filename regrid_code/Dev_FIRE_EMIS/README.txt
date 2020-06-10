This code produces regridded FINN emissions for CESM/CAM-Chem Spectral Element (SE)
and the Model for Prediction Across Scales (MPAS) models.

/////////////////////////////////////////////////////////////////////////////////////

Section 1: Directory Strucuture

The src directory stores all model code
The tst directory is where the code is run and includes example namelists

/////////////////////////////////////////////////////////////////////////////////////

Section 2: To Run the Code

Prior to compiling the code, you will need to add several lines to your .bashrc file. Please adapt
the following code to your system and shell environment. The following example is for a bash shell
and NCAR's modeling1.acom.ucar linux machine. 

export FC=pgf95
export CC=pgcc
export NETCDF_DIR=/usr/local/netcdf-4.7.0

Now to compile the code, go to the src directory.
Run
./make_fire_emis

If the code has compiled correctly, you should see the following lines:

++++++++++++++++++
fire_emis build Ok
++++++++++++++++++

Now to run the code, go to the tst directory
Update the namelist as needed.

Some example namelists are provided:
example.camse_ne30.inp -> example namelist to produce FINN emissions on the Spectral Element ne30 grid for TS1 chemical mechanism

example.camse_ne0CONUS.inp -> example namelist to produce FINN emissions on the Spectral Element ne0CONUS grid for TS1 chemical mechanism

example.mpas.inp -> example namelist to produce FINN emissions on the MPAS grid for the TS1 chemical mechanism

The namelist parameters are described briefly here:

fire_directory   || The directory storing the raw FINN data files.
fire_filename(1) || The file name of the raw FINN data (must be .txt format). You can optionally add up to 5 files
                    here by adding a new line fire_filename(2), etc. (See namelist example for SE ne30 resolution).
mdlFilenm        || The file describing the grid (must be .nc format), see examples provided. 
mdlDir           || The directory storing the grid data files. Default is the current working directory ('./').
start_date       || The start date in YYYY-MM-DD format.
end_date         || The end date in YYYY-MM-DD format.
output_timing    || The output timing. Only 'daily' is currently available.
Model            || Model dycore: for CAMSE this must be 'CAMSE' and for MPAS this must be 'MPAS'.
resol            || The resolution of the model (only used for naming the output files).
EmisType         || The emission type (only used for naming the output files).
FinnVers         || The Finn Version (only used for naming the output files).
defaultUnits     || The default units (molecules/cm^2/s). This namelist value is only used for
                    adding unit attributes to the ouput files. If you want to convert to different units,
		    you must update this name and use the glb2fire_map namelist parameter with the
		    appropriate scalling factors for unit conversion. All unit conversions beyond
		    molecules/cm^2/s are done in this namelist, and not by the code itself (see examples).
		    If you want to change units for only certain species, use the
		    optional_non_default_unit_attribute in the glb2fire_map namelist parameter, as
		    explained below.

glb2fire_map     || The mapping description (See example files and Section 3 below for more detail).

                    Generally the format for gases is
                    'OUTPUT_SPECIES_NAME(optional_non_default_unit_attribute)->SCALING_FACTOR*INPUT_SPECIES_NAME'
		    Note that for gases the units in FINN are moles/day and the code will
		    automatically convert to molecules/cm^2/s

                    Generally the format for aerosols is 'OUTPUT_SPECIES_NAME(optional_non_default_unit_attribute)->SCALING_FACTOR*INPUT_SPECIES_NAME;aerosol'
		    
		    Note that for aerosols, the units in FINN are Kg/day and the code will
		    automatically convert to molecules/cm^2/s as long as you add ';aerosol',
		    as shown above.

Then once you are satisfied with the namelist parameters. To run the code,
./fire_emis < example.camse_ne30.inp > example.camse_ne30.out

Where example.camse_ne30.inp is your input namelist and example.camse_ne30.out stores the output
description.

Once completed, at the end of example.camse_ne30.out you should see:

 =================================
 fire_emis: Completed successfully
 =================================

/////////////////////////////////////////////////////////////////////////////////////

Section 3: Example and description of glb2fire_map namelist variable for MOZART-TS1 mechanism for CESM/CAM-Chem SE

glb2fire_map = 'CO -> CO',
               'NO -> NO',
	       'NO2 -> NO2',
	       'SO2 -> SO2',
	       'NH3 -> NH3',
	       'BIGALD -> BIGALD',
	       'BIGALK -> BIGALK',
	       'BIGENE -> BIGENE',
	       'MTERP -> C10H16',
	       'C2H4 -> C2H4',
	       'C2H5OH -> C2H5OH',
	       'C2H6 -> C2H6',
	       'C3H6 -> C3H6',
	       'C3H8 -> C3H8',
	       'CH2O -> CH2O',
	       'CH3CHO -> CH3CHO',
	       'CH3COCH3 -> CH3COCH3',
	       'CH3COCHO -> CH3COCHO',
	       'CH3COOH -> CH3COOH',
	       'CH3OH -> CH3OH',
	       'CRESOL -> CRESOL',
	       'GLYALD -> GLYALD',
	       'HYAC -> HYAC',
	       'ISOP -> ISOP',
	       'MACR -> MACR',
	       'MEK -> MEK',
	       'MVK->MVK',
	       'HCN->HCN',
	       'CH3CN -> CH3CN',
	       'TOLUENE->0.33333333*TOLUENE',
	       'BENZENE->0.33333333*TOLUENE',
	       'XYLENES->0.33333333*TOLUENE',
	       'pom_a4->1.4*OC;aerosol',
	       'bc_a4->BC;aerosol',
	       'HCOOH->HCOOH',
	       'C2H2->C2H2',
	       'num_bc_a4((particles/cm2/s)(molecules/mole)(g/kg))->5.60298303e18*BC;aerosol',
	       'num_pom_a4((particles/cm2/s)(molecules/mole)(g/kg))->1.33350996e19*OC;aerosol',
	       'svoc->0.03251613*OC;aerosol',
               'ivoc->0.04565217*C3H6+0.04782609*C3H8+0.03260870*C2H6+0.03043478*C2H4+0.06086957*BIGENE+0.0782609*BIGALK+0.06304348*CH3COCH3+0.0782609*MEK+0.04782609*CH3CHO+0.03260870*CH2O+0.1*TOLUENE'

Note the aerosol particle number variables (num_bc_a4 and num_pom_a4)
have different units ((particles/cm2/s)(molecules/mole)(g/kg)) and the
formulas for their calculation are provided below:

General formula: emis_num = emis_mol * mw / (rho * (PI/6) * (diam)^3)

For num_bc_a4
emis_num_BC = emis_mol_BC * 12. / (1700. * (PI/6) * (0.134e-6)^3)
emis_num_BC = emis_mol_BC * 5.60298303e18

For num_pom_a4
emis_num_OC = emis_mol_OC * 1.4 * 12. / (1000. * (PI/6) * (0.134e-6)^3)
emis_num_OC = emis_mol_OC * 1.33350996e19

The formulas to calculate SVOC and IVOC are provided below:


IVOC = 0.2*42./184.*C3H6 + 0.2*44./184.*C3H8 + 0.2*30./184.*C2H6 + 0.2*28./184.*C2H4
       + 0.2*56./184.*BIGENE + 0.2*72./184.*BIGALK + 0.2*58./184.*CH3COCH3
       + 0.2*72./184.*MEK + 0.2*44./184.*CH3CHO + 0.2*30./184.*CH2O + 0.2*78./184.*0.33333333*TOLUENE
       + 0.2*92./184.*0.3333*TOLUENE + 0.2*106./184.*0.3333*TOLUENE

IVOC = 0.04565217*C3H6 + 0.04782609*C3H8 + 0.03260870*C2H6 + 0.03043478*C2H4
       + 0.06086957*BIGENE + 0.0782609*BIGALK + 0.06304348*CH3COCH3
       + 0.0782609*MEK + 0.04782609*CH3CHO + 0.03260870*CH2O + 0.1*TOLUENE

Here for IVOC,
BIGALK = pentane,
BENZENE (0.33333333*TOLUENE) = benzene,
TOLUENE (0.33333333*TOLUENE) = toluene,
XYLENE (0.33333333*TOLUENE) = xylene


SVOC = 0.6*12./310.*1.4*OC;aerosol

SVOC = 0.03251613*OC;aerosol
