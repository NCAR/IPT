# IPT
Input Processing Tools for MUSICA

A collection of scripts and programs to process the emissions, meterological data and initial conditions files required for running CAM-chem-SE-RR (regionally refined grid of the Spectral Element configuration of CESM), also known as MUSICA-V0.

Further explanations of these tools are at https://wiki.ucar.edu/display/MUSICA/MUSICA+Home

1. Emissions
  a. FINN fire emissions: grid_fire_emis_finn_v2000 grids text files of emissions for each fire
  
  b. QFED fire emissions: regrid_fv2se_qfed_bb.ncl regrids global regular grid files to SE or SE-RR
  
  c. EPA U.S. anthropogenic emissions: anthro_epa_camse (and anthro_epa_camfv)
  
  d. CAMS global anthropogenic emissions: ncl scripts
  
  e. CMIP6 global anthropogenic and biomass burning emissions: ncl scripts
  
2. Meteorological data for nudging

3. Initial conditions
 
