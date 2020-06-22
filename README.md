# IPT
Input Processing Tools for MUSICA

Detailed discussion will be provided on the MUSICA wiki. [https://wiki.ucar.edu/display/MUSICA/MUSICA+Home]

This collection of scripts and code are for the production of input fields that are required to run MUSICA-V0.  They target emissions, meterological reanalyses for nudging, and initial conditions.


## Emissions

### Fire (Biomass burning)

1. FINN (Fire INventory from NCAR)

A tool for regridding FINN v1.5 and v2.0 \*.txt files to FV, SE and SE-RR are available. The tool includes a README file, which describes in detail how to use the code, and several example namelists files (\*.inp), which can be easily adapted as needed.
FINN fire emissions are available between 2002 and 2020 from http://bai.acom.ucar.edu/Data/fire/

2. QFED (Quick Fire Emissions Dataset)

QFED fire emissions are available between 2000 and 2020 from NASA Goddard/GMAO  https://gmao.gsfc.nasa.gov/research/science_snapshots/global_fire_emissions.php

### Anthropogenic

1. EPA (Environmental Protection Agency)

NEI emissions for the U.S. are only available for the year 20XX, but scripts are available to correct them for other years.
Source data can be found [somewhere over the rainbow](http://somewhere_over_the_rainbow)

2. CAMS (Copernicus Atmosphere Monitoring Service)

Global CAMS anthropogenic emissions are available between 2000 and 2020 from [somewhere over the rainbow](http://somewhere_over_the_rainbow)

Interpolate original CAMS emissions to the desired grid: regrid\_fv2se\_cams\_anthro.ncl

Rename and calculated regridded CAMS emissions to be able to run in CESM2: rename\_cams\_anthro\_se.ncl

### CMIP6 (Climate Model Intercomparison Project)

CMIP6 standard historical anthropogenic and biomass burning emissions are available for 1750-2014. 
Interpolate original CMIP6 emissions to the desired grid: regrid\_fv2se\_cmip6\_main.ncl

Rename and calculate MOZART species to be able to run in CESM2: 

Source data is available [somewhere over the rainbow](http://somewhere_over_the_rainbow)

## Meteorological Data Regridding

MERRA2 original data script to regrid to the desired resolution

MERRA2 reanalysis between 1975 and 2020 on the original horizontal and vertical grid are available [here](https://rda.ucar.edu/datasets/ds313.3/) or local on cheyenne

## Initial Conditions Regridding

Atmospheric initial condition files are regridded using a script that allows us to regrid from the standard spectral element 1 degree resolution to the new grid. 
/glade/u/home/tilmes/ncl/SE/regrid\_all\_spectral\_data.ncl

To run with CESM, the resulting file has to be converted to a cdf5 format:

`nccopy -k cdf5 oldfile newfile`

Source data can be found [somewhere over the rainbow](http://somewhere_over_the_rainbow)
