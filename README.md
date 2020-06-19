# IPT
Input Processing Tools for MUSICA
Detailed discussion is in the [wiki](https://github.com/NCAR/EMIT/wiki)

This collection of scripts supports the production of input fields that are required to run MUSICA V0, with spectral element horizontal grids that have not been produced within the release. 

Surface and external emissions provided by the model have to be interpolated to the horizontal grid that is used for the simulations before the model can run. A conservative regridding within CESM2 is not available at this point. This collection of scripts include offline emission processing. Different types of emissions are available for regridding for different periods, including global CAMS anthropogenic emissions are available between 2000 and 2010, NEI emissions are only available only for the year 20XX, but scripts are available to correct them for  other years, and CMIP6 standard historical emissions are available between 1750 and 2014. Both Finn and Qfed fire emissions are available between xx and 2020.

Meteorological data have to be also regridded in both horizontal and vertical resolution of the model simulations. For this we provide MERRA2 reanalysis between 1975 and 2020 on the original horizontal and vertical grid  (https://rda.ucar.edu/datasets/ds313.3/) or local on cheyenne, as well as scripts that interpolate those to the vertical and horizontal resolution of the model.

Finally, atmospheric initial condition files have to be regridded to the required model resolution. We provide an ncl script that allows us to regrid from the standard spectral element 1 degree resolution to the new grid. 

1.Emissions
1.1. Fire (Biomass burning)
1.1.1. FINN 
A tool for regridding FINN v1.5 and v2.0 (Fire INventory from NCAR) \*.txt files to FV, SE and SE-RR and MPAS are available. The tool includes a README file, which describes in detail how to use the code, and several example namelists files (\*.inp), which can be easily adapted as needed.
Produced emission data can be found here:
/glade/p/acom/acom-nsc/musicav0/emissions/FINN\_2020.
- [ ] source data location here
 
1.1.2. Qfed (Biomass burning)
Qfed biomass burning emissions have been produced on a regular fv 1 degree resolution. A tool is available to regrid these emissions to spectral element resolution. Higher resolution Qfed emissions will be provided later.
- [ ] source data location here

1.2. Anthropogenic
1.2.1. EPA
- [ ] source data location here
1.2.2. CAMS
- [ ] source data location here
Interpolate original CAMS emissions to the desired grid: regrid\_fv2se\_cams\_anthro.ncl
Rename and calculated regridded CAMS emissions to be able to run in CESM2: rename\_cams\_anthro\_se.ncl

1.3. CMIP6
Interpolate original CMIP6 emissions to the desired grid: regrid\_fv2se\_cmip6\_main.ncl
Rename and calculate regridded to be able to run in CESM2: 
- [ ] source data location here

2. Meteorological Data Regridding
2.1 MERRA2 original data script to regrid to the desired resolution
- [ ] source data location here

3. Initial Conditions Regridding
/glade/u/home/tilmes/ncl/SE/regrid\_all\_spectral\_data.ncl
To run with CESM, this has to be converted to a cdf5 format:
nccopy -k cdf5 oldfile newfile
- [ ] source data location here

