# IPT
Input Processing Tool for MUSICA version 0

A collection of scripts and code to convert emissions to CAM-SE regionally refined grids

### Four Components
1. Regridding of point sources, biomass burning emissions and NEI emissions 
2. Combining NEI emissions with global emissions
3. Regridding CMIP6 emissions using ncl and ESMF tools
4. Regridding CAMS anthropogenic emissions using ncl and ESMF tools

The scientific description can be found in the [wiki](https://github.com/NCAR/EMIT/wiki)

# 1. Regrid point sources, biomass burning emissions and NEI emissions

## Where to get netcdf files for testing (And where to place them)
- [ ] source data location here
## Compile and run tests
```
cd regrid_code/Dev_FIRE_EMIS/src
gmake Makefile
cd ../tst
./run.exe
```

# 2. Combining NEI emissions with global emissions
- [ ] Should something be here?

# 3. Regrid CMIP6 emissions to different grids using ncl and ESMF

### 3.1. Regrid CMIP6 original resolution emissions to desired resolution
Script directory: regrid/FV2SE \
Batch script (for cheyenne) for conversion is in regrid\_main \
Script to specify grids: Regrid\_fv2se\_cmip6\_main.ncl

### 3.2 Emission name remapping
Script directory: regrid/emissions\_hist \
Batch submission script (slurm) script\_rename \
Script to map names: rename\_and\_convert\_cmip\_hist.ncl 

# Regridded emission data
* CSLAM: /glade/p/acom/acom-climate/cmip6inputs/historical\_ne30pg3
* ne30\_ne30: /glade/p/acom/acom-climate/cmip6inputs/historical\_ne30
* FV3: /glade/p/acom/acom-climate/cmip6inputs/historical\_C96

# 4. Regrid CAMS anthropogenic emissions 
Script directory: regrid/CAMS \
Source data /glade/p/acom/acom-climate/tilmes/emis/download \
Regridding the origin emissions to the desired SE mesh: Regrid\_fv2se\_cams\_anthro.ncl \
Regridding CAMS anthropogenic:  emissions\_cams/regrid\_cams\_anthro\_se.ncl \
