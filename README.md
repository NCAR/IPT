# EMIT
Emissions Model Interface Tool

A collection of scripts and code to convert emissions to CAM-SE regionally refined grids

## Four Components
1. Regridding of point sources, biomass burning emissions and NEI emissions 
2. Combining NEI emissions with global emissions
3. Regridding CMIP6 emissions using ncl and ESMF tools
4. Regridding CAMS anthropogenic emissions using ncl and ESMF tools

## Regridding of point sources, biomass burning emissions and NEI emissions

###  checkout
git clone https://github.com/NCAR/EMIT
### where to get netcdf files for testing (And where to place them)
* source data location here *
### Compile and run tests
```
cd regrid/Dev_FIRE_EMIS/src
gmake Makefile
cd ../tst
./run.exe
```


