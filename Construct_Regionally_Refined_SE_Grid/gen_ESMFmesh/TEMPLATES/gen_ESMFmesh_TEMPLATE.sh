#!/bin/bash
#

#
#  Load modules for EMSF build for 1 CPU (non-MPI) using netcdf libraries
#
#  If these current values are out of date, run:
#
#    casper% module spider esmf-VERSION-ncdfio-uni-O 
#
#  for a list of the modules that need to be loaded.
#
#-----------------------------------------------------------------------------
module purge
module load intel/17.0.1
module load esmflibs/7.1.0r
module load esmf-7.1.0r-ncdfio-uni-O


# Set the VR grid name, the location for related files
#------------------------------------------------------
VRrepoPath="RRRRRRRRRR"
VRgridName="GGGGGGGGGG"
VRscripFile="SSSSSSSSSS"
VRmeshFile="MMMMMMMMMM"
VRgridPath="${VRrepoPath}/${VRgridName}/grids/"


ESMF_Scrip2Unstruct  ${VRgridPath}/${VRscripFile}  ${VRgridPath}/${VRmeshFile} 0

