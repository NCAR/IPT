#!/bin/bash
#
# Colin Zarzycki
# Created: 9/15/2015
# Last updated: 12/25/2018
#
# Batch script to generate domain files for VR grids.
# VRrepoPath: the low-level directory where files will be stored
# VRscripfile: Absolute path to SCRIP descriptor for VR grid
# VRname: Long VR name
# VRshort: Short VR name
#
# NOTES
# Folder will be generated at $VRrepoPath/$VRname for storage
# Resultant output will be (1 each) *.lnd.* and *.ocn.* domain files
# Use tx0.1v2 mask unless using an exceptionally low-res grid (per Mariana)
#
# If you are getting errors with gen_domain binary, make sure you follow the README
# instructions for building, as the tool now configures mach specific
#

##=======================================================================
#PBS -N sub_gendomain_AAAAAAAAAA
#PBS -A Pxxxxxxxx 
#PBS -l walltime=00:19:00
#PBS -q premium
#PBS -j oe
#PBS -M user@ucar.edu
#PBS -l select=1:ncpus=36
################################################################


# Set the path to the tag being used
#-------------------------------------
CESMroot="TTTTTTTTTT"

# Set the path for the VRM files, the name of the grid as used 
# by CESM, and the path to the VRM SCRIP file. 
#---------------------------------------------------------------
VRrepoPath="RRRRRRRRRR"
VRgridName="GGGGGGGGGG"
VRscripFile="${VRrepoPath}/${VRgridName}/grids/SSSSSSSSSS"
cdate="DDDDDDDDDD"

# Path to CIME mapping tools
# Note: write access to this directory is necessary, so either 
# copy the exec or checkout and build in your own dir
#-----------------------------------------------------------
PATH_TO_MAPPING="${CESMroot}/cime/tools/mapping/"

# Set OCN grid to use.
# May need to change these settings, but safe to use t12 for mask 
# for anything >ne30
#----------------------------------------------------------------
#ocnName="tx0.1v2"
#ocnGridName="/glade/p/cesmdata/cseg/mapping/grids/tx0.1v2_090127.nc"
ocnName="gx1v7"
ocnGridName="/glade/p/cesmdata/cseg/mapping/grids/gx1v7_151008.nc"

# Set the map file needed
#-------------------------------------------
wgtFileDir="${VRrepoPath}/${VRgridName}/maps"
CurDate=`date +%y%m%d`
aaveMap="${wgtFileDir}/map_${ocnName}_TO_${VRgridName}_aave.${cdate}.nc"
echo ${aaveMap}


#----------------------------------------------------------------------
# CREATE MAPPING FILE
#----------------------------------------------------------------------
cd $PATH_TO_MAPPING/gen_domain_files

# do ATM2OCN_FMAPNAME (aave)
#-----------------------------
#ESMFBIN_PATH="/glade/u/apps/ch/opt/esmf/7.0.0-ncdfio-mpi/intel/17.0.1/bin/binO/Linux.intel.64.mpi.default"
#ESMFBIN_PATH="/glade/u/apps/ch/opt/esmf/7.0.0-ncdfio/intel/17.0.1/bin/bing/Linux.intel.64.mpiuni.default"
interp_method="conserve"   # bilinear, patch, conserve
ESMF_RegridWeightGen --ignore_unmapped -m ${interp_method} -w ${aaveMap} -s ${ocnGridName} -d ${VRscripFile}
##${ESMFBIN_PATH}/ESMF_RegridWeightGen --ignore_unmapped -m ${interp_method} -w ${aaveMap} -s ${ocnGridName} -d ${VRscripFile}

#----------------------------------------------------------------------
# CREATE DOMAIN FILES
#----------------------------------------------------------------------

set +e
./gen_domain -m ${aaveMap} -o gx1v7 -l ${VRgridName}

#----------------------------------------------------------------------
# MOVING FILES + CLEANUP 
#----------------------------------------------------------------------
# Move domain files to VRrepoPath dir
mkdir -p ${VRrepoPath}/${VRgridName}/domains
## mv domain*${VRgridName}*${cdate}*nc ${VRrepoPath}/${VRgridName}/domains
mv domain.lnd.${VRgridName}_${ocnName}.${CurDate}.nc ${VRrepoPath}/${VRgridName}/domains/domain.lnd.${VRgridName}_${ocnName}.${cdate}.nc
mv domain.ocn.${VRgridName}_${ocnName}.${CurDate}.nc ${VRrepoPath}/${VRgridName}/domains/domain.ocn.${VRgridName}_${ocnName}.${cdate}.nc



# Remove mapping files since they are large and we really only needed aave for domains anyway
###rm map_*.nc


