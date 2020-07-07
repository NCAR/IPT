#!/bin/bash

# CESMROOT = top-level directory from CESM source tree
# LEVELFILE = CAM ncdata file containing the vertical coordinate you want to use (any dycore)
# INTERPFILE = CAM ncdata file to interpolate from (FV preferred)
#
# VRlatlonFile = SE "latlon" file containing the ncol points on the VR grid
# VRinicFile = CAM ncdata file for use with VR grid
#
# Notes
# A.) This script requires interpic_new to be build. See below.
#======================================================================================

# Set the tag location for the interpic program, a level file to provide 
# vertical level information, and a startup/inic source file 
#------------------------------------------------------------------------
CESMROOT="TTTTTTTTTT"
INTERPIC="${CESMROOT}/components/cam/tools/interpic_new/interpic"
LEVELFILE="/glade/p/cesmdata/inputdata/atm/cam/inic/cam_vcoords_L32_c180105.nc"
INTERPFILE="/glade/p/cesmdata/inputdata/atm/cam/inic/fv/cami-mam4_0000-01-01_0.9x1.25_L32_c150403.nc"

# Set the VR grid name, the location for related files
#------------------------------------------------------
VRdate="DDDDDDDDDD"
VRgridName="GGGGGGGGGG"
VRrepoPath="RRRRRRRRRR"
VRgridPath="${VRrepoPath}/${VRgridName}"
VRinicPath="${VRrepoPath}/${VRgridName}/inic"

# Set the LATLON file for the grid, and the output filename
#-----------------------------------------------------------
VRlatlonFile="${VRgridPath}/grids/LLLLLLLLLL"
VRinicFile="${VRinicPath}/cami-mam4_0000-01-01_${VRgridName}_L32_c${VRdate}.nc"


# Check the status of interpic
#----------------------------------------
if [ ! -f ${INTERPIC} ]; then
  echo "interpic binary doesn't exist. Either re/build interpic or correct path"
  echo "BINARY: ${INTERPIC}"
  echo "... aborting"
  # To build...
  # $USER> cd ${CESMROOT}/components/cam/tools/interpic_new/
  # $USER> make
  exit
fi

# Ensure that the destination path exists
#-----------------------------------------
mkdir -p ${VRinicPath}

# Strip lat/lon arrays from SE template file
#---------------------------------------------
ncks -v lat,lon ${VRlatlonFile} SE_template.nc

# Attach vertical coordinates to template file
#-----------------------------------------------
ncks -A -v hyai,hyam,hybi,hybm ${LEVELFILE} SE_template.nc

# Run interpic to interpolate to SE grid.
#------------------------------------------
${INTERPIC} -t SE_template.nc ${INTERPFILE} ${VRinicFile}

# Rename US and VS (assuming FV inic) to U,V for SE
#---------------------------------------------------
ncrename -v US,U -v VS,V ${VRinicFile}

# Cleanup unneeded files
#-------------------------
rm SE_template.nc
