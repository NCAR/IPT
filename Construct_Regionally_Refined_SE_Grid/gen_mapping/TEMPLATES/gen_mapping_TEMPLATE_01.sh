################################################################
#PBS -N ncl-maps
#PBS -A Pxxxxxxxx 
#PBS -l walltime=01:00:00
#PBS -q regular
#PBS -k oe
#PBS -m a 
#PBS -M user@ucar.edu
#PBS -l select=1:ncpus=36
################################################################

#----------------------------------------------------------------------
#  PURPOSE: This script optionally generates mapping weights between 
#           various grids atm/ocn/lnd/rof/glc/wav that users may need. 
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#  Set the path where VRM files are stored, the name of the VRMgrid 
#  as used by CESM, and the Name of the SCRIP file for the grid.
#
#  It is expected that the EXODUS/SCRIP/LATLON grid files generated 
#  by the VRM_Editor and Gen_Control_Volumes programs are stored 
#  in the directory:
#                    ${VRrepoPath}/${VRgridName}/grids/
#
#  This Script will create the directory containing mapping weights
#  between atm/ocn/lnd/rof/glc/wav grids:
#
#                    ${VRrepoPath}/${VRgridName}/maps/
#--------------------------------------------------------------------
CESMroot="TTTTTTTTTT"
VRrepoPath="RRRRRRRRRR"
VRgridName="GGGGGGGGGG"
VRscripName="SSSSSSSSSS"
cdate="DDDDDDDDDD"

##mapDir="${CESMroot}/cime/tools/mapping/gen_mapping_files"
mapDir="${PWD}/../mapping/gen_mapping_files"

#--------------------------------------------------------------------
# For Typical runs the ocn/lnd/atm are on the same grid. 
# The specified grids are either VRM scrip files or grids 
# specified from the CESM inputdata repository. When 
# the same grid is specified for 2 components, map creation 
# is skipped.
#--------------------------------------------------------------------
atmName=${VRgridName}
atmGridName="${VRrepoPath}/${atmName}/grids/${VRscripName}"
lndName=${VRgridName}
lndGridName="${VRrepoPath}/${lndName}/grids/${VRscripName}"
### lndName="1x1"
### lndGridName="/glade/p/cesm/cseg/mapping/grids/1x1d.nc"
ocnName=${VRgridName}
ocnGridName="${VRrepoPath}/${ocnName}/grids/${VRscripName}"
### ocnName="gx1v7"
### ocnGridName="/glade/p/cesmdata/cseg/mapping/grids/gx1v7_151008.nc"
rofName="r05"
rofGridName="/glade/p/cesmdata/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_0.5x0.5_nomask_c110308.nc"
glcName="gland4km"
glcGridName="/glade/p/cesm/cseg/inputdata/glc/cism/griddata/SCRIPgrid_greenland_4km_epsg3413_c170414.nc"
wavName="ww3a"
wavGridName="/glade/p/cesm/cseg/mapping/grids/ww3a_120222.nc"

#wgtFileDir="${VRrepoPath}/${VRgridName}/maps.${cdate}/"
wgtFileDir="${VRrepoPath}/${VRgridName}/maps/"
mkdir -p $wgtFileDir

echo " "
echo "==============================================================="
echo " Create mapping weight files for atm/ocn/lnd/rof/glc/wav grids:"
echo " "
echo " Atmosphere:"
echo " -----------"
echo " Grid Name= "${atmName}
echo " Grid File= "${atmGridName}
echo " "
echo " Land:"
echo " -----------"
echo " Grid Name= "${lndName}
echo " Grid File= "${lndGridName}
echo " "
echo " Ocean:"
echo " -----------"
echo " Grid Name= "${ocnName}
echo " Grid File= "${ocnGridName}
echo " "
echo " ROF Grid:"
echo " -----------"
echo " Grid Name= "${rofName}
echo " Grid File= "${rofGridName}
echo " "
echo " GLC Grid:"
echo " -----------"
echo " Grid Name= "${glcName}
echo " Grid File= "${glcGridName}
echo " "
echo " WAV Grid:"
echo " -----------"
echo " Grid Name= "${wavName}
echo " Grid File= "${wavGridName}
echo " "
echo " Mapping files will be created in the directory:"
echo " ----------------------------------------------- "
echo ${wgtFileDir}
echo " "
echo "==============================================================="
echo " "


#-------------------------------------
# It's all automatic from here on out.
#-------------------------------------
cd ${wgtFileDir}

if [ "$atmName" != "$lndName" ]; then
  echo "Generating ATM <-> [LND,OCN,ROF,GLC] maps..... "

  ${mapDir}/gen_cesm_maps.sh -fatm ${atmGridName} -natm ${atmName} -flnd ${lndGridName} -nlnd ${lndName} -focn ${ocnGridName} -nocn ${ocnName} -frtm ${rofGridName} -nrtm ${rofName} -fglc ${glcGridName} -nglc ${glcName} -idate ${cdate}

else

  echo "Generating ATM <-> [OCN,ROF,GLC] maps..... "
  ${mapDir}/gen_cesm_maps.sh -fatm ${atmGridName} -natm ${atmName} -focn ${ocnGridName} -nocn ${ocnName} -frtm ${rofGridName} -nrtm ${rofName} -fglc ${glcGridName} -nglc ${glcName} -idate ${cdate}

fi

############################# WAV <-> ATM ########################################

if [ "$wavName" != "$atmName" ]; then
  echo "Generating ATM <-> WAV maps..... "

  interp_method="blin"
  ${mapDir}/gen_ESMF_mapping_file/create_ESMF_map.sh -fsrc ${atmGridName} -nsrc ${atmName} -fdst ${wavGridName} -ndst ${wavName} -map ${interp_method} -idate ${cdate}

  interp_method="blin"
  ${mapDir}/gen_ESMF_mapping_file/create_ESMF_map.sh -fsrc ${wavGridName} -nsrc ${wavName} -fdst ${atmGridName} -ndst ${atmName} -map ${interp_method} -idate ${cdate}

fi

