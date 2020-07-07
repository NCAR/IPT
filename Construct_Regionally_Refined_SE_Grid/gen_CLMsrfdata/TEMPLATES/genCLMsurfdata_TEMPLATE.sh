#!/bin/bash
#
# Colin Zarzycki 9/17/2015
#
# This script will generate the CLM fsurdat and transient land use files for 1850/2000 runs
# 
# NOTES:
# User needs to have write access to $CESMROOT (files are generated within this substructure
# User needs to have write access to $VRrepoPath
# Output files are written to ${VRrepoPath}/${VRgridName}/clm_surfdata_${CLMVERSION}

##=======================================================================
#PBS -N sub_genfsurdat_AAAAAAAAAA
#PBS -A Pxxxxxxxx                     <-- SET PROJECT  HERE
#PBS -l walltime=3:59:00
#PBS -q premium
#PBS -j oe
#PBS -l select=4:ncpus=2:mpiprocs=2:mem=109GB
################################################################

set +e

module load mpt

CESMROOT="TTTTTTTTTT"
VRrepoPath="RRRRRRRRRR"
VRgridName="GGGGGGGGGG"
VRscripFile="${VRrepoPath}/${VRgridName}/grids/SSSSSSSSSS"
MYcdate="DDDDDDDDDD"

VRshort=${VRgridName}
## TMPDIRBASE="/glade/scratch/USERNAME/"
TMPDIRBASE=${VRrepoPath}/${VRgridName}

#FIX  ESMFBIN_PATH="/glade/u/apps/ch/opt/esmf/7.0.0-ncdfio-mpi/intel/17.0.1/bin/binO/Linux.intel.64.mpi.default"
ESMFBIN_PATH="/glade/u/apps/ch/opt/esmf/7.1.0r-ncdfio/intel/17.0.1/bin/binO/Linux.intel.64.mpi.default"
CLMVERSION="5_0" # options are 4_0 or 5_0
DO_SP_ONLY=true   # true (only create SP surdats) or false (create full crop surdats)

#----------------------------------------------------------------------
# First, we need to generate the mapping files
# This may next up to 12 hours to finish with very refined grids
# NOTE: May need to change {mpitype} in ESMF binary path to mpiuni (mkmapdata.sh)
# probably since we are spawning in serial.
#----------------------------------------------------------------------
cdate=`date +%y%m%d` # Get data in YYMMDD format

# Create TMPDIR
### TMPDIR=${TMPDIRBASE}/tmp.clmsurfdata.${cdate}/
TMPDIR=${TMPDIRBASE}/maps_clm/
mkdir -p ${TMPDIR}

# Use for CESM2.0xx
MKMAPDATADIR=${CESMROOT}/components/clm/tools/mkmapdata/

cd ${TMPDIR}
regrid_num_proc=8
time env ESMFBIN_PATH=${ESMFBIN_PATH} REGRID_PROC=$regrid_num_proc ${MKMAPDATADIR}/mkmapdata.sh -b -v --gridfile ${VRscripFile} --res ${VRgridName} --gridtype global

cd ${CESMROOT}/components/clm/tools/mksurfdata_map/

if ($DO_SP_ONLY); then
  CROPSTRING="-no-crop"
else
  CROPSTRING=""
fi
./mksurfdata.pl -years 1850-2000,1850,2000 ${CROPSTRING} -res usrspec -usr_gname ${VRgridName} -usr_gdate ${cdate} -usr_mapdir ${TMPDIR}

## Move the surface datasets
mkdir -p ${VRrepoPath}/${VRgridName}/clm_surfdata_${CLMVERSION}
mv landuse*${VRgridName}*nc surfdata_${VRgridName}_*.nc ${VRrepoPath}/${VRgridName}/clm_surfdata_${CLMVERSION}

cd ${VRrepoPath}/${VRgridName}/clm_surfdata_${CLMVERSION}
mv landuse.timeseries_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr1850-2015_c${cdate}.nc  landuse.timeseries_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr1850-2015_c${MYcdate}.nc
mv surfdata_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr1850_c${cdate}.nc surfdata_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr1850_c${MYcdate}.nc 
mv surfdata_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr2000_c${cdate}.nc surfdata_${VRgridName}_hist_16pfts_Irrig_CMIP6_simyr2000_c${MYcdate}.nc 

# Delete mapping files. Comment this lie out if you need the map_clm files for some reason.
rm -rf ${TMPDIR}

