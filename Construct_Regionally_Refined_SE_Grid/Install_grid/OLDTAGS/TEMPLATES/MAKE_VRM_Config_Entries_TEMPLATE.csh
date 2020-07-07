#! /bin/csh -f

##--------------------------------------------------------
##  set CESM_TAG    = "/path/to/cesm/tag/"
##  set REPO_PATH   = "/path/to/repo/"  
##  set GRID_TAG    = "NEWGRIDNAME_ne30x4"
##  set GRID_NAME   = "ne0np4.NEWGRIDNAME.ne30x4"
##  set DATE        = "200303"
##
##  set SCRIP_FILE  = $GRID_TAG"_np4_SCRIP.nc"
##  set LATLON_FILE = $GRID_TAG"_np4_LATLON.nc"
##  set EXODUS_FILE = $GRID_TAG"_np4_EXODUS.nc"
##--------------------------------------------------------

  set GRID_TAG   = "AXAXAXAXAX"
  set GRIDNAME   = "GXGXGXGXGX"
  set GRIDDESC   = "ADD A DESCRIPTION OF THE GRID HERE"

  set REPO_PATH  = "RXRXRXRXRX"
  set DATE       = "DXDXDXDXDX"
  set NCOL       = "NXNXNXNXNX"

  set EXODUS_FILE = $GRID_TAG"_np4_EXODUS.nc"
  set RESNAME_tx  = $GRIDNAME'_mt12'
  set RESNAME_gx  = $GRIDNAME'_mg17'
  set FULLPATH    = $REPO_PATH'/'$GRIDNAME'/'
  set OUTFILE     = 'VRMODS_config_grids.xml'

  # model_grid Entry:
  #---------------------
  echo ' '                                         >  $OUTFILE
  echo ' ******************************'           >> $OUTFILE
  echo ' ADD VarMesh: model_grid entry:'           >> $OUTFILE
  echo ' ******************************'           >> $OUTFILE
  echo ' '                                         >> $OUTFILE
  echo '<model_grid alias="'$RESNAME_gx'">'        >> $OUTFILE
  echo '   <grid name="atm">'$GRIDNAME'</grid>'    >> $OUTFILE
  echo '   <grid name="lnd">'$GRIDNAME'</grid>'    >> $OUTFILE
  echo '   <grid name="ocnice">'$GRIDNAME'</grid>' >> $OUTFILE
  echo '   <mask>gx1v7</mask>'                     >> $OUTFILE
  echo '</model_grid>'                             >> $OUTFILE
  echo ' '                                         >> $OUTFILE
  echo '<model_grid alias="'$RESNAME_tx'">'        >> $OUTFILE
  echo '   <grid name="atm">'$GRIDNAME'</grid>'    >> $OUTFILE
  echo '   <grid name="lnd">'$GRIDNAME'</grid>'    >> $OUTFILE
  echo '   <grid name="ocnice">'$GRIDNAME'</grid>' >> $OUTFILE
  echo '   <mask>tx0.1v2</mask>'                   >> $OUTFILE
  echo '</model_grid>'                             >> $OUTFILE
  echo ' '                                         >> $OUTFILE
  echo ' '                                         >> $OUTFILE

  # domain Entry:
  #---------------------
  echo ' '                                         >> $OUTFILE
  echo ' **************************'               >> $OUTFILE
  echo ' ADD VarMesh: domain entry:'               >> $OUTFILE
  echo ' **************************'               >> $OUTFILE
  echo ' '                                         >> $OUTFILE
  echo '<domain name="'$GRIDNAME'">'                                                                                >> $OUTFILE
  echo '   <nx>'$NCOL'</nx> <ny>1</ny>'                                                                             >> $OUTFILE
  echo '   <file grid="atm|lnd" mask="gx1v7"  >'$FULLPATH'domains/domain.lnd.'$GRIDNAME'_gx1v7.'$DATE'.nc</file>'   >> $OUTFILE
  echo '   <file grid="ocnice"  mask="gx1v7"  >'$FULLPATH'domains/domain.ocn.'$GRIDNAME'_gx1v7.'$DATE'.nc</file>'   >> $OUTFILE
  echo '   <file grid="atm|lnd" mask="tx0.1v2">'$FULLPATH'domains/domain.lnd.'$GRIDNAME'_tx0.1v2.'$DATE'.nc</file>' >> $OUTFILE
  echo '   <file grid="ocnice"  mask="tx0.1v2">'$FULLPATH'domains/domain.ocn.'$GRIDNAME'_tx0.1v2.'$DATE'.nc</file>' >> $OUTFILE
  echo '   <desc>'$GRIDNAME' '$GRIDDESC'</desc>'                                                                    >> $OUTFILE
  echo '   <support>Test support only</support>'                                                                    >> $OUTFILE
  echo '</domain>'                                                                                                  >> $OUTFILE
  echo ' '                                                                                                          >> $OUTFILE
  echo ' '                                                                                                          >> $OUTFILE

  # gridmap Entries:
  #---------------------
  echo ' '                                         >> $OUTFILE
  echo ' *****************************'            >> $OUTFILE
  echo ' ADD VarMesh: gridmap entries:'            >> $OUTFILE
  echo ' *****************************'            >> $OUTFILE
  echo ' '                                         >> $OUTFILE
  echo '<gridmap atm_grid="'$GRIDNAME'" ocn_grid="gx1v7">'                                                 >> $OUTFILE
  echo '   <map name="ATM2OCN_FMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_gx1v7_aave.'$DATE'.nc</map>'    >> $OUTFILE
  echo '   <map name="ATM2OCN_SMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_gx1v7_blin.'$DATE'.nc</map>'    >> $OUTFILE
  echo '   <map name="ATM2OCN_VMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_gx1v7_blin.'$DATE'.nc</map>'    >> $OUTFILE
  echo '   <map name="OCN2ATM_FMAPNAME">'$FULLPATH'maps/map_gx1v7_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>'    >> $OUTFILE
  echo '   <map name="OCN2ATM_SMAPNAME">'$FULLPATH'maps/map_gx1v7_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>'    >> $OUTFILE
  echo '</gridmap>'                                                                                        >> $OUTFILE
  echo '<gridmap atm_grid="'$GRIDNAME'" ocn_grid="tx0.1v2">'                                               >> $OUTFILE
  echo '   <map name="ATM2OCN_FMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_tx0.1v2_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '   <map name="ATM2OCN_SMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_tx0.1v2_blin.'$DATE'.nc</map>'  >> $OUTFILE
  echo '   <map name="ATM2OCN_VMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_tx0.1v2_blin.'$DATE'.nc</map>'  >> $OUTFILE
  echo '   <map name="OCN2ATM_FMAPNAME">'$FULLPATH'maps/map_tx0.1v2_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '   <map name="OCN2ATM_SMAPNAME">'$FULLPATH'maps/map_tx0.1v2_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '</gridmap>'                                                                                        >> $OUTFILE
  echo '<gridmap lnd_grid="'$GRIDNAME'" glc_grid="gland4" >'                                               >> $OUTFILE
  echo '   <map name="LND2GLC_FMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_gland4km_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '   <map name="LND2GLC_SMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_gland4km_blin.'$DATE'.nc</map>' >> $OUTFILE
  echo '   <map name="GLC2LND_FMAPNAME">'$FULLPATH'maps/map_gland4km_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '   <map name="GLC2LND_SMAPNAME">'$FULLPATH'maps/map_gland4km_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '</gridmap>'                                                                                        >> $OUTFILE
  echo '<gridmap lnd_grid="'$GRIDNAME'" rof_grid="r05">'                                                   >> $OUTFILE
  echo '   <map name="LND2ROF_FMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_r05_aave.'$DATE'.nc</map>'      >> $OUTFILE
  echo '   <map name="ROF2LND_FMAPNAME">'$FULLPATH'maps/map_r05_TO_'$GRIDNAME'_aave.'$DATE'.nc</map>'      >> $OUTFILE
  echo '</gridmap>'                                                                                        >> $OUTFILE
  echo '<gridmap atm_grid="'$GRIDNAME'" wav_grid="ww3a">'                                                  >> $OUTFILE
  echo '   <map name="ATM2WAV_SMAPNAME">'$FULLPATH'maps/map_'$GRIDNAME'_TO_ww3a_blin.'$DATE'.nc</map>'     >> $OUTFILE
  echo '</gridmap>'                                                                                        >> $OUTFILE
  echo ' '                                                                                                 >> $OUTFILE
  echo ' '                                                                                                 >> $OUTFILE


  set OUTFILE  = 'VRMODS_config_component_cesm.xml'

  # ADD ATM_NCPL value for VRM grid:
  #---------------------------------
  echo ' '                                                        >  $OUTFILE
  echo ' ******************************'                          >> $OUTFILE
  echo ' ADD VarMesh: ATM_NCPL entry:'                            >> $OUTFILE
  echo ' ******************************'                          >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo ' Add an ATM_NCPL entry for the VR grid:'                  >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo '<entry id="ATM_NCPL">'                                    >> $OUTFILE
  echo '  <type>integer</type>'                                   >> $OUTFILE
  echo '  <default_value>48</default_value>'                      >> $OUTFILE
  echo '  <values match="last">'                                  >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo ' Add one of the following:'                               >> $OUTFILE
  echo '  To set a 900 sec timestep:'                             >> $OUTFILE
  echo '    <value compset=".+" grid="a%'$GRIDNAME'">96</value>'  >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo '  To set a 600 sec timestep:'                             >> $OUTFILE
  echo '    <value compset=".+" grid="a%'$GRIDNAME'">144</value>' >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo '  To set a 450 sec timestep:'                             >> $OUTFILE
  echo '    <value compset=".+" grid="a%'$GRIDNAME'">192</value>' >> $OUTFILE
  echo ' '                                                        >> $OUTFILE
  echo ' '                                                        >> $OUTFILE

  set OUTFILE  = 'VRMODS_horiz_grid.xml'

  # ADD horiz_grid Entry:
  #---------------------------------
  echo ' '                                                                                      >  $OUTFILE
  echo ' ******************************'                                                        >> $OUTFILE
  echo ' ADD VarMesh: horiz_grid entry:'                                                        >> $OUTFILE
  echo ' ******************************'                                                        >> $OUTFILE
  echo ' '                                                                                      >> $OUTFILE
  echo '  <horiz_grid dyn="se" hgrid="'$GRIDNAME'" ncol="'$NCOL'" csne="0" csnp="4" npg="0" />' >> $OUTFILE
  echo ' '                                                                                      >> $OUTFILE
  echo ' '                                                                                      >> $OUTFILE

  set OUTFILE  = 'VRMODS_namelist_defaults_cam.xml'

  # ADD namelist default values for CAM:
  #--------------------------------------
  echo ' '                                                                                                                             >  $OUTFILE
  echo ' ********************************************'                                                                                 >> $OUTFILE
  echo ' ADD VarMesh: namelist_defaults for VRM grid:'                                                                                 >> $OUTFILE
  echo ' ********************************************'                                                                                 >> $OUTFILE
  echo ' '                                                                                                                             >> $OUTFILE
  echo ' These are the default namelist settings for the new VR grid:'                                                                 >> $OUTFILE
  echo '     NOTE:  <dtime>      must be consistent with the ATM_NCPL value.'                                                          >> $OUTFILE
  echo '            <se_fine_ne> should be set to the hi resolution value for the VR grid.'                                            >> $OUTFILE
  echo ' '                                                                                                                             >> $OUTFILE
  echo ' '                                                                                                                             >> $OUTFILE
  echo '    <dtime dyn="se"    hgrid="'$GRIDNAME'">600</dtime>'                                                                        >> $OUTFILE
  echo '    <ncdata dyn="se" hgrid="'$GRIDNAME'"  nlev="32"  ic_ymd="101" >'$FULLPATH'inic/cami-mam4_0000-01-01_'$GRIDNAME'_L32.nc</ncdata>' >> $OUTFILE
  echo '    <bnd_topo hgrid="'$GRIDNAME'"  >'$FULLPATH'topo/topo_'$GRIDNAME'_blin_'$DATE'.nc</bnd_topo>'                               >> $OUTFILE
  echo '    <se_mesh_file hgrid="'$GRIDNAME'">'$FULLPATH'grids/'$EXODUS_FILE'</se_mesh_file>'                                          >> $OUTFILE
  echo '    <drydep_srf_file hgrid="'$GRIDNAME'">'$FULLPATH'atmsrf/atmsrf_'$GRIDNAME'_'$DATE'.nc</drydep_srf_file>'                    >> $OUTFILE
  echo '    <se_fine_ne hgrid="'$GRIDNAME'" > 240 </se_fine_ne>'                                                                       >> $OUTFILE
  echo '    <se_hypervis_subcycle  hgrid="'$GRIDNAME'"> 10 </se_hypervis_subcycle>'                                                    >> $OUTFILE
  echo '    <se_nsplit hgrid="'$GRIDNAME'"> 6 </se_nsplit>'                                                                            >> $OUTFILE
  echo '    <se_nu_div se_refined_mesh="1" hypervis_type="scalar" hgrid="'$GRIDNAME'" > 1.0000e13 </se_nu_div>'                        >> $OUTFILE
  echo '    <se_nu_p se_refined_mesh="1" hypervis_type="scalar" hgrid="'$GRIDNAME'" > 1.0000e13 </se_nu_p>'                            >> $OUTFILE
  echo '    <se_nu_top hgrid="'$GRIDNAME'"> 1.0e5 </se_nu_top>'                                                                        >> $OUTFILE
  echo '    <se_refined_mesh hgrid="'$GRIDNAME'" > .true. </se_refined_mesh>'                                                          >> $OUTFILE
  echo '    <se_rsplit hgrid="'$GRIDNAME'"> 4 </se_rsplit>'                                                                            >> $OUTFILE
  echo '    <focndomain hgrid="'$GRIDNAME'" mask="gx1v7"   >'$FULLPATH'domains/domain.ocn.'$GRIDNAME'_gx1v7.'$DATE'.nc</focndomain>'   >> $OUTFILE
  echo '    <focndomain hgrid="'$GRIDNAME'" mask="tx0.1v2" >'$FULLPATH'domains/domain.ocn.'$GRIDNAME'_tx0.1v2.'$DATE'.nc</focndomain>' >> $OUTFILE
  echo ' '                                                                                                                             >> $OUTFILE
  echo ' '                                                                                                                             >> $OUTFILE

  set OUTFILE  = 'VRMODS_namelist_defaults_clm4_5.xml'

  # ADD namelist default values for CLM:
  #--------------------------------------
  echo ' '                                                                                                                                 >  $OUTFILE
  echo ' ********************************************'                                                                                     >> $OUTFILE
  echo ' ADD VarMesh: namelist_defaults for VRM grid:'                                                                                     >> $OUTFILE
  echo ' ********************************************'                                                                                     >> $OUTFILE
  echo ' '                                                                                                                                 >> $OUTFILE
  echo ' <fsurdat hgrid="'$GRIDNAME'"    sim_year="2000" use_crop=".false." irrigate=".true.">'                                            >> $OUTFILE
  echo $FULLPATH'clm_surfdata_5_0/surfdata_'$GRIDNAME'_hist_16pfts_Irrig_CMIP6_simyr2000_c'$DATE'.nc</fsurdat>'                            >> $OUTFILE
  echo ' '                                                                                                                                 >> $OUTFILE
  echo ' <fsurdat hgrid="'$GRIDNAME'"    sim_year="1850" use_crop=".false." irrigate=".true.">'                                            >> $OUTFILE
  echo $FULLPATH'clm_surfdata_5_0/surfdata_'$GRIDNAME'_hist_16pfts_Irrig_CMIP6_simyr1850_c'$DATE'.nc</fsurdat>'                            >> $OUTFILE
  echo ' '                                                                                                                                 >> $OUTFILE
  echo ' <flanduse_timeseries hgrid="'$GRIDNAME'"      sim_year_range="1850-2000"  irrigate=".true." use_crop=".false."  >'                >> $OUTFILE
  echo $FULLPATH'clm_surfdata_5_0/landuse.timeseries_'$GRIDNAME'_hist_16pfts_Irrig_CMIP6_simyr1850-2015_c'$DATE'.nc</flanduse_timeseries>' >> $OUTFILE
  echo ' '                                                                                                                                 >> $OUTFILE
  echo ' '                                                                                                                                 >> $OUTFILE

