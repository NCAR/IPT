#! /bin/csh -f

  # EDIT THESE VALUES FOR THE NEW GRID
  #----------------------------------------
##--------------------------------------------------------
##  set REPO_PATH   = "/path/to/repo/" 
##  set GRID_TAG    = "NEWGRIDNAME_ne30x4"
##  set GRID_NAME   = "ne0np4.NEWGRIDNAME.ne30x4"
##  set GRID_DESC   = "Description for the VRM file BLAH, BLAH, BLAH....."
##  set NCOL        = 'NNNNNN'
##  set DATE        = "YYMMDD"
##
##  set EXODUS_FILE = $GRID_TAG"_EXODUS.nc"
##  set FULLPATH    = $REPO_PATH'/'$GRID_NAME'/'
##  set OUTFILE     = $FULLPATH"config_grids.xml"
##--------------------------------------------------------

  set GRID_TAG   = "AXAXAXAXAX"
  set GRID_NAME  = "GXGXGXGXGX"
  set GRIDDESC   = "ADD A DESCRIPTION OF THE GRID HERE"

  set REPO_PATH  = "RXRXRXRXRX"
  set DATE       = "DXDXDXDXDX"
  set NCOL       = "NXNXNXNXNX"

  #-------------------------------------------------------------------------------
  # NO NEED TO EDIT PAST HERE
  #-------------------------------------------------------------------------------
  set EXODUS_FILE = $GRID_TAG"_EXODUS.nc"
  set MESH_FILE   = $GRID_TAG"_MESH.nc"
  set FULLPATH    = $REPO_PATH'/'$GRID_NAME'/'
  set OUTFILE     = $FULLPATH"config_grids.xml"

  cp ./CONFIG-TEMPLATES/shell_commands-TEMPLATE $FULLPATH/shell_commands


  set SEDcmdG="s:GGGGGGGGGG:"$GRID_NAME":g"
  set SEDcmdP="s:PPPPPPPPPP:"$FULLPATH":g"
  set SEDcmdE="s:EEEEEEEEEE:"$EXODUS_FILE":g"
  set SEDcmdD="s:DDDDDDDDDD:"$DATE":g"
  sed -e $SEDcmdG -e $SEDcmdP -e $SEDcmdE -e $SEDcmdD < ./CONFIG-TEMPLATES/user_nl_cam-TEMPLATE  >  $FULLPATH/user_nl_cam
  sed -e $SEDcmdG -e $SEDcmdP -e $SEDcmdD             < ./CONFIG-TEMPLATES/user_nl_clm-TEMPLATE  >  $FULLPATH/user_nl_clm



  cat ./CONFIG-TEMPLATES/config_grids.xml-TEMPLATE_PT1                            >  $OUTFILE
  echo '   <!--- VR-CESM grids with CAM-SE -->'                                   >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <model_grid alias="'$GRID_NAME'_g17"> '                               >> $OUTFILE
  echo '      <grid name="atm">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="lnd">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="ocnice">gx1v7</grid>'                                   >> $OUTFILE
  echo '      <mask>gx1v7</mask>'                                                 >> $OUTFILE
  echo '    </model_grid>'                                                        >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <model_grid alias="'$GRID_NAME'_mg17" not_compset="_POP">'            >> $OUTFILE
  echo '      <grid name="atm">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="lnd">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="ocnice">'$GRID_NAME'</grid>'                            >> $OUTFILE
  echo '      <mask>gx1v7</mask>'                                                 >> $OUTFILE
  echo '    </model_grid>'                                                        >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <model_grid alias="'$GRID_NAME'_t12">'                                >> $OUTFILE
  echo '      <grid name="atm">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="lnd">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="ocnice">tx0.1v2</grid>'                                 >> $OUTFILE
  echo '      <mask>tx0.1v2</mask>'                                               >> $OUTFILE
  echo '    </model_grid>'                                                        >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <model_grid alias="'$GRID_NAME'_mt12" not_compset="_POP">'            >> $OUTFILE
  echo '      <grid name="atm">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="lnd">'$GRID_NAME'</grid>'                               >> $OUTFILE
  echo '      <grid name="ocnice">'$GRID_NAME'</grid>'                            >> $OUTFILE
  echo '      <mask>tx0.1v2</mask>'                                               >> $OUTFILE
  echo '    </model_grid>'                                                        >> $OUTFILE
  cat ./CONFIG-TEMPLATES/config_grids.xml-TEMPLATE_PT2                            >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo '    <!--- ATM/LND domains global -->'                                     >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <domain name="'$GRID_NAME'">'                                         >> $OUTFILE
  echo '      <nx>'$NCOL'</nx> <ny>1</ny>'                                        >> $OUTFILE
  echo '      <file grid="atm|lnd" mask="gx1v7"  >'$FULLPATH'domains/domain.lnd.'$GRID_NAME'_gx1v7.'$DATE'.nc</file>'   >> $OUTFILE
  echo '      <file grid="ocnice"  mask="gx1v7"  >'$FULLPATH'domains/domain.ocn.'$GRID_NAME'_gx1v7.'$DATE'.nc</file>'   >> $OUTFILE
  echo '      <file grid="atm|lnd" mask="tx0.1v2">'$FULLPATH'domains/domain.lnd.'$GRID_NAME'_tx0.1v2.'$DATE'.nc</file>' >> $OUTFILE
  echo '      <file grid="ocnice"  mask="tx0.1v2">'$FULLPATH'domains/domain.ocn.'$GRID_NAME'_tx0.1v2.'$DATE'.nc</file>' >> $OUTFILE
  echo '      <mesh driver="nuopc">'$FULLPATH'grids/'$MESH_FILE'</mesh>' >> $OUTFILE
  echo '      <desc>'$GRIDDESC'</desc>'                                           >> $OUTFILE
  echo '      <support>Test support only</support>'                               >> $OUTFILE
  echo '    </domain>'                                                            >> $OUTFILE
  cat ./CONFIG-TEMPLATES/config_grids.xml-TEMPLATE_PT3                            >> $OUTFILE
  echo '    <!--- VRM grid mappings -->'                                          >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo '    <!--- gridS: lnd to glc and glc to lnd mapping                 -->'   >> $OUTFILE
  echo '    <!---                                                          -->'   >> $OUTFILE
  echo '    <!--- Note that we use aave maps for GLC2LND_SMAPNAME, because -->'   >> $OUTFILE
  echo '    <!--- glc is typically much higher resolution than lnd, which  -->'   >> $OUTFILE
  echo '    <!--- means that area-conservative remapping is preferable.    -->'   >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- GLC: gland4           -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <gridmap lnd_grid="'$GRID_NAME'" glc_grid="gland4" >'                  >> $OUTFILE
  echo '      <map name="LND2GLC_FMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_gland4km_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="LND2GLC_SMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_gland4km_blin.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="GLC2LND_FMAPNAME">'$FULLPATH'maps/map_gland4km_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="GLC2LND_SMAPNAME">'$FULLPATH'maps/map_gland4km_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '    </gridmap>'                                                           >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- GLC: gland5UM         -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- GLC: gland20          -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo '    <!---- river to land and land to river mapping files           -->'   >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- ROF: r05              -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <gridmap lnd_grid="'$GRID_NAME'" rof_grid="r05">'                     >> $OUTFILE
  echo '      <map name="LND2ROF_FMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_r05_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="ROF2LND_FMAPNAME">'$FULLPATH'maps/map_r05_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '    </gridmap>'                                                           >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- ROF: r01              -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo '    <!--- ROF: JRA025           -->'                                      >> $OUTFILE
  echo '    <!--- ********************  -->'                                      >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '<!---  ********************************************************************************** -->' >> $OUTFILE
  echo '<!---  ************    config_grids_mct.xml   ******************************************* -->' >> $OUTFILE
  echo '<!---  ********************************************************************************** -->' >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo '    <!---- atm to ocean and ocean to atm mapping files             -->'   >> $OUTFILE
  echo '    <!--- ======================================================== -->'   >> $OUTFILE
  echo ''                                                                         >> $OUTFILE
  echo '    <gridmap atm_grid="'$GRID_NAME'" ocn_grid="gx1v7">'                   >> $OUTFILE
  echo '      <map name="ATM2OCN_FMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_gx1v7_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '      <map name="ATM2OCN_SMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_gx1v7_blin.'$DATE'.nc</map>'  >> $OUTFILE
  echo '      <map name="ATM2OCN_VMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_gx1v7_blin.'$DATE'.nc</map>'  >> $OUTFILE
  echo '      <map name="OCN2ATM_FMAPNAME">'$FULLPATH'maps/map_gx1v7_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '      <map name="OCN2ATM_SMAPNAME">'$FULLPATH'maps/map_gx1v7_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>'  >> $OUTFILE
  echo '    </gridmap>'                                                                                      >> $OUTFILE
  echo '    <gridmap atm_grid="'$GRID_NAME'" ocn_grid="tx0.1v2">'                                            >> $OUTFILE
  echo '      <map name="ATM2OCN_FMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_tx0.1v2_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="ATM2OCN_SMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_tx0.1v2_blin.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="ATM2OCN_VMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_tx0.1v2_blin.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="OCN2ATM_FMAPNAME">'$FULLPATH'maps/map_tx0.1v2_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '      <map name="OCN2ATM_SMAPNAME">'$FULLPATH'maps/map_tx0.1v2_TO_'$GRID_NAME'_aave.'$DATE'.nc</map>' >> $OUTFILE
  echo '    </gridmap>'                                                                                      >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo '    <!--- ======================================================== -->'                              >> $OUTFILE
  echo '    <!---- atm to wav mapping files                                -->'                              >> $OUTFILE
  echo '    <!--- ======================================================== -->'                              >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo '    <gridmap atm_grid="'$GRID_NAME'" wav_grid="ww3a">'                                               >> $OUTFILE
  echo '      <map name="ATM2WAV_SMAPNAME">'$FULLPATH'maps/map_'$GRID_NAME'_TO_ww3a_blin.200131.nc</map>'    >> $OUTFILE
  echo '    </gridmap>'                                                                                      >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo '  </gridmaps>'                                                                                       >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo ''                                                                                                    >> $OUTFILE
  echo '</grid_data>'                                                                                        >> $OUTFILE
