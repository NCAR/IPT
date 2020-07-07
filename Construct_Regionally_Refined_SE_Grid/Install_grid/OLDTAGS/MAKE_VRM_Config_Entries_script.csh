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

  set REPO_PATH   = $VRM_REPO_PATH
  set GRID_TAG    = $VRM_GRID_TAG
  set GRID_NAME   = $VRM_GRID_NAME
  set DATE        = $VRM_DATE
  set NCOL        = $VRM_GRID_NCOL

  set SEDcmdA="s:AXAXAXAXAX:"$GRID_TAG":g"
  set SEDcmdG="s:GXGXGXGXGX:"$GRID_NAME":g"
  set SEDcmdR="s:RXRXRXRXRX:"$REPO_PATH":g"
  set SEDcmdD="s:DXDXDXDXDX:"$DATE":g"
  set SEDcmdN="s:NXNXNXNXNX:"$NCOL":g"

  echo " "
  echo "  #--------------------------------------------------------- "
  echo "  # Creating MAKE_VRM_Config_Entries script in directory: "  $GRID_NAME
  echo "  # "
  echo "  #  Review the variable settings at the begining of the "
  echo "  #  script and ADD a grid description for the variable GRIDDESC"
  echo "  # "
  echo "  #  Then run the script to generate XML file entries:"
  echo "  # "
  echo "  #   cheyenne% chmod 744 MAKE_VRM_Config_Entries_"$GRID_TAG".csh"
  echo "  #   cheyenne% ./MAKE_VRM_Config_Entries_"$GRID_TAG".csh"
  echo "  # "
  echo "  #  Then follow the instructions to 'mouse' the entries into the XML files"
  echo "  # "
  echo "  #--------------------------------------------------------- "
  echo " "
  mkdir -p $GRID_NAME
  sed -e $SEDcmdA -e $SEDcmdG -e $SEDcmdR -e $SEDcmdD -e $SEDcmdN < ./TEMPLATES/MAKE_VRM_Config_Entries_TEMPLATE.csh >  ./$GRID_NAME/MAKE_VRM_Config_Entries_$GRID_TAG.csh 

