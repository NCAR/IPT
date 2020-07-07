#! /bin/csh -f

##--------------------------------------------------------
##  set CESM_TAG    = "/path/to/cesm/tag/"
##  set REPO_PATH   = "/path/to/repo/"  
##  set GRID_TAG    = "NEWGRIDNAME_ne30x4"
##  set GRID_NAME   = "ne0np4.NEWGRIDNAME.ne30x4"
##  set DATE        = "YYMMDD"
##
##  set SCRIP_FILE  = $GRID_TAG"_np4_SCRIP.nc"
##  set LATLON_FILE = $GRID_TAG"_np4_LATLON.nc"
##  set EXODUS_FILE = $GRID_TAG"_np4_EXODUS.nc"
##--------------------------------------------------------

  set CESM_TAG    = $VRM_CESM_TAG
  set REPO_PATH   = $VRM_REPO_PATH
  set GRID_TAG    = $VRM_GRID_TAG
  set GRID_NAME   = $VRM_GRID_NAME
  set DATE        = $VRM_DATE

  set SCRIP_FILE  = $GRID_TAG"_np4_SCRIP.nc"
  set LATLON_FILE = $GRID_TAG"_np4_LATLON.nc"
  set EXODUS_FILE = $GRID_TAG"_np4_EXODUS.nc"

  set SEDcmdA="s:AAAAAAAAAA:"$GRID_TAG":g"
  set SEDcmdT="s:TTTTTTTTTT:"$CESM_TAG":g"
  set SEDcmdR="s:RRRRRRRRRR:"$REPO_PATH":g"
  set SEDcmdG="s:GGGGGGGGGG:"$GRID_NAME":g"
  set SEDcmdS="s:SSSSSSSSSS:"$SCRIP_FILE":g"
  set SEDcmdL="s:LLLLLLLLLL:"$LATLON_FILE":g"
  set SEDcmdE="s:EEEEEEEEEE:"$EXODUS_FILE":g"
  set SEDcmdD="s:DDDDDDDDDD:"$DATE":g"


  echo " "
  echo "  #--------------------------------------------------------- "
  echo "  # Creating gen_CAMncdata scripts in directory: "  $GRID_NAME
  echo "  # "
  echo "  #  Review the variable settings at the begining of the "
  echo "  #  script and then run it:"
  echo "  # "
  echo "  #   cheyenne% module load ncl"
  echo "  #   cheyenne% module load nco"
  echo "  #   cheyenne% qcmd -- 'sh interpic_script_"$GRID_TAG".sh >& LOG_"$GRID_TAG"'"
  echo "  # "
  echo "  #--------------------------------------------------------- "
  echo " "
  mkdir -p $GRID_NAME
  sed -e $SEDcmdT -e $SEDcmdG -e $SEDcmdR -e $SEDcmdL -e $SEDcmdD < ./TEMPLATES/interpic_script_TEMPLATE.sh >  ./$GRID_NAME/interpic_script_$GRID_TAG.sh 

