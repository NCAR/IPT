load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ; These need to be modified by the user
  ;--------------------------------------------
  cdate        = "DDDDDDDDDD"
  VRgridName   = "GGGGGGGGGG"
  VRgridPath   = "RRRRRRRRRR"+VRgridName
  VRscripFile  =  VRgridPath+"/grids/"+"SSSSSSSSSS"
  VRlatlonFile =  VRgridPath+"/grids/"+"LLLLLLLLLL"
  ATMSRF_File  = "atmsrf_"+VRgridName+"_"+cdate+".nc"

  ;============================================
  ; default settings only modify if necessary
  ;============================================
  InterpMethod="patch"   ;bilinear, patch, conserve
  infile="/glade/p/cesmdata/inputdata/atm/cam/chem/trop_mam/atmsrf_ne120np4_181018.nc"
  srcName="ne120np4"
  srcGridDir="/glade/p/cesmdata/cseg/mapping/grids/"
  srcGridFile="ne120np4_pentagons_100310.nc"
  
  if (InterpMethod .eq. "patch") then
    interpString="patc"
  end if
  if (InterpMethod .eq. "bilinear") then
    interpString="blin"
  end if
  if (InterpMethod .eq. "conserve") then
    interpString="aave"
  end if


  ;=======================================
  ; Now generate weights file
  ;=======================================
  srcGridName = srcGridDir+"/"+srcGridFile
  dstGridName = VRscripFile
  wgtFileName = VRgridPath+"/maps/map_"+srcName+"_TO_"+VRgridName+"_"+interpString+"."+cdate+".nc"
  
  Opt                = True
  Opt@InterpMethod   = InterpMethod
  Opt@ForceOverwrite = True
  Opt@PrintTimings   = True
  ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)
  delete(Opt)


  ;==================================
  ; Read in values from source file
  ;==================================
  fin = addfile(infile,"r")
  fraction_landuse_in = fin->fraction_landuse
  soilw_in = fin->soilw
  
  Opt = True
  fraction_landuse_out = ESMF_regrid_with_weights(fraction_landuse_in,wgtFileName,Opt)
  soilw_out            = ESMF_regrid_with_weights(soilw_in           ,wgtFileName,Opt)
  
  if(any(ismissing(soilw_out))) then
    print("soilw_out contains some missing values. Exiting.")
    exit
  end if
    if(any(ismissing(fraction_landuse_out))) then
    print("fraction_landuse_out contains some missing values. Exiting.")
    exit
  end if
  
  soilw_out!1            = "ncol"
  fraction_landuse_out!1 = "ncol"
  
  delete_VarAtts(soilw_out           ,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))
  delete_VarAtts(fraction_landuse_out,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

  ; Get lat,lon,area from VRlatlonFile
  ;-----------------------------------
  cin             = addfile(VRlatlonFile,"r")
  lon             = cin->lon
  lat             = cin->lat
  area            = cin->area
  element_corners = cin->element_corners
  delete(cin)


  ;===============================
  ; Create the output file
  ;===============================
  system("mkdir -p "  +VRgridPath+"/atmsrf/"            )  ; ensure that the path exists
  system("/bin/rm -f "+VRgridPath+"/atmsrf/"+ATMSRF_File)  ; remove any pre-existing file
  ncdf = addfile(VRgridPath+"/atmsrf/"+ATMSRF_File ,"c")   ; open output netCDF file
  
  fAtt               = True            ; assign file attributes
  fAtt@source_file   = srcGridDir+"/"+srcGridFile
  fAtt@regrid_file   = wgtFileName
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes

  ncdf->lon              = lon
  ncdf->lat              = lat
  ncdf->area             = area
  ncdf->element_corners  = element_corners
  ncdf->fraction_landuse = fraction_landuse_out
  ncdf->soilw            = soilw_out
  
  ; cleanup
  ;----------
  print("cleanin up...")
  ;system("rm "+wgtFileName)

end
