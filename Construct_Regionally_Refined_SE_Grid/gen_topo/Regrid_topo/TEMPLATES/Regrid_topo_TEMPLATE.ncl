load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ; These need to be modified by the user
  ;--------------------------------------------
  VRgridName   = "GGGGGGGGGG"
  VRgridPath   = "RRRRRRRRRR"+VRgridName 
  VRscripFile  =  VRgridPath+"/grids/"+"SSSSSSSSSS"
  VRlatlonFile =  VRgridPath+"/grids/"+"LLLLLLLLLL"
  cdate        = "DDDDDDDDDD"

  ;============================================
  ; default settings only modify if necessary
  ;============================================
  InterpMethod="bilinear"   ;bilinear, patch, conserve
  infile="/glade/p/cesmdata/inputdata/atm/cam/topo/ne120np4_nc3000_Nsw010_Nrs002_Co015_Fi001_20180926.nc"
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
;;;  delete(Opt)


  ;==================================
  ; Read in values from source file
  ;==================================
  fin = addfile(infile,"r")
    PHIS_in         = fin->PHIS
    LANDFRAC_in     = fin->LANDFRAC
    SGH_in          = fin->SGH
    SGH30_in        = fin->SGH30
    LANDM_COSLAT_in = fin->LANDM_COSLAT
    SGH_UF_in       = fin->SGH_UF
    TERR_UF_in      = fin->TERR_UF
    GBXAR_in        = fin->GBXAR
    MXDIS_in        = fin->MXDIS
    RISEQ_in        = fin->RISEQ
    FALLQ_in        = fin->FALLQ
    MXVRX_in        = fin->MXVRX
    MXVRY_in        = fin->MXVRY
    ANGLL_in        = fin->ANGLL
    ANGLX_in        = fin->ANGLX
    ANISO_in        = fin->ANISO
    ANIXY_in        = fin->ANIXY
    HWDTH_in        = fin->HWDTH
    WGHTS_in        = fin->WGHTS
    CWGHT_in        = fin->CWGHT
    CLNGT_in        = fin->CLNGT
    COUNT_in        = fin->COUNT
  delete(fin)

  ; Get lat,lon,area from VRlatlonFile
  ;-----------------------------------
  cin             = addfile(VRlatlonFile,"r")
  lon             = cin->lon
  lat             = cin->lat
  area            = cin->area
  element_corners = cin->element_corners
  delete(cin)

  ;===============================
  ; Interpolate Values 
  ;===============================
;;;  Opt = True

     ; PHIS
     ;----------
     PHIS = ESMF_regrid_with_weights(PHIS_in,wgtFileName,Opt)
     copy_VarAtts(PHIS_in,PHIS)
     PHIS!0 = "ncol"
     if(any(ismissing(PHIS))) then
       print("PHIS contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(PHIS,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; LANDFRAC
     ;----------
     LANDFRAC = ESMF_regrid_with_weights(LANDFRAC_in,wgtFileName,Opt)
     copy_VarAtts(LANDFRAC_in,LANDFRAC)
     LANDFRAC!0 = "ncol"
     if(any(ismissing(LANDFRAC))) then
       print("LANDFRAC contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(LANDFRAC,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; SGH
     ;----------
     SGH = ESMF_regrid_with_weights(SGH_in,wgtFileName,Opt)
     copy_VarAtts(SGH_in,SGH)
     SGH!0 = "ncol"
     if(any(ismissing(SGH))) then
       print("SGH contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(SGH,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; SGH30
     ;----------
     SGH30 = ESMF_regrid_with_weights(SGH30_in,wgtFileName,Opt)
     copy_VarAtts(SGH30_in,SGH30)
     SGH30!0 = "ncol"
     if(any(ismissing(SGH30))) then
       print("SGH30 contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(SGH30,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; LANDM_COSLAT
     ;----------
     LANDM_COSLAT = ESMF_regrid_with_weights(LANDM_COSLAT_in,wgtFileName,Opt)
     copy_VarAtts(LANDM_COSLAT_in,LANDM_COSLAT)
     LANDM_COSLAT!0 = "ncol"
     if(any(ismissing(LANDM_COSLAT))) then
       print("LANDM_COSLAT contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(LANDM_COSLAT,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; SGH_UF
     ;----------
     SGH_UF = ESMF_regrid_with_weights(SGH_UF_in,wgtFileName,Opt)
     copy_VarAtts(SGH_UF_in,SGH_UF)
     SGH_UF!0 = "ncol"
     if(any(ismissing(SGH_UF))) then
       print("SGH_UF contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(SGH_UF,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; TERR_UF
     ;----------
     TERR_UF = ESMF_regrid_with_weights(TERR_UF_in,wgtFileName,Opt)
     copy_VarAtts(TERR_UF_in,TERR_UF)
     TERR_UF!0 = "ncol"
     if(any(ismissing(TERR_UF))) then
       print("TERR_UF contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(TERR_UF,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; GBXAR
     ;----------
     GBXAR = ESMF_regrid_with_weights(GBXAR_in,wgtFileName,Opt)
     copy_VarAtts(GBXAR_in,GBXAR)
     GBXAR!0 = "ncol"
     if(any(ismissing(GBXAR))) then
       print("GBXAR contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(GBXAR,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; MXDIS
     ;----------
     MXDIS = ESMF_regrid_with_weights(MXDIS_in,wgtFileName,Opt)
     copy_VarAtts(MXDIS_in,MXDIS)
     MXDIS!0 = "nrdg"
     MXDIS!1 = "ncol"
     if(any(ismissing(MXDIS))) then
       print("MXDIS contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(MXDIS,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; RISEQ
     ;----------
     RISEQ = ESMF_regrid_with_weights(RISEQ_in,wgtFileName,Opt)
     copy_VarAtts(RISEQ_in,RISEQ)
     RISEQ!0 = "nrdg"
     RISEQ!1 = "ncol"
     if(any(ismissing(RISEQ))) then
       print("RISEQ contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(RISEQ,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; FALLQ
     ;----------
     FALLQ = ESMF_regrid_with_weights(FALLQ_in,wgtFileName,Opt)
     copy_VarAtts(FALLQ_in,FALLQ)
     FALLQ!0 = "nrdg"
     FALLQ!1 = "ncol"
     if(any(ismissing(FALLQ))) then
       print("FALLQ contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(FALLQ,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; MXVRX
     ;----------
     MXVRX = ESMF_regrid_with_weights(MXVRX_in,wgtFileName,Opt)
     copy_VarAtts(MXVRX_in,MXVRX)
     MXVRX!0 = "nrdg"
     MXVRX!1 = "ncol"
     if(any(ismissing(MXVRX))) then
       print("MXVRX contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(MXVRX,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; MXVRY
     ;----------
     MXVRY = ESMF_regrid_with_weights(MXVRY_in,wgtFileName,Opt)
     copy_VarAtts(MXVRY_in,MXVRY)
     MXVRY!0 = "nrdg"
     MXVRY!1 = "ncol"
     if(any(ismissing(MXVRY))) then
       print("MXVRY contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(MXVRY,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; ANGLL
     ;----------
     ANGLL = ESMF_regrid_with_weights(ANGLL_in,wgtFileName,Opt)
     copy_VarAtts(ANGLL_in,ANGLL)
     ANGLL!0 = "nrdg"
     ANGLL!1 = "ncol"
     if(any(ismissing(ANGLL))) then
       print("ANGLL contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(ANGLL,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; ANGLX
     ;----------
     ANGLX = ESMF_regrid_with_weights(ANGLX_in,wgtFileName,Opt)
     copy_VarAtts(ANGLX_in,ANGLX)
     ANGLX!0 = "nrdg"
     ANGLX!1 = "ncol"
     if(any(ismissing(ANGLX))) then
       print("ANGLX contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(ANGLX,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; ANISO
     ;----------
     ANISO = ESMF_regrid_with_weights(ANISO_in,wgtFileName,Opt)
     copy_VarAtts(ANISO_in,ANISO)
     ANISO!0 = "nrdg"
     ANISO!1 = "ncol"
     if(any(ismissing(ANISO))) then
       print("ANISO contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(ANISO,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; ANIXY
     ;----------
     ANIXY = ESMF_regrid_with_weights(ANIXY_in,wgtFileName,Opt)
     copy_VarAtts(ANIXY_in,ANIXY)
     ANIXY!0 = "nrdg"
     ANIXY!1 = "ncol"
     if(any(ismissing(ANIXY))) then
       print("ANIXY contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(ANIXY,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; HWDTH
     ;----------
     HWDTH = ESMF_regrid_with_weights(HWDTH_in,wgtFileName,Opt)
     copy_VarAtts(HWDTH_in,HWDTH)
     HWDTH!0 = "nrdg"
     HWDTH!1 = "ncol"
     if(any(ismissing(HWDTH))) then
       print("HWDTH contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(HWDTH,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; WGHTS
     ;----------
     WGHTS = ESMF_regrid_with_weights(WGHTS_in,wgtFileName,Opt)
     copy_VarAtts(WGHTS_in,WGHTS)
     WGHTS!0 = "nrdg"
     WGHTS!1 = "ncol"
     if(any(ismissing(WGHTS))) then
       print("WGHTS contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(WGHTS,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; CWGHT
     ;----------
     CWGHT = ESMF_regrid_with_weights(CWGHT_in,wgtFileName,Opt)
     copy_VarAtts(CWGHT_in,CWGHT)
     CWGHT!0 = "nrdg"
     CWGHT!1 = "ncol"
     if(any(ismissing(CWGHT))) then
       print("CWGHT contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(CWGHT,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; CLNGT
     ;----------
     CLNGT = ESMF_regrid_with_weights(CLNGT_in,wgtFileName,Opt)
     copy_VarAtts(CLNGT_in,CLNGT)
     CLNGT!0 = "nrdg"
     CLNGT!1 = "ncol"
     if(any(ismissing(CLNGT))) then
       print("CLNGT contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(CLNGT,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))

     ; COUNT
     ;----------
     COUNT = ESMF_regrid_with_weights(COUNT_in,wgtFileName,Opt)
     copy_VarAtts(COUNT_in,COUNT)
     COUNT!0 = "nrdg"
     COUNT!1 = "ncol"
     if(any(ismissing(COUNT))) then
       print("COUNT contains some missing values. Exiting.")
       exit
     end if
     delete_VarAtts(COUNT,(/"lat1d","lon1d","remap","_FillValue","missing_value"/))


  ;===============================
  ; Create the output file
  ;===============================
  TOPO_File   = "topo_"+VRgridName+"_"+interpString+"_"+cdate+".nc"
  system("mkdir -p "  +VRgridPath+"/topo/"          )  ; ensure that the path exists
  system("/bin/rm -f "+VRgridPath+"/topo/"+TOPO_File)  ; remove any pre-existing file
  ncdf = addfile(VRgridPath+"/topo/"+TOPO_File ,"c")   ; open output netCDF file
  
  fAtt               = True            ; assign file attributes
  fAtt@source_file   = srcGridDir+"/"+srcGridFile
  fAtt@regrid_file   = wgtFileName
  fAtt@creation_date = systemfunc ("date")
  fileattdef( ncdf, fAtt )            ; copy file attributes
  
  ncdf->lon             = lon
  ncdf->lat             = lat
  ncdf->area            = area
  ncdf->element_corners = element_corners
  ncdf->PHIS            = PHIS
  ncdf->LANDFRAC        = LANDFRAC
  ncdf->SGH             = SGH
  ncdf->SGH30           = SGH30
  ncdf->LANDM_COSLAT    = LANDM_COSLAT
  ncdf->SGH_UF          = SGH_UF
  ncdf->TERR_UF         = TERR_UF
  ncdf->GBXAR           = GBXAR
  ncdf->MXDIS           = MXDIS
  ncdf->RISEQ           = RISEQ
  ncdf->FALLQ           = FALLQ
  ncdf->MXVRX           = MXVRX
  ncdf->MXVRY           = MXVRY
  ncdf->ANGLL           = ANGLL
  ncdf->ANGLX           = ANGLX
  ncdf->ANISO           = ANISO
  ncdf->ANIXY           = ANIXY
  ncdf->HWDTH           = HWDTH
  ncdf->WGHTS           = WGHTS
  ncdf->CWGHT           = CWGHT
  ncdf->CLNGT           = CLNGT
  ncdf->COUNT           = COUNT
  
  ; cleanup
  ;----------
  print("cleanin up...")
  ;system("rm "+wgtFileName)

end
