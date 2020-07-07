#include "ReadWrite.h"


//=====================================================================
void Write_Exodus_File(std::string     ExodusName,
                       VRM_Param &     VRM,
                       RefinementMap & OutputRefineMap,
                       CubeGrid &      OutputGrid)
{
    const int ParamFour      = 4;
    const int ParamLenString = 33;

    // Output to a NetCDF Exodus file
    //--------------------------------
    NcFile ncOut(ExodusName.c_str(), NcFile::FileMode::Replace);

    // Random Exodus dimensions
    //--------------------------
    NcDim *dimLenString = ncOut.add_dim("len_string", ParamLenString);
    NcDim *dimLenLine   = ncOut.add_dim("len_line"  , 81            );
    NcDim *dimFour      = ncOut.add_dim("four"      , ParamFour     );
    NcDim *dimTime      = ncOut.add_dim("time_step"                 );
    NcDim *dimDimension = ncOut.add_dim("num_dim"   , 3             );

    // Number of nodes
    //-------------------
    int nNodeCount =  OutputGrid.Nodes.size();
    NcDim *dimNodes = ncOut.add_dim("num_nodes", nNodeCount);

    // Number of elements
    //---------------------
    int nElementCount = OutputGrid.Faces.size();


    NcDim *dimElements         = ncOut.add_dim("num_elem"       , nElementCount);
    NcDim *dimNumElementBlocks = ncOut.add_dim("num_el_blk"     , 1            );
    NcDim *dimNumQARec         = ncOut.add_dim("num_qa_rec"     , 1            );
    NcDim *dimElementBlock1    = ncOut.add_dim("num_el_in_blk1" , nElementCount);
    NcDim *dimNodesPerElement  = ncOut.add_dim("num_nod_per_el1", 4            );
    NcDim *dimAttBlock1        = ncOut.add_dim("num_att_in_blk1", 1            );

    // Global attributes
    //--------------------
    ncOut.add_att("api_version"             , 4.98f);
    ncOut.add_att("version"                 , 4.98f);
    ncOut.add_att("floating_point_word_size", 8    );
    ncOut.add_att("file_size"               , 0    );

    char szTitle[512];
    sprintf(szTitle, "VRM_Editor(%s) ", ExodusName.c_str());
    ncOut.add_att("title", szTitle);

    // Time_whole (unused)
    //---------------------
    ncOut.add_var("time_whole", ncDouble, dimTime);

    // QA records
    //-------------
    char szQARecord[ParamFour][ParamLenString] = {"CUBIT","13.0","01/01/2013","00:00:00"};

    NcVar *varQARecords = ncOut.add_var("qa_records",ncChar,dimNumQARec,dimFour,dimLenString);
    varQARecords->set_cur(0,0,0);
    varQARecords->put(&(szQARecord[0][0]), 1, 4, ParamLenString);

    // Coordinate names
    //---------------------
    char szCoordNames[3][ParamLenString] = {"x", "y", "z"};

    NcVar *varCoordNames = ncOut.add_var("coor_names", ncChar, dimDimension, dimLenString);
    varCoordNames->set_cur(0,0,0);
    varCoordNames->put(&(szCoordNames[0][0]), 3, ParamLenString);

    // Element block names
    //-----------------------
    NcVar *varElementBlockNames = ncOut.add_var("eb_names",ncChar,dimNumElementBlocks,dimLenString);

    // Element map
    //---------------
    int *nElementMap = new int[nElementCount];
    for (int i = 0; i < nElementCount; i++) { nElementMap[i] = i+1; }

    NcVar * varElementMap = ncOut.add_var("elem_map", ncInt, dimElements);
    varElementMap->put(nElementMap, nElementCount);
    delete[] nElementMap;

    // Element block status
    //-----------------------
    int nOne = 1;

    NcVar *varElementBlockStatus = ncOut.add_var("eb_status", ncInt, dimNumElementBlocks);
    varElementBlockStatus->put(&nOne, 1);

    NcVar *varElementProperty = ncOut.add_var("eb_prop1", ncInt, dimNumElementBlocks);
    varElementProperty->put    (&nOne, 1);
    varElementProperty->add_att("name", "ID");

    // Attributes
    //--------------
    double *dAttrib1 = new double[nElementCount];
    for (int i = 0; i < nElementCount; i++) { dAttrib1[i] = 1.0; }

    NcVar *varAttrib1 = ncOut.add_var("attrib1", ncDouble, dimElementBlock1, dimAttBlock1);
    varAttrib1->put(dAttrib1, nElementCount, 1);
    delete[] dAttrib1;

    // Face nodes (1-indexed)
    //-------------------------
    NcVar *varFaces = ncOut.add_var("connect1", ncInt, dimElementBlock1, dimNodesPerElement);

    varFaces->add_att("elem_type", "SHELL4");

    int nNodesPerElement = 4;
    int *nConnect        = new int[nNodesPerElement];
    for (int i = 0; i < nElementCount   ; i++) {
      for (int k = 0; k < nNodesPerElement; k++) {
        nConnect[k] = OutputGrid.Faces[i][k] + 1;
      }
      varFaces->set_cur(i,0);
      varFaces->put(nConnect, 1, nNodesPerElement);
    }
    delete[] nConnect;

    // Node list
    //-------------
    NcVar  *varNodes = ncOut.add_var("coord", ncDouble, dimDimension, dimNodes);
    double *dCoord   = new double[nNodeCount];

    for (int i = 0; i < nNodeCount; i++) { dCoord[i] = OutputGrid.Nodes[i].x; }
    varNodes->set_cur(0,0);
    varNodes->put(dCoord, 1, nNodeCount);

    for (int i = 0; i < nNodeCount; i++) { dCoord[i] = OutputGrid.Nodes[i].y; }
    varNodes->set_cur(1,0);
    varNodes->put(dCoord, 1, nNodeCount);

    for (int i = 0; i < nNodeCount; i++) { dCoord[i] = OutputGrid.Nodes[i].z; }
    varNodes->set_cur(2,0);
    varNodes->put(dCoord, 1, nNodeCount);

    delete[] dCoord;

    // Add dimensions/variables for refinement map
    //--------------------------------------------
    int NLON = OutputRefineMap.nRefLon;
    int NLAT = OutputRefineMap.nRefLat;
    double lat[NLAT];
    double lon[NLON];
    double rval[NLAT][NLON];
    NcDim *dimLat = ncOut.add_dim("ref_lat",NLAT);
    NcDim *dimLon = ncOut.add_dim("ref_lon",NLON);
    NcVar *varLat = ncOut.add_var("ref_lat",ncDouble,dimLat);
    NcVar *varLon = ncOut.add_var("ref_lon",ncDouble,dimLon);
    NcVar *varMap = ncOut.add_var("refineMap",ncDouble,dimLat,dimLon);
    varMap->add_att("long_name","refinement map");
    varMap->add_att("description","Lat,Lon map of [0.=low res,1.=hi res] refinement values");

    for(int ii=0; ii<NLON; ii++) {
        lon[ii] = OutputRefineMap.RefLon[ii][0];
    }
    for(int ii=0; ii<NLAT; ii++) {
        lat[ii] = OutputRefineMap.RefLat[ii][0];
    }
    for(int jj=0; jj<NLAT; jj++) {
    for(int ii=0; ii<NLON; ii++) {
        rval[jj][ii] = OutputRefineMap.val[ii][jj];
    }
    }

    varLat->put(lat,NLAT);
    varLon->put(lon,NLON);
    varMap->put(&rval[0][0],NLAT,NLON);

    // Add VRM_Editor processing parameters
    //   --when it comes to string data, netcdf is *really* stupid....
    //------------------------------------------------------------------
    int   MyStrLen  = 12;
    char *OutString = new char[MyStrLen];
    char *OutBlank  = new char[MyStrLen];
    strcpy(OutBlank,"           ");
    NcDim *dimString             = ncOut.add_dim("StrLen",MyStrLen);
    NcVar *VRMresolution         = ncOut.add_var("VRM_Resolution"        ,ncInt   );
    NcVar *VRMrefinement         = ncOut.add_var("VRM_Refinement"        ,ncInt   );
    NcVar *VRMtessellations      = ncOut.add_var("VRM_Tessellations"     ,ncInt   );
    NcVar *VRMsubcellresolution  = ncOut.add_var("VRM_SubCellResolution" ,ncInt   );
    NcVar *VRMtranssmoothdist    = ncOut.add_var("VRM_TransSmoothDist"   ,ncInt   );
    NcVar *VRMsmoothiterations   = ncOut.add_var("VRM_SmoothIterations"  ,ncInt   );
    NcVar *VRMgridxrotate        = ncOut.add_var("VRM_GridXRotate"       ,ncDouble);
    NcVar *VRMgridyrotate        = ncOut.add_var("VRM_GridYRotate"       ,ncDouble);
    NcVar *VRMgridlonshift       = ncOut.add_var("VRM_GridShiftLon"      ,ncDouble);
    NcVar *VRMreverseorientation = ncOut.add_var("VRM_ReverseOrientation",ncChar  );
    NcVar *VRMrefinetype         = ncOut.add_var("VRM_RefineType",ncChar,dimString);
    NcVar *VRMgridtype           = ncOut.add_var("VRM_GridType"  ,ncChar,dimLenString);
    NcVar *VRMsmoothtype         = ncOut.add_var("VRM_SmoothType",ncChar,dimLenString);

    VRMresolution        ->put(&VRM.Resolution);
    VRMrefinement        ->put(&VRM.RefinementLevel);
    VRMtessellations     ->put(&VRM.Tessellations);
    VRMsubcellresolution ->put(&VRM.SubCellResolution);
    VRMtranssmoothdist   ->put(&VRM.TransSmoothDist);
    VRMsmoothiterations  ->put(&VRM.SmoothIterations);
    VRMgridxrotate       ->put(&VRM.GridXRotate);
    VRMgridyrotate       ->put(&VRM.GridYRotate);
    VRMgridlonshift      ->put(&VRM.GridLonShift);
    if(VRM.ReverseOrientation) {
        VRMreverseorientation->put("T");
    } else {
        VRMreverseorientation->put("F");
    }
    strcpy(OutString,OutBlank);
    VRM.RefineType.copy(OutString,VRM.RefineType.length(),0);
    VRMrefinetype->put(OutString,MyStrLen);
    strcpy(OutString,OutBlank);
    VRM.GridType.copy(OutString,VRM.GridType.length(),0);
    VRMgridtype->put(OutString,MyStrLen);
    strcpy(OutString,OutBlank);
    VRM.SmoothType.copy(OutString,VRM.SmoothType.length(),0);
    VRMsmoothtype->put(OutString,MyStrLen);
    delete [] OutBlank;
    delete [] OutString;

    // Add information about grid
    //----------------------------
    NcVar *GRIDnumclockwise = ncOut.add_var("GRID_Num_Clockwise",ncInt   );
    NcVar *GRIDnumcounter   = ncOut.add_var("GRID_Num_Counter"  ,ncInt   );
    NcVar *GRIDtotalarea    = ncOut.add_var("GRID_Total_Area"   ,ncDouble);
    NcVar *GRIDarearatio    = ncOut.add_var("GRID_Area_Ratio"   ,ncDouble);
    NcVar *GRIDminarea      = ncOut.add_var("GRID_MinArea"      ,ncDouble);
    NcVar *GRIDminarealon   = ncOut.add_var("GRID_MinArea_Lon"  ,ncDouble);
    NcVar *GRIDminarealat   = ncOut.add_var("GRID_MinArea_Lat"  ,ncDouble);
    NcVar *GRIDmaxarea      = ncOut.add_var("GRID_MaxArea"      ,ncDouble);
    NcVar *GRIDmaxarealon   = ncOut.add_var("GRID_MaxArea_Lon"  ,ncDouble);
    NcVar *GRIDmaxarealat   = ncOut.add_var("GRID_MaxArea_Lat"  ,ncDouble);
    NcVar *GRIDqarc         = ncOut.add_var("GRID_Qarc"         ,ncDouble);
    NcVar *GRIDminarc       = ncOut.add_var("GRID_MinArc"       ,ncDouble);
    NcVar *GRIDminarclon    = ncOut.add_var("GRID_MinArc_Lon"   ,ncDouble);
    NcVar *GRIDminarclat    = ncOut.add_var("GRID_MinArc_Lat"   ,ncDouble);
    NcVar *GRIDmaxarc       = ncOut.add_var("GRID_MaxArc"       ,ncDouble);
    NcVar *GRIDmaxarclon    = ncOut.add_var("GRID_MaxArc_Lon"   ,ncDouble);
    NcVar *GRIDmaxarclat    = ncOut.add_var("GRID_MaxArc_Lat"   ,ncDouble);
    NcVar *GRIDqang         = ncOut.add_var("GRID_Qang"         ,ncDouble);
    NcVar *GRIDminang       = ncOut.add_var("GRID_MinAng"       ,ncDouble);
    NcVar *GRIDminanglon    = ncOut.add_var("GRID_MinAng_Lon"   ,ncDouble);
    NcVar *GRIDminanglat    = ncOut.add_var("GRID_MinAng_Lat"   ,ncDouble);
    NcVar *GRIDmaxang       = ncOut.add_var("GRID_MaxAng"       ,ncDouble);
    NcVar *GRIDmaxanglon    = ncOut.add_var("GRID_MaxAng_Lon"   ,ncDouble);
    NcVar *GRIDmaxanglat    = ncOut.add_var("GRID_MaxAng_Lat"   ,ncDouble);

    GRIDnumclockwise ->put(&OutputGrid.Num_Clockwise);
    GRIDnumcounter   ->put(&OutputGrid.Num_Counter);
    GRIDtotalarea    ->put(&OutputGrid.Total_Area);
    GRIDarearatio    ->put(&OutputGrid.Area_Ratio);
    GRIDminarea      ->put(&OutputGrid.MinArea);
    GRIDminarealon   ->put(&OutputGrid.MinArea_Lon);
    GRIDminarealat   ->put(&OutputGrid.MinArea_Lat);
    GRIDmaxarea      ->put(&OutputGrid.MaxArea);
    GRIDmaxarealon   ->put(&OutputGrid.MaxArea_Lon);
    GRIDmaxarealat   ->put(&OutputGrid.MaxArea_Lat);
    GRIDqarc         ->put(&OutputGrid.Qarc);
    GRIDminarc       ->put(&OutputGrid.MinArc);
    GRIDminarclon    ->put(&OutputGrid.MinArc_Lon);
    GRIDminarclat    ->put(&OutputGrid.MinArc_Lat);
    GRIDmaxarc       ->put(&OutputGrid.MaxArc);
    GRIDmaxarclon    ->put(&OutputGrid.MaxArc_Lon);
    GRIDmaxarclat    ->put(&OutputGrid.MaxArc_Lat);
    GRIDqang         ->put(&OutputGrid.Qang);
    GRIDminang       ->put(&OutputGrid.MinAng);
    GRIDminanglon    ->put(&OutputGrid.MinAng_Lon);
    GRIDminanglat    ->put(&OutputGrid.MinAng_Lat);
    GRIDmaxang       ->put(&OutputGrid.MaxAng);
    GRIDmaxanglon    ->put(&OutputGrid.MaxAng_Lon);
    GRIDmaxanglat    ->put(&OutputGrid.MaxAng_Lat);

    // End Function OutputNetCDFFile()
    //-----------------------------------
    ncOut.close();
}
//=====================================================================


//=====================================================================
void Read_Refinement_Map( std::string     RefMapName, 
                          RefinementMap & InputRefineMap)
{
    int    NLAT;
    int    NLON;

    // Open a netCDF file for input
    //---------------------------------
    NcFile ncRef(RefMapName.c_str(), NcFile::ReadOnly);

    // inquire if there are lat,lon dimensions,
    // read values if present.
    //-----------------------------------------
    NcDim *latDim, *lonDim;
    if(!(latDim = ncRef.get_dim("latitude"))) {
        std::cout << "Read RefineMap ERROR: NO LAT DIMENSION" << std::endl;
    } else {
        NLAT = latDim->size();
    }
    if(!(lonDim = ncRef.get_dim("longitude"))) {
        std::cout << "Read RefineMap ERROR: NO LON DIMENSION" << std::endl;
    } else {
        NLON = lonDim->size();
    }

    // check for lat,lon coordinate variables,
    // read values if present
    //----------------------------------------
    double lat[NLAT];
    NcVar *latVar;
    if(!(latVar = ncRef.get_var("latitude"))) {
        std::cout << "Read RefineMap ERROR: error getting latitude variable" << std::endl;
    } else {
        if(!(latVar->get(lat,NLAT))) {
            std::cout << "Read RefineMap ERROR: error reading latirude values" << std::endl;
        }
    }

    double lon[NLON];
    NcVar *lonVar = ncRef.get_var("longitude");
    if(!(lonVar = ncRef.get_var("longitude"))) {
        std::cout << "Read RefineMap ERROR: error getting longitude variable" << std::endl;
    } else {
        if(!(lonVar->get(lon,NLON))) {
            std::cout << "Read RefineMap ERROR: error reading longitude values" << std::endl;
        }
    }

    // check that the lat,lon ranges conform to refeineMap requirements
    //------------------------------------------------------------------
    double latMin = lat[0];
    double latMax = lat[0];
    for(int jj=0; jj < NLAT; jj++) {
        if(latMin > lat[jj]) {latMin = lat[jj];}
        if(latMax < lat[jj]) {latMax = lat[jj];}
    }

    double lonMin = lon[0];
    double lonMax = lon[0];
    for(int ii=0; ii < NLON; ii++) {
        if(lonMin > lon[ii]) {lonMin = lon[ii];}
        if(lonMax < lon[ii]) {lonMax = lon[ii];}
    }

    if((latMin < -90.)||(latMax > 90.)||(lonMin < 0.)||(lonMax > 360.)) {
        std::cout << "Read RefineMap ERROR: Lat/Lon values out of range for refineMap file" << std::endl;
        std::cout << " LatMin= " << latMin << " LatMaX=" << latMax << std::endl;
        std::cout << " LonMin= " << lonMin << " LonMaX=" << lonMax << std::endl;
    }

    // check for refineMap variable,
    // read it in of present.
    //--------------------------------
    double rval[NLAT][NLON];
    NcVar *refMap;
    if(!(refMap = ncRef.get_var("refineMap"))) {
        std::cout << "Read RefineMap ERROR: error getting refineMap variable" << std::endl;
    } else {
        if(!(refMap->get(&rval[0][0],NLAT,NLON))) {
            std::cout << "Read RefineMap ERROR: error reading refineMap values" << std::endl;
        }
    }

    // check the range of the refine values
    //--------------------------------------
    double rvalMin = rval[0][0];
    double rvalMax = rval[0][0];
    for(int jj=0; jj < NLAT; jj++) {
        for(int ii=0; ii < NLON; ii++) {
            if(rvalMin > rval[jj][ii]) {rvalMin = rval[jj][ii];}
            if(rvalMax < rval[jj][ii]) {rvalMax = rval[jj][ii];}
        }
    }

    if((rvalMin < 0.)||(rvalMax > 1.)) {
        std::cout << "Read refineMap ERROR: refineMap values are outside of [0,1]" << std::endl;
    }
    if((rvalMin > 0.)||(rvalMax < 1.)) {
        std::cout << "Read refineMap ERROR: refineMap values do not span the range [0,1]" << std::endl;
    }

    // close the input file
    //-----------------------
    ncRef.close();

    // Transfer the refineMap values to the global data structure
    //------------------------------------------------------------
    std::cout << " Read in refineMap: NLAT=" << NLAT << " NLON=" << NLON << std::endl;
    InputRefineMap.resize(NLON,NLAT);

    for(int ii=0; ii < NLON; ii++) {
    for(int jj=0; jj < NLAT; jj++) {
        InputRefineMap.val[ii][jj] = rval[jj][ii];
    }
    }

}
//=====================================================================


//=====================================================================
void Write_Refinement_Map( std::string     RefMapName, 
                           RefinementMap & OutputRefineMap)
{
    // Copy refineMap values into local arrays
    //----------------------------------------
    int    NLAT = OutputRefineMap.nRefLat;
    int    NLON = OutputRefineMap.nRefLon;
    double lat[NLAT];
    double lon[NLON];
    double rval[NLAT][NLON];

    for(int ii=0; ii < NLON; ii++) {
        lon[ii] = OutputRefineMap.RefLon[ii][0];
    }
    for(int jj=0; jj < NLAT; jj++) {
        lat[jj] = OutputRefineMap.RefLat[jj][0];
    }
    for(int jj=0; jj < NLAT; jj++) {
        for(int ii=0; ii < NLON; ii++) {
            rval[jj][ii] = OutputRefineMap.val[ii][jj];
        }
    }

    // Open a netCDF file for output
    //---------------------------------
    NcFile ncRef(RefMapName.c_str(), NcFile::FileMode::Replace);

    // Add dimensions to output
    //---------------------------
    NcDim *latDim = ncRef.add_dim("latitude" , NLAT);
    NcDim *lonDim = ncRef.add_dim("longitude", NLON);

    // Add variable declarations to output
    //-------------------------------------
    NcVar *latVar = ncRef.add_var("latitude" , ncDouble, latDim);
    NcVar *lonVar = ncRef.add_var("longitude", ncDouble, lonDim);
    NcVar *refMap = ncRef.add_var("refineMap", ncDouble, latDim, lonDim);

    // Write values to output
    //-----------------------
    latVar->put(lat,NLAT);
    lonVar->put(lon,NLON);
    refMap->put( &rval[0][0], NLAT, NLON);

    // Close up shop
    //----------------
    ncRef.close();
}
//=====================================================================


//=====================================================================
void Read_Reference_Map( std::string     RefMapName,
                         RefinementMap & OutputRefineMapRef)
{
    int    NLAT;
    int    NLON;

    // Open a netCDF file for input
    //---------------------------------
    NcFile ncRef(RefMapName.c_str(), NcFile::ReadOnly);

    // inquire if there are lat,lon dimensions,
    // read values if present.
    //-----------------------------------------
    NcDim *latDim, *lonDim;
    if(!(latDim = ncRef.get_dim("latitude"))) {
        std::cout << "Read ReferenceMap ERROR: NO LAT DIMENSION" << std::endl;
    } else {
        NLAT = latDim->size();
    }
    if(!(lonDim = ncRef.get_dim("longitude"))) {
        std::cout << "Read ReferenceMap ERROR: NO LON DIMENSION" << std::endl;
    } else {
        NLON = lonDim->size();
    }

    // check for lat,lon coordinate variables,
    // read values if present
    //----------------------------------------
    double lat[NLAT];
    NcVar *latVar;
    if(!(latVar = ncRef.get_var("latitude"))) {
        std::cout << "Read ReferenceMap ERROR: error getting latitude variable" << std::endl;
    } else {
        if(!(latVar->get(lat,NLAT))) {
            std::cout << "Read ReferenceMap ERROR: error reading latirude values" << std::endl;
        }
    }

    double lon[NLON];
    NcVar *lonVar = ncRef.get_var("longitude");
    if(!(lonVar = ncRef.get_var("longitude"))) {
        std::cout << "Read ReferenceMap ERROR: error getting longitude variable" << std::endl;
    } else {
        if(!(lonVar->get(lon,NLON))) {
            std::cout << "Read ReferenceMap ERROR: error reading longitude values" << std::endl;
        }
    }

    // check that the lat,lon ranges conform to referenceMap requirements
    //------------------------------------------------------------------
    double latMin = lat[0];
    double latMax = lat[0];
    for(int jj=0; jj < NLAT; jj++) {
        if(latMin > lat[jj]) {latMin = lat[jj];}
        if(latMax < lat[jj]) {latMax = lat[jj];}
    }

    double lonMin = lon[0];
    double lonMax = lon[0];
    for(int ii=0; ii < NLON; ii++) {
        if(lonMin > lon[ii]) {lonMin = lon[ii];}
        if(lonMax < lon[ii]) {lonMax = lon[ii];}
    }

    if((latMin < -90.)||(latMax > 90.)||(lonMin < 0.)||(lonMax > 360.)) {
        std::cout << "Read ReferenceMap ERROR: Lat/Lon values out of range for refineMap file" << std::endl;
        std::cout << " LatMin= " << latMin << " LatMaX=" << latMax << std::endl;
        std::cout << " LonMin= " << lonMin << " LonMaX=" << lonMax << std::endl;
    }

    // check for refineMap variable,
    // read it in of present.
    //--------------------------------
    double rval[NLAT][NLON];
    NcVar *refMap;
    if(!(refMap = ncRef.get_var("refineMap"))) {
        std::cout << "Read ReferenceMap ERROR: error getting refineMap variable" << std::endl;
    } else {
        if(!(refMap->get(&rval[0][0],NLAT,NLON))) {
            std::cout << "Read ReferenceMap ERROR: error reading refineMap values" << std::endl;
        }
    }

    // calc the Max/Min refineMap values
    //--------------------------------------
    double rvalMin = rval[0][0];
    double rvalMax = rval[0][0];
    for(int jj=0; jj < NLAT; jj++) {
        for(int ii=0; ii < NLON; ii++) {
            if(rvalMin > rval[jj][ii]) {rvalMin = rval[jj][ii];}
            if(rvalMax < rval[jj][ii]) {rvalMax = rval[jj][ii];}
        }
    }

    // close the input file
    //-----------------------
    ncRef.close();

    // Transfer the refineMap values to the global data structure
    //------------------------------------------------------------
    std::cout << " Read in referenceMap: NLAT=" << NLAT << " NLON=" << NLON << std::endl;
    std::cout << "       Min=" << rvalMin << " Max=" << rvalMax << std::endl;
    OutputRefineMapRef.resize(NLON,NLAT);
    OutputRefineMapRef.Min = rvalMin;
    OutputRefineMapRef.Max = rvalMax;

    for(int ii=0; ii < NLON; ii++) {
    for(int jj=0; jj < NLAT; jj++) {
        OutputRefineMapRef.val[ii][jj] = rval[jj][ii];
    }
    }
//=====================================================================

}
