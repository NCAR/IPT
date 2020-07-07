#include <cmath>
#include "SQuadGen/GridElements.h"
#include "SQuadGen/CubedSphereGrid.h"
#include "SQuadGen/IcosahedralFlagGrid.h"
#include "SQuadGen/RefinementTemplateCUBIT.h"
#include "SQuadGen/RefinementTemplateLOWCONN.h"
#include "SQuadGen/RefinementTemplateLOWCONNOLD.h"
#include "SQuadGen/Tessellate.h"
#include "SQuadGen/SpringDynamics.h"
#include "SQuadGen/MathHelper.h"
#include "SQuadGen/Exception.h"

#include "CommandLine.h"
#include "ReadWrite.h"
#include "RefinementMap.h"
#include "RefinementCube.h"
#include "CubeGrid.h"
#include "RefineCubeGrid.h"
#include "NeighborAdjustments.h"


int main(int argc, char** argv)
//
// Create_VRMgrid: For the current refinement map and processing options,
//               calcualte the ...
//
//=====================================================================
{

    int         MyResolution         = 30;
    std::string MyRefineType         = "CUBIT";
    std::string MyGridType           = "CubeSquared";
    double      MyGridXRotate        = 0.0;
    double      MyGridYRotate        = 0.0;
    double      MyGridLonShift       = 0.0;
    int         MyRefinementLevel    = 0;
    int         MyTessellations      = 0;
    int         MySubCellResolution  = 0;
    bool        MyReverseOrientation = false;
    std::string MySmoothType         = "NONE";
    int         MyTransSmoothDist    = 0;
    int         MySmoothIterations   = 0;
    std::string MyRefineFile         = "NULL";
    std::string MyOutputFile         = "NULL";

    RefinementMap  refineMap;
    RefinementCube refineCube;
    CubeGrid       refineGrid;
    CubeGrid       baseGrid;
    NodeVector     vecNodes;
    FaceVector     vecFaces;
    

    // Process command line and create a grid
    //=========================================
try {

    // Parse the command line
    //------------------------
    std::cout << "  " << std::endl;
    BeginCommandLine()
       CommandLineStringD(MyRefineType        ,"refine_type"   ,"CUBIT", 
                                               "(Options: CUBIT | LOWCONN | LOWCONNold)");
       CommandLineStringD(MyGridType          ,"grid_type"     ,"CubeSquared", 
                                               "(Options: Icosahedral | CubeSquared)");
       CommandLineStringD(MySmoothType        ,"smooth_type"   ,"NONE", 
                                               "(Options: NONE | SPRING | PRESSURE)");
       CommandLineInt    (MyResolution        ,"resolution"    ,30);
       CommandLineInt    (MyRefinementLevel   ,"refine_level"  ,0);
       CommandLineInt    (MyTessellations     ,"tessellate"    ,0);
       CommandLineInt    (MySubCellResolution ,"subcells"      ,0);
       CommandLineInt    (MyTransSmoothDist   ,"smooth_dist"   ,0);
       CommandLineInt    (MySmoothIterations  ,"smooth_iter"   ,0);
       CommandLineBool   (MyReverseOrientation,"reverse_orient"  );
       CommandLineDouble (MyGridXRotate       ,"x_rotate"      ,0.0);
       CommandLineDouble (MyGridYRotate       ,"x_rotate"      ,0.0);
       CommandLineDouble (MyGridLonShift      ,"lon_shift"     ,0.0);
       CommandLineString (MyRefineFile        ,"refine_file"   ,"NULL");
       CommandLineString (MyOutputFile        ,"output"        ,"NULL");

       ParseCommandLine(argc, argv);
    EndCommandLine(argv)
    std::cout << "  " << std::endl;

    // Check input values
    //------------------------
    if(MyOutputFile == "NULL") {
      std::cout << argv[0] << ": No output file specified" << std::endl;
      return(-1);
    }

    if(MyRefineFile == "NULL") {
      std::cout << argv[0] << ": No refinement file specified" << std::endl;
      return(-1);
    }

    // Read in Refinement Map
    //------------------------
    Read_Refinement_Map(MyRefineFile, refineMap);

    // Check for odd resolution for LOWCONN refile types
    //----------------------------------------------------
    int nResolution = MyResolution;
    if ((MyRefineType == "LOWCONN") || (MyRefineType == "LOWCONNold") ) {
        if (nResolution % 2 == 1) {
            std::cout << " ERROR:: RefineType = "          << MyRefineType << std::endl;
            std::cout << " ERROR:: nResolution is ODD! = " << nResolution  << std::endl;
            return(-1);
        } else {
            nResolution /= 2;
        }
    }

    // Generate Base Grid
    //----------------------------
    if(MyGridType == "CubeSquared") {
//        std::cout << " GridType = CubeSquared" << std::endl;
//        std::cout << " nResolution = "<< nResolution << std::endl;
        GenerateCubedSphere( nResolution, vecNodes, vecFaces);
    } else if (MyGridType == "Icosahedral") {
//        std::cout << " GridType = Icosahedral" << std::endl;
//        std::cout << " nResolution = "<< nResolution << std::endl;
        GenerateIcosahedralQuadGrid( nResolution, vecNodes, vecFaces);
    } else {
        std::cout << " ERROR:: GridType = UNKNOWN" << std::endl;
    }

    // Rotate that Grid around the X-Axis
    //------------------------------------
    if( MyGridXRotate != 0.0) {
        double dCosTheta = cos(MyGridXRotate*M_PI/180.0);
        double dSinTheta = sin(MyGridXRotate*M_PI/180.0);
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            double dTempY = vecNodes[i].y;
            double dTempZ = vecNodes[i].z;
            vecNodes[i].y = dCosTheta*dTempY - dSinTheta*dTempZ;
            vecNodes[i].z = dSinTheta*dTempY + dCosTheta*dTempZ;
        }
    }

    // Rotate that Grid around the Y-Axis
    //------------------------------------
    if( MyGridYRotate != 0.0) {
        double dCosTheta = cos(MyGridYRotate*M_PI/180.0);
        double dSinTheta = sin(MyGridYRotate*M_PI/180.0);
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            double dTempX = vecNodes[i].x;
            double dTempZ = vecNodes[i].z;
            vecNodes[i].x = dCosTheta*dTempX - dSinTheta*dTempZ;
            vecNodes[i].z = dSinTheta*dTempX + dCosTheta*dTempZ;
        }
    }

    // Add Longitude Shift to that Grid
    //-----------------------------------
    if (MyGridLonShift != 0.0) {
        double dLon;
        double dLat;
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            dLat = asin(vecNodes[i].z);
            dLon = atan2(vecNodes[i].y,vecNodes[i].x);
            dLon += MyGridLonShift * M_PI/180.0;
            vecNodes[i].x = cos(dLat)*cos(dLon);
            vecNodes[i].y = cos(dLat)*sin(dLon);
        }
    }

    // Store Base grid values
    //-------------------------
    baseGrid.loadGrid(vecNodes, vecFaces);
//    std::cout << " baseGrid nodes set N= "<< vecNodes.size() << " : " << baseGrid.Nodes.size() << std::endl;
//    std::cout << " baseGrid faces set N= "<< vecFaces.size() << " : " << baseGrid.Faces.size() << std::endl;
//    std::cout << " baseGrid grid  set N= " << baseGrid.Grid.size() << std::endl;

    // Initialize the reference Cube values
    //---------------------------------------
    refineCube.resize(nResolution, MyRefinementLevel, MyGridLonShift, MyGridXRotate, MyGridYRotate);
    refineCube.loadCubeVals(refineMap);
    refineCube.Normalize();

    // Perform Refinement if requested
    //----------------------------------
    if (MyRefinementLevel > 0) {
        if(MyRefineType == "CUBIT") {
            RefinementTemplateCUBIT RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else if( MyRefineType == "LOWCONN") {
            RefinementTemplateLOWCONN RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else if ( MyRefineType == "LOWCONNold") {
            RefinementTemplateLOWCONNOLD RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else {
            // *** ERROR ***
        }
    }

    // Tessellate that Grid
    //------------------------
    if ((MyTessellations < 0) || (MyTessellations > 100)) {
        std::cout << " ERROR:: MyTesselations out of Range = "<< MyTessellations << std::endl;
    } else {
        for (int n=0; n<MyTessellations; n++) {
            Tessellate( vecNodes, vecFaces);
        }
    }

    // Add Sub-Cell Resolution to that Grid
    //---------------------------------------
    if (MySubCellResolution < 0) {
        std::cout << " ERROR:: MySubCellResolution out of Range = "<< MySubCellResolution << std::endl;
    } else if (MySubCellResolution == 0){
        // Do Nothing
    } else {
        RefineEverything( vecNodes, vecFaces, MySubCellResolution+1 );
    }

    // Reverse Orientation of Faces
    //-------------------------------
    if ((MyReverseOrientation) && (MyGridType != "CubeSquared")) {
        ReverseFaceOrientation( vecFaces);
    } else if ((!MyReverseOrientation) && (MyGridType == "CubeSquared")) {
        ReverseFaceOrientation( vecFaces);
    }


    // Smooth that Grid
    //-------------------
    if (MySmoothType == "NONE") {
        // No Smoothing
        //--------------
    } else if (MySmoothType == "SPRING") {
        SpringDynamics( vecNodes, vecFaces, MyTransSmoothDist, MySmoothIterations);
    } else if (MySmoothType == "PRESSURE") {
        PressureDynamics( vecNodes, vecFaces, MyTransSmoothDist, MySmoothIterations);
    } else if (MySmoothType == "DISTANCE") {
//        int NE_lo = MyResolution;
//        int NE_hi = MyResolution*IntPow(2,MyRefinementLevel);
//        std::cout << " NE_lo= "<< NE_lo << " NE_hi= " << NE_hi << std::endl;
//        EquializeDistances(refineMap, vecNodes, vecFaces, MyTransSmoothDist, 
//                                                          MySmoothIterations, NE_lo, NE_hi);
    } else {
        std::cout << " ERROR:: Unknown MySmoothType = "<< MySmoothType << std::endl;
    }

    // Store final Refined grid values
    //---------------------------------
    refineGrid.loadGrid( vecNodes, vecFaces);

    // Write out Grid Quality parameters:
    //------------------------------------
    std::cout << "  " << std::endl;
    std::cout << "Grid Quality Parameters: " << std::endl;
    std::cout << "-------------------------  " << std::endl;
    std::cout << " Number of Nodes= "<< refineGrid.Nodes.size() << std::endl;
    std::cout << " Number of Faces= "<< refineGrid.Faces.size() << std::endl;
    std::cout << "         Clockwise Faces= "<< refineGrid.Num_Clockwise << std::endl;
    std::cout << " Counter-Clockwise Faces= "<< refineGrid.Num_Counter   << std::endl;
    std::cout << "              Total Area= "<< refineGrid.Total_Area    << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "                Min Area= "<< refineGrid.MinArea     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MinArea_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MinArea_Lat << std::endl;
    std::cout << "                Max Area= "<< refineGrid.MaxArea     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MaxArea_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MaxArea_Lat << std::endl;
    std::cout << "              Area Ratio= "<< refineGrid.Area_Ratio  << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "            Arc Length Q= "<< refineGrid.Qarc       << std::endl;
    std::cout << "          Min Arc Length= "<< refineGrid.MinArc     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MinArc_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MinArc_Lat << std::endl;
    std::cout << "          Max Arc Length= "<< refineGrid.MaxArc     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MaxArc_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MaxArc_Lat << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "                 Angle Q= "<< refineGrid.Qang       << std::endl;
    std::cout << "               Min Angle= "<< refineGrid.MinAng     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MinAng_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MinAng_Lat << std::endl;
    std::cout << "               Max Angle= "<< refineGrid.MaxAng     << std::endl;
    std::cout << "                    @Lon= "<< refineGrid.MaxAng_Lon << std::endl;
    std::cout << "                    @Lat= "<< refineGrid.MaxAng_Lat << std::endl;
    std::cout << "------------------------------------------------  " << std::endl;
    std::cout << "  " << std::endl;

    double           Qarc, MinArc, MaxArc;
    double           MinArc_Lon,MinArc_Lat,MaxArc_Lon,MaxArc_Lat;
    double           Qang, MinAng, MaxAng;
    double           MinAng_Lon,MinAng_Lat,MaxAng_Lon,MaxAng_Lat;


    // The Read/Write functions need to be independent of Qt(Strings) 
    // so they can be shared by the standalone program.
    // So copy values into the VRM structure passed to Write_Exodus_File()
    //-----------------------------------------------------------------
    VRM_Param VRM;
    VRM.Resolution         = MyResolution;
    VRM.RefineType         = MyRefineType;
    VRM.GridType           = MyGridType;
    VRM.GridXRotate        = MyGridXRotate;
    VRM.GridYRotate        = MyGridYRotate;
    VRM.GridLonShift       = MyGridLonShift;
    VRM.RefinementLevel    = MyRefinementLevel;
    VRM.Tessellations      = MyTessellations;
    VRM.SubCellResolution  = MySubCellResolution;
    VRM.ReverseOrientation = MyReverseOrientation;
    VRM.SmoothType         = MySmoothType;
    VRM.TransSmoothDist    = MyTransSmoothDist;
    VRM.SmoothIterations   = MySmoothIterations;

    // Now write out the EXODUS file
    //-------------------------------
    Write_Exodus_File(MyOutputFile,VRM, refineMap, refineGrid);

    vecNodes.clear();
    vecFaces.clear();
    return 0;

} catch (Exception e) {
    std::cout << e.ToString() << std::endl;
    std::cout << "Create_VRMgrid failed" << std::endl;
}

}
