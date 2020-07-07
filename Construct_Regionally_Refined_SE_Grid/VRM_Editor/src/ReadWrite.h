#ifndef READWRITE_H
#define READWRITE_H

#include "SQuadGen/netcdfcpp.h"
#include <iostream>

#include "CubeGrid.h"
#include "RefinementMap.h"

//=====================================================================
struct VRM_Param {
    int         Resolution;
    std::string RefineType;
    std::string GridType;
    double      GridXRotate;
    double      GridYRotate;
    double      GridLonShift;
    int         RefinementLevel;
    int         Tessellations;
    int         SubCellResolution;
    bool        ReverseOrientation;
    std::string SmoothType;
    int         TransSmoothDist;
    int         SmoothIterations;
};
//=====================================================================


//=====================================================================
void Write_Exodus_File(std::string     ExodusFile,
                       VRM_Param &     VRM,
                       RefinementMap & OutputRefineMap,
                       CubeGrid &      OutputGrid);
//=====================================================================


//=====================================================================
void Read_Refinement_Map( std::string     RefMapName, 
                          RefinementMap & InputRefineMap);
//=====================================================================


//=====================================================================
void Write_Refinement_Map(std::string     RefMapName,
                           RefinementMap & OutputRefineMap);
//=====================================================================


//=====================================================================
void Read_Reference_Map( std::string     RefMapName,
                         RefinementMap & OutputRefineMapRef);
//=====================================================================

#endif // READWRITE_H
