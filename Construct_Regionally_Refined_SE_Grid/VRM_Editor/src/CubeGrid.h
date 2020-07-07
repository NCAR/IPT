#ifndef CUBEGRID_H
#define CUBEGRID_H

#include <cmath>
#include "SQuadGen/GridElements.h"

class CubeGrid
{
public:
    CubeGrid();
    CubeGrid( const CubeGrid &srcGrid);
    CubeGrid &operator =(const CubeGrid srcGrid);
    ~CubeGrid();

    void loadGrid(const NodeVector &srcNodes, const FaceVector &srcFaces);

private:
    void computeGridQuality();

public:
    NodeVector       Nodes;
    LonLatNodeVector Grid;
    FaceVector       Faces;
    int              Num_Clockwise;
    int              Num_Counter;
    double           Total_Area;
    double           Qarc, MinArc, MaxArc;
    double           MinArc_Lon,MinArc_Lat,MaxArc_Lon,MaxArc_Lat;
    double           Qang, MinAng, MaxAng;
    double           MinAng_Lon,MinAng_Lat,MaxAng_Lon,MaxAng_Lat;
    double           Area_Ratio, MinArea, MaxArea;
    double           MinArea_Lon,MinArea_Lat,MaxArea_Lon,MaxArea_Lat;
};

#endif // CUBEGRID_H
