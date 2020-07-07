#ifndef GRIDINFOITEM_H
#define GRIDINFOITEM_H

#include "CubeGrid.h"

class gridInfoItem
{
public:
    gridInfoItem();

    void Set_GridInfo(const CubeGrid & I_CubeGrid);
    void Set_GridInfo(int    I_Num_Clockwise,
                      int    I_Num_Counter,
                      double I_Qarc,
                      double I_MinArc,
                      double I_MinArc_Lon,
                      double I_MinArc_Lat,
                      double I_MaxArc,
                      double I_MaxArc_Lon,
                      double I_MaxArc_Lat,
                      double I_Qang,
                      double I_MinAng,
                      double I_MinAng_Lon,
                      double I_MinAng_Lat,
                      double I_MaxAng,
                      double I_MaxAng_Lon,
                      double I_MaxAng_Lat,
                      double I_Total_Area,
                      double I_Area_Ratio,
                      double I_MinArea,
                      double I_MaxArea,
                      double I_MinArea_Lon,
                      double I_MinArea_Lat,
                      double I_MaxArea_Lon,
                      double I_MaxArea_Lat);
    void Set_Enabled(bool I_Enabled);

    bool   Get_Enabled();
    int    Get_Num_Clockwise();
    int    Get_Num_Counter();

    double Get_Qarc();
    double Get_MinArc();
    double Get_MinArc_Lon();
    double Get_MinArc_Lat();
    double Get_MaxArc();
    double Get_MaxArc_Lon();
    double Get_MaxArc_Lat();

    double Get_Qang();
    double Get_MinAng();
    double Get_MinAng_Lon();
    double Get_MinAng_Lat();
    double Get_MaxAng();
    double Get_MaxAng_Lon();
    double Get_MaxAng_Lat();

    double Get_Total_Area();
    double Get_Area_Ratio();
    double Get_MinArea();
    double Get_MinArea_Lon();
    double Get_MinArea_Lat();
    double Get_MaxArea();
    double Get_MaxArea_Lon();
    double Get_MaxArea_Lat();

private:
    bool             Enabled;
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

#endif // GRIDINFOITEM_H
