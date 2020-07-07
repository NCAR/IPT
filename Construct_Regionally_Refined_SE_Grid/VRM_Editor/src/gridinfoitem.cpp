#include "gridinfoitem.h"

gridInfoItem::gridInfoItem()
{
    // Initially disable Grid info
    //-----------------------------
    Enabled = false;

    // Initialize Values
    //--------------------
    Num_Clockwise = 0;
    Num_Counter   = 0;

    Qarc          = 1.0;
    MinArc        =  999.;
    MinArc_Lon    = -999.;
    MinArc_Lat    = -999.;
    MaxArc        = -999.;
    MaxArc_Lon    = -999.;
    MaxArc_Lat    = -999.;

    Qang          = 1.0;
    MinAng        =  999.;
    MinAng_Lon    = -999.;
    MinAng_Lat    = -99.;
    MaxAng        = -999.;
    MaxAng_Lon    = -999.;
    MaxAng_Lat    = -99.;

    Total_Area    = 1.0;
    Area_Ratio    = 1.0;
    MinArea       =  999.;
    MinArea_Lon   = -999.;
    MinArea_Lat   = -99.;
    MaxArea       = -999.;
    MaxArea_Lon   = -999.;
    MaxArea_Lat   = -99.;
}

void gridInfoItem::Set_GridInfo(const CubeGrid &I_CubeGrid)
{
    Num_Clockwise = I_CubeGrid.Num_Clockwise;
    Num_Counter   = I_CubeGrid.Num_Counter;
    Qarc          = I_CubeGrid.Qarc;
    MinArc        = I_CubeGrid.MinArc;
    MinArc_Lon    = I_CubeGrid.MinArc_Lon;
    MinArc_Lat    = I_CubeGrid.MinArc_Lat;
    MaxArc        = I_CubeGrid.MaxArc;
    MaxArc_Lon    = I_CubeGrid.MaxArc_Lon;
    MaxArc_Lat    = I_CubeGrid.MaxArc_Lat;
    Qang          = I_CubeGrid.Qang;
    MinAng        = I_CubeGrid.MinAng;
    MinAng_Lon    = I_CubeGrid.MinAng_Lon;
    MinAng_Lat    = I_CubeGrid.MinAng_Lat;
    MaxAng        = I_CubeGrid.MaxAng;
    MaxAng_Lon    = I_CubeGrid.MaxAng_Lon;
    MaxAng_Lat    = I_CubeGrid.MaxAng_Lat;
    Total_Area    = I_CubeGrid.Total_Area;
    Area_Ratio    = I_CubeGrid.Area_Ratio;
    MinArea       = I_CubeGrid.MinArea;
    MaxArea       = I_CubeGrid.MaxArea;
    MinArea_Lon   = I_CubeGrid.MinArea_Lon;
    MinArea_Lat   = I_CubeGrid.MinArea_Lat;
    MaxArea_Lon   = I_CubeGrid.MaxArea_Lon;
    MaxArea_Lat   = I_CubeGrid.MaxArea_Lat;
}

void gridInfoItem::Set_GridInfo(int    I_Num_Clockwise,
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
                                double I_MaxArea_Lat)
{
    Num_Clockwise = I_Num_Clockwise;
    Num_Counter   = I_Num_Counter;
    Qarc          = I_Qarc;
    MinArc        = I_MinArc;
    MinArc_Lon    = I_MinArc_Lon;
    MinArc_Lat    = I_MinArc_Lat;
    MaxArc        = I_MaxArc;
    MaxArc_Lon    = I_MaxArc_Lon;
    MaxArc_Lat    = I_MaxArc_Lat;
    Qang          = I_Qang;
    MinAng        = I_MinAng;
    MinAng_Lon    = I_MinAng_Lon;
    MinAng_Lat    = I_MinAng_Lat;
    MaxAng        = I_MaxAng;
    MaxAng_Lon    = I_MaxAng_Lon;
    MaxAng_Lat    = I_MaxAng_Lat;
    Total_Area    = I_Total_Area;
    Area_Ratio    = I_Area_Ratio;
    MinArea       = I_MinArea;
    MaxArea       = I_MaxArea;
    MinArea_Lon   = I_MinArea_Lon;
    MinArea_Lat   = I_MinArea_Lat;
    MaxArea_Lon   = I_MaxArea_Lon;
    MaxArea_Lat   = I_MaxArea_Lat;
}

void gridInfoItem::Set_Enabled(bool I_Enabled)
{
    Enabled = I_Enabled;
}

bool gridInfoItem::Get_Enabled()
{
    return Enabled;
}

int gridInfoItem::Get_Num_Clockwise()
{
    return Num_Clockwise;
}

int gridInfoItem::Get_Num_Counter()
{
    return Num_Counter;
}

double gridInfoItem::Get_Qarc()
{
    return Qarc;
}

double gridInfoItem::Get_MinArc()
{
    return MinArc;
}

double gridInfoItem::Get_MinArc_Lon()
{
    return MinArc_Lon;
}

double gridInfoItem::Get_MinArc_Lat()
{
    return MinArc_Lat;
}

double gridInfoItem::Get_MaxArc()
{
    return MaxArc;
}

double gridInfoItem::Get_MaxArc_Lon()
{
    return MaxArc_Lon;
}

double gridInfoItem::Get_MaxArc_Lat()
{
    return MaxArc_Lat;
}

double gridInfoItem::Get_Qang()
{
    return Qang;
}

double gridInfoItem::Get_MinAng()
{
    return MinAng;
}

double gridInfoItem::Get_MinAng_Lon()
{
    return MinAng_Lon;
}

double gridInfoItem::Get_MinAng_Lat()
{
    return MinAng_Lat;
}

double gridInfoItem::Get_MaxAng()
{
    return MaxAng;
}

double gridInfoItem::Get_MaxAng_Lon()
{
    return MaxAng_Lon;
}

double gridInfoItem::Get_MaxAng_Lat()
{
    return MaxAng_Lat;
}

double gridInfoItem::Get_Total_Area()
{
    return Total_Area;
}

double gridInfoItem::Get_Area_Ratio()
{
    return Area_Ratio;
}

double gridInfoItem::Get_MinArea()
{
    return MinArea;
}

double gridInfoItem::Get_MinArea_Lon()
{
    return MinArea_Lon;
}

double gridInfoItem::Get_MinArea_Lat()
{
    return MinArea_Lat;
}

double gridInfoItem::Get_MaxArea()
{
    return MaxArea;
}

double gridInfoItem::Get_MaxArea_Lon()
{
    return MaxArea_Lon;
}

double gridInfoItem::Get_MaxArea_Lat()
{
    return MaxArea_Lat;
}
