#include "vrmview.h"
#include <iostream>

VrmView::VrmView(QWidget *parent)
{

}

VrmView::VrmView(int I_HRES, int I_VRES, QWidget *parent) : QGraphicsView(parent)
{
    HRES   = I_HRES;
    VRES   = I_VRES;
    BRES   = VRES/10;
    DltLon = 360.;
    DltLat = 180.;

    Lon0   =  180.;
    Lat0   =    0.;
    LonMin = -180.;
    LonMax =  540.;
    LatMin =  -90.;
    LatMax =   90.;

    ZoomScl    = 1.0;
    ZoomSclMax = 100.0 ; // 6.0;
    ZoomSclMin = 1.0;
    ZoomSclInc = 1.10;
    SclLon =  (HRES/DltLon);
    SclLat = -(VRES/DltLat);

    LonC1 = Lon0 - DltLon/2.;
    LatC1 = Lat0 - DltLat/2.;
    LonC2 = LonC1 + DltLon;
    LatC2 = LatC1 + DltLat;

    if(false) {   // DIAG
        ZoomScl = 1.;
        LonC1 = Lon0 - DltLon/(2.*ZoomScl);
        LatC1 = Lat0 - DltLat/(2.*ZoomScl);
        LonC2 = LonC1 + DltLon/ZoomScl;
        LatC2 = LatC1 + DltLat/ZoomScl;
    }
}

void VrmView::initialize()
{
    int HSIZE   = HRES+2;
    int VSIZE   = VRES+2;
    int LONSIZE = BRES+2;
    int LATSIZE = BRES+2;
    double slon = SclLon*ZoomScl;
    double slat = SclLat*ZoomScl;
    emit setViewWindows(HSIZE,  VSIZE,  LATSIZE, LONSIZE, Lon0, Lat0, slon, slat);
}

double VrmView::getLonMin()
{
    return LonMin;
}

double VrmView::getLonMax()
{
    return LonMax;
}

double VrmView::getLon0()
{
    return Lon0;
}

void VrmView::setLon0(double I_Lon0, double I_LonC1, double I_LonC2)
{
    Lon0  = I_Lon0;
    LonC1 = I_LonC1;
    LonC2 = I_LonC2;
}

double VrmView::getZoomScl()
{
    return ZoomScl;
}

double VrmView::getZoomSclMin()
{
    return ZoomSclMin;
}

double VrmView::getZoomSclMax()
{
    return ZoomSclMax;
}

double VrmView::getZoomSclInc()
{
    return ZoomSclInc;
}

void VrmView::setZoomScl(double I_ZoomScl)
{
    ZoomScl = I_ZoomScl;
}

double VrmView::getDltLon()
{
    return DltLon;
}

double VrmView::getSclLon()
{
    return SclLon;
}

double VrmView::getLatMin()
{
    return LatMin;
}

double VrmView::getLatMax()
{
    return LatMax;
}

double VrmView::getLat0()
{
    return Lat0;
}

void VrmView::setLat0(double I_Lat0, double I_LatC1, double I_LatC2)
{
    Lat0  = I_Lat0;
    LatC1 = I_LatC1;
    LatC2 = I_LatC2;
}

double VrmView::getDltLat()
{
    return DltLat;
}

double VrmView::getSclLat()
{
    return SclLat;
}

int VrmView::getBarHeight()
{
    return BRES;
}

void VrmView::wheelEvent(QWheelEvent *event)
{

}

