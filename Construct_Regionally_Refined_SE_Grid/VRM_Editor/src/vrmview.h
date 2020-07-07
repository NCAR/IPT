#ifndef VRMVIEW_H
#define VRMVIEW_H

#include <QGraphicsView>
#include <QtWidgets>
#include <QWidget>
#include <QtGui>
#include <QKeyEvent>

class VrmView : public QGraphicsView
{
    Q_OBJECT

public:
    VrmView(QWidget *parent = 0);
    VrmView(int I_HRES, int I_VRES, QWidget *parent = 0);
    void   initialize();
    int    getBarHeight();

    double getLonMin();
    double getLonMax();
    double getLon0();
    double getDltLon();
    double getSclLon();
    void   setLon0(double I_Lon0, double I_LonC1, double I_LonC2);

    double getZoomScl();
    double getZoomSclMin();
    double getZoomSclMax();
    double getZoomSclInc();
    void   setZoomScl(double I_ZoomScl);

    double getLatMin();
    double getLatMax();
    double getLat0();
    double getDltLat();
    double getSclLat();
    void   setLat0(double I_Lat0, double I_LatC1, double I_LatC2);


    virtual void wheelEvent(QWheelEvent *event);

private:
    int    HRES;
    int    VRES;
    int    BRES;
    double SclLon;
    double SclLat;
    double DltLon;
    double DltLat;
    double LonMin;
    double LonMax;
    double LatMin;
    double LatMax;
    double ZoomScl;
    double ZoomSclMin;
    double ZoomSclMax;
    double ZoomSclInc;
    double Lon0;
    double Lat0;
    double LonC1;
    double LatC1;
    double LonC2;
    double LatC2;

signals:
    void setViewWindows(int HSIZE, int VSIZE, int LATSIZE, int LONSIZE,
                        double Lon0, double Lat0, double SclLon, double SclLat);
    void scrollViews(double Lon0, double Lat0);
};


#endif // VRMVIEW_H
