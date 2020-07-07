#ifndef REFINECUBEITEM_H
#define REFINECUBEITEM_H
#include <QGraphicsItem>
#include <QPainter>
#include "SQuadGen/DataMatrix3D.h"
#include "RefinementCube.h"


class RefineCubeItem : public QGraphicsItem
{
public:
    RefineCubeItem();
    RefineCubeItem(double I_LonMin, double I_LatMin, double I_LonMax, double I_LatMax);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

private:
    double LonMin;
    double LonMax;
    double LatMin;
    double LatMax;

public:
    int                  nMaxResolution;
    int                  MaxVal;
    DataMatrix3D<int>    val;
    DataMatrix3D<double> CubeLon;
    DataMatrix3D<double> CubeLat;

};

#endif // REFINECUBEITEM_H
