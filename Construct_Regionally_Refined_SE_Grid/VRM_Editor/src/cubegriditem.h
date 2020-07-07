#ifndef CUBEGRIDITEM_H
#define CUBEGRIDITEM_H

#include <QGraphicsItem>
#include <QPainter>
#include "CubeGrid.h"
#include "SQuadGen/DataMatrix.h"

class CubeGridItem : public QGraphicsItem
{
public:
    CubeGridItem();

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void loadGrid( const CubeGrid &I_CubeGrid);

    bool Initialized;

private:
    double LonMin;
    double LonMax;
    double LatMin;
    double LatMax;

    int                NumFaces;
    DataMatrix<double> FaceLon;
    DataMatrix<double> FaceLat;
};

#endif // CUBEGRIDITEM_H
