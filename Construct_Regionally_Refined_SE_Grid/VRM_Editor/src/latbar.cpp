#include "latbar.h"

LatBar::LatBar(double I_LatMin, double I_LatMax, int I_BarWidth)
{
    LatMin   = I_LatMin;
    LatMax   = I_LatMax;
    LatWidth = static_cast<double>(I_BarWidth);
}

QRectF LatBar::boundingRect() const
{
    return QRectF(0., LatMin, LatWidth, (LatMax-LatMin));
}

void LatBar::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    double hgt05   = 0.90*LatWidth;
    double hgt10   = 0.80*LatWidth;
    double hgt30   = 0.70*LatWidth;
    double hgt90   = 0.60*LatWidth;

    QPen   pen;
    pen.setWidthF(0.);

    for( double Lat=LatMin; Lat<LatMax; Lat += 5.) {
        pen.setColor(Qt::black);
        painter->setPen(pen);
        painter->drawLine(hgt05,Lat,LatWidth,Lat);
    }

    for( double Lat=LatMin; Lat<LatMax; Lat += 10.) {
        pen.setColor(Qt::black);
        painter->setPen(pen);
        painter->drawLine(hgt10,Lat,LatWidth,Lat);
    }

    for( double Lat=LatMin; Lat<LatMax; Lat += 30.) {
        pen.setColor(Qt::black);
        painter->setPen(pen);
        painter->drawLine(hgt30,Lat,LatWidth,Lat);
    }

    for( double Lat=LatMin; Lat<LatMax; Lat += 90.) {
        pen.setColor(Qt::black);
        painter->setPen(pen);
        painter->drawLine(hgt90,Lat,LatWidth,Lat);
    }
}

void LatBar::setLatBar(double I_LatMin, double I_LatMax, double I_LatWidth)
{
    LatMin   = I_LatMin;
    LatMax   = I_LatMax;
    LatWidth = I_LatWidth;
}
