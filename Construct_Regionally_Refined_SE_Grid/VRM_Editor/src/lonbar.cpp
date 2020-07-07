#include "lonbar.h"
#include <math.h>
#include <iostream>

LonBar::LonBar(double I_LonMin, double I_LonMax, int I_BarHeight)
{
    LonMin    = I_LonMin;
    LonMax    = I_LonMax;
    LonHeight = static_cast<double>(I_BarHeight);
}

QRectF LonBar::boundingRect() const
{
    return QRectF(LonMin, 0., (LonMax-LonMin), LonHeight);
}

void LonBar::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    double hgt05   = 0.10*LonHeight;
    double hgt10   = 0.20*LonHeight;
    double hgt30   = 0.30*LonHeight;
    double hgt90   = 0.40*LonHeight;
    double hgt180  = 0.50*LonHeight;

    QPen   pen;
    pen.setWidthF(0.);

    for( double Lon=LonMin; Lon<LonMax; Lon += 5.) {
        if((Lon >= 0.) && (Lon <= 360.)) {
            pen.setColor(Qt::black);
        } else {
            pen.setColor(Qt::red);
        }
        painter->setPen(pen);
        painter->drawLine(Lon,0.,Lon,hgt05);
    }

    for( double Lon=LonMin; Lon<LonMax; Lon += 10.) {
        if((Lon >= 0.) && (Lon <= 360.)) {
            pen.setColor(Qt::black);
        } else {
            pen.setColor(Qt::red);
        }
        painter->setPen(pen);
        painter->drawLine(Lon,0.,Lon,hgt10);
    }

    for( double Lon=LonMin; Lon<LonMax; Lon += 30.) {
        if((Lon >= 0.) && (Lon <= 360.)) {
            pen.setColor(Qt::black);
        } else {
            pen.setColor(Qt::red);
        }
        painter->setPen(pen);
        painter->drawLine(Lon,0.,Lon,hgt30);
    }

    for( double Lon=LonMin; Lon<LonMax; Lon += 90.) {
        if((Lon >= 0.) && (Lon <= 360.)) {
            pen.setColor(Qt::black);
        } else {
            pen.setColor(Qt::red);
        }
        painter->setPen(pen);
        painter->drawLine(Lon,0.,Lon,hgt90);
    }

    pen.setColor(Qt::black);
    painter->setPen(pen);
    painter->drawLine(  0.,0.,  0.,hgt180);
    painter->drawLine(180.,0.,180.,hgt180);
    painter->drawLine(360.,0.,360.,hgt180);
}

void LonBar::setLonBar(double I_LonMin, double I_LonMax, double I_LonHeight)
{
    LonMin    = I_LonMin;
    LonMax    = I_LonMax;
    LonHeight = I_LonHeight;
}

