#include "refinecubeitem.h"


RefineCubeItem::RefineCubeItem()
{

}

RefineCubeItem::RefineCubeItem(double I_LonMin, double I_LatMin, double I_LonMax, double I_LatMax)
{
    LonMin = I_LonMin;
    LonMax = I_LonMax;
    LatMin = I_LatMin;
    LatMax = I_LatMax;
}

QRectF RefineCubeItem::boundingRect() const
{
    return QRectF(LonMin, LatMin, (LonMax-LonMin), (LatMax-LatMin));
}

void RefineCubeItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    double r1=1.0;
    double r2=0.0;
    double b1=1.0;
    double b2=0.0;
    double g1=1.0;
    double g2=1.0;
    double alpha = 0.4; //0.25;
    double red;
    double green;
    double blue;
    double rval;
    QColor MyColor;

    QPen MyPen = QPen(Qt::black);
    MyPen.setWidthF(0.);
//    painter->setPen(MyPen);
//    painter->setBrush(QBrush(QColor(0,255,0,128)));
    painter->setRenderHint(QPainter::Antialiasing,true);

    for( int iP=0; iP < 6             ; iP++) {
    for( int iA=0; iA < nMaxResolution; iA++) {
    for( int iB=0; iB < nMaxResolution; iB++) {
        double Lon1 = CubeLon[iP][iA  ][iB  ];
        double Lon2 = CubeLon[iP][iA+1][iB  ];
        double Lon3 = CubeLon[iP][iA+1][iB+1];
        double Lon4 = CubeLon[iP][iA  ][iB+1];
        double AvgLon = (Lon1 + Lon2 + Lon3 + Lon4)/4.0;
        double dLon1 = abs(Lon1 - AvgLon);
        double dLon2 = abs(Lon2 - AvgLon);
        double dLon3 = abs(Lon3 - AvgLon);
        double dLon4 = abs(Lon4 - AvgLon);
        double mLon1;
        double mLon2;
        double mLon3;
        double mLon4;

        if( (dLon1 > 135.)||(dLon2 > 135.)||(dLon3 > 135.)||(dLon4 >135.)) {
            Lon1 = (180./M_PI)*atan( tan ( Lon1 *M_PI/180.));
            Lon2 = (180./M_PI)*atan( tan ( Lon2 *M_PI/180.));
            Lon3 = (180./M_PI)*atan( tan ( Lon3 *M_PI/180.));
            Lon4 = (180./M_PI)*atan( tan ( Lon4 *M_PI/180.));
            mLon1 = Lon1 + 360.;
            mLon2 = Lon2 + 360.;
            mLon3 = Lon3 + 360.;
            mLon4 = Lon4 + 360.;
        } else if (AvgLon < 180.) {
            mLon1 = Lon1 +360.;
            mLon2 = Lon2 +360.;
            mLon3 = Lon3 +360.;
            mLon4 = Lon4 +360.;
        } else {
            mLon1 = Lon1 -360.;
            mLon2 = Lon2 -360.;
            mLon3 = Lon3 -360.;
            mLon4 = Lon4 -360.;
        }

        if(MaxVal > 0) {
            rval = static_cast<double>( val[iP][iA][iB]) / static_cast<double>(MaxVal);
        } else {
            rval = 0.0;
        }

        red   = r1 + rval*(r2-r1);
        green = g1 + rval*(g2-g1);
        blue  = b1 + rval*(b2-b1);
        MyColor.setRgbF(red,green,blue,alpha);
   //     painter->setPen(Qt::NoPen);
        painter->setPen(MyPen);
        painter->setBrush(QBrush(MyColor));
        QPolygonF poly1;
        QPolygonF poly2;
        poly1 << QPointF( Lon1,CubeLat[iP][iA  ][iB  ])
              << QPointF( Lon2,CubeLat[iP][iA+1][iB  ])
              << QPointF( Lon3,CubeLat[iP][iA+1][iB+1])
              << QPointF( Lon4,CubeLat[iP][iA  ][iB+1]) ;
        poly2 << QPointF(mLon1,CubeLat[iP][iA  ][iB  ])
              << QPointF(mLon2,CubeLat[iP][iA+1][iB  ])
              << QPointF(mLon3,CubeLat[iP][iA+1][iB+1])
              << QPointF(mLon4,CubeLat[iP][iA  ][iB+1]) ;
   //     painter->drawPolygon(poly1);
   //     painter->drawPolygon(poly2);
        painter->drawPolyline(poly1);
        painter->drawPolyline(poly2);
    }
    }
    }

}

