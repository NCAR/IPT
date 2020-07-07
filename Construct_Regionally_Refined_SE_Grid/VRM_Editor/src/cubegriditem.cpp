#include "cubegriditem.h"

CubeGridItem::CubeGridItem()
{
    Initialized = false;
    LonMin = -180.;
    LonMax =  540.;
    LatMin = -90.;
    LatMax =  90.;
}

QRectF CubeGridItem::boundingRect() const
{
    return QRectF(LonMin, LatMin, (LonMax-LonMin), (LatMax-LatMin));
}

void CubeGridItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    double  Lat1,  Lat2,  Lat3,  Lat4;
    double  Lon1,  Lon2,  Lon3,  Lon4;
    double dLon1, dLon2, dLon3, dLon4;
    double mLon1, mLon2, mLon3, mLon4;
    double AvgLon;

    if(!Initialized) { return; }

    QPen MyPen = QPen(Qt::black);
    MyPen.setWidthF(0.);

    painter->setPen(MyPen);
    painter->setBrush(QBrush(Qt::NoBrush));
    painter->setRenderHint(QPainter::Antialiasing,true);

    for( int ii=0; ii < NumFaces; ii++) {
        Lat1 = FaceLat[ii][0];
        Lon1 = FaceLon[ii][0];
        Lat2 = FaceLat[ii][1];
        Lon2 = FaceLon[ii][1];
        Lat3 = FaceLat[ii][2];
        Lon3 = FaceLon[ii][2];
        Lat4 = FaceLat[ii][3];
        Lon4 = FaceLon[ii][3];

        AvgLon = (Lon1 + Lon2 + Lon3 + Lon4)/4.0;
        dLon1  = abs(Lon1 - AvgLon);
        dLon2  = abs(Lon2 - AvgLon);
        dLon3  = abs(Lon3 - AvgLon);
        dLon4  = abs(Lon4 - AvgLon);

        if((dLon1 > 135.)||(dLon2 > 135.)||(dLon3 > 135.)||(dLon4 >135.)) {
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

        QPolygonF poly1;
        QPolygonF poly2;
        poly1  << QPointF( Lon1,Lat1) << QPointF( Lon2,Lat2)
               << QPointF( Lon3,Lat3) << QPointF( Lon4,Lat4) << QPointF( Lon1,Lat1);
        poly2  << QPointF(mLon1,Lat1) << QPointF(mLon2,Lat2)
               << QPointF(mLon3,Lat3) << QPointF(mLon4,Lat4) << QPointF(mLon1,Lat1);
        painter->drawPolyline(poly1);
        painter->drawPolyline(poly2);
    }
}

void CubeGridItem::loadGrid(const CubeGrid &I_CubeGrid)
{
    NumFaces = I_CubeGrid.Faces.size();
    if(FaceLon.IsInitialized()) { FaceLon.Deinitialize();}
    if(FaceLat.IsInitialized()) { FaceLat.Deinitialize();}
    FaceLon.Initialize(NumFaces,4);
    FaceLat.Initialize(NumFaces,4);

    for( unsigned int ii=0; ii < I_CubeGrid.Faces.size(); ii++) {
        int ixNode0 = I_CubeGrid.Faces[ii].ixNode[0];
        int ixNode1 = I_CubeGrid.Faces[ii].ixNode[1];
        int ixNode2 = I_CubeGrid.Faces[ii].ixNode[2];
        int ixNode3 = I_CubeGrid.Faces[ii].ixNode[3];

        FaceLat[ii][0] = I_CubeGrid.Grid[ixNode0].lat;
        FaceLat[ii][1] = I_CubeGrid.Grid[ixNode1].lat;
        FaceLat[ii][2] = I_CubeGrid.Grid[ixNode2].lat;
        FaceLat[ii][3] = I_CubeGrid.Grid[ixNode3].lat;

        FaceLon[ii][0] = I_CubeGrid.Grid[ixNode0].lon;
        FaceLon[ii][1] = I_CubeGrid.Grid[ixNode1].lon;
        FaceLon[ii][2] = I_CubeGrid.Grid[ixNode2].lon;
        FaceLon[ii][3] = I_CubeGrid.Grid[ixNode3].lon;
    }

    Initialized = true;
}
