#include "rectedititem.h"
#include <iostream>

RectEditItem::RectEditItem(QGraphicsItem *parent) : QGraphicsRectItem(parent)
{
 //   setFlag(QGraphicsItem::ItemIsMovable,true);
    setFlag(QGraphicsItem::ItemSendsGeometryChanges,true);
//    setFlag(QGraphicsItem::ItemIsSelectable,true);
    setRect(MyLonMin,MyLatMin,MyLonWidth,MyLatWidth);
}

void RectEditItem::ResetRect(double I_Lon0, double I_Lat0, double I_LonWidth, double I_LatWidth)
{
    MyLon0     = I_Lon0;
    MyLat0     = I_Lat0;
    MyLonWidth = I_LonWidth;
    MyLatWidth = I_LatWidth;
    MyDltLon   = 0.;
    MyDltLat   = 0.;
    MyLonMin   = MyLon0 - (MyLonWidth/2.);
    MyLatMin   = MyLat0 - (MyLatWidth/2.);
    MyLonMax   = MyLon0 + (MyLonWidth/2.);
    MyLatMax   = MyLat0 + (MyLatWidth/2.);
    MyInsideRectVisible  = false;
    MyOutsideRectVisible = false;
}

QString RectEditItem::Get_FillOpt()
{
    return MyFillOpt;
}

void RectEditItem::Set_FillOpt(const QString &I_FillOpt)
{
    MyFillOpt = I_FillOpt;
}

double RectEditItem::Get_FillVal()
{
    return MyVal;
}

void RectEditItem::Set_FillVal(double I_val)
{
    MyVal = I_val;
}

bool RectEditItem::Is_OutsideRectVisible()
{
    return MyOutsideRectVisible;
}

void RectEditItem::Set_OutsideRectVisible(bool I_val)
{
    MyOutsideRectVisible = I_val;
}

QRectF RectEditItem::Get_OutsideRect()
{
    double MinLon = MyLonMin - MyDltLon/2.;
    double MaxLon = MyLonMax + MyDltLon/2.;
    double MinLat = MyLatMin - MyDltLat/2.;
    double MaxLat = MyLatMax + MyDltLat/2.;
    return QRectF(MinLon,MinLat,(MaxLon-MinLon),(MaxLat-MinLat));
}

bool RectEditItem::Is_InsideRectVisible()
{
    return MyInsideRectVisible;
}

void RectEditItem::Set_InsideRectVisible(bool I_val)
{
    MyInsideRectVisible = I_val;
}

QRectF RectEditItem::Get_InsideRect()
{
    double MinLon = MyLonMin + MyDltLon/2.;
    double MaxLon = MyLonMax - MyDltLon/2.;
    double MinLat = MyLatMin + MyDltLat/2.;
    double MaxLat = MyLatMax - MyDltLat/2.;
    return QRectF(MinLon,MinLat,(MaxLon-MinLon),(MaxLat-MinLat));
}

QRectF RectEditItem::Get_Rect()
{
    return QRectF(MyLonMin,MyLatMin,MyLonWidth,MyLatWidth);
}

void RectEditItem::Set_Rect(double I_LonMin, double I_LatMin, double I_LonWidth, double I_LatWidth)
{
    MyLonMin = I_LonMin;
    MyLatMin = I_LatMin;
    MyLonWidth = I_LonWidth;
    MyLatWidth = I_LatWidth;

    MyLonMax = MyLonMin + MyLonWidth;
    MyLatMax = MyLatMin + MyLatWidth;
    MyLon0   = (MyLonMin + MyLonMax)/2.;
    MyLat0   = (MyLatMin + MyLatMax)/2.;

    // Check if the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }

}

double RectEditItem::Get_DltLon()
{
    return MyDltLon;
}

void RectEditItem::Set_DltLon(double I_val)
{
    MyDltLon = I_val;

    // Check if the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_DltLat()
{
    return MyDltLat;
}

void RectEditItem::Set_DltLat(double I_val)
{
    MyDltLat = I_val;

    // Check of the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_LonMax()
{
    return MyLonMax;
}

void RectEditItem::Set_LonMax(double I_val)
{
    if( I_val > MyLonMin + 1.) {
        MyLonMax = I_val;
    } else {
        MyLonMax = MyLonMin + 1.;
    }
    MyLon0     = (MyLonMax + MyLonMin)/2.;
    MyLonWidth =  MyLonMax - MyLonMin;

    // Check of the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_LonMin()
{
    return MyLonMin;
}

void RectEditItem::Set_LonMin(double I_val)
{
    if( I_val < MyLonMax - 1.) {
        MyLonMin = I_val;
    } else {
        MyLonMin = MyLonMax - 1.;
    }
    MyLon0     = (MyLonMax + MyLonMin)/2.;
    MyLonWidth =  MyLonMax - MyLonMin;

    // Check of the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_LatMax()
{
    return MyLatMax;
}

void RectEditItem::Set_LatMax(double I_val)
{
    if( I_val > MyLatMin + 1.) {
        MyLatMax = I_val;
    } else {
        MyLatMax = MyLatMin + 1.;
    }
    MyLat0     = (MyLatMax + MyLatMin)/2.;
    MyLatWidth =  MyLatMax - MyLatMin;

    // Check of the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_LatMin()
{
    return MyLatMin;
}

void RectEditItem::Set_LatMin(double I_val)
{
    if( I_val < MyLatMax - 1.) {
        MyLatMin = I_val;
    } else {
        MyLatMin = MyLatMax - 1.;
    }
    MyLat0     = (MyLatMax + MyLatMin)/2.;
    MyLatWidth =  MyLatMax - MyLatMin;

    // Check of the inside rect needs to be switched on/off
    //------------------------------------------------------------
    if((MyDltLon > MyLonWidth/2.) || (MyDltLat > MyLatWidth/2.)) {
        MyInsideRectVisible = false;
    } else if((MyDltLon ==0.) && (MyDltLat ==0.)) {
        MyInsideRectVisible  = false;
        MyOutsideRectVisible = false;
    } else {
        MyInsideRectVisible  = true;
        MyOutsideRectVisible = true;
    }
}

double RectEditItem::Get_Lon0()
{
    return MyLon0;
}

double RectEditItem::Get_Lat0()
{
    return MyLat0;
}

double RectEditItem::Get_LonWidth()
{
    return MyLonWidth;
}

double RectEditItem::Get_LatWidth()
{
    return MyLatWidth;
}

QVariant RectEditItem::itemChange(QGraphicsItem::GraphicsItemChange change, const QVariant &value)
{
    if(change == ItemPositionHasChanged) {
        QRectF MyRect;
        MyRect = this->rect();
        QPolygonF MyPoly = this->mapToScene(MyRect);
        for(int ii=0; ii<MyPoly.size(); ii++) {
//            std::cout << " Poly: "<< MyPoly[ii].x() << " " << MyPoly[ii].y() << std::endl;
        }
//        std::cout << " Change: L:"<< MyRect.left() << " B:"<< MyRect.bottom()<< " W:"<<MyRect.width()<< " H:"<<MyRect.height() << std::endl;

        this->Set_Rect(MyPoly[0].x(),MyPoly[0].y(),(MyPoly[1].x()-MyPoly[0].x()),(MyPoly[2].y()-MyPoly[1].y()));
        MyRect = this->rect();
//        std::cout << " Change: L:"<< MyRect.left() << " B:"<< MyRect.bottom()<< " W:"<<MyRect.width()<< " H:"<<MyRect.height() << std::endl;
    //    emit do_UIupdate();
        emit do_SceneUpdate();
    }
    return QGraphicsItem::itemChange(change,value);
}

void RectEditItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  //  if(event->
   //     emit do_UIupdate();
  //      emit do_SceneUpdate();
  //      event->accept();
  //  }
}
