#include "polyedititem.h"
#include <iostream>

PolyEditItem::PolyEditItem(QGraphicsItem *parent): QGraphicsPolygonItem(parent)
{
    setFlag(QGraphicsItem::ItemIsMovable,true);
  //  setFlag(QGraphicsItem::ItemIsSelectable,true);
    setFlag(QGraphicsItem::ItemSendsGeometryChanges,true);
}

double PolyEditItem::Val()
{
    return MyVal;
}

void PolyEditItem::setVal(double I_val)
{
    MyVal = I_val;
}

void PolyEditItem::UpdatePolygon(QPolygonF I_poly)
{
    QPolygonF MyPoly = I_poly;
    QPointF   MyPos  = this->pos();

    for(int ii=0; ii < MyPoly.size(); ii++) {
       MyPoly[ii].setX(MyPoly[ii].x() - MyPos.x());
       MyPoly[ii].setY(MyPoly[ii].y() - MyPos.y());
    }
    this->setPolygon(MyPoly);
}

//void PolyEditItem::resetPosition()
//{
//    QPolygonF MyPoly = this->mapToScene(this->polygon());
//    this->setPos(0.,0.);
//    this->setPolygon(MyPoly);
//}


QVariant PolyEditItem::itemChange(QGraphicsItem::GraphicsItemChange change, const QVariant &value)
{
    if(change == ItemPositionHasChanged) {
        QPolygonF MyPoly;
        MyPoly = this->mapToScene(this->polygon());

 //       std::cout << " emitting do_UpdateNodes " << std::endl;
//        for(int ii=0; ii < MyPoly.size(); ii++) {
 //           std::cout << " Poly: =" << MyPoly[ii].x()  <<" " << MyPoly[ii].y() << std::endl;
 //       }
 //       std::cout << " " << std::endl;



        emit do_UpdateNodes(MyPoly);

 //       QPolygonF MyPoly = this->polygon();
 //       QPointF   MyPos  = this->pos();

 //       QPoint NewPoint;
 //       QPolygonF MyPolyS = this->mapToScene(this->polygon());
 //       for(int ii=0; ii < MyPoly.size(); ii++) {
 //          std::cout << " Poly: =" << MyPoly[ii].x()  <<" " << MyPoly[ii].y() <<" : "
 //                                  << MyPolyS[ii].x() <<" " << MyPolyS[ii].y()<< std::endl;
 //          MyPolyS[ii].setX(MyPolyS[ii].x() - MyPos.x());
 //          MyPolyS[ii].setY(MyPolyS[ii].y() - MyPos.y());
 //       }

    }
    return QGraphicsItem::itemChange(change, value);
}

void PolyEditItem::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if(event->button() == Qt::RightButton) {
        QPolygonF *MyPoly;
        MyPoly = new QPolygonF(this->mapToScene(this->polygon()));
        emit do_SetPolyVal( MyPoly, MyVal);
        emit do_SceneUpdate();
        event->accept();

    //    resetPosition();
/*
        QPolygonF MyPolyS = this->mapToScene(this->polygon());
        QPointF   MyPos = this->pos();
        for(int ii=0; ii < MyPolyS.size(); ii++) {
           MyPolyS[ii].setX(MyPolyS[ii].x() - MyPos.x());
           MyPolyS[ii].setY(MyPolyS[ii].y() - MyPos.y());
        }
        MyPolyS[1].setX(MyPolyS[1].x() + 5.);
        this->setPolygon(MyPolyS);
*/
    }
}
