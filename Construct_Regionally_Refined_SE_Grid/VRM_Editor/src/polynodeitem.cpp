#include "polynodeitem.h"
#include <iostream>

PolyNodeItem::PolyNodeItem(QPointF I_Point, QGraphicsItem *parent)
{
    setFlag(QGraphicsItem::ItemIsMovable,true);
  //  setFlag(QGraphicsItem::ItemIsSelectable,true);
    setFlag(QGraphicsItem::ItemSendsGeometryChanges,true);

    MyPoint = I_Point;
    this->setRect(I_Point.x() - size_x/2.,I_Point.y() - size_y/2.,size_x,size_y);
    QColor MyColor;
    QPen   MyPen = QPen(Qt::black);
    MyPen.setWidthF(0.);
    this->setPen(MyPen);
    MyColor.setRgbF(0.7,0.7,0.7,0.4);
    QBrush MyBrush = QBrush(MyColor);
    this->setBrush(MyBrush);
    this->setZValue(0.75);
}

void PolyNodeItem::updateNode(QPointF I_point)
{
//    std::cout << " updateNode: "<< I_point.x() << " " << I_point.y() << std::endl;
    MyPoint = I_point;
    QPointF   MyPos  = this->pos();
    this->setRect(MyPoint.x() - size_x/2. - MyPos.x(),MyPoint.y() - size_y/2. - MyPos.y(),size_x,size_y);
}

QVariant PolyNodeItem::itemChange(QGraphicsItem::GraphicsItemChange change, const QVariant &value)
{
    if(change == ItemPositionHasChanged) {
      //  MyPoint = this->scenePos();
        QRectF MyRect;
        MyPoint = this->mapToScene(this->rect().center());
//        std::cout << " Node Changed: Point: " << MyPoint.x() << " " << MyPoint.y() << std::endl;
//        std::cout << "   Pos: " << this->pos().x() << " " << this->pos().y() << std::endl;
        emit do_UpdateNodePoly();

    }
    return QGraphicsItem::itemChange(change, value);
}
