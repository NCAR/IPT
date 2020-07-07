#ifndef POLYNODEITEM_H
#define POLYNODEITEM_H

#include <QGraphicsEllipseItem>
#include <QPen>
#include <QCollator>

class PolyNodeItem : public QObject, public QGraphicsEllipseItem
{
    Q_OBJECT

public:
    PolyNodeItem(QPointF I_Point, QGraphicsItem *parent=0);
    QPointF MyPoint;

public slots:
    void updateNode(QPointF I_point);
//    resetPosition();

signals:
    void do_UpdateNodePoly();

protected:
    QVariant itemChange(QGraphicsItem::GraphicsItemChange change, const QVariant &value);
//    virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);

private:

    double  size_x = 2.0;
    double  size_y = 2.0;
};

#endif // POLYNODEITEM_H
