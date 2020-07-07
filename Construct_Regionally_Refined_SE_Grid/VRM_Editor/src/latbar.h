#ifndef LATBAR_H
#define LATBAR_H

#include <QPainter>
#include <QGraphicsItem>
#include <QPen>

class LatBar : public QGraphicsItem
{
public:
    LatBar(double I_LatMin, double I_LatMax, int I_BarWidth);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void setLatBar(double I_LatMin, double I_LatMax, double I_LatWidth);

private:
    double LatMin;
    double LatMax;
    double LatWidth;

};

#endif // LATBAR_H
