#ifndef LONBAR_H
#define LONBAR_H

#include <QPainter>
#include <QGraphicsItem>
#include <QPen>
#include <QtGui>
#include <QKeyEvent>
#include <QtWidgets>
#include <QWidget>

class LonBar : public QGraphicsItem
{
public:
    LonBar(double I_LonMin, double I_LonMax, int I_BarHeight);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    void setLonBar(double I_LonMin, double I_LonMax, double I_LonHeight);

private:
    double LonMin;
    double LonMax;
    double LonHeight;

};

#endif // LONBAR_H
