#ifndef RECTEDITITEM_H
#define RECTEDITITEM_H

#include <QString>
#include <QtGui>
#include <QGraphicsRectItem>
#include <QGraphicsSceneMouseEvent>

class RectEditItem : public QObject, public QGraphicsRectItem
{
    Q_OBJECT

public:
    RectEditItem(QGraphicsItem *parent=0);
    void    ResetRect(double I_Lon0, double I_Lat0, double I_LonWidth, double I_LatWidth);
    double  Get_FillVal();
    void    Set_FillVal(double I_val);
    QString Get_FillOpt();
    void    Set_FillOpt(const QString & I_FillOpt);
    bool    Is_OutsideRectVisible();
    void    Set_OutsideRectVisible(bool I_val);
    QRectF  Get_OutsideRect();
    bool    Is_InsideRectVisible();
    void    Set_InsideRectVisible(bool I_val);
    QRectF  Get_InsideRect();
    QRectF  Get_Rect();
    void    Set_Rect(double I_LonMin, double I_LatMin, double I_LonWidth, double I_LatWidth);
    double  Get_DltLon();
    void    Set_DltLon(double I_val);
    double  Get_DltLat();
    void    Set_DltLat(double I_val);
    double  Get_LonMax();
    void    Set_LonMax(double I_val);
    double  Get_LonMin();
    void    Set_LonMin(double I_val);
    double  Get_LatMax();
    void    Set_LatMax(double I_val);
    double  Get_LatMin();
    void    Set_LatMin(double I_val);   
    double  Get_Lon0();
    double  Get_Lat0();
    double  Get_LonWidth();
    double  Get_LatWidth();

public slots:

signals:
    void do_SceneUpdate();
    void do_UIupdate();

protected:
    QVariant itemChange(QGraphicsItem::GraphicsItemChange change, const QVariant &value);
    virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);

private:

    QString MyFillOpt  = "Fill All";
    double  MyVal      = 1.;
    double  MyLon0     = 180.;
    double  MyLat0     = 0.;
    double  MyLonWidth = 30.;
    double  MyLatWidth = 15.;
    double  MyDltLon   = 0.;
    double  MyDltLat   = 0.;
    double  MyLonMin   = MyLon0 - (MyLonWidth/2.);
    double  MyLatMin   = MyLat0 - (MyLatWidth/2.);
    double  MyLonMax   = MyLon0 + (MyLonWidth/2.);
    double  MyLatMax   = MyLat0 + (MyLatWidth/2.);

    bool    MyInsideRectVisible  = false;
    bool    MyOutsideRectVisible = false;

};

#endif // RECTEDITITEM_H
