#ifndef VRMWIN_H
#define VRMWIN_H
#include "vrmgrid.h"
#include "refinecubeitem.h"
#include "cubegriditem.h"
#include "polyedititem.h"
#include "polynodeitem.h"
#include "rectedititem.h"
#include "gridinfoitem.h"

#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QString>
#include <QImage>


class VrmWin : public QGraphicsScene
{
    Q_OBJECT

public:
    VrmWin(double I_LonMin, double I_LatMin, double I_LonMax, double I_LatMax);

    QRectF boundingRect() const;
    void initialize();
    void update_scene(VrmGrid *vrmGrid);

    double BackgroundFilterAlpha();
    void   setBackgroundFilterAlpha( double I_Alpha);
    QString BackgroundImgVal();
    void    setBackgroundImgVal(QString I_Image);

    bool    ReferenceMapVisible();
    void    setReferenceMapVisible(bool I_Visible);
    double  ReferenceMapAlpha();
    void    setReferenceMapAlpha( double I_Alpha);

    bool    RefineMapEditVisible();
    void    setRefineMapEditVisible(bool I_Visible);
    bool    RefineMapVisible();
    void    setRefineMapVisible(bool I_Visible);
    double  RefineMapAlpha();
    void    setRefineMapAlpha( double I_Alpha);

    bool    RefineCubeVisible();
    void    setRefineCubeVisible(bool I_Visible);
    double  RefineCubeAlpha();
    void    setRefineCubeAlpha( double I_Alpha);
    void    updateRefineCubeImg( VrmGrid *I_vrmGrid);
    void    adjustRefineCubeImgAlpha(double I_Alpha);

    bool    DisplayGridVisible();
    void    setDisplayGridVisible(bool I_Visible);
    QString DisplayGridName();
    void    setDisplayGridName( QString I_GridName);

    void  setSmoothKernel     (int    I_val);
    void  setSmoothPersistence(double I_val);
    void  setSmoothIterNum    (int    I_val);
    void  applySmoothing( VrmGrid *I_vrmGrid);

signals:
    void do_gridInfoUpdate();

public slots:
    void UpdateNodes(QPolygonF MyPoly);
    void UpdateNodePoly();

private:
    double LonMin;
    double LonMax;
    double LatMin;
    double LatMax;

    int PYSIZE = 1*1024; //2048;
    int PXSIZE = 4*PYSIZE;
//    int GYSIZE = 1*1024;
//    int GXSIZE = 4*GYSIZE;

    double   backgroundFilter_Alpha;
    QString  backgroundImg;
    QImage  *backgroundImg_Terrain;
    QImage  *backgroundImg_LandWater;
    QImage  *backgroundImg_None;
    double   backgroundImg_TerrainScl;
    double   backgroundImg_LandWaterScl;
    double   backgroundImg_NoneScl;

    bool     refineMapRefVisible;
    double   refineMapRefAlpha;
    bool     refineMapEditVisible;
    bool     refineMapVisible;
    double   refineMapAlpha;



    bool     refineCubeVisible;
    double   refineCubeAlpha;
    QImage  *refineCubeImg;
    bool     refineCubeImgUpdate;

    bool     displayGridVisible;
    QString  displayGridName;
    bool     baseGridVisible;
    bool     refineGridVisible;

    int      smoothKernel      = 3;
    double   smoothPersistence = 1.0;
    int      smoothIterNum     = 1;


public:
    QGraphicsScene *Scene;

    PolyEditItem   *polyEdit_Item;
    PolyNodeItem   *Node1;
    PolyNodeItem   *Node2;
    PolyNodeItem   *Node3;
    PolyNodeItem   *Node4;
    PolyNodeItem   *Node5;
    PolyNodeItem   *Node6;
    PolyNodeItem   *Node7;
    PolyNodeItem   *Node8;

    RectEditItem      *rectEdit_Item;
    QGraphicsRectItem *rectEdit_Inside;
    QGraphicsRectItem *rectEdit_Outside;

private:
    QPixmap             *backgroundImage_Pix;
    QGraphicsPixmapItem *backgroundImage_Item;

    QPixmap             *backgroundFilter_Pix;
    QGraphicsPixmapItem *backgroundFilter_Item;

    QGraphicsPixmapItem *refineMap_Item;
    QGraphicsPixmapItem *refineMapRef_Item;
    QGraphicsPixmapItem *refineMapEdit_Item;

    QPixmap             *refineCube_Pix;
    QGraphicsPixmapItem *refineCube_Item;

    CubeGridItem        *baseGrid_Item;
    CubeGridItem        *refineGrid_Item;


};

#endif // VRMWIN_H
