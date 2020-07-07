#include "vrmwin.h"
#include "vrmgrid.h"
#include <QGraphicsRectItem>
#include <QImage>
#include <iostream>
#include <QPixmap>
#include <QBrush>
#include <QPen>

VrmWin::VrmWin(double I_LonMin, double I_LatMin, double I_LonMax, double I_LatMax)
{
    LonMin = I_LonMin;
    LonMax = I_LonMax;
    LatMin = I_LatMin;
    LatMax = I_LatMax;
    Scene  = new QGraphicsScene(this);
    Scene->setSceneRect(LonMin,LatMin,(LonMax-LonMin),(LatMax-LatMin));
}

QRectF VrmWin::boundingRect() const
{
    return QRectF(LonMin, LatMin, (LonMax-LonMin), (LatMax-LatMin));
}

void VrmWin::initialize()
{
    QColor MyColor;
    QPen   MyPen;
    int    ii_size;
    int    jj_size;

    // Read in Land/Water map image
    //-------------------------------
    QImage landwater_in = QImage(":/images/earthspec2k.jpg").mirrored(false,true);
    ii_size   = landwater_in.width();
    jj_size   = landwater_in.height();
    backgroundImg_LandWater = new QImage(2*ii_size, jj_size, landwater_in.format());
    for( int ii=0; ii < ii_size; ii++) {
    for( int jj=0; jj < jj_size; jj++) {
        backgroundImg_LandWater->setPixel(ii          ,jj,landwater_in.pixel(ii,jj));
        backgroundImg_LandWater->setPixel(ii + ii_size,jj,landwater_in.pixel(ii,jj));
    }
    }
    backgroundImg_LandWaterScl = 360./static_cast<double>(ii_size);

    // Read in Terrain map image
    //----------------------------
    QImage terrain_in   = QImage( ":/images/earthmap2k.jpg").mirrored(false,true);
    ii_size   = terrain_in.width();
    jj_size   = terrain_in.height();
    backgroundImg_Terrain = new QImage(2*ii_size, jj_size, terrain_in.format());
    for( int ii=0; ii < ii_size; ii++) {
    for( int jj=0; jj < jj_size; jj++) {
        backgroundImg_Terrain->setPixel(ii          ,jj,terrain_in.pixel(ii,jj));
        backgroundImg_Terrain->setPixel(ii + ii_size,jj,terrain_in.pixel(ii,jj));
    }
    }
    backgroundImg_TerrainScl = 360./static_cast<double>(ii_size);

    // Set NULL image based on Terrain size,etc..
    //----------------------------------------------
    backgroundImg_None = new QImage(2*ii_size, jj_size, terrain_in.format());
    MyColor.setRgbF(0.1,0.1,0.5,1.);
    backgroundImg_None->fill(MyColor);
    backgroundImg_NoneScl = 360./static_cast<double>(ii_size);

    // Initialize Background Image Item
    //-----------------------------------
    backgroundImage_Pix = new QPixmap(2*terrain_in.width(), terrain_in.height());
    backgroundImage_Pix->convertFromImage(*backgroundImg_Terrain);
    backgroundImage_Item = new QGraphicsPixmapItem();
    backgroundImage_Item->setPixmap(*backgroundImage_Pix);
    backgroundImage_Item->setPos(-180.,-90.);
    backgroundImage_Item->setScale(backgroundImg_TerrainScl);
    backgroundImage_Item->setZValue(-3.);
    Scene->addItem(backgroundImage_Item);

    // Initialize Background Filter Item
    //------------------------------------
    backgroundFilter_Alpha = 0.2;
    MyColor.setRgbF(1.,1.,1.,backgroundFilter_Alpha);
    backgroundFilter_Pix   = new QPixmap(PXSIZE,PYSIZE);
    backgroundFilter_Pix->fill(MyColor);
    backgroundFilter_Item = new QGraphicsPixmapItem();
    backgroundFilter_Item->setPixmap( *backgroundFilter_Pix);
    backgroundFilter_Item->setPos(-180.,-90.);
    backgroundFilter_Item->setScale(180./PYSIZE);
    backgroundFilter_Item->setZValue(-2.);
    Scene->addItem(backgroundFilter_Item);

    // Initialize Reference Map Item
    //----------------------------------
    refineMapRefVisible = false;
    refineMapRefAlpha   = 0.5;
    refineMapRef_Item = new QGraphicsPixmapItem();
    refineMapRef_Item->setVisible(refineMapRefVisible);
    refineMapRef_Item->setZValue(-1.);
    Scene->addItem(refineMapRef_Item);

    // Initialize Refine Map Item
    //----------------------------
    refineMapVisible = false;
    refineMapAlpha   = 0.5;
    refineMap_Item    = new QGraphicsPixmapItem();
    refineMap_Item->setVisible(refineMapVisible);
    refineMap_Item->setZValue(0.);
    Scene->addItem(refineMap_Item);

    // Initialize Refine Map Edit Item
    //---------------------------------
    refineMapEditVisible = false;
    refineMapEdit_Item    = new QGraphicsPixmapItem();
    refineMapEdit_Item->setVisible(refineMapEditVisible);
    refineMapEdit_Item->setZValue(0.5);
    Scene->addItem(refineMapEdit_Item);

    rectEdit_Item    = new RectEditItem();
    rectEdit_Inside  = new QGraphicsRectItem();
    rectEdit_Inside->setRect(rectEdit_Item->Get_InsideRect());
    rectEdit_Outside = new QGraphicsRectItem();
    rectEdit_Outside->setRect(rectEdit_Item->Get_OutsideRect());

    MyPen = QPen(Qt::black);
    MyPen.setWidthF(0.);
    rectEdit_Item->setPen(MyPen);
    rectEdit_Item->setBrush(Qt::NoBrush);
    MyPen = QPen(Qt::green);
    MyPen.setWidthF(0.);
    rectEdit_Inside->setPen(MyPen);
    rectEdit_Inside->setBrush(Qt::NoBrush);
    rectEdit_Outside->setPen(MyPen);
    rectEdit_Outside->setBrush(Qt::NoBrush);
    rectEdit_Item->setZValue(0.6);
    rectEdit_Item->setVisible(false);
    Scene->addItem(rectEdit_Item);
    rectEdit_Inside->setZValue(0.575);
    rectEdit_Inside->setVisible(false);
    Scene->addItem(rectEdit_Inside);
    rectEdit_Outside->setZValue(0.575);
    rectEdit_Outside->setVisible(false);
    Scene->addItem(rectEdit_Outside);


    double x0 = (LonMax + LonMin)/2.;
    double y0 = (LatMin + LatMax)/2.;
    double dx = 0.05*(LonMax - LonMin);
    double dy = 0.05*(LatMax - LatMin);
    QPointF MyPoint1(x0-dx,y0-dy);
    QPointF MyPoint2(x0   ,y0-dy);
    QPointF MyPoint3(x0+dx,y0-dy);
    QPointF MyPoint4(x0+dx,y0   );
    QPointF MyPoint5(x0+dx,y0+dy);
    QPointF MyPoint6(x0   ,y0+dy);
    QPointF MyPoint7(x0-dx,y0+dy);
    QPointF MyPoint8(x0-dx,y0   );
    QPolygonF MyPoly;
    MyPoly << MyPoint1 << MyPoint2 << MyPoint3 << MyPoint4
           << MyPoint5 << MyPoint6 << MyPoint7 << MyPoint8 << MyPoint1;

    Node1 = new PolyNodeItem(MyPoint1);
    Node2 = new PolyNodeItem(MyPoint2);
    Node3 = new PolyNodeItem(MyPoint3);
    Node4 = new PolyNodeItem(MyPoint4);
    Node5 = new PolyNodeItem(MyPoint5);
    Node6 = new PolyNodeItem(MyPoint6);
    Node7 = new PolyNodeItem(MyPoint7);
    Node8 = new PolyNodeItem(MyPoint8);
    Node1->setVisible(false);
    Node2->setVisible(false);
    Node3->setVisible(false);
    Node4->setVisible(false);
    Node5->setVisible(false);
    Node6->setVisible(false);
    Node7->setVisible(false);
    Node8->setVisible(false);
    Scene->addItem(Node1);
    Scene->addItem(Node2);
    Scene->addItem(Node3);
    Scene->addItem(Node4);
    Scene->addItem(Node5);
    Scene->addItem(Node6);
    Scene->addItem(Node7);
    Scene->addItem(Node8);

    polyEdit_Item = new PolyEditItem();
    polyEdit_Item->setPolygon(MyPoly);
    MyPen = QPen(Qt::green);
    MyPen.setWidthF(0.);
    polyEdit_Item->setPen(MyPen);
    MyColor.setRgbF(0.,1.,0.,0.1);
    QBrush MyBrush = QBrush(MyColor);
    polyEdit_Item->setBrush(MyBrush);
    polyEdit_Item->setVal(1.0);
    polyEdit_Item->setZValue(0.7);
    polyEdit_Item->setVisible(false);
    Scene->addItem(polyEdit_Item);

    // Initialize Refine Cube Item
    //----------------------------
    refineCubeVisible   = false;
    refineCubeAlpha     = 0.5;
    refineCubeImgUpdate = false;
    refineCubeImg   = new QImage(PXSIZE,PYSIZE, QImage::Format_ARGB32);
    MyColor.setRgbF(1.,1.,1.,0.0);
    refineCubeImg->fill(MyColor);
    refineCube_Pix  = new QPixmap(PXSIZE,PYSIZE);
    refineCube_Pix->convertFromImage( *refineCubeImg);
    refineCube_Item = new QGraphicsPixmapItem();
    refineCube_Item->setPixmap( *refineCube_Pix);
    refineCube_Item->setPos(-180.,-90.);
    refineCube_Item->setScale(180./PYSIZE);
    refineCube_Item->setVisible(refineCubeVisible);
    refineCube_Item->setZValue(1.);
    Scene->addItem(refineCube_Item);

    // Initialize Base/Refine Grid Items
    //------------------------------------
    displayGridVisible = true;
    displayGridName = "Variable Resolution";

    baseGridVisible = false;
    baseGrid_Item = new CubeGridItem();
    baseGrid_Item->setVisible(baseGridVisible);
    baseGrid_Item->setZValue(2.);
    Scene->addItem(baseGrid_Item);

    refineGridVisible = true;
    refineGrid_Item = new CubeGridItem();
    refineGrid_Item->setVisible(refineGridVisible);
    refineGrid_Item->setZValue(3.);
    Scene->addItem(refineGrid_Item);

    // Initialize Lat-Lon Grid Item
    //------------------------------
//  latlonGrid_Item = new LatLonGridItem();
//  latlonGrid_ItemX->setVisible(false);
//  latlonGrid_ItemX->setZValue(4.);
//  Scene->addItem(latlonGrid_ItemX);




}

void VrmWin::update_scene(VrmGrid *vrmGrid)
{
    // Branch to a VRM scene or a edit refineMap scene
    //===================================================






    // If selected/active add the reference refineMap
    //------------------------------------------------
    if(refineMapRefVisible) {
        int ii_size = vrmGrid->refineMapRef.nRefLon;
        int jj_size = vrmGrid->refineMapRef.nRefLat;
        int ii_offset = ii_size/2;

        double r1 = 1.0;
        double g1 = 1.0;
        double b1 = 0.9;

        double r2 = 1.0;
        double g2 = 1.0;
        double b2 = 0.0;

        double alpha;
        double red;
        double green;
        double blue;
        double rval;
        QColor MyColor;

        // create a Qimage to work with and set pixels from refineMap values
        //------------------------------------------------------------------
        QImage MyImg(2*ii_size, jj_size, QImage::Format_ARGB32);
        for( int ii=0; ii < ii_size; ii++) {
            for( int jj=0; jj < jj_size; jj++) {
                rval  = (vrmGrid->refineMapRef.val[ii][jj] - vrmGrid->refineMapRef.Min)
                       /(vrmGrid->refineMapRef.Max         - vrmGrid->refineMapRef.Min);
                if(rval > 0.01) {
                    alpha = refineMapRefAlpha;
                } else {
                    alpha = 0.0;
                }
                red   = r1 + rval*(r2-r1);
                green = g1 + rval*(g2-g1);
                blue  = b1 + rval*(b2-b1);
                MyColor.setRgbF(red,green,blue,alpha);
                MyImg.setPixelColor(ii + ii_offset,jj, MyColor);
            }
            if( ii < ii_offset) {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = (vrmGrid->refineMapRef.val[ii][jj] - vrmGrid->refineMapRef.Min)
                           /(vrmGrid->refineMapRef.Max         - vrmGrid->refineMapRef.Min);
                    if(rval > 0.01) {
                        alpha = refineMapRefAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha);
                    MyImg.setPixelColor(ii + ii_size + ii_offset,jj, MyColor);
                }
            } else {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = (vrmGrid->refineMapRef.val[ii][jj] - vrmGrid->refineMapRef.Min)
                           /(vrmGrid->refineMapRef.Max         - vrmGrid->refineMapRef.Min);
                    if(rval > 0.01) {
                        alpha = refineMapRefAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha);
                    MyImg.setPixelColor(ii - ii_offset,jj, MyColor);
                }
            }
        }

        // Update refineMap_Item and make sure it is visible
        //--------------------------------------------------
        QPixmap MyPmapR;
        double dltLat = (180./(jj_size));
        MyPmapR.convertFromImage(MyImg);
        refineMapRef_Item->setPixmap(MyPmapR);
        refineMapRef_Item->setPos(-180.,-90. );
        refineMapRef_Item->setScale(dltLat);
        refineMapRef_Item->setVisible(true);
        refineMapRef_Item->update();
    } else {
        QPixmap MyPmap;
        MyPmap.fill(QColor(0,0,0,0));
        refineMapRef_Item->setPixmap(MyPmap);
        refineMapRef_Item->setVisible(false);
    }

    // If selected/active add the current refineMap to the Scene
    //------------------------------------------------------------
    if(refineMapVisible) {
        int ii_size = vrmGrid->refineMap.nRefLon;
        int jj_size = vrmGrid->refineMap.nRefLat;
        int ii_offset = ii_size/2;

        double r1 = 1.0;
        double g1 = 0.9;  // DIAG 1.0;
        double b1 = 0.9 ; // DIAG 1.0;

        double r2 = 1.0;
        double g2 = 0.0;
        double b2 = 0.0;

        double alpha;
        double red;
        double green;
        double blue;
        double rval;
        QColor MyColor;

        // create a Qimage to work with and set pixels from refineMap values
        //------------------------------------------------------------------
        QImage MyImg(2*ii_size, jj_size, QImage::Format_ARGB32);
        for( int ii=0; ii < ii_size; ii++) {
            for( int jj=0; jj < jj_size; jj++) {
                rval  = vrmGrid->refineMap.val[ii][jj];
                if(rval > 0.01) {
                    alpha = refineMapAlpha;
                } else {
                    alpha = 0.0;
                }
                red   = r1 + rval*(r2-r1);
                green = g1 + rval*(g2-g1);
                blue  = b1 + rval*(b2-b1);
                MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                MyImg.setPixelColor(ii + ii_offset,jj, MyColor);
            }
            if( ii < ii_offset) {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = vrmGrid->refineMap.val[ii][jj];
                    if(rval > 0.01) {
                        alpha = refineMapAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                    MyImg.setPixelColor(ii + ii_size + ii_offset,jj, MyColor);
                }
            } else {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = vrmGrid->refineMap.val[ii][jj];
                    if(rval > 0.01) {
                        alpha = refineMapAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                    MyImg.setPixelColor(ii - ii_offset,jj, MyColor);
                }
            }
        }

        // Update refineMap_Item and make sure it is visible
        //--------------------------------------------------
        QPixmap MyPmap;
        double dltLat = (180./(jj_size));
        MyPmap.convertFromImage(MyImg);
        refineMap_Item->setPixmap(MyPmap);
        refineMap_Item->setPos(-180.,-90. );
        refineMap_Item->setScale(dltLat);
        refineMap_Item->setVisible(true);
        refineMap_Item->update();
    } else {
        QPixmap MyPmap;
        MyPmap.fill(QColor(0,0,0,0));
        refineMap_Item->setPixmap(MyPmap);
        refineMap_Item->setVisible(false);
    }

    // If selected/active add the current refineMapEdit to the Scene
    // the Edit copy of the refineMap uses the display settings
    // (color/alpha) for refineMap since it is deactivated in edit mode.
    //------------------------------------------------------------
    if(refineMapEditVisible) {
        int ii_size = vrmGrid->refineMapEdit.nRefLon;
        int jj_size = vrmGrid->refineMapEdit.nRefLat;
        int ii_offset = ii_size/2;

        double r1 = 1.0;
        double g1 = 0.9;  // DIAG 1.0;
        double b1 = 0.9 ; // DIAG 1.0;

        double r2 = 1.0;
        double g2 = 0.0;
        double b2 = 0.0;

        double alpha;
        double red;
        double green;
        double blue;
        double rval;
        QColor MyColor;

        // create a Qimage to work with and set pixels from refineMap values
        //------------------------------------------------------------------
        QImage MyImg(2*ii_size, jj_size, QImage::Format_ARGB32);
        for( int ii=0; ii < ii_size; ii++) {
            for( int jj=0; jj < jj_size; jj++) {
                rval  = vrmGrid->refineMapEdit.val[ii][jj];
                if(rval > 0.01) {
                    alpha = refineMapAlpha;
                } else {
                    alpha = 0.0;
                }
                red   = r1 + rval*(r2-r1);
                green = g1 + rval*(g2-g1);
                blue  = b1 + rval*(b2-b1);
                MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                MyImg.setPixelColor(ii + ii_offset,jj, MyColor);
            }
            if( ii < ii_offset) {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = vrmGrid->refineMapEdit.val[ii][jj];
                    if(rval > 0.01) {
                        alpha = refineMapAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                    MyImg.setPixelColor(ii + ii_size + ii_offset,jj, MyColor);
                }
            } else {
                for( int jj=0; jj < jj_size; jj++) {
                    rval  = vrmGrid->refineMapEdit.val[ii][jj];
                    if(rval > 0.01) {
                        alpha = refineMapAlpha;
                    } else {
                        alpha = 0.0;
                    }
                    red   = r1 + rval*(r2-r1);
                    green = g1 + rval*(g2-g1);
                    blue  = b1 + rval*(b2-b1);
                    MyColor.setRgbF(red,green,blue,alpha); // rval*alpha); // DIAG alpha);
                    MyImg.setPixelColor(ii - ii_offset,jj, MyColor);
                }
            }
        }

        // Update refineMap_Item and make sure it is visible
        //--------------------------------------------------
        QPixmap MyPmap;
        double dltLat = (180./(jj_size));
        MyPmap.convertFromImage(MyImg);
        refineMapEdit_Item->setPixmap(MyPmap);
        refineMapEdit_Item->setPos(-180.,-90. );
        refineMapEdit_Item->setScale(dltLat);
        refineMapEdit_Item->setVisible(true);
        refineMapEdit_Item->update();
    } else {
        QPixmap MyPmap;
        MyPmap.fill(QColor(0,0,0,0));
        refineMapEdit_Item->setPixmap(MyPmap);
        refineMapEdit_Item->setVisible(false);
    }

    // refineCube_Item
    //-------------------
    if(refineCubeVisible) {
        if(refineCubeImgUpdate) {
            this->updateRefineCubeImg( vrmGrid);
        }
        refineCube_Pix->convertFromImage( *refineCubeImg);
        refineCube_Item->setPixmap( *refineCube_Pix);
        refineCube_Item->setVisible(true);
    } else {
        refineCube_Item->setVisible(false);
    }

    // DisplayGrid value determine the visibility of grids
    //----------------------------------------------------
    if(displayGridVisible) {
        if(displayGridName == "Base Resolution") {
            baseGridVisible   = true;
            refineGridVisible = false;
        } else if (displayGridName == "Variable Resolution") {
            baseGridVisible   = false;
            refineGridVisible = true;
        }
    } else {
        baseGridVisible   = false;
        refineGridVisible = false;
    }

    //  baseGrid_Item
    //---------------------
    if(baseGridVisible) {
        baseGrid_Item->loadGrid(vrmGrid->baseGrid);
        baseGrid_Item->setVisible(true);
        baseGrid_Item->update();
    } else {
        baseGrid_Item->setVisible(false);
    }

    // refineGrid_Item
    //---------------------
    if(refineGridVisible) {
        refineGrid_Item->loadGrid(vrmGrid->refineGrid);
        refineGrid_Item->setVisible(true);
        refineGrid_Item->update();
    } else {
        refineGrid_Item->setVisible(false);
    }

    // Update Displayed Grid Info
    //----------------------------
    emit do_gridInfoUpdate();

    // boundary rectangles for rect edit
    //-----------------------------------
    if(rectEdit_Item->Is_InsideRectVisible()) {
        rectEdit_Inside->setRect(rectEdit_Item->Get_InsideRect());
        rectEdit_Inside->setVisible(true);
    } else {
        rectEdit_Inside->setRect(rectEdit_Item->Get_InsideRect());
        rectEdit_Inside->setVisible(false);
    }
    if(rectEdit_Item->Is_OutsideRectVisible()) {
        rectEdit_Outside->setRect(rectEdit_Item->Get_OutsideRect());
        rectEdit_Outside->setVisible(true);
    } else {
        rectEdit_Outside->setRect(rectEdit_Item->Get_OutsideRect());
        rectEdit_Outside->setVisible(false);
    }

    // latlonGrid_Item
    //------------------
    update();

}


//=======================================================
// User Interface routines to Change/Update DISPLAY values
//=======================================================
double VrmWin::BackgroundFilterAlpha()
{
    return(backgroundFilter_Alpha);
}

void VrmWin::setBackgroundFilterAlpha(double I_Alpha)
{
    backgroundFilter_Alpha = I_Alpha;
    QColor MyColor;
    MyColor.setRgbF(1.,1.,1.,backgroundFilter_Alpha);
    backgroundFilter_Pix->fill(MyColor);
    backgroundFilter_Item->setPixmap( *backgroundFilter_Pix);
}

QString VrmWin::BackgroundImgVal()
{
    return(backgroundImg);
}

void VrmWin::setBackgroundImgVal(QString I_Image)
{
    backgroundImg = I_Image;
    if(backgroundImg == "TERRAIN") {
        backgroundImage_Pix->convertFromImage(*backgroundImg_Terrain);
        backgroundImage_Item->setPixmap( *backgroundImage_Pix);
        backgroundImage_Item->setScale(backgroundImg_TerrainScl);
    } else if (backgroundImg == "LAND/WATER") {
        backgroundImage_Pix->convertFromImage(*backgroundImg_LandWater);
        backgroundImage_Item->setPixmap( *backgroundImage_Pix);
        backgroundImage_Item->setScale(backgroundImg_LandWaterScl);
    } else if (backgroundImg == "NONE") {
        backgroundImage_Pix->convertFromImage(*backgroundImg_None);
        backgroundImage_Item->setPixmap( *backgroundImage_Pix);
        backgroundImage_Item->setScale(backgroundImg_NoneScl);
    } else {
        std::cout << " ADD ERROR MESSAGE HERE: Unknown background Image" << std::endl;
    }

}

bool VrmWin::ReferenceMapVisible()
{
    return(refineMapRefVisible);
}

void VrmWin::setReferenceMapVisible(bool I_Visible)
{
    refineMapRefVisible = I_Visible;
}

double VrmWin::ReferenceMapAlpha()
{
    return(refineMapRefAlpha);
}

void VrmWin::setReferenceMapAlpha(double I_Alpha)
{
    refineMapRefAlpha = I_Alpha;
}

bool VrmWin::RefineMapEditVisible()
{
    return(refineMapEditVisible);
}

void VrmWin::setRefineMapEditVisible(bool I_Visible)
{
    refineMapEditVisible = I_Visible;
}

bool VrmWin::RefineMapVisible()
{
    return(refineMapVisible);
}

void VrmWin::setRefineMapVisible(bool I_Visible)
{
    refineMapVisible = I_Visible;
}

double VrmWin::RefineMapAlpha()
{
    return(refineMapAlpha);
}

void VrmWin::setRefineMapAlpha(double I_Alpha)
{
    refineMapAlpha = I_Alpha;
}

bool VrmWin::RefineCubeVisible()
{
    return(refineCubeVisible);
}

void VrmWin::setRefineCubeVisible(bool I_Visible)
{
    refineCubeVisible = I_Visible;
}

double VrmWin::RefineCubeAlpha()
{
    return(refineCubeAlpha);
}

void VrmWin::setRefineCubeAlpha(double I_Alpha)
{
    refineCubeAlpha = I_Alpha;
}

void VrmWin::updateRefineCubeImg(VrmGrid *I_vrmGrid)
{
    double dLon;
    double dLat;
    int    iA,iB,iP;

    double r1=1.0;
    double r2=0.0;
    double b1=1.0;
    double b2=0.0;
    double g1=1.0;
    double g2=1.0;

    double red;
    double green;
    double blue;
    double alpha;
    double rval;
    QColor MyColor;
    MyColor.setRgbF(0.,0.,0.,0.);

    int nRes = I_vrmGrid->refineCube.GetMaxResolution();
    int MaxVal = 0;
    for( int iP=0; iP < 6   ; iP++) {
    for( int iA=0; iA < nRes; iA++) {
    for( int iB=0; iB < nRes; iB++) {
        if( I_vrmGrid->refineCube.val[iP][iA][iB] > MaxVal ) {
            MaxVal = I_vrmGrid->refineCube.val[iP][iA][iB];
        }
    }
    }
    }

//    std::cout << " loading Image...  MaxVal=" << MaxVal << std::endl;
    double dltLon = 720./static_cast<double>(PXSIZE);
    for( int ii=0; ii < PXSIZE; ii++) {
        dLon = dltLon*(ii) - 180. - .5*dltLon;
        for( int jj=0; jj < PYSIZE; jj++) {
            dLat = -1.*(-180.*(jj/(static_cast<double>(PYSIZE)-1.)) + 90.);
            I_vrmGrid->refineCube.get_CubeIndices(dLon,dLat, iP, iA, iB);
            if(MaxVal > 0) {
                rval = static_cast<double>( I_vrmGrid->refineCube.val[iP][iA][iB]) / static_cast<double>(MaxVal);
            } else {
                rval = 0.0;
            }
            if(rval > 0.01) {
                alpha = refineCubeAlpha;
            } else {
                alpha = 0.0;
            }
            red   = r1 + rval*(r2-r1);
            green = g1 + rval*(g2-g1);
            blue  = b1 + rval*(b2-b1);
            MyColor.setRgbF(red,green,blue,alpha);
            refineCubeImg->setPixelColor(ii,jj, MyColor);
        }
    }
    refineCube_Pix->convertFromImage( *refineCubeImg);
    refineCube_Item->setPixmap( *refineCube_Pix);  // IS THIS NEEDED ???
//    std::cout << "  .....Done." << std::endl;
}

void VrmWin::adjustRefineCubeImgAlpha(double I_Alpha)
{
    QColor MyColor;

    for( int ii=0; ii < PXSIZE; ii++) {
        for( int jj=0; jj < PYSIZE; jj++) {
            MyColor = refineCubeImg->pixelColor(ii,jj);
            if(MyColor.alphaF() > 0.) {
                MyColor.setAlphaF(I_Alpha);
                refineCubeImg->setPixelColor(ii,jj,MyColor);
            }
        }
    }
    refineCube_Pix->convertFromImage( *refineCubeImg);
    refineCube_Item->setPixmap( *refineCube_Pix);
}

bool VrmWin::DisplayGridVisible()
{
    return(displayGridVisible);
}

void VrmWin::setDisplayGridVisible(bool I_Visible)
{
    displayGridVisible = I_Visible;
}

QString VrmWin::DisplayGridName()
{
    return(displayGridName);
}

void VrmWin::setDisplayGridName(QString I_GridName)
{
    displayGridName = I_GridName;
}

void VrmWin::setSmoothKernel(int I_val)
{
    smoothKernel = I_val;
}

void VrmWin::setSmoothPersistence(double I_val)
{
    smoothPersistence = I_val;
}

void VrmWin::setSmoothIterNum(int I_val)
{
    smoothIterNum = I_val;
}

void VrmWin::applySmoothing(VrmGrid *I_vrmGrid)
{
    // Apply the desired smoothing to refineMapEdit
    //------------------------------------------------
//    std::cout << " Apply Smoothing " << std::endl;

    int NLON = I_vrmGrid->refineMapEdit.nRefLon;
    int NLAT = I_vrmGrid->refineMapEdit.nRefLat;

    double sval[NLON][NLAT];
    double rval[NLON][NLAT];

    // For now just use a simple averageing kernel
    //--------------------------------------------
    int    KDIM = smoothKernel;
    int    KLEN  =((KDIM-1)/2);

    double wval[KDIM][KDIM];
    for(int ii=0; ii < KDIM; ii++) {
    for(int jj=0; jj < KDIM; jj++) {
        wval [ii][jj] = 1.0/static_cast<double>(KDIM*KDIM);
    }
    }

    for(int ii=0; ii < NLON; ii++) {
    for(int jj=0; jj < NLAT; jj++) {
        rval[ii][jj] = I_vrmGrid->refineMapEdit.val[ii][jj];
    }
    }

    // Loop over the Number of Iterations
    //-------------------------------------
    for(int iter=0; iter < smoothIterNum; iter++) {

        // Loop over refineMapEdit gridpoints
        //-------------------------------------
        for(int ii=0; ii < NLON; ii++) {
        for(int jj=0; jj < NLAT; jj++) {

            sval[ii][jj] = 0.0;

            // Apply smoothing Kernel
            //------------------------
            for(int ik=0; ik < KDIM; ik++) {
            for(int jk=0; jk < KDIM; jk++) {
                int ix = (ii + ik - KLEN);
                int jx = (jj + jk - KLEN);
                if(ix <     0) { ix += NLON;}
                if(ix >= NLON) { ix -= NLON;}
                if(jx < 0    ) {
                    jx = -jx;
                    ix += NLON;
                    if(ix >= NLON) {ix -= NLON;}
                }
                if(jx >= NLAT) {
                    jx = 2*NLAT - jx;
                    ix += NLON;
                    if(ix >= NLON) {ix -= NLON;}
                }

                sval[ii][jj] += wval[ik][jk] * rval[ix][jx];
            }
            }
        }
        }

        // Copy smooth values to refineMapEdit
        //----------------------------------------
        double MyMax = -9999.;
        double MyMin =  9999.;
        for(int ii=0; ii < NLON; ii++) {
        for(int jj=0; jj < NLAT; jj++) {
            if(sval[ii][jj] < 0.) { sval[ii][jj] = 0.;}
            if(sval[ii][jj] > 1.) { sval[ii][jj] = 1.;}

            if(I_vrmGrid->refineMapEdit.val[ii][jj] > 0.999) {
                rval[ii][jj] =      smoothPersistence *I_vrmGrid->refineMapEdit.val[ii][jj]
                              + (1.-smoothPersistence)*sval[ii][jj];
            } else {
                rval[ii][jj] = sval[ii][jj];
            }
            if(sval[ii][jj] > MyMax) { MyMax = sval[ii][jj];}
            if(sval[ii][jj] < MyMin) { MyMin = sval[ii][jj];}
        }
        }
//        std::cout << " MyMin=" << MyMin << " MyMax=" << MyMax << std::endl;


    } // iter-LOOP

    for(int ii=0; ii < NLON; ii++) {
    for(int jj=0; jj < NLAT; jj++) {
        I_vrmGrid->refineMapEdit.val[ii][jj] = rval[ii][jj];
    }
    }

}

void VrmWin::UpdateNodes(QPolygonF MyPoly)
{
//    std::cout << " UpdateNodes Called" << std::endl;
//    for( int ii=0; ii < 7; ii++) {
//        std::cout << "Poly:" << MyPoly[ii].x() << " " << MyPoly[ii].y() << std::endl;
//    }
//    std::cout << " " << std::endl;

    Node1->updateNode(MyPoly[0]);
    Node2->updateNode(MyPoly[1]);
    Node3->updateNode(MyPoly[2]);
    Node4->updateNode(MyPoly[3]);
    Node5->updateNode(MyPoly[4]);
    Node6->updateNode(MyPoly[5]);
    Node7->updateNode(MyPoly[6]);
    Node8->updateNode(MyPoly[7]);
}

void VrmWin::UpdateNodePoly()
{
    QPolygonF MyPoly;
    MyPoly << Node1->MyPoint << Node2->MyPoint << Node3->MyPoint << Node4->MyPoint
           << Node5->MyPoint << Node6->MyPoint << Node7->MyPoint << Node8->MyPoint << Node1->MyPoint;
    polyEdit_Item->UpdatePolygon(MyPoly);
}

