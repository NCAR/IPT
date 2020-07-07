#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QGraphicsItem>
#include <QGraphicsRectItem>
#include <QImage>
#include <QPixmap>
#include <QPen>
#include <QBrush>
#include <iostream>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    int  HRES = 1024;
    int  VRES = 512;
    int  TRES = 250;

    ui->setupUi(this);
    ui->tabWidget->setMinimumWidth(TRES);
    ui->tabWidget->setMaximumWidth(TRES);
    ui->VrmOpts->setEnabled(true);
    ui->editActionsGroup->setEnabled(false);

    vrmView = new VrmView(HRES,VRES, parent);
    vrmWin  = new VrmWin(vrmView->getLonMin(),vrmView->getLatMin(),vrmView->getLonMax(),vrmView->getLatMax());
    lonbar  = new LonBar(vrmView->getLonMin(), vrmView->getLonMax(), vrmView->getBarHeight());
    latbar  = new LatBar(vrmView->getLatMin(), vrmView->getLatMax(), vrmView->getBarHeight());
    initLonScene(vrmView->getBarHeight());
    initLatScene(vrmView->getBarHeight());
    vrmWin->initialize();
    connect(vrmView, SIGNAL(   setViewWindows(int,int,int,int,double,double,double,double)),
               this, SLOT  (on_setViewWindows(int,int,int,int,double,double,double,double)));
    vrmView->initialize();

    // Initialize VRM grid parameters from UI defaults.
    //--------------------------------------------------
    vrmGrid.Set_RefineType        ( ui->refineTypeVal->currentText()    );
    vrmGrid.Set_SmoothType        ( ui->smoothTypeVal->currentText()    );
    vrmGrid.Set_GridType          ( ui->refineGridVal->currentText()    );
    vrmGrid.Set_Resolution        ( ui->refineBaseVal->value()          );
    vrmGrid.Set_RefinementLevel   ( ui->refineLevelVal->value()         );
    vrmGrid.Set_SmoothIterations  ( ui->smoothIterVal->value()          );
    vrmGrid.Set_TransSmoothDist   ( ui->smoothDistVal->value()          );
    vrmGrid.Set_GridLonShift      ( ui->lonShiftVal->value()            );
    vrmGrid.Set_GridXRotate       ( ui->XrotateVal->value()             );
    vrmGrid.Set_GridYRotate       ( ui->YrotateVal->value()             );
    vrmGrid.Set_Tessellations     ( ui->tessellationsVal->value()       );
    vrmGrid.Set_SubCellResolution ( ui->subCellResolutionVal->value()   );
    vrmGrid.Set_ReverseOrientation( ui->reverseOrientation->isChecked() );

    // Initialize Refinement map/cube
    //--------------------------------
//TEST     vrmGrid.refineMap.initialize(360,180);
//TEST     vrmGrid.refineMapRef.initialize(360,180);
//TEST     vrmGrid.refineMapEdit.initialize(360,180);
    vrmGrid.refineMap.initialize(720,360);
    vrmGrid.refineMapRef.initialize(360,180);
    vrmGrid.refineMapEdit.initialize(720,360);
    vrmGrid.refineCube.resize(vrmGrid.Get_Resolution(),vrmGrid.Get_RefinementLevel()); 

    QPalette mypal;
    QColor   mycol;
    mycol.setRgbF(1.,1.,0.);
    mypal.setColor(QPalette::Background,mycol);
    ui->referenceMapColorVal->setPalette(mypal);
    mycol.setRgbF(1.,0.,0.);
    mypal.setColor(QPalette::Background,mycol);
    ui->refineMapColorVal->setPalette(mypal);
    mycol.setRgbF(0.,1.,0.);
    mypal.setColor(QPalette::Background,mycol);
    ui->cubeMapColorVal->setPalette(mypal);

    connect( vrmWin->polyEdit_Item, SIGNAL(do_SetPolyVal(QPolygonF*,double)),
            &vrmGrid              , SLOT  (   SetPolyVal(QPolygonF*,double)));
    connect( vrmWin->polyEdit_Item, SIGNAL(do_SceneUpdate()),
             this                 , SLOT  (do_SceneUpdate()));
    connect( vrmWin->rectEdit_Item, SIGNAL(do_SceneUpdate()),
             this                 , SLOT  (do_SceneUpdate()));
    connect( vrmWin->rectEdit_Item, SIGNAL(do_UIupdate()),
             this                 , SLOT  (do_rectUIupdate()));
    connect( vrmWin->polyEdit_Item, SIGNAL(do_UpdateNodes(QPolygonF)),
             vrmWin               , SLOT  (   UpdateNodes(QPolygonF)));
    connect( vrmWin->Node1        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node2        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node3        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node4        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node5        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node6        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node7        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin->Node8        , SIGNAL(do_UpdateNodePoly()),
             vrmWin               , SLOT  (   UpdateNodePoly()));
    connect( vrmWin               , SIGNAL(do_gridInfoUpdate()),
             this                 , SLOT  (   gridInfoUpdate()));

    // Initialize Grid Info values
    //-----------------------------
    gridInfoUpdate();

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initLonScene(int I_htext)
{
    double htext = 0.55 * static_cast<double> (I_htext);

    LonScene = new QGraphicsScene(this);
    LonScene->setSceneRect(lonbar->boundingRect());
    LonScene->addItem(lonbar);

    for( double Lon = 0.; Lon < 540.; Lon += 30.) {
        double mLon = fmod(Lon,360);
        QString StrLon = QString::number(mLon);
        QGraphicsTextItem *LonLab = new QGraphicsTextItem(StrLon);
        LonLab->setPos(Lon,htext);
        LonLab->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        LonLab->document()->setDocumentMargin(0);
        LonLab->document()->setDefaultTextOption(QTextOption(Qt::AlignRight | Qt::AlignVCenter));
        LonScene->addItem(LonLab);
    }
    for( double Lon = -30.; Lon > -180.; Lon -= 30.) {
        double mLon = fmod(Lon,360);
        QString StrLon = QString::number(mLon);
        QGraphicsTextItem *LonLab = new QGraphicsTextItem(StrLon);
        LonLab->setPos(Lon,htext);
        LonLab->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        LonLab->document()->setDocumentMargin(0);
        LonLab->document()->setDefaultTextOption(QTextOption(Qt::AlignJustify));
        LonScene->addItem(LonLab);
    }
}

void MainWindow::initLatScene(int I_htext)
{
    double htext = 0.15 * static_cast<double> (I_htext);

    LatScene = new QGraphicsScene(this);
    LatScene->setSceneRect(latbar->boundingRect());
    LatScene->addItem(latbar);
    for( double Lat = -90.; Lat < 90.; Lat += 30.) {
        QString StrLat = QString::number(Lat);
        QGraphicsTextItem *LatLab = new QGraphicsTextItem(StrLat);
        LatLab->setPos(htext,Lat);
        LatLab->setFlag(QGraphicsItem::ItemIgnoresTransformations);
        LatLab->document()->setDocumentMargin(0);
        LatLab->document()->setDefaultTextOption(QTextOption(Qt::AlignCenter));
        LatScene->addItem(LatLab);
    }
}

bool MainWindow::eventFilter(QObject *watched, QEvent *event)
{
    if(watched == ui->LatScale) {
//        std::cout << " LatScale Filter Event" << std::endl;
        if(event->type() == QEvent::GraphicsSceneWheel) {
//            std::cout << " LatScale Wheel Event" << std::endl;
        }
        event->accept();
        return true;
    } else {
        return false;
    }

}

void MainWindow::keyPressEvent(QKeyEvent *event)
{
    double step_factor = 30.;
    double lat_step;
    double new_lat0;
    double lon_step;
    double new_lon0;
    double Lat0;
    double Lon0;
    double DltLat;
    double DltLon;
    double LatMin;
    double LatMax;
    double LonMin;
    double LonMax;
    double LatC1;
    double LatC2;
    double LonC1;
    double LonC2;
    double ZoomScl;
    double ZoomSclMin;
    double SclLon;
    double SclLat;

    switch(event->key()) {
    case Qt::Key_Escape:
        ZoomSclMin = vrmView->getZoomSclMin();
        SclLon     = vrmView->getSclLon();
        SclLat     = vrmView->getSclLat();
        DltLon     = vrmView->getDltLon();
        DltLat     = vrmView->getDltLat();
        LatMin     = vrmView->getLatMin();
        LatMax     = vrmView->getLatMax();
        LonMin     = vrmView->getLonMin();
        LonMax     = vrmView->getLonMax();

        ZoomScl = ZoomSclMin;
        Lon0 = (LonMax + LonMin)/2.;
        Lat0 = (LatMax + LatMin)/2.;
        LonC1 = Lon0 - DltLon/2.;
        LatC1 = Lat0 - DltLat/2.;
        LonC2 = LonC1 + DltLon;
        LatC2 = LatC1 + DltLat;
        vrmView->setLat0(Lat0,LatC1,LatC2);
        vrmView->setLon0(Lon0,LonC1,LonC2);
        vrmView->setZoomScl(ZoomScl);

        ui->graphicsView->resetTransform();
        ui->graphicsView->scale   ( SclLon, SclLat);
        ui->graphicsView->centerOn( Lon0, Lat0);
        ui->LonScale->resetTransform();
        ui->LonScale->scale   ( SclLon,  1.);
        ui->LonScale->centerOn( Lon0,  0.);
        ui->LatScale->resetTransform();
        ui->LatScale->scale   ( 1., SclLat);
        ui->LatScale->centerOn( 0., Lat0);
        event->accept();
        update();
        break;
    case Qt::Key_Up:
        DltLat = vrmView->getDltLat()/vrmView->getZoomScl();
        LatMax = vrmView->getLatMax();
        LatMin = vrmView->getLatMin();
        Lat0   = vrmView->getLat0();
        Lon0   = vrmView->getLon0();
        lat_step = DltLat/step_factor;
        new_lat0 = Lat0 + lat_step;
        if((new_lat0 + DltLat/2.) > LatMax) {
            new_lat0 = LatMax - DltLat/2.;
        } else if ((new_lat0 - DltLat/2.) < LatMin) {
            new_lat0 = LatMin + DltLat/2.;
        }
        Lat0 = new_lat0;
        LatC1 = Lat0 - DltLat/2.;
        LatC2 = LatC1 + DltLat;
        vrmView->setLat0(Lat0,LatC1,LatC2);
        ui->graphicsView->centerOn(Lon0, Lat0);
        ui->LonScale->centerOn    (Lon0,   0.);
        ui->LatScale->centerOn    (  0., Lat0);
        event->accept();
        update();
        break;
    case Qt::Key_Down:
        DltLat = vrmView->getDltLat()/vrmView->getZoomScl();
        LatMax = vrmView->getLatMax();
        LatMin = vrmView->getLatMin();
        Lat0   = vrmView->getLat0();
        Lon0   = vrmView->getLon0();
        lat_step = DltLat/step_factor;
        new_lat0 = Lat0 - lat_step;
        if((new_lat0 + DltLat/2.) > LatMax) {
            new_lat0 = LatMax - DltLat/2.;
        } else if ((new_lat0 - DltLat/2.) < LatMin) {
            new_lat0 = LatMin + DltLat/2.;
        }
        Lat0 = new_lat0;
        LatC1 = Lat0 - DltLat/2.;
        LatC2 = LatC1 + DltLat;
        vrmView->setLat0(Lat0,LatC1,LatC2);
        ui->graphicsView->centerOn(Lon0, Lat0);
        ui->LonScale->centerOn    (Lon0,   0.);
        ui->LatScale->centerOn    (  0., Lat0);
        event->accept();
        update();
        break;
    case Qt::Key_Right:
        DltLon = vrmView->getDltLon()/vrmView->getZoomScl();
        LonMax = vrmView->getLonMax();
        LonMin = vrmView->getLonMin();
        Lat0   = vrmView->getLat0();
        Lon0   = vrmView->getLon0();
        lon_step = DltLon/step_factor;
        new_lon0 = Lon0 + lon_step;
        if((new_lon0 + DltLon/2.) > LonMax) {
            new_lon0 = LonMax - DltLon/2.;
        } else if ((new_lon0 - DltLon/2.) < LonMin) {
            new_lon0 = LonMin + DltLon/2.;
        }
        Lon0 = new_lon0;
        LonC1 = Lon0 - DltLon/2.;
        LonC2 = LonC1 + DltLon;
        vrmView->setLon0(Lon0,LonC1,LonC2);
        ui->graphicsView->centerOn(Lon0, Lat0);
        ui->LonScale->centerOn    (Lon0,   0.);
        ui->LatScale->centerOn    (  0., Lat0);
        event->accept();
        update();
        break;
    case Qt::Key_Left:
        DltLon = vrmView->getDltLon()/vrmView->getZoomScl();
        LonMax = vrmView->getLonMax();
        LonMin = vrmView->getLonMin();
        Lat0   = vrmView->getLat0();
        Lon0   = vrmView->getLon0();
        lon_step = DltLon/step_factor;
        new_lon0 = Lon0 - lon_step;
        if((new_lon0 + DltLon/2.) > LonMax) {
            new_lon0 = LonMax - DltLon/2.;
        } else if ((new_lon0 - DltLon/2.) < LonMin) {
            new_lon0 = LonMin + DltLon/2.;
        }
        Lon0 = new_lon0;
        LonC1 = Lon0 - DltLon/2.;
        LonC2 = LonC1 + DltLon;
        vrmView->setLon0(Lon0,LonC1,LonC2);
        ui->graphicsView->centerOn(Lon0, Lat0);
        ui->LonScale->centerOn    (Lon0,   0.);
        ui->LatScale->centerOn    (  0., Lat0);
        event->accept();
        update();
        break;
    }
}

void MainWindow::wheelEvent(QWheelEvent *event)
{
    double ZoomScl    = vrmView->getZoomScl();
    double ZoomSclMax = vrmView->getZoomSclMax();
    double ZoomSclMin = vrmView->getZoomSclMin();
    double ZoomSclInc = vrmView->getZoomSclInc();
    double SclLon     = vrmView->getSclLon();
    double SclLat     = vrmView->getSclLat();
    double DltLon     = vrmView->getDltLon();
    double DltLat     = vrmView->getDltLat();
    double Lon0       = vrmView->getLon0();
    double Lat0       = vrmView->getLat0();
    double LatMin     = vrmView->getLatMin();
    double LatMax     = vrmView->getLatMax();
    double LonMin     = vrmView->getLonMin();
    double LonMax     = vrmView->getLonMax();

    if(event->delta() > 0) {
        ZoomScl *= ZoomSclInc;
        if(ZoomScl > ZoomSclMax) {ZoomScl = ZoomSclMax;}

        double slon = SclLon*ZoomScl;
        double slat = SclLat*ZoomScl;
        double dlon = DltLon/ZoomScl;
        double dlat = DltLat/ZoomScl;

        double LonC1 = Lon0 - dlon/2.;
        double LatC1 = Lat0 - dlat/2.;
        double LonC2 = LonC1 + dlon;
        double LatC2 = LatC1 + dlat;
        if(LatC1 < LatMin) {
            LatC1 = LatMin;
            Lat0  = LatMin + dlat/2.;
            LatC2 = LatMin + dlat;
        } else if (LatC2 > LatMax) {
            LatC1 = LatMax - dlat;
            Lat0  = LatMax - dlat/2.;
            LatC2 = LatMax;
        }
        if(LonC1 < LonMin) {
            LonC1 = LonMin;
            Lon0  = LonMin + dlon/2.;
            LonC2 = LonMin + dlon;
        } else if (LonC2 > LonMax) {
            LonC1 = LonMax - dlon;
            Lon0  = LonMax - dlon/2.;
            LonC2 = LonMax;
        }
        vrmView->setLat0(Lat0,LatC1,LatC2);
        vrmView->setLon0(Lon0,LonC1,LonC2);
        vrmView->setZoomScl(ZoomScl);

        ui->graphicsView->resetTransform();
        ui->graphicsView->scale   ( slon, slat);
        ui->graphicsView->centerOn( Lon0, Lat0);
        ui->LonScale->resetTransform();
        ui->LonScale->scale   ( slon,  1.);
        ui->LonScale->centerOn( Lon0,  0.);
        ui->LatScale->resetTransform();
        ui->LatScale->scale   ( 1., slat);
        ui->LatScale->centerOn( 0., Lat0);
        update();

    } else {
        ZoomScl /= ZoomSclInc;
        if(ZoomScl < ZoomSclMin) {ZoomScl = ZoomSclMin;}

        double slon = SclLon*ZoomScl;
        double slat = SclLat*ZoomScl;
        double dlon = DltLon/ZoomScl;
        double dlat = DltLat/ZoomScl;

        double LonC1 = Lon0 - dlon/2.;
        double LatC1 = Lat0 - dlat/2.;
        double LonC2 = LonC1 + dlon;
        double LatC2 = LatC1 + dlat;
        if(LatC1 < LatMin) {
            LatC1 = LatMin;
            Lat0  = LatMin + dlat/2.;
            LatC2 = LatMin + dlat;
        } else if (LatC2 > LatMax) {
            LatC1 = LatMax - dlat;
            Lat0  = LatMax - dlat/2.;
            LatC2 = LatMax;
        }
        if(LonC1 < LonMin) {
            LonC1 = LonMin;
            Lon0  = LonMin + dlon/2.;
            LonC2 = LonMin + dlon;
        } else if (LonC2 > LonMax) {
            LonC1 = LonMax - dlon;
            Lon0  = LonMax - dlon/2.;
            LonC2 = LonMax;
        }

        vrmView->setLat0(Lat0,LatC1,LatC2);
        vrmView->setLon0(Lon0,LonC1,LonC2);
        vrmView->setZoomScl(ZoomScl);

        ui->graphicsView->resetTransform();
        ui->graphicsView->scale   ( slon, slat);
        ui->graphicsView->centerOn( Lon0, Lat0);
        ui->LonScale->resetTransform();
        ui->LonScale->scale   ( slon,  1.);
        ui->LonScale->centerOn( Lon0,  0.);
        ui->LatScale->resetTransform();
        ui->LatScale->scale   ( 1., slat);
        ui->LatScale->centerOn( 0., Lat0);
        update();
    }
    event->accept();
}

void MainWindow::wheelEvent(QGraphicsSceneEvent *event)
{
    event->accept();
}

void MainWindow::do_SceneUpdate()
{
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::do_rectUIupdate()
{
    double LonMin = vrmWin->rectEdit_Item->Get_LonMin();
    double LonMax = vrmWin->rectEdit_Item->Get_LonMax();
    double LatMin = vrmWin->rectEdit_Item->Get_LatMin();
    double LatMax = vrmWin->rectEdit_Item->Get_LatMax();

    if(ui->rectLonMinVal->value() != LonMin) {ui->rectLonMinVal->setValue(LonMin);}
    if(ui->rectLonMaxVal->value() != LonMax) {ui->rectLonMaxVal->setValue(LonMax);}
    if(ui->rectLatMinVal->value() != LatMin) {ui->rectLonMinVal->setValue(LatMin);}
    if(ui->rectLatMaxVal->value() != LatMax) {ui->rectLonMaxVal->setValue(LatMax);}
}

void MainWindow::gridInfoUpdate()
{
    // Update grid info values for the current visible grid
    //-------------------------------------------------------
    QString GridName = vrmWin->DisplayGridName();
    if(GridName == "Base Resolution") {
        gridInfo.Set_GridInfo(vrmGrid.baseGrid);
    } else if (GridName == "Variable Resolution") {
        gridInfo.Set_GridInfo(vrmGrid.refineGrid);
    }

    // Update Grid Info values
    //-----------------------------
    QString S1;
    QString S2;
    ui->gridInformation ->setEnabled(gridInfo.Get_Enabled());
    ui->faceClockwiseVal->setText(    S1.setNum(gridInfo.Get_Num_Clockwise()));
    ui->faceCounterVal  ->setText(    S1.setNum(gridInfo.Get_Num_Counter  ()));
    ui->faceTotalAreaVal->setText(    S1.setNum(gridInfo.Get_Total_Area() ,'g',4));
    ui->areaRatioVal    ->setText(    S1.setNum(gridInfo.Get_Area_Ratio() ,'g',4));
    ui->areaMaxVal      ->setText(    S1.setNum(gridInfo.Get_MaxArea()    ,'g',4));
    ui->areaMaxLonLat   ->setText("("+S1.setNum(gridInfo.Get_MaxArea_Lon(),'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MaxArea_Lat(),'f',1)+")");
    ui->areaMinVal      ->setText(    S1.setNum(gridInfo.Get_MinArea()    ,'g',4));
    ui->areaMinLonLat   ->setText("("+S1.setNum(gridInfo.Get_MinArea_Lon(),'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MinArea_Lat(),'f',1)+")");
    ui->arcQualityVal   ->setText(    S1.setNum(gridInfo.Get_Qarc()       ,'g',4));
    ui->arcMaxVal       ->setText(    S1.setNum(gridInfo.Get_MaxArc()     ,'g',4));
    ui->arcMaxLonLat    ->setText("("+S1.setNum(gridInfo.Get_MaxArc_Lon() ,'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MaxArc_Lat() ,'f',1)+")");
    ui->arcMinVal       ->setText(    S1.setNum(gridInfo.Get_MinArc()     ,'g',4));
    ui->arcMinLonLat    ->setText("("+S1.setNum(gridInfo.Get_MinArc_Lon() ,'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MinArc_Lat() ,'f',1)+")");
    ui->angQualityVal   ->setText(    S1.setNum(gridInfo.Get_Qang()       ,'g',4));
    ui->angMaxVal       ->setText(    S1.setNum(gridInfo.Get_MaxAng()     ,'g',4));
    ui->angMaxLonLat    ->setText("("+S1.setNum(gridInfo.Get_MaxAng_Lon() ,'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MaxAng_Lat() ,'f',1)+")");
    ui->angMinVal       ->setText(    S1.setNum(gridInfo.Get_MinAng()     ,'g',4));
    ui->angMinLonLat    ->setText("("+S1.setNum(gridInfo.Get_MinAng_Lon() ,'f',1)+" , "
                                     +S2.setNum(gridInfo.Get_MinAng_Lat() ,'f',1)+")");
}

void MainWindow::on_setViewWindows(int HSIZE, int VSIZE, int LATSIZE, int LONSIZE,
                                   double lon0, double lat0, double slon, double slat)
{
    ui->graphicsView->setMinimumHeight(VSIZE);
    ui->graphicsView->setMinimumWidth (HSIZE);
    ui->graphicsView->setMaximumHeight(VSIZE);
    ui->graphicsView->setMaximumWidth (HSIZE);
    ui->graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    ui->graphicsView->setVerticalScrollBarPolicy  (Qt::ScrollBarAlwaysOff);
    ui->graphicsView->scale   ( slon, slat);
    ui->graphicsView->centerOn( lon0, lat0);

    ui->LonScale->setMinimumHeight(LATSIZE);
    ui->LonScale->setMinimumWidth (HSIZE);
    ui->LonScale->setMaximumHeight(LATSIZE);
    ui->LonScale->setMaximumWidth (HSIZE);
    ui->LonScale->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    ui->LonScale->setVerticalScrollBarPolicy  (Qt::ScrollBarAlwaysOff);
    ui->LonScale->setStyleSheet("background: transparent");
    ui->LonScale->scale   ( slon,  1.);
    ui->LonScale->centerOn( lon0,  0.);
    ui->LonScale->setFocusPolicy(Qt::NoFocus);
    ui->LonScale->setInteractive(false);

    ui->LatScale->setMinimumHeight(VSIZE);
    ui->LatScale->setMinimumWidth (LONSIZE);
    ui->LatScale->setMaximumHeight(VSIZE);
    ui->LatScale->setMaximumWidth (LONSIZE);
    ui->LatScale->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    ui->LatScale->setVerticalScrollBarPolicy  (Qt::ScrollBarAlwaysOff);
    ui->LatScale->setStyleSheet("background: transparent");
    ui->LatScale->scale   ( 1., slat);
    ui->LatScale->centerOn( 0., lat0);
    ui->LatScale->setFocusPolicy(Qt::NoFocus);
    ui->LatScale->setInteractive(false);

    ui->LonScale->setScene(LonScene);
    ui->LatScale->setScene(LatScene);
    ui->graphicsView->setScene(vrmWin->Scene);
}

//============================================================
// User Interface for selected ACTIONS on menu bar
//============================================================
void MainWindow::on_actionSave_VRM_File_triggered()
{

    // If the application is in EDIT mode, this action must be ignored
    //-----------------------------------------------------------------
    if( ui->beginEditMode->isChecked()) {
        std::cout << "Can't Write a VRM File while in EDIT mode." << std::endl;
        return;
    } else {
        std::cout << " Not Checked?" << std::endl;
    }
}

void MainWindow::on_actionRead_VRM_File_triggered()
{

    // If the application is in EDIT mode, this action must be ignored
    //-----------------------------------------------------------------
    if( ui->beginEditMode->isChecked()) {
        std::cout << "Can't Read a VRM File while in EDIT mode." << std::endl;
        return;
    } else {
        std::cout << " Not Checked?" << std::endl;
    }
}

void MainWindow::on_actionWrite_Exodus_File_triggered()
{

    // If the application is in EDIT mode, this action must be ignored
    //-----------------------------------------------------------------
    if( ui->beginEditMode->isChecked()) {
        std::cout << "Can't Write an Exodus File while in EDIT mode." << std::endl;
        return;
    } else {
        std::cout << " Not Checked?" << std::endl;
    }

    //######################################################################################
    // NEED A TEST OF refineGrid STATUS - calc_VarMesh() need to be called at least 1 time.
    //######################################################################################

    // Get the desired file from the user.
    //-------------------------------------
    QString ExodusFile = QFileDialog::getSaveFileName(this,"Write VRM mesh to Exodus File",
                                                      QDir::currentPath(),"All files (*.*)");

    // return if no file was selected
    //--------------------------------
    if(ExodusFile.isEmpty()) { return; }

    // The Read/Write functions need to be independent of Qt(Strings) 
    // so they can be shared by the standalone program.
    // So copy values into the VRM structure passed to Write_Exodus_File()
    //-----------------------------------------------------------------
    VRM_Param VRM;
    VRM.Resolution         = vrmGrid.Get_Resolution();
    VRM.RefineType         = vrmGrid.Get_RefineType().toStdString();
    VRM.GridType           = vrmGrid.Get_GridType().toStdString();
    VRM.GridXRotate        = vrmGrid.Get_GridXRotate();
    VRM.GridYRotate        = vrmGrid.Get_GridYRotate();
    VRM.GridLonShift       = vrmGrid.Get_GridLonShift();
    VRM.RefinementLevel    = vrmGrid.Get_RefinementLevel();
    VRM.Tessellations      = vrmGrid.Get_Tessellations();
    VRM.SubCellResolution  = vrmGrid.Get_SubCellResolution();
    VRM.ReverseOrientation = vrmGrid.Get_ReverseOrientation();
    VRM.SmoothType         = vrmGrid.Get_SmoothType().toStdString();
    VRM.TransSmoothDist    = vrmGrid.Get_TransSmoothDist();
    VRM.SmoothIterations   = vrmGrid.Get_SmoothIterations();

    // Output to a NetCDF Exodus file
    //--------------------------------
    std::string ExodusName =  ExodusFile.toStdString();
    Write_Exodus_File(ExodusName , VRM, vrmGrid.refineMap, vrmGrid.refineGrid);

}

void MainWindow::on_actionRead_Refinement_Map_triggered()
{

    // If the application is in EDIT mode, this action must be ignored
    //-----------------------------------------------------------------
    if( ui->beginEditMode->isChecked()) {
        std::cout << "Can't Read a refineMap while in EDIT mode." << std::endl;
        return;
    } else {
        std::cout << " Not Checked?" << std::endl;
    }

    QString RefMapFile = QFileDialog::getOpenFileName(this,"Open Refinement Map",QDir::currentPath(),"All files (*.*)");

    // return if no file was selected
    //--------------------------------
    if(RefMapFile.isEmpty()) { return; }

    //##########################################################################
    // TODO: Check input values span the range [0,1], that the longitudes
    //       are (0,360), and the latitudes are (-90,90).
    //     OPTIONS? if values exceed [0,1], if the min !=0 and max !=1. --> rescale?
    //    What if it is not a refineMap file?
    //###########################################################################

    // Open a netCDF file for input
    //---------------------------------
    std::string RefMapName =  RefMapFile.toStdString();
    Read_Refinement_Map( RefMapName, vrmGrid.refineMap);

    // update vrmWin
    //---------------
    vrmWin->update_scene(&vrmGrid);
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}


void MainWindow::on_actionSave_Refinement_Map_triggered()
{
    QString RefMapFile = QFileDialog::getSaveFileName(this,"Save Refinement Map",QDir::currentPath(),"All files (*.*)");

    // return if no file was selected
    //--------------------------------
    if(RefMapFile.isEmpty()) { return; }

    // Open a netCDF file for output
    //---------------------------------
    std::string RefMapName =  RefMapFile.toStdString();
    Write_Refinement_Map( RefMapName, vrmGrid.refineMap);

}

void MainWindow::on_actionRead_Reference_Map_triggered()
{

    QString RefMapFile = QFileDialog::getOpenFileName(this,"Open Reference Map",QDir::currentPath(),"All files (*.*)");

    // return if no file was selected
    //--------------------------------
    if(RefMapFile.isEmpty()) { return; }

    // Open a netCDF file for input
    //---------------------------------
    std::string RefMapName =  RefMapFile.toStdString();
    Read_Reference_Map( RefMapName, vrmGrid.refineMapRef);

    // update vrmWin
    //---------------
    vrmWin->update_scene(&vrmGrid);
}

void MainWindow::on_actionQuit_triggered()
{
    close();
}

//============================================================
// User Interface for changes to VRM settings
//============================================================
void MainWindow::on_refineGridVal_currentTextChanged(const QString &arg1)
{
    vrmGrid.Set_GridType( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_refineTypeVal_currentTextChanged(const QString &arg1)
{
    vrmGrid.Set_RefineType( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_refineBaseVal_valueChanged(int arg1)
{
    vrmGrid.Set_Resolution( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_refineLevelVal_valueChanged(int arg1)
{
    vrmGrid.Set_RefinementLevel( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_smoothTypeVal_currentTextChanged(const QString &arg1)
{
    vrmGrid.Set_SmoothType( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_smoothIterVal_valueChanged(int arg1)
{
    vrmGrid.Set_SmoothIterations( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_smoothDistVal_valueChanged(int arg1)
{
    vrmGrid.Set_TransSmoothDist( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_lonShiftVal_valueChanged(double arg1)
{
    vrmGrid.Set_GridLonShift( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_XrotateVal_valueChanged(double arg1)
{
    vrmGrid.Set_GridXRotate( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_YrotateVal_valueChanged(double arg1)
{
    vrmGrid.Set_GridYRotate( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_tessellationsVal_valueChanged(int arg1)
{
    vrmGrid.Set_Tessellations( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_subCellResolutionVal_valueChanged(int arg1)
{
    vrmGrid.Set_SubCellResolution( arg1 );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_reverseOrientation_clicked(bool checked)
{
    vrmGrid.Set_ReverseOrientation( checked );
    ui->calcVarMesh->setEnabled(true);
    gridInfo.Set_Enabled(false);
    gridInfoUpdate();
}

void MainWindow::on_calcVarMesh_clicked()
{
    bool Status = vrmGrid.calc_VarMesh();

    if( Status ) {
        ui->calcVarMesh->setEnabled(false);
        gridInfo.Set_Enabled(true);
        gridInfoUpdate();
    } else {
        ui->calcVarMesh->setEnabled(true);
        gridInfo.Set_Enabled(false);
        gridInfoUpdate();
    }

    vrmWin->updateRefineCubeImg( &vrmGrid);
    vrmWin->update_scene( &vrmGrid);
    update();
}

//============================================================
// User Interface for changes to DISPLAY settings
//============================================================
void MainWindow::on_backgroundFadeVal_valueChanged(int value)
{
    double fval = static_cast<double>(value)/100.;
    if(fval != vrmWin->BackgroundFilterAlpha()) {
        vrmWin->setBackgroundFilterAlpha(fval);
    }
}

void MainWindow::on_backgroundImgVal_currentTextChanged(const QString &arg1)
{
    if( arg1 != vrmWin->BackgroundImgVal()) {
        vrmWin->setBackgroundImgVal( arg1);
    }
}

void MainWindow::on_referenceMapOnOff_clicked(bool checked)
{
    if( checked != vrmWin->ReferenceMapVisible()) {
        vrmWin->setReferenceMapVisible(checked);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_referenceMapFadeVal_valueChanged(int value)
{
    double fval = static_cast<double>(value)/100.;
    if(fval != vrmWin->ReferenceMapAlpha()) {
        vrmWin->setReferenceMapAlpha(fval);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_referenceMapColorMod_clicked()
{

}

void MainWindow::on_refineMapOnOff_clicked(bool checked)
{
    if(ui->editActionsGroup->isEnabled()) {
        if( checked != vrmWin->RefineMapEditVisible()) {
            vrmWin->setRefineMapEditVisible(checked);
            vrmWin->update_scene( &vrmGrid);
        }
    } else {
        if( checked != vrmWin->RefineMapVisible()) {
            vrmWin->setRefineMapVisible(checked);
            vrmWin->update_scene( &vrmGrid);
        }
    }
}

void MainWindow::on_refineMapFadeVal_valueChanged(int value)
{
    double fval = static_cast<double>(value)/100.;
    if(fval != vrmWin->RefineMapAlpha()) {
        vrmWin->setRefineMapAlpha(fval);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_refineMapColorMod_clicked(bool checked)
{

}

void MainWindow::on_cubeMapOnOff_clicked(bool checked)
{
    if( checked != vrmWin->RefineCubeVisible()) {
        vrmWin->setRefineCubeVisible(checked);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_cubeMapFadeVal_valueChanged(int value)
{
    double fval = static_cast<double>(value)/100.;
    if(fval != vrmWin->RefineCubeAlpha()) {
        vrmWin->setRefineCubeAlpha(fval);
        vrmWin->adjustRefineCubeImgAlpha(fval);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_cubeGridOnOff_clicked(bool checked)
{
    if( checked != vrmWin->DisplayGridVisible()) {
        vrmWin->setDisplayGridVisible(checked);
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_cubeGridVal_currentTextChanged(const QString &arg1)
{
    if( arg1 != vrmWin->DisplayGridName()) {
        vrmWin->setDisplayGridName( arg1);
        vrmWin->update_scene( &vrmGrid);
    }
}

//============================================================
// User Interface for changes to EDIT settings
//============================================================
void MainWindow::on_beginEditMode_toggled(bool checked)
{
    if(checked) {

        // Disable VRM options while in edit mode
        //----------------------------------------
        ui->VrmOpts->setEnabled(false);
        ui->editActionsGroup->setEnabled(true);

        // Copy current refineMap to refineMapEdit
        //-----------------------------------------
        vrmGrid.refineMapEdit.resize(vrmGrid.refineMap.nRefLon, vrmGrid.refineMap.nRefLat);
        vrmGrid.refineMapEdit.Min = vrmGrid.refineMap.Min;
        vrmGrid.refineMapEdit.Max = vrmGrid.refineMap.Max;

        for(int ii=0; ii < vrmGrid.refineMapEdit.nRefLon; ii++) {
        for(int jj=0; jj < vrmGrid.refineMapEdit.nRefLat; jj++) {
            vrmGrid.refineMapEdit.val[ii][jj] = vrmGrid.refineMap.val[ii][jj];
        }
        }

        // Save current display options, make refineMapEdit Active
        //---------------------------------------------------------
        Edit_refineMapOnOff = ui->refineMapOnOff->isChecked();
        vrmWin->setRefineMapVisible    (false);
        vrmWin->setRefineMapEditVisible(true);
        ui->refineMapOnOff->setChecked (true);

        // Set Editor Items, depending on which tool box tab
        // was last active.
        //---------------------------------------------------
        int index = ui->editActionsGroup->currentIndex();
        on_editActionsGroup_currentChanged(index);

        // update display.
        //----------------------------------------------
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_endEditMode_toggled(bool checked)
{
    if(checked) {

        // The user needs to decide whether to keep the edited map or discard it.
        //-----------------------------------------------------------------------
        QMessageBox::StandardButton keepEdits;
        keepEdits = QMessageBox::question(this, "Exiting EDIT Mode", " Do you wish to keep the \n"
                                                         "modified refineMap values?",
                                                  QMessageBox::Yes | QMessageBox::No);


        // Copy refineMapEdit to refineMap if we are going to keep the edits.
        // Enable the calcVarMesh button.
        //-----------------------------------------------------------------------
        if(keepEdits == QMessageBox::Yes) {
            vrmGrid.refineMap.resize(vrmGrid.refineMapEdit.nRefLon, vrmGrid.refineMapEdit.nRefLat);
            vrmGrid.refineMap.Min = vrmGrid.refineMapEdit.Min;
            vrmGrid.refineMap.Max = vrmGrid.refineMapEdit.Max;

            for(int ii=0; ii < vrmGrid.refineMap.nRefLon; ii++) {
            for(int jj=0; jj < vrmGrid.refineMap.nRefLat; jj++) {
                vrmGrid.refineMap.val[ii][jj] = vrmGrid.refineMapEdit.val[ii][jj];
            }
            }
            ui->calcVarMesh->setEnabled(true);
            gridInfo.Set_Enabled(false);
            gridInfoUpdate();
        }

        // Deactivate Editor Items
        //--------------------------
        vrmWin->rectEdit_Item->setVisible(false);
        vrmWin->rectEdit_Inside->setVisible(false);
        vrmWin->rectEdit_Outside->setVisible(false);
        vrmWin->rectEdit_Item->Set_InsideRectVisible(false);
        vrmWin->rectEdit_Item->Set_OutsideRectVisible(false);

        vrmWin->Node1->setVisible(false);
        vrmWin->Node2->setVisible(false);
        vrmWin->Node3->setVisible(false);
        vrmWin->Node4->setVisible(false);
        vrmWin->Node5->setVisible(false);
        vrmWin->Node6->setVisible(false);
        vrmWin->Node7->setVisible(false);
        vrmWin->Node8->setVisible(false);
        vrmWin->polyEdit_Item->setVisible(false);
     // vrmWin->ellipseEdit_Item(false);

        // Deactivate refineMapEdit, restore refineMap display settings
        //--------------------------------------------------------------
        vrmWin->setRefineMapEditVisible(false);
        ui->refineMapOnOff->setChecked(Edit_refineMapOnOff);
        vrmWin->setRefineMapVisible   (Edit_refineMapOnOff);

        // Enable VRM options tab
        //---------------------------
        ui->VrmOpts->setEnabled(true);
        ui->editActionsGroup->setEnabled(false);

        // update display.
        //----------------------------------------------
        vrmWin->update_scene( &vrmGrid);
    }
}

void MainWindow::on_editActionsGroup_currentChanged(int index)
{
    switch(index) {
    case 0: // Grid Rest Active

        // Make sure everthing else is deactivated
        //------------------------------------------
        vrmWin->rectEdit_Item->setVisible(false);
        vrmWin->rectEdit_Inside->setVisible(false);
        vrmWin->rectEdit_Outside->setVisible(false);
        vrmWin->rectEdit_Item->Set_InsideRectVisible(false);
        vrmWin->rectEdit_Item->Set_OutsideRectVisible(false);

        vrmWin->Node1->setVisible(false);
        vrmWin->Node2->setVisible(false);
        vrmWin->Node3->setVisible(false);
        vrmWin->Node4->setVisible(false);
        vrmWin->Node5->setVisible(false);
        vrmWin->Node6->setVisible(false);
        vrmWin->Node7->setVisible(false);
        vrmWin->Node8->setVisible(false);
        vrmWin->polyEdit_Item->setVisible(false);

     // vrmWin->ellipseEdit_Item(false);

        update();
        break;
    case 1: // Grid Smoothing Active

        // Make sure everthing else is deactivated
        //------------------------------------------
        vrmWin->rectEdit_Item->setVisible(false);
        vrmWin->rectEdit_Inside->setVisible(false);
        vrmWin->rectEdit_Outside->setVisible(false);
        vrmWin->rectEdit_Item->Set_InsideRectVisible(false);
        vrmWin->rectEdit_Item->Set_OutsideRectVisible(false);

        vrmWin->Node1->setVisible(false);
        vrmWin->Node2->setVisible(false);
        vrmWin->Node3->setVisible(false);
        vrmWin->Node4->setVisible(false);
        vrmWin->Node5->setVisible(false);
        vrmWin->Node6->setVisible(false);
        vrmWin->Node7->setVisible(false);
        vrmWin->Node8->setVisible(false);
        vrmWin->polyEdit_Item->setVisible(false);

     // vrmWin->ellipseEdit_Item(false);

        update();
        break;
    case 2: // Polygon Edit Active

        // Activate an initialize Polygon Editor
        //-------------------------------------------
        vrmWin->Node1->setVisible(true);
        vrmWin->Node2->setVisible(true);
        vrmWin->Node3->setVisible(true);
        vrmWin->Node4->setVisible(true);
        vrmWin->Node5->setVisible(true);
        vrmWin->Node6->setVisible(true);
        vrmWin->Node7->setVisible(true);
        vrmWin->Node8->setVisible(true);
        vrmWin->polyEdit_Item->setVisible(true);

        // Make sure everthing else is deactivated
        //------------------------------------------
        vrmWin->rectEdit_Item->setVisible(false);
        vrmWin->rectEdit_Inside->setVisible(false);
        vrmWin->rectEdit_Outside->setVisible(false);
        vrmWin->rectEdit_Item->Set_InsideRectVisible(false);
        vrmWin->rectEdit_Item->Set_OutsideRectVisible(false);


     // vrmWin->ellipseEdit_Item(false);

        update();
        break;
    case 3: // Rectangle Edit Active

        // Activate an initialize Reactangle Editor
        //-------------------------------------------
        vrmWin->rectEdit_Item->setVisible(true);

        if(vrmWin->rectEdit_Item->Is_InsideRectVisible()) {
            vrmWin->rectEdit_Inside->setRect(vrmWin->rectEdit_Item->Get_InsideRect());
            vrmWin->rectEdit_Inside->setVisible(true);
        }

        if(vrmWin->rectEdit_Item->Is_OutsideRectVisible()) {
            vrmWin->rectEdit_Outside->setRect(vrmWin->rectEdit_Item->Get_OutsideRect());
            vrmWin->rectEdit_Outside->setVisible(true);
        }


        // Make sure everthing else is deactivated
        //------------------------------------------
        vrmWin->Node1->setVisible(false);
        vrmWin->Node2->setVisible(false);
        vrmWin->Node3->setVisible(false);
        vrmWin->Node4->setVisible(false);
        vrmWin->Node5->setVisible(false);
        vrmWin->Node6->setVisible(false);
        vrmWin->Node7->setVisible(false);
        vrmWin->Node8->setVisible(false);
        vrmWin->polyEdit_Item->setVisible(false);

     // vrmWin->ellipseEdit_Item(false);

        update();
        break;
    case 4: // Ellipse Edit Avtive

        // Activate an initialize Ellipse Editor
        //-------------------------------------------
     // vrmWin->ellipseEdit_Item(true);

        // Make sure everthing else is deactivated
        //------------------------------------------
        vrmWin->rectEdit_Item->setVisible(false);
        vrmWin->rectEdit_Inside->setVisible(false);
        vrmWin->rectEdit_Outside->setVisible(false);
        vrmWin->rectEdit_Item->Set_InsideRectVisible(false);
        vrmWin->rectEdit_Item->Set_OutsideRectVisible(false);


        vrmWin->Node1->setVisible(false);
        vrmWin->Node2->setVisible(false);
        vrmWin->Node3->setVisible(false);
        vrmWin->Node4->setVisible(false);
        vrmWin->Node5->setVisible(false);
        vrmWin->Node6->setVisible(false);
        vrmWin->Node7->setVisible(false);
        vrmWin->Node8->setVisible(false);
        vrmWin->polyEdit_Item->setVisible(false);

        update();
        break;

    default:

        update();
    }

}


  //===========================================
  // Reset Tab Interface
  //===========================================
void MainWindow::on_resetEditMap_clicked()
{
    // TODO:: Need to add check for EDIT state active.
    //---------------------------------------------------

    // Zero the refineMapEdit values, then update the scene
    //-------------------------------------------------------
    for(int ii=0; ii < vrmGrid.refineMapEdit.nRefLon; ii++) {
    for(int jj=0; jj < vrmGrid.refineMapEdit.nRefLat; jj++) {
        vrmGrid.refineMapEdit.val[ii][jj] = 0.0;
    }
    }
    vrmWin->update_scene( &vrmGrid);
}


  //===========================================
  // Smoothing Tab Interface
  //===========================================
void MainWindow::on_kernalSizeVal_valueChanged(int arg1)
{
    vrmWin->setSmoothKernel(arg1);
}

void MainWindow::on_persistenceVal_valueChanged(double arg1)
{
    vrmWin->setSmoothPersistence(arg1);
}

void MainWindow::on_iterationsVal_valueChanged(int arg1)
{
    vrmWin->setSmoothIterNum(arg1);
}

void MainWindow::on_applySmoothing_clicked()
{
    vrmWin->applySmoothing( &vrmGrid);
    vrmWin->update_scene( &vrmGrid);
}


  //===========================================
  // Polygon Edit Tab Interface
  //===========================================
void MainWindow::on_polyEditVal_valueChanged(double arg1)
{
     vrmWin->polyEdit_Item->setVal(arg1);
}

void MainWindow::on_resetPolyButton_clicked()
{
    double x0 = vrmView->getLon0();
    double y0 = vrmView->getLat0();
    double dx = 0.10*(vrmView->getDltLon()/vrmView->getZoomScl());
    double dy = 0.05*(vrmView->getDltLat()/vrmView->getZoomScl());

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

    vrmWin->Node1->updateNode(MyPoly[0]);
    vrmWin->Node2->updateNode(MyPoly[1]);
    vrmWin->Node3->updateNode(MyPoly[2]);
    vrmWin->Node4->updateNode(MyPoly[3]);
    vrmWin->Node5->updateNode(MyPoly[4]);
    vrmWin->Node6->updateNode(MyPoly[5]);
    vrmWin->Node7->updateNode(MyPoly[6]);
    vrmWin->Node8->updateNode(MyPoly[7]);
    vrmWin->polyEdit_Item->UpdatePolygon(MyPoly);
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_applyPolyButton_clicked()
{
    QPolygonF *MyPoly;
    MyPoly = new QPolygonF(vrmWin->polyEdit_Item->mapToScene(vrmWin->polyEdit_Item->polygon()));
    double MyVal = vrmWin->polyEdit_Item->Val();
    vrmGrid.SetPolyVal(MyPoly,MyVal);
    vrmWin->update_scene( &vrmGrid);
}


  //===========================================
  // Rect Edit Tab Interface
  //===========================================
void MainWindow::on_rectEditVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_FillVal( arg1 );
}

void MainWindow::on_rectFillMode_currentTextChanged(const QString &arg1)
{
    vrmWin->rectEdit_Item->Set_FillOpt( arg1 );
}

void MainWindow::on_resetRectButton_clicked()
{
    double x0 = vrmView->getLon0();
    double y0 = vrmView->getLat0();
    double dx = 0.20*(vrmView->getDltLon()/vrmView->getZoomScl());
    double dy = 0.10*(vrmView->getDltLat()/vrmView->getZoomScl());

    vrmWin->rectEdit_Item->setPos(0.,0.);
    vrmWin->rectEdit_Item->ResetRect(x0, y0, dx, dy);
    vrmWin->rectEdit_Item->setRect(vrmWin->rectEdit_Item->Get_Rect());
    ui->rectDltLatVal->setValue(vrmWin->rectEdit_Item->Get_DltLat());
    ui->rectDltLonVal->setValue(vrmWin->rectEdit_Item->Get_DltLon());
    ui->rectLonMinVal->setValue(vrmWin->rectEdit_Item->Get_LonMin());
    ui->rectLonMaxVal->setValue(vrmWin->rectEdit_Item->Get_LonMax());
    ui->rectLatMinVal->setValue(vrmWin->rectEdit_Item->Get_LatMin());
    ui->rectLatMaxVal->setValue(vrmWin->rectEdit_Item->Get_LatMax());
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectDltLonVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_DltLon( arg1 );
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectDltLatVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_DltLat( arg1 );
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectLonMinVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_LonMin( arg1 );
    vrmWin->rectEdit_Item->setRect(vrmWin->rectEdit_Item->Get_Rect());
 //   ui->rectLonMinVal->setValue(vrmWin->rectEdit_Item->Get_LonMin());
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectLonMaxVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_LonMax( arg1 );
    vrmWin->rectEdit_Item->setRect(vrmWin->rectEdit_Item->Get_Rect());
//    ui->rectLonMaxVal->setValue(vrmWin->rectEdit_Item->Get_LonMax());
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectLatMinVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_LatMin( arg1 );
    vrmWin->rectEdit_Item->setRect(vrmWin->rectEdit_Item->Get_Rect());
 //   ui->rectLatMinVal->setValue(vrmWin->rectEdit_Item->Get_LatMin());
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_rectLatMaxVal_valueChanged(double arg1)
{
    vrmWin->rectEdit_Item->Set_LatMax( arg1 );
    vrmWin->rectEdit_Item->setRect(vrmWin->rectEdit_Item->Get_Rect());
 //   ui->rectLatMaxVal->setValue(vrmWin->rectEdit_Item->Get_LatMax());
    vrmWin->update_scene( &vrmGrid);
}

void MainWindow::on_applyRectButton_clicked()
{
    vrmGrid.SetRectVal(vrmWin->rectEdit_Item->Get_Lon0()    ,vrmWin->rectEdit_Item->Get_Lat0()    ,
                       vrmWin->rectEdit_Item->Get_LonWidth(),vrmWin->rectEdit_Item->Get_LatWidth(),
                       vrmWin->rectEdit_Item->Get_DltLon()  ,vrmWin->rectEdit_Item->Get_DltLat()  ,
                       vrmWin->rectEdit_Item->Get_FillVal() ,vrmWin->rectEdit_Item->Get_FillOpt() );
    vrmWin->update_scene( &vrmGrid);
}
