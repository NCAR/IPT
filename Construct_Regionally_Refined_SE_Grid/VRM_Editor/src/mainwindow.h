#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ReadWrite.h"

#include <QMainWindow>
#include <QGraphicsScene>
#include "lonbar.h"
#include "latbar.h"
#include "vrmview.h"
#include "vrmwin.h"
#include "vrmgrid.h"
#include "polyedititem.h"
#include "rectedititem.h"
#include "gridinfoitem.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void initLonScene(int I_htext);
    void initLatScene(int I_htext);

    VrmView        *vrmView;
    LonBar         *lonbar;
    LatBar         *latbar;

    VrmWin         *vrmWin;
    QGraphicsScene *LonScene;
    QGraphicsScene *LatScene;
    gridInfoItem    gridInfo;

    VrmGrid vrmGrid;

protected:
    bool eventFilter(QObject *watched, QEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);
    virtual void wheelEvent( QGraphicsSceneEvent *event);

public slots:
    void do_SceneUpdate();
    void do_rectUIupdate();
    void gridInfoUpdate();

private slots:
    void on_setViewWindows(int HSIZE, int VSIZE, int LATSIZE, int LONSIZE,
                           double lon0, double lat0, double slon, double slat);

private slots:
    void on_actionSave_VRM_File_triggered();

    void on_actionRead_VRM_File_triggered();

    void on_actionWrite_Exodus_File_triggered();

    void on_actionRead_Refinement_Map_triggered();

    void on_actionSave_Refinement_Map_triggered();

    void on_actionRead_Reference_Map_triggered();

     void on_actionQuit_triggered();



    void on_refineGridVal_currentTextChanged(const QString &arg1);

    void on_refineTypeVal_currentTextChanged(const QString &arg1);

    void on_refineBaseVal_valueChanged(int arg1);

    void on_refineLevelVal_valueChanged(int arg1);

    void on_smoothTypeVal_currentTextChanged(const QString &arg1);

    void on_smoothIterVal_valueChanged(int arg1);

    void on_smoothDistVal_valueChanged(int arg1);

    void on_lonShiftVal_valueChanged(double arg1);

    void on_XrotateVal_valueChanged(double arg1);

    void on_YrotateVal_valueChanged(double arg1);

    void on_tessellationsVal_valueChanged(int arg1);

    void on_subCellResolutionVal_valueChanged(int arg1);

    void on_reverseOrientation_clicked(bool checked);

    void on_calcVarMesh_clicked();



    void on_backgroundFadeVal_valueChanged(int value);

    void on_backgroundImgVal_currentTextChanged(const QString &arg1);

    void on_referenceMapOnOff_clicked(bool checked);

    void on_referenceMapFadeVal_valueChanged(int value);

    void on_referenceMapColorMod_clicked();

    void on_refineMapOnOff_clicked(bool checked);

    void on_refineMapFadeVal_valueChanged(int value);

    void on_refineMapColorMod_clicked(bool checked);

    void on_cubeMapOnOff_clicked(bool checked);

    void on_cubeMapFadeVal_valueChanged(int value);

    void on_cubeGridOnOff_clicked(bool checked);

    void on_cubeGridVal_currentTextChanged(const QString &arg1);



    void on_beginEditMode_toggled(bool checked);

    void on_endEditMode_toggled(bool checked);

    void on_editActionsGroup_currentChanged(int index);


    void on_resetEditMap_clicked();


    void on_kernalSizeVal_valueChanged(int arg1);

    void on_persistenceVal_valueChanged(double arg1);

    void on_iterationsVal_valueChanged(int arg1);

    void on_applySmoothing_clicked();    


    void on_polyEditVal_valueChanged(double arg1);

    void on_resetPolyButton_clicked();

    void on_applyPolyButton_clicked();


    void on_rectEditVal_valueChanged(double arg1);

    void on_rectFillMode_currentTextChanged(const QString &arg1);

    void on_resetRectButton_clicked();

    void on_rectDltLonVal_valueChanged(double arg1);

    void on_rectDltLatVal_valueChanged(double arg1);

    void on_rectLonMinVal_valueChanged(double arg1);

    void on_rectLonMaxVal_valueChanged(double arg1);

    void on_rectLatMinVal_valueChanged(double arg1);

    void on_rectLatMaxVal_valueChanged(double arg1);

    void on_applyRectButton_clicked();

private:
    Ui::MainWindow *ui;


    bool    Edit_refineMapOnOff;

};

#endif // MAINWINDOW_H
