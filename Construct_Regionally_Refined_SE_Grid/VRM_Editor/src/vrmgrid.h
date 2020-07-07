#ifndef VRMGRID_H
#define VRMGRID_H

#include <QObject>
#include <QString>
#include <QPolygonF>
#include <cmath>

#include "SQuadGen/GridElements.h"
#include "SQuadGen/CubedSphereGrid.h"
#include "SQuadGen/IcosahedralFlagGrid.h"
#include "SQuadGen/Tessellate.h"
#include "SQuadGen/SpringDynamics.h"
#include "SQuadGen/RefinementTemplateCUBIT.h"
#include "SQuadGen/RefinementTemplateLOWCONN.h"
#include "SQuadGen/RefinementTemplateLOWCONNOLD.h"
#include "SQuadGen/MathHelper.h"

#include "RefinementMap.h"
#include "RefinementCube.h"
#include "CubeGrid.h"
#include "RefineCubeGrid.h"
#include "NeighborAdjustments.h"


class VrmGrid : public QObject
{
    Q_OBJECT
public:
    explicit VrmGrid(QObject *parent = nullptr);
    bool    calc_VarMesh();

    // Functions to Get/Set Grid Processin gOptions
    //-----------------------------------------------
    QString Get_RefineType();
    void    Set_RefineType(const QString & I_RefineType);
    QString Get_SmoothType();
    void    Set_SmoothType(const QString & I_SmoothType);
    QString Get_GridType();
    void    Set_GridType( const QString & I_GridType);
    int     Get_Resolution();
    void    Set_Resolution(int I_Resolution);
    int     Get_RefinementLevel();
    void    Set_RefinementLevel(int I_RefinementLevel);
    int     Get_SmoothIterations();
    void    Set_SmoothIterations(int I_SmoothIterations);
    int     Get_TransSmoothDist();
    void    Set_TransSmoothDist(int I_TransSmoothDist);
    double  Get_GridLonShift();
    void    Set_GridLonShift(double I_GridLonShift);
    double  Get_GridXRotate();
    void    Set_GridXRotate(double I_GridXRotate);
    double  Get_GridYRotate();
    void    Set_GridYRotate(double I_GridYRotate);
    int     Get_Tessellations();
    void    Set_Tessellations(int I_Tessellations);
    int     Get_SubCellResolution();
    void    Set_SubCellResolution(int I_SubCellResolution);
    bool    Get_ReverseOrientation();
    void    Set_ReverseOrientation(bool I_ReverseOrientation);

    void SetRectVal(double I_Lon0    , double I_Lat0    ,
                    double I_LonWidth, double I_LatWidth,
                    double I_DltLon  , double I_DltLat  , double I_MyVal, const QString &I_FillMode);

public slots:
    void SetPolyVal( QPolygonF *I_poly, double I_val);


public:

    RefinementMap  refineMap;
    RefinementMap  refineMapRef;
    RefinementMap  refineMapEdit;
    RefinementCube refineCube;
    CubeGrid       refineGrid;
    CubeGrid       refineGrid_copy;
    CubeGrid       baseGrid;

private:
    // Grid Processing Options
    //--------------------------
    QString MyRefineType;
    QString MySmoothType;
    QString MyGridType;
    int     MyResolution;
    int     MyRefinementLevel;
    int     MySmoothIterations;
    int     MyTransSmoothDist;
    double  MyGridLonShift;
    double  MyGridXRotate;
    double  MyGridYRotate;
    int     MyTessellations;
    int     MySubCellResolution;
    bool    MyReverseOrientation;

};

#endif // VRMGRID_H
