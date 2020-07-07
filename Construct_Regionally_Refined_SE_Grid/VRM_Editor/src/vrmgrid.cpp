#include "vrmgrid.h"

VrmGrid::VrmGrid(QObject *parent)
{

}

bool VrmGrid::calc_VarMesh()
//
// calc_VarMesh: For the current refinement map and processing options,
//               calcualte the ...
//
//=====================================================================
{

    NodeVector vecNodes;
    FaceVector vecFaces;
    int nResolution = MyResolution;


    // Check for odd resolution for LOWCONN refile types
    //----------------------------------------------------
    if ((MyRefineType == "LOWCONN") || (MyRefineType == "LOWCONNold") ) {
        if (nResolution % 2 == 1) {
            std::cout << " ERROR:: RefineType = "          << MyRefineType.toStdString() << std::endl;
            std::cout << " ERROR:: nResolution is ODD! = " << nResolution                << std::endl;
        } else {
            nResolution /= 2;
        }
    }

 //DIAG    std::vector<Node>().swap(vecNodes);
 //DIAG    std::vector<Face>().swap(vecFaces);

    // Generate Base Grid
    //----------------------------
    if(MyGridType == "CubeSquared") {
//        std::cout << " GridType = CubeSquared" << std::endl;
//        std::cout << " nResolution = "<< nResolution << std::endl;
        GenerateCubedSphere( nResolution, vecNodes, vecFaces);
    } else if (MyGridType == "Icosahedral") {
//        std::cout << " GridType = Icosahedral" << std::endl;
//        std::cout << " nResolution = "<< nResolution << std::endl;
        GenerateIcosahedralQuadGrid( nResolution, vecNodes, vecFaces);
    } else {
        std::cout << " ERROR:: GridType = UNKNOWN" << std::endl;
    }

    // Rotate that Grid around the X-Axis
    //------------------------------------
    if( MyGridXRotate != 0.0) {
        double dCosTheta = cos(MyGridXRotate*M_PI/180.0);
        double dSinTheta = sin(MyGridXRotate*M_PI/180.0);
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            double dTempY = vecNodes[i].y;
            double dTempZ = vecNodes[i].z;
            vecNodes[i].y = dCosTheta*dTempY - dSinTheta*dTempZ;
            vecNodes[i].z = dSinTheta*dTempY + dCosTheta*dTempZ;
        }
    }

    // Rotate that Grid around the Y-Axis
    //------------------------------------
    if( MyGridYRotate != 0.0) {
        double dCosTheta = cos(MyGridYRotate*M_PI/180.0);
        double dSinTheta = sin(MyGridYRotate*M_PI/180.0);
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            double dTempX = vecNodes[i].x;
            double dTempZ = vecNodes[i].z;
            vecNodes[i].x = dCosTheta*dTempX - dSinTheta*dTempZ;
            vecNodes[i].z = dSinTheta*dTempX + dCosTheta*dTempZ;
        }
    }

    // Add Longitude Shift to that Grid
    //-----------------------------------
    if (MyGridLonShift != 0.0) {
        double dLon;
        double dLat;
        for( unsigned int i=0; i < vecNodes.size(); i++) {
            dLat = asin(vecNodes[i].z);
            dLon = atan2(vecNodes[i].y,vecNodes[i].x);
            dLon += MyGridLonShift * M_PI/180.0;
            vecNodes[i].x = cos(dLat)*cos(dLon);
            vecNodes[i].y = cos(dLat)*sin(dLon);
        }
    }

    // Store Base grid values
    //-------------------------
    baseGrid.loadGrid(vecNodes, vecFaces);
//    std::cout << " baseGrid nodes set N= "<< vecNodes.size() << " : " << baseGrid.Nodes.size() << std::endl;
//    std::cout << " baseGrid faces set N= "<< vecFaces.size() << " : " << baseGrid.Faces.size() << std::endl;
//    std::cout << " baseGrid grid  set N= " << baseGrid.Grid.size() << std::endl;

    // Initialize the reference Cube values
    //---------------------------------------
    refineCube.resize(nResolution, MyRefinementLevel, MyGridLonShift, MyGridXRotate, MyGridYRotate);
    refineCube.loadCubeVals(refineMap);
    refineCube.Normalize();
//    refineCube.write( "refine_map.dat");

    // Perform Refinement if requested
    //----------------------------------
    if (MyRefinementLevel > 0) {
        if(MyRefineType == "CUBIT") {
            RefinementTemplateCUBIT RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else if( MyRefineType == "LOWCONN") {
            RefinementTemplateLOWCONN RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else if ( MyRefineType == "LOWCONNold") {
            RefinementTemplateLOWCONNOLD RefTemp;
            RefineGrid( vecNodes, vecFaces, RefTemp, refineCube);
        } else {
            // *** ERROR ***
        }
    }

    // Tessellate that Grid
    //------------------------
    if ((MyTessellations < 0) || (MyTessellations > 100)) {
        std::cout << " ERROR:: MyTesselations out of Range = "<< MyTessellations << std::endl;
    } else {
        for (int n=0; n<MyTessellations; n++) {
            Tessellate( vecNodes, vecFaces);
        }
    }

    // Add Sub-Cell Resolution to that Grid
    //---------------------------------------
    if (MySubCellResolution < 0) {
        std::cout << " ERROR:: MySubCellResolution out of Range = "<< MySubCellResolution << std::endl;
    } else if (MySubCellResolution == 0){
        // Do Nothing
    } else {
        RefineEverything( vecNodes, vecFaces, MySubCellResolution+1 );
    }

    // Reverse Orientation of Faces
    //-------------------------------
    if ((MyReverseOrientation) && (MyGridType != "CubeSquared")) {
        ReverseFaceOrientation( vecFaces);
    } else if ((!MyReverseOrientation) && (MyGridType == "CubeSquared")) {
        ReverseFaceOrientation( vecFaces);
    }


    // Smooth that Grid
    //-------------------
    if (MySmoothType == "NONE") {
        // No Smoothing
        //--------------
    } else if (MySmoothType == "SPRING") {
        SpringDynamics( vecNodes, vecFaces, MyTransSmoothDist, MySmoothIterations);
    } else if (MySmoothType == "PRESSURE") {
        PressureDynamics( vecNodes, vecFaces, MyTransSmoothDist, MySmoothIterations);
    } else if (MySmoothType == "DISTANCE") {
//        int NE_lo = MyResolution;
//        int NE_hi = MyResolution*IntPow(2,MyRefinementLevel);
//        std::cout << " NE_lo= "<< NE_lo << " NE_hi= " << NE_hi << std::endl;
//        EquializeDistances(refineMap, vecNodes, vecFaces, MyTransSmoothDist,
//                                                          MySmoothIterations, NE_lo, NE_hi);
    } else {
        std::cout << " ERROR:: Unknown MySmoothType = "<< MySmoothType.toStdString() << std::endl;
    }

    // Compute the Mesh Quality
    //--------------------------


    // Store final Refined grid values
    //---------------------------------
    refineGrid.loadGrid( vecNodes, vecFaces);
    refineGrid_copy = refineGrid;
//    std::cout << " refineGrid nodes set N= "<< vecNodes.size() << " : " << refineGrid.Nodes.size() << std::endl;
//    std::cout << " refineGrid faces set N= "<< vecFaces.size() << " : " << refineGrid.Faces.size() << std::endl;
//    std::cout << " refineGrid grid  set N= " << refineGrid.Grid.size() << std::endl;


    vecNodes.clear();
    vecFaces.clear();
    return true;

}

// Functions to Get/Set Grid Processing Options
//----------------------------------------------
QString VrmGrid::Get_RefineType()
{
    return MyRefineType;
}

void VrmGrid::Set_RefineType(const QString &I_RefineType)
{
    MyRefineType = I_RefineType;
}

QString VrmGrid::Get_SmoothType()
{
    return MySmoothType;
}

void VrmGrid::Set_SmoothType(const QString &I_SmoothType)
{
    MySmoothType = I_SmoothType;
}

QString VrmGrid::Get_GridType()
{
    return MyGridType;
}

void VrmGrid::Set_GridType(const QString &I_GridType)
{
    MyGridType = I_GridType;
}

int VrmGrid::Get_Resolution()
{
    return MyResolution;
}

void VrmGrid::Set_Resolution(int I_Resolution)
{
    MyResolution = I_Resolution;
}

int VrmGrid::Get_RefinementLevel()
{
    return MyRefinementLevel;
}

void VrmGrid::Set_RefinementLevel(int I_RefinementLevel)
{
    MyRefinementLevel = I_RefinementLevel;
}

int VrmGrid::Get_SmoothIterations()
{
    return MySmoothIterations;
}

void VrmGrid::Set_SmoothIterations(int I_SmoothIterations)
{
    MySmoothIterations = I_SmoothIterations;
}

int VrmGrid::Get_TransSmoothDist()
{
    return MyTransSmoothDist;
}

void VrmGrid::Set_TransSmoothDist(int I_TransSmoothDist)
{
    MyTransSmoothDist = I_TransSmoothDist;
}

double VrmGrid::Get_GridLonShift()
{
    return MyGridLonShift;
}

void VrmGrid::Set_GridLonShift(double I_GridLonShift)
{
    MyGridLonShift = I_GridLonShift;
}

double VrmGrid::Get_GridXRotate()
{
    return MyGridXRotate;
}

void VrmGrid::Set_GridXRotate(double I_GridXRotate)
{
    MyGridXRotate = I_GridXRotate;
}

double VrmGrid::Get_GridYRotate()
{
    return MyGridYRotate;
}

void VrmGrid::Set_GridYRotate(double I_GridYRotate)
{
    MyGridYRotate = I_GridYRotate;
}

int VrmGrid::Get_Tessellations()
{
    return MyTessellations;
}

void VrmGrid::Set_Tessellations(int I_Tessellations)
{
    MyTessellations = I_Tessellations;
}

int VrmGrid::Get_SubCellResolution()
{
    return MySubCellResolution;
}

void VrmGrid::Set_SubCellResolution(int I_SubCellResolution)
{
    MySubCellResolution = I_SubCellResolution;
}

bool VrmGrid::Get_ReverseOrientation()
{
    return MyReverseOrientation;
}

void VrmGrid::Set_ReverseOrientation(bool I_ReverseOrientation)
{
    MyReverseOrientation = I_ReverseOrientation;
}

void VrmGrid::SetPolyVal(QPolygonF *I_poly, double I_val)
{
    int     iLon0;
    int     iLat0;
    QRectF  bounds = I_poly->boundingRect();
    QPointF center = bounds.center();
    double  dLon0 = center.x();
    double  dLat0 = center.y();
    double  dLon;
    double  dLat;

    refineMapEdit.getInd(center.x(),center.y(), iLon0, iLat0);

    int nLon2 = ((1 + floor(bounds.width() /refineMapEdit.dltLon))/2);
    int nLat2 = ((1 + floor(bounds.height()/refineMapEdit.dltLat))/2);

    for(int ii=(-nLon2); ii <= nLon2; ii++) {
        for(int jj=(-nLat2); jj <= nLat2; jj++) {
            dLon = dLon0 + ii*refineMapEdit.dltLon;
            dLat = dLat0 + jj*refineMapEdit.dltLat;
            QPointF MyPoint(dLon,dLat);
            if(I_poly->containsPoint(MyPoint,Qt::WindingFill)) {
                int Lon = iLon0 + 1 + ii;
                int Lat = iLat0 + 1 + jj;
                if(Lat <             0            ) {Lat = 0;}
                if(Lat > (refineMapEdit.nRefLat-1)) {Lat = refineMapEdit.nRefLat-1;}
                if(Lon <             0            ) {Lon += refineMapEdit.nRefLon; }
                if(Lon > (refineMapEdit.nRefLon-1)) {Lon -= refineMapEdit.nRefLon; }

                if(refineMapEdit.val[Lon][Lat] < I_val) {   // DIAG for DEMO
                refineMapEdit.val[Lon][Lat] = I_val;
                } // DIAG
            }
        }
    }
}

void VrmGrid::SetRectVal(double I_Lon0    , double I_Lat0    ,
                         double I_LonWidth, double I_LatWidth,
                         double I_DltLon  , double I_DltLat  , double I_MyVal, const QString &I_FillMode)
{
    double Lon0 = I_Lon0;
    double Lat0 = I_Lat0;
    double LonMin;
    double LonMax;
    double LatMin;
    double LatMax;

    // Shift the window into 0-360.
    //----------------------------
    if(Lon0 <    0.) {Lon0 = Lon0 + 360.;}
    if(Lon0 >= 360.) {Lon0 = Lon0 - 360.;}

    LonMin =      - (I_LonWidth/2.);
    LonMax =        (I_LonWidth/2.);
    LatMin = Lat0 - (I_LatWidth/2.);
    LatMax = Lat0 + (I_LatWidth/2.);

    // Calculate min/max of RAW window function
    //------------------------------------------
    double Val1_p = (1. + tanh(( 180. - LonMin)/I_DltLon))/2.;
    double Val1_0 = (1. + tanh((   0. - LonMin)/I_DltLon))/2.;
    double Val1_n = (1. + tanh((-180. - LonMin)/I_DltLon))/2.;

    double Val2_p = (1. + tanh((LonMax - 180.)/I_DltLon))/2.;
    double Val2_0 = (1. + tanh((LonMax -   0.)/I_DltLon))/2.;
    double Val2_n = (1. + tanh((LonMax + 180.)/I_DltLon))/2.;

    double Val3_p = (1. + tanh(( 90. - LatMin)/I_DltLat))/2.;
    double Val3_0 = (1. + tanh((Lat0 - LatMin)/I_DltLat))/2.;
    double Val3_n = (1. + tanh((-90. - LatMin)/I_DltLat))/2.;

    double Val4_p = (1. + tanh((LatMax -  90.)/I_DltLat))/2.;
    double Val4_0 = (1. + tanh((LatMax - Lat0)/I_DltLat))/2.;
    double Val4_n = (1. + tanh((LatMax +  90.)/I_DltLat))/2.;

    double Wmax = Val1_0 * Val2_0 * Val3_0 * Val4_0;
    double Wmin = Val1_p * Val2_p * Val3_n * Val4_n;
    if( Wmin >Val1_p * Val2_p * Val3_p * Val4_p ) {Wmin = Val1_p * Val2_p * Val3_p * Val4_p;}
    if( Wmin >Val1_n * Val2_n * Val3_p * Val4_p ) {Wmin = Val1_n * Val2_n * Val3_p * Val4_p;}
    if( Wmin >Val1_n * Val2_n * Val3_n * Val4_n ) {Wmin = Val1_n * Val2_n * Val3_n * Val4_n;}

//    std::cout << " Wmin = " << Wmin << " Wmax= " << Wmax << std::endl;
//    std::cout << " I_Lon0=" << I_Lon0 << "LonMax=" << LonMax << " LonMin="<< LonMin << std::endl;
//    std::cout << " Lat0=" << Lat0 << " Lon0=" << Lon0 << "LonWidth=" << I_LonWidth << " LatWidth="<< I_LatWidth << std::endl;
//    std::cout << " Val1_0= " << Val1_0 << " Val2_0=" << Val2_0 << " Val3_0= " << Val3_0 << " Val4_0=" << Val4_0 << std::endl;

    // [min-max] values are mapped to [0-MyVal] range.
    //------------------------------------------------
    for(int iLon=0; iLon < refineMapEdit.nRefLon; iLon++) {
    for(int iLat=0; iLat < refineMapEdit.nRefLat; iLat++) {

        double rLon = refineMapEdit.RefLon[iLon][0];
        double rLat = refineMapEdit.RefLat[iLat][0];

        rLon = rLon - Lon0;
        if(rLon < -180.) {rLon = rLon + 360.;}
        if(rLon >  180.) {rLon = rLon - 360.;}

        double Wval =  ((1. + tanh((rLon-LonMin)/I_DltLon))/2.)
                      *((1. + tanh((LonMax-rLon)/I_DltLon))/2.)
                      *((1. + tanh((rLat-LatMin)/I_DltLat))/2.)
                      *((1. + tanh((LatMax-rLat)/I_DltLat))/2.);

        if(Wmax <= Wmin) {
            Wval = 1.0*I_MyVal;
        } else {
            Wval = I_MyVal*(Wval-Wmin)/(Wmax-Wmin);
        }

        if(I_FillMode == "Fill All") {
            refineMapEdit.val[iLon][iLat] = Wval;
        } else { // "Fill Max"
            if(refineMapEdit.val[iLon][iLat] < Wval) {
               refineMapEdit.val[iLon][iLat] = Wval;
            }
        }
    }
    }


}
