#include "CubeGrid.h"
#include <iostream>

CubeGrid::CubeGrid()
{
    NodeVector       ().swap(Nodes);
    FaceVector       ().swap(Faces);
    LonLatNodeVector ().swap(Grid);

    // Initialize Values
    //--------------------
    Num_Clockwise = 0;
    Num_Counter   = 0;
    Total_Area    = 1.0;
    Qarc          = 1.0;
    MinArc        =  999.;
    MaxArc        = -999.;
    MinArc_Lon    = -999.;
    MinArc_Lat    = -999.;
    MaxArc_Lon    = -999.;
    MaxArc_Lat    = -999.;
    Qang          = 1.0;
    MinAng        =  999.;
    MaxAng        = -999.;
    MinAng_Lon    = -999.;
    MinAng_Lat    = -999.;
    MaxAng_Lon    = -999.;
    MaxAng_Lat    = -999.;
    Area_Ratio    = 1.0;
    MinArea       =  999.;
    MaxArea       = -999.;
    MinArea_Lon   = -999.;
    MinArea_Lat   = -999.;
    MaxArea_Lon   =  999.;
    MaxArea_Lat   = -999.;
}

CubeGrid::CubeGrid(const CubeGrid &srcGrid)
{
    NodeVector       ().swap(Nodes);
    FaceVector       ().swap(Faces);
    LonLatNodeVector ().swap(Grid);

    for(unsigned int nn=0; nn < srcGrid.Nodes.size(); nn++) {
        Nodes.push_back(srcGrid.Nodes[nn]);
        Grid.push_back(srcGrid.Grid[nn]);
    }
    for(unsigned int nn=0; nn < srcGrid.Faces.size(); nn++) {
        Faces.push_back(srcGrid.Faces[nn]);
    }
    // Copy Values
    //--------------------
    Num_Clockwise = srcGrid.Num_Clockwise;
    Num_Counter   = srcGrid.Num_Counter;
    Total_Area    = srcGrid.Total_Area;
    Qarc          = srcGrid.Qarc;
    MinArc        = srcGrid.MinArc;
    MaxArc        = srcGrid.MaxArc;
    MinArc_Lon    = srcGrid.MinArc_Lon;
    MinArc_Lat    = srcGrid.MinArc_Lat;
    MaxArc_Lon    = srcGrid.MaxArc_Lon;
    MaxArc_Lat    = srcGrid.MaxArc_Lat;
    Qang          = srcGrid.Qang;
    MinAng        = srcGrid.MinAng;
    MaxAng        = srcGrid.MaxAng;
    MinAng_Lon    = srcGrid.MinAng_Lon;
    MinAng_Lat    = srcGrid.MinAng_Lat;
    MaxAng_Lon    = srcGrid.MaxAng_Lon;
    MaxAng_Lat    = srcGrid.MaxAng_Lat;
    Area_Ratio    = srcGrid.Area_Ratio;
    MinArea       = srcGrid.MinArea;
    MaxArea       = srcGrid.MaxArea;
    MinArea_Lon   = srcGrid.MinArea_Lon;
    MinArea_Lat   = srcGrid.MinArea_Lat;
    MaxArea_Lon   = srcGrid.MaxArea_Lon;
    MaxArea_Lat   = srcGrid.MaxArea_Lat;
}

CubeGrid &CubeGrid::operator =(const CubeGrid srcGrid)
{
    NodeVector       ().swap(this->Nodes);
    FaceVector       ().swap(this->Faces);
    LonLatNodeVector ().swap(this->Grid);

    for(unsigned int nn=0; nn < srcGrid.Nodes.size(); nn++){
        this->Nodes.push_back(srcGrid.Nodes[nn]);
        this->Grid.push_back(srcGrid.Grid[nn]);
    }
    for(unsigned int nn=0; nn < srcGrid.Faces.size(); nn++){
        this->Faces.push_back(srcGrid.Faces[nn]);
    }
    // Copy Values
    //--------------------
    this->Num_Clockwise = srcGrid.Num_Clockwise;
    this->Num_Counter   = srcGrid.Num_Counter;
    this->Total_Area    = srcGrid.Total_Area;
    this->Qarc          = srcGrid.Qarc;
    this->MinArc        = srcGrid.MinArc;
    this->MaxArc        = srcGrid.MaxArc;
    this->MinArc_Lon    = srcGrid.MinArc_Lon;
    this->MinArc_Lat    = srcGrid.MinArc_Lat;
    this->MaxArc_Lon    = srcGrid.MaxArc_Lon;
    this->MaxArc_Lat    = srcGrid.MaxArc_Lat;
    this->Qang          = srcGrid.Qang;
    this->MinAng        = srcGrid.MinAng;
    this->MaxAng        = srcGrid.MaxAng;
    this->MinAng_Lon    = srcGrid.MinAng_Lon;
    this->MinAng_Lat    = srcGrid.MinAng_Lat;
    this->MaxAng_Lon    = srcGrid.MaxAng_Lon;
    this->MaxAng_Lat    = srcGrid.MaxAng_Lat;
    this->Area_Ratio    = srcGrid.Area_Ratio;
    this->MinArea       = srcGrid.MinArea;
    this->MaxArea       = srcGrid.MaxArea;
    this->MinArea_Lon   = srcGrid.MinArea_Lon;
    this->MinArea_Lat   = srcGrid.MinArea_Lat;
    this->MaxArea_Lon   = srcGrid.MaxArea_Lon;
    this->MaxArea_Lat   = srcGrid.MaxArea_Lat;
    return *this;
}

CubeGrid::~CubeGrid()
{
    NodeVector       ().swap(Nodes);
    FaceVector       ().swap(Faces);
    LonLatNodeVector ().swap(Grid);
}

void CubeGrid::loadGrid(const NodeVector &srcNodes, const FaceVector &srcFaces)
{
    NodeVector       ().swap(Nodes);
    FaceVector       ().swap(Faces);
    LonLatNodeVector ().swap(Grid);

    // For Each Source Node, copy the node into the CubeGrid nodes.
    // Calculate the Lon,Lat values for the node and store the
    // gridpoints in cubeGrid.
    //--------------------------------------------------------
    for(unsigned int nn=0; nn < srcNodes.size(); nn++) {
        Nodes.push_back(srcNodes[nn]);

        double dLat;
        double dLon;
        dLat = (180./M_PI)*asin ( srcNodes[nn].z );
        if( (srcNodes[nn].x ==0.) && (srcNodes[nn].y ==0. )) {
            dLon = 0.;
        } else {
            dLon = (180./M_PI)*atan2( srcNodes[nn].y,srcNodes[nn].x );
        }
        if(dLon < 0.) { dLon += 360.; }
        Grid.push_back(LonLatNode(dLon,dLat));
    }

    // For each Source Face, copy the face into the CubeGris faces.
    //-------------------------------------------------------------
    for(unsigned int nn=0; nn < srcFaces.size(); nn++) {
        Faces.push_back(srcFaces[nn]);
    }

    // Calculate Grid Quality Parameters
    //-----------------------------------
    computeGridQuality();
}

void CubeGrid::computeGridQuality()
{

    // Initialize Values
    //--------------------
    Num_Clockwise = 0;
    Num_Counter   = 0;
    Total_Area    = 0.0;
    Qarc          = 0.0;
    MinArc        =  999.;
    MaxArc        = -999.;
    MinArc_Lon    = -999.;
    MinArc_Lat    = -999.;
    MaxArc_Lon    = -999.;
    MaxArc_Lat    = -999.;
    Qang          = 0.0;
    MinAng        =  999.;
    MaxAng        = -999.;
    MinAng_Lon    = -999.;
    MinAng_Lat    = -999.;
    MaxAng_Lon    = -999.;
    MaxAng_Lat    = -999.;
    Area_Ratio    = 0.0; 
    MinArea       =  999.;
    MaxArea       = -999.;
    MinArea_Lon   = -999.;
    MinArea_Lat   = -999.;
    MaxArea_Lon   =  999.;
    MaxArea_Lat   = -999.;

    // Loop over the number of Faces
    //------------------------------
    int nFaces = Faces.size();
    for(unsigned int nn=0; nn < nFaces; nn++) {
      double arc_min = +10.;
      double arc_max = -10.;
      double ang_min = +10.;
      double ang_max = -10.;
      double area    =  0.0;
      double ang_sum =  0.0;
      double x0      =  0.0;
      double y0      =  0.0;
      double DTR     = M_PI/180.;
      double sin0    = sin(DTR*Grid[Faces[nn][0]].lat); 
      double cos0    = cos(DTR*Grid[Faces[nn][0]].lat); 
      bool   pole_panel;

      if(fabs(Nodes[Faces[nn][0]].z) > 0.5) {
        pole_panel = false;
      } else {
        pole_panel = true;
      }

      // Loop over face nodes
      //-----------------------
      for(int k0=0; k0 < 4; k0++) {
        int k1 = (k0+1)%4;
        int k2 = (k0+2)%4;

        // Calculate arc length of this edge
        //-----------------------------------
        double dX0 = Nodes[Faces[nn][k0]].x;
        double dY0 = Nodes[Faces[nn][k0]].y;
        double dZ0 = Nodes[Faces[nn][k0]].z;
        double dX1 = Nodes[Faces[nn][k1]].x;
        double dY1 = Nodes[Faces[nn][k1]].y;
        double dZ1 = Nodes[Faces[nn][k1]].z;
        double dDeltaX = dX1 - dX0;
        double dDeltaY = dY1 - dY0;
        double dDeltaZ = dZ1 - dZ0;
        double dCartLength = sqrt(dDeltaX*dDeltaX + dDeltaY*dDeltaY + dDeltaZ*dDeltaZ);
        double dArcLength  = 2.0 * asin(0.5 * dCartLength);

        // keep min/max arc length for face, for Qarc calc
        //-------------------------------------------------
        if(dArcLength < arc_min) { arc_min = dArcLength; }
        if(dArcLength > arc_max) { arc_max = dArcLength; }

        // Keep track of Global min/max arc length and locations
        //-------------------------------------------------------
        if(dArcLength < MinArc) { 
          MinArc     = dArcLength; 
          MinArc_Lon = Grid[Faces[nn][k0]].lon; 
          MinArc_Lat = Grid[Faces[nn][k0]].lat; 
         
        }
        if(dArcLength > MaxArc) { 
          MaxArc     = dArcLength; 
          MaxArc_Lon = Grid[Faces[nn][k0]].lon; 
          MaxArc_Lat = Grid[Faces[nn][k0]].lat; 
        }

        // Calculate the angle at corner K1 using normals to the face
        //-----------------------------------------------------------
        double dDeltaX1 = + Nodes[Faces[nn][k0]].y * Nodes[Faces[nn][k1]].z
                          - Nodes[Faces[nn][k0]].z * Nodes[Faces[nn][k1]].y;
        double dDeltaY1 = + Nodes[Faces[nn][k0]].z * Nodes[Faces[nn][k1]].x
                          - Nodes[Faces[nn][k0]].x * Nodes[Faces[nn][k1]].z;
        double dDeltaZ1 = + Nodes[Faces[nn][k0]].x * Nodes[Faces[nn][k1]].y
                          - Nodes[Faces[nn][k0]].y * Nodes[Faces[nn][k1]].x;
        double dDeltaX2 = + Nodes[Faces[nn][k2]].y * Nodes[Faces[nn][k1]].z
                          - Nodes[Faces[nn][k2]].z * Nodes[Faces[nn][k1]].y;
        double dDeltaY2 = + Nodes[Faces[nn][k2]].z * Nodes[Faces[nn][k1]].x
                          - Nodes[Faces[nn][k2]].x * Nodes[Faces[nn][k1]].z;
        double dDeltaZ2 = + Nodes[Faces[nn][k2]].x * Nodes[Faces[nn][k1]].y
                          - Nodes[Faces[nn][k2]].y * Nodes[Faces[nn][k1]].x;
        double dDotProd =  dDeltaX1 * dDeltaX2 
                         + dDeltaY1 * dDeltaY2 
                         + dDeltaZ1 * dDeltaZ2;
        double dNorm1 = sqrt( dDeltaX1 * dDeltaX1 
                            + dDeltaY1 * dDeltaY1 
                            + dDeltaZ1 * dDeltaZ1); 
        double dNorm2 = sqrt( dDeltaX2 * dDeltaX2 
                            + dDeltaY2 * dDeltaY2 
                            + dDeltaZ2 * dDeltaZ2);
        dDotProd /= (dNorm1 * dNorm2);

        double dCurrentAngle = acos(dDotProd) - 0.5 * M_PI;

        // keep min/max angle for face, for Qang calc
        //-------------------------------------------------
        if(dCurrentAngle < ang_min) { ang_min = dCurrentAngle; }
        if(dCurrentAngle > ang_max) { ang_max = dCurrentAngle; }

        // Keep track of Global min/max angle and locations
        //-------------------------------------------------------
        if((180./M_PI)*fabs(dCurrentAngle + 0.5 *M_PI) < MinAng) { 
          MinAng     = (180./M_PI)*fabs(dCurrentAngle + 0.5 *M_PI); 
          MinAng_Lon = Grid[Faces[nn][k1]].lon; 
          MinAng_Lat = Grid[Faces[nn][k1]].lat; 
        }
        if((180./M_PI)*fabs(dCurrentAngle + 0.5 *M_PI) > MaxAng) { 
          MaxAng     = (180./M_PI)*fabs(dCurrentAngle + 0.5 *M_PI); 
          MaxAng_Lon = Grid[Faces[nn][k1]].lon; 
          MaxAng_Lat = Grid[Faces[nn][k1]].lat; 
        }

        // Unit sphere element area
        //--------------------------
        area += dCurrentAngle;

        // Determine if the face is counter-clockwise by summing Lat/Lon 
        // around the face. if the sum is + then the nodes have a 
        // counter-clockwise orientation. Avoid the poles by rotating 
        // those faces before summing.
        //--------------------------------------------------------------
//DIAG        double dLat1;
//DIAG        double dLon1;
//DIAG        double dLat2;
//DIAG        double dLon2;
//DIAG        double dLat;
//DIAG        double dLon;
//DIAG
//DIAG        if(pole_panel) {
//DIAG          dLat1 = -asin(Nodes[Faces[nn][k0]].x);
//DIAG          dLon1 = atan2(Nodes[Faces[nn][k0]].y,Nodes[Faces[nn][k0]].z);
//DIAG          dLat2 = -asin(Nodes[Faces[nn][k1]].x);
//DIAG          dLon2 = atan2(Nodes[Faces[nn][k1]].y,Nodes[Faces[nn][k1]].z);
//DIAG        } else {
//DIAG          dLat1 =  asin(Nodes[Faces[nn][k0]].z);
//DIAG          dLon1 = atan2(Nodes[Faces[nn][k0]].y,Nodes[Faces[nn][k0]].x);
//DIAG          dLat2 =  asin(Nodes[Faces[nn][k1]].z);
//DIAG          dLon2 = atan2(Nodes[Faces[nn][k1]].y,Nodes[Faces[nn][k1]].x);
//DIAG        }
//DIAG        if((dLon2 < -0.5 * M_PI) && (dLon1 >  0.5 * M_PI)) {
//DIAG          dLon2 += 2.0 * M_PI;
//DIAG        }
//DIAG        if ((dLon2 >  0.5 * M_PI) && (dLon1 < -0.5 * M_PI)) {
//DIAG          dLon1 += 2.0 * M_PI;
//DIAG        }
//DIAG         ang_sum += (dLon2 - dLon1) * (dLat2 + dLat1);
        double sin1  = sin(DTR*Grid[Faces[nn][k1]].lat); 
        double cos1  = cos(DTR*Grid[Faces[nn][k1]].lat); 
        double sinD  = sin(DTR*(Grid[Faces[nn][k1]].lon - Grid[Faces[nn][0]].lon)); 
        double cosD  = cos(DTR*(Grid[Faces[nn][k1]].lon - Grid[Faces[nn][0]].lon)); 
        double Denom = sin0*sin1 + cos1*cos1*cosD;
        double x1    = (cos1*sinD)/Denom;
        double y1    = (cos0*sin1 - sin0*cos1*cosD)/Denom;

        ang_sum += (x0*y1) - (y0*x1);
        x0 = x1;
        y0 = y1;

      } // for( int k0=0; k0 < 4; k0++) 

      // Quality using edge length metric
      //----------------------------------
      Qarc += (1.0 - arc_min/arc_max);

      // Quality using element angles
      //----------------------------------
      Qang += (1.0 - ang_min/ang_max);

      // Summ up total area
      //--------------------
      Total_Area += area;

      // Keep track of Global min/max face area and locations
      //  (tag the face with face node 0).
      //-------------------------------------------------------
      if(area < MinArea) { 
        MinArea     = area; 
        MinArea_Lon = Grid[Faces[nn][0]].lon; 
        MinArea_Lat = Grid[Faces[nn][0]].lat; 
      }
      if(area > MaxArea) { 
        MaxArea     = area; 
        MaxArea_Lon = Grid[Faces[nn][0]].lon; 
        MaxArea_Lat = Grid[Faces[nn][0]].lat; 
      }

      // increment the Clockwise - CounterClockwise counts
      //--------------------------------------------------
      if(ang_sum <= 0.) {
        Num_Clockwise += 1;
      } else {
        Num_Counter += 1;
      }

    }  // for(unsigned int nn=0; nn < nFaces; nn++) 

    // Calc Area Ratio and scale the Q-values by the number of faces.
    //---------------------------------------------------------------
    double MeshSize = 4.0 * static_cast<double>(nFaces);
    Qang       = Qang / MeshSize;
    Qarc       = Qarc / MeshSize;
    Total_Area = Total_Area /(4.0*M_PI);
    Area_Ratio = MaxArea/MinArea; 

}
