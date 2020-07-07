#ifndef _NEIGHBORADJUSTMENTS_H_
#define _NEIGHBORADJUSTMENTS_H_

#include "SQuadGen/GridElements.h"
#include "RefinementMap.h"
#include <cmath>

//=====================================================================
//bool PRESERVE_VERTICIES  = true;

enum NeighborType {NEAR,DISTANT};

struct Neighbor {
  int          ixNode;
  NeighborType type;
  double       Alpha;
  double       CosAlpha;
  double       SinAlpha;
  double       Ugc;
  double       Vgc;
  double       Beta;
  

  // Constructor
  //==============
  Neighbor( int I_ixNode, NeighborType I_type) 
  {
    ixNode   = I_ixNode;
    type     = I_type;
    Alpha    = 0.;
    CosAlpha = 0.;
    SinAlpha = 0.;
    Ugc      = 0.;
    Vgc      = 0.;
    Beta     = 0.;
  }
};

typedef std::vector< std::vector<Neighbor> > NeighborList;

struct NodeData {
  double       Lat;
  double       Lon;
  double       CosLat;
  double       SinLat;
  double       CosLon;
  double       SinLon;
  double       Rval;
  bool         HasMoved;

  // Constructor
  //==============
  NodeData( double I_Lat   , double I_Lon   ,
            double I_CosLat, double I_SinLat,
            double I_CosLon, double I_SinLon, 
            double I_Rval  , bool   I_HasMoved)
  {
    Lat      = I_Lat;
    Lon      = I_Lon;
    CosLat   = I_CosLat;
    SinLat   = I_SinLat;
    CosLon   = I_CosLon;
    SinLon   = I_SinLon;
    Rval     = I_Rval;
    HasMoved = I_HasMoved;
  }
};
//=====================================================================


//=====================================================================
void EquializeDistances(RefinementMap & I_refMap,
                        NodeVector    & vecNodes, 
                        FaceVector    & vecFaces,
                        int             nPercentStepsize,
                        int             nSmoothIterations,
                        int             NE_lo,
                        int             NE_hi);
//=====================================================================


//=====================================================================
void RefmapSmooth(RefinementMap & I_refMap,
                  NodeVector    & vecNodes, 
                  FaceVector    & vecFaces,
                  int             nPercentStepsize,
                  int             nSmoothIterations,
                  int             NE_lo,
                  int             NE_hi);
//=====================================================================


//=====================================================================
void RefmapDistances(NodeVector & vecNodes, 
                     FaceVector & vecFaces,
                     int          nPercentStepsize,
                     int          nSmoothIterations);
//=====================================================================

#endif
