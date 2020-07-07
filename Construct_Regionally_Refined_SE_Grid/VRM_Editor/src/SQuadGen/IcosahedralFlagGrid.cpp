///////////////////////////////////////////////////////////////////////////////
///
///	\file    IcosahedralFlagGrid.cpp
///	\author  Paul Ullrich
///	\version February 21, 2012
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "IcosahedralFlagGrid.h"

#include <cmath>


//=============================================================================
void GenerateFacesFromTriangle( int          nRefineLevel,
                                const Edge & edge0       ,
                                const Edge & edge1       ,
                                const Edge & edge2       ,
                                NodeVector & vecNodes    ,
                                FaceVector & vecFaces    ) 
{
    int  i;
    int  j;
    int  k;
    int  ixEndNode;
    int  ixInt;
    Edge edgeBot;      // Edges
    Edge edgeMid;
    Edge edgeTop;
    
    // Initial bottom edge
    //----------------------
    edgeBot.push_back(edge0[0]);
    
    // Loop over all refined faces
    //------------------------------
    for (j = 0; j < nRefineLevel; j++) {
    
      // Generate mid level vertices
      //-------------------------------
      GenerateEdgeVertices(2*j+1, edge0[2*j+1], edge1[2*j+1], vecNodes, edgeMid);
    
      // Generate top level vertices
      //-------------------------------
      if(j == nRefineLevel-1) {
        edgeTop = edge2;
      } else {
        GenerateEdgeVertices(2*j+2, edge0[2*j+2], edge1[2*j+2], vecNodes, edgeTop);
      }
    
      // Generate faces
      //---------------------
      for (i = 0; i < 2*j+1; i++) {
        if (i % 2 == 0) {
          // Downward pointing faces
          //--------------------------
          ixInt = InsertTriFaceCentroidNode(edgeMid[i  ], 
                                            edgeMid[i+1], 
                                            edgeTop[i+1],vecNodes);
          vecFaces.push_back(Face(edgeBot[i  ],edgeMid[i  ],ixInt,edgeMid[i+1]));
          vecFaces.push_back(Face(edgeTop[i  ],edgeTop[i+1],ixInt,edgeMid[i  ]));
          vecFaces.push_back(Face(edgeTop[i+2],edgeMid[i+1],ixInt,edgeTop[i+1]));
        } else {
          // Upward pointing faces
          //--------------------------
          ixInt = InsertTriFaceCentroidNode(edgeMid[i  ], 
                                            edgeMid[i+1], 
                                            edgeBot[i  ],vecNodes);
          vecFaces.push_back(Face(edgeBot[i-1],edgeMid[i  ],ixInt,edgeBot[i  ]));
          vecFaces.push_back(Face(edgeTop[i+1],edgeMid[i+1],ixInt,edgeMid[i  ]));
          vecFaces.push_back(Face(edgeBot[i+1],edgeBot[i  ],ixInt,edgeMid[i+1]));
        }
      }
    
      // New bottom edge
      //-----------------
      edgeBot = edgeTop;
    }

} // GenerateFacesFromTriangle
//=============================================================================


//=============================================================================
void ConvertFromLonLatToCartesian( const LonLatNodeVector & vecLonLatNodes,
                                         NodeVector       & vecNodes      ) 
{
    unsigned int i;
    vecNodes.resize(vecLonLatNodes.size());

    // Loop over all nodes
    //----------------------
    for (i = 0; i < vecLonLatNodes.size(); i++) {
      vecNodes[i].x = sin(vecLonLatNodes[i].lon)*cos(vecLonLatNodes[i].lat);
      vecNodes[i].y = cos(vecLonLatNodes[i].lon)*cos(vecLonLatNodes[i].lat);
      vecNodes[i].z = sin(vecLonLatNodes[i].lat);
    }
} // ConvertFromLonLatToCartesian
//=============================================================================


//=============================================================================
void GenerateIcosahedralQuadGrid( int          nRefineLevel,
                                  NodeVector & vecNodes    ,
                                  FaceVector & vecFaces    ) 
{
    // Latitude of nodes (Northern Hemisphere)
    //-------------------------------------------
    const double NodeLat = atan(0.5);
    
    // Store all icosahedral nodes
    //-------------------------------
    LonLatNodeVector vecLonLatNodes;
    
    vecLonLatNodes.push_back(LonLatNode(0.0         , -0.5*M_PI));
    vecLonLatNodes.push_back(LonLatNode(0.0         , -NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.2, -NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.4, -NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.6, -NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.8, -NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.1, +NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.3, +NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.5, +NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.7, +NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.9, +NodeLat ));
    vecLonLatNodes.push_back(LonLatNode(0.0         , +0.5*M_PI));
    
    // Convert icosahedral nodes to Cartesian geometry
    //-------------------------------------------------
    ConvertFromLonLatToCartesian(vecLonLatNodes, vecNodes);
    
    // Vector of edges
    //-----------------
    EdgeVector vecEdges;
    vecEdges.resize(30);
    
    // Generate vertices along edges
    //-------------------------------
    for (int i = 0; i < 5; i++) {
      GenerateEdgeVertices(2*nRefineLevel, 0, i+1, vecNodes, vecEdges[i]);
    }
    
    for (int i = 0; i < 5; i++) {
      GenerateEdgeVertices(2*nRefineLevel, i+1, ((i+1)%5)+1, vecNodes, vecEdges[i+5]);
    }
    
    GenerateEdgeVertices(2*nRefineLevel,  1,  6, vecNodes, vecEdges[10]);
    GenerateEdgeVertices(2*nRefineLevel,  6,  2, vecNodes, vecEdges[11]);
    GenerateEdgeVertices(2*nRefineLevel,  2,  7, vecNodes, vecEdges[12]);
    GenerateEdgeVertices(2*nRefineLevel,  7,  3, vecNodes, vecEdges[13]);
    GenerateEdgeVertices(2*nRefineLevel,  3,  8, vecNodes, vecEdges[14]);
    GenerateEdgeVertices(2*nRefineLevel,  8,  4, vecNodes, vecEdges[15]);
    GenerateEdgeVertices(2*nRefineLevel,  4,  9, vecNodes, vecEdges[16]);
    GenerateEdgeVertices(2*nRefineLevel,  9,  5, vecNodes, vecEdges[17]);
    GenerateEdgeVertices(2*nRefineLevel,  5, 10, vecNodes, vecEdges[18]);
    GenerateEdgeVertices(2*nRefineLevel, 10,  1, vecNodes, vecEdges[19]);
    
    for (int i = 0; i < 5; i++) {
      GenerateEdgeVertices(2*nRefineLevel, i+6, ((i+1)%5)+6, vecNodes, vecEdges[i+20]);
    }
    
    for (int i = 0; i < 5; i++) {
      GenerateEdgeVertices(2*nRefineLevel, i+6, 11, vecNodes, vecEdges[i+25]);
    }
    
    // Generate south polar faces
    //-----------------------------
    for (int i = 0; i < 5; i++) {
      GenerateFacesFromTriangle(nRefineLevel, vecEdges[i      ], 
                                              vecEdges[(i+1)%5], 
                                              vecEdges[ i+5   ], vecNodes, vecFaces);
    }
    
    // Generate south equatorial faces
    //----------------------------------
    for (int i = 0; i < 5; i++) {
      GenerateFacesFromTriangle(nRefineLevel, vecEdges[2*i+10], 
                                              vecEdges[  i+5 ], 
                                              vecEdges[2*i+11], vecNodes, vecFaces);
    }
    
    // Generate north equatorial faces
    //---------------------------------
    for (int i = 0; i < 5; i++) {
      GenerateFacesFromTriangle(nRefineLevel, vecEdges[  i+20        ], 
                                              vecEdges[2*i+11        ], 
                                              vecEdges[2*((i+1)%5)+10].Flip(), vecNodes, vecFaces);
    }

    // Generate north polar faces
    //--------------------------------
    for (int i = 0; i < 5; i++) {
      GenerateFacesFromTriangle(nRefineLevel, vecEdges[i+25        ], 
                                              vecEdges[i+20        ], 
                                              vecEdges[((i+1)%5)+25].Flip(), vecNodes, vecFaces);
    }

} // GenerateIcosahedralQuadGrid
//=============================================================================


