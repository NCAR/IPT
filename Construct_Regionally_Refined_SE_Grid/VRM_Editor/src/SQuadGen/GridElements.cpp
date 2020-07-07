///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.cpp
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

#include "GridElements.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

#define EQUIANGULAR_SUB_NODES


//=============================================================================
void ReverseFaceOrientation( FaceVector & vecFaces) 
{
    for (unsigned int i = 0; i < vecFaces.size(); i++) {
      int     ixNode = vecFaces[i][0];
      vecFaces[i][0] = vecFaces[i][2];
      vecFaces[i][2] = ixNode;
    }
} // ReverseFaceOrientation
//=============================================================================


//=============================================================================
int InsertSubNode( int ix0, int ix1, double alpha, NodeVector & vecNodes) 
{
#ifdef EQUIANGULAR_SUB_NODES
    //alpha = tan((2.0*alpha-1.0)*(M_PI/6.0)) / tan(M_PI/6.0);
    //alpha = 0.5*(alpha+1.0);
    //dArcLength = 2.0 * asin(0.5 * dCartLength);
    double     dDeltaX = (vecNodes[ix1].x - vecNodes[ix0].x);
    double     dDeltaY = (vecNodes[ix1].y - vecNodes[ix0].y);
    double     dDeltaZ = (vecNodes[ix1].z - vecNodes[ix0].z);
    double dCartLength = sqrt(dDeltaX*dDeltaX + dDeltaY*dDeltaY + dDeltaZ*dDeltaZ);
    //double dArcLength = 2.0 * asin(0.5 * dCartLength);
    
    double      dGamma = acos(0.5 * dCartLength);
    double      dTheta = acos(1.0 - 0.5 * dCartLength * dCartLength);
    double dAlphaTheta = alpha * dTheta;
    double       dBeta = M_PI - dGamma - dAlphaTheta;

    alpha = sin(dAlphaTheta) / sin(dBeta) / dCartLength;
#endif

    double dX = vecNodes[ix0].x + (vecNodes[ix1].x - vecNodes[ix0].x)*alpha;
    double dY = vecNodes[ix0].y + (vecNodes[ix1].y - vecNodes[ix0].y)*alpha;
    double dZ = vecNodes[ix0].z + (vecNodes[ix1].z - vecNodes[ix0].z)*alpha;
    
    // Project to sphere
    //-----------------------
    double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);
    dX /= dRadius;
    dY /= dRadius;
    dZ /= dRadius;
    
    // Index
    //-------
    int ix = vecNodes.size();
    
    // Insert node
    //------------
    vecNodes.push_back(Node(dX, dY, dZ));
    
    return ix;

} // InsertSubNode
//=============================================================================


//=============================================================================
int InsertTriFaceCentroidNode( int ix0, int ix1, int ix2, NodeVector & vecNodes)
{
    double dX = (vecNodes[ix0].x + vecNodes[ix1].x + vecNodes[ix2].x) / 3.0;
    double dY = (vecNodes[ix0].y + vecNodes[ix1].y + vecNodes[ix2].y) / 3.0;
    double dZ = (vecNodes[ix0].z + vecNodes[ix1].z + vecNodes[ix2].z) / 3.0;

    // Project to sphere
    //---------------------
    double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);
    
    dX /= dRadius;
    dY /= dRadius;
    dZ /= dRadius;
    
    // Index
    //-------
    int ix = vecNodes.size();
    
    // Insert node
    //-------------
    vecNodes.push_back(Node(dX, dY, dZ));
    
    return ix;

} // InsertTriFaceCentroidNode
//=============================================================================


//=============================================================================
int InsertQuadNodeCenter( NodeVector & vecNodes, const Node & node0, const Node & node1, 
                                                 const Node & node2, const Node & node3) 
{
    double dX = 0.25*(node0.x + node1.x + node2.x + node3.x);
    double dY = 0.25*(node0.y + node1.y + node2.y + node3.y);
    double dZ = 0.25*(node0.z + node1.z + node2.z + node3.z);
    
    double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);
    
    dX /= dRadius;
    dY /= dRadius;
    dZ /= dRadius;
    
    int ixNode = vecNodes.size();
    vecNodes.push_back(Node(dX, dY, dZ));

    return ixNode;

} // InsertQuadNodeCenter
//=============================================================================


//=============================================================================
int InsertQuadNode( NodeVector & vecNodes, const Node & node0, const Node & node1, 
                                           const Node & node2, const Node & node3, 
                                           double dAlpha     , double dBeta      ) 
{
    double dX0 = (1.0 - dBeta)*node0.x + dBeta*node3.x;
    double dY0 = (1.0 - dBeta)*node0.y + dBeta*node3.y;
    double dZ0 = (1.0 - dBeta)*node0.z + dBeta*node3.z;
    
    double dX1 = (1.0 - dBeta)*node1.x + dBeta*node2.x;
    double dY1 = (1.0 - dBeta)*node1.y + dBeta*node2.y;
    double dZ1 = (1.0 - dBeta)*node1.z + dBeta*node2.z;
    
    double dX = (1.0 - dAlpha)*dX0 + dAlpha*dX1;
    double dY = (1.0 - dAlpha)*dY0 + dAlpha*dY1;
    double dZ = (1.0 - dAlpha)*dZ0 + dAlpha*dZ1;
    
    double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);
    
    dX /= dRadius;
    dY /= dRadius;
    dZ /= dRadius;
    
    int ixNode = vecNodes.size();
    vecNodes.push_back(Node(dX, dY, dZ));
    
    return ixNode;

} // InsertQuadNode
//=============================================================================


//=============================================================================
void GenerateEdgeVertices( int nRefineLevel, int ix0, int ix1, 
                           NodeVector & vecNodes, Edge & edge)
{
    edge.clear();
    edge.push_back(ix0);
    
    for (int i = 1; i < nRefineLevel; i++) {
    
      // Nodes along line in Cartesian geometry
      //----------------------------------------
      double alpha = static_cast<double>(i) / static_cast<double>(nRefineLevel);
    
      // Insert node along edge
      //-------------------------
      int ixNode = InsertSubNode(ix0, ix1, alpha, vecNodes);
    
      // Add node to edge
      //-------------------
      edge.push_back(ixNode);
    }
    edge.push_back(ix1);

} // GenerateEdgeVertices
//=============================================================================


//=============================================================================
void ConstructSegmentMap( FaceVector & vecFaces, SegmentMap & mapSegments, int nRefineLevel)
{
    unsigned int i;
    int k;
    
    // Clear the segment list from all faces
    //-----------------------------------------
    for (unsigned int i = 0; i < vecFaces.size(); i++) {
      vecFaces[i].nSegments = 0;
    }
    
    // Construct the segment map
    //---------------------------
    mapSegments.clear();
    for (i = 0; i < vecFaces.size(); i++) {
      if(vecFaces[i].nRefineLevel >= nRefineLevel) {
        for (k = 0; k < 4; k++) {
          Segment            seg(vecFaces[i][k],vecFaces[i][(k+1)%4]);
          FacePair           facepair;
          SegmentMapIterator iter = mapSegments.insert(SegmentMapPair(seg,facepair)).first;
          iter->second.AddFace(i);
          if(vecFaces[i].ixSplitMidpoints[k] != InvalidNode) {
            iter->second.nMidpoints        = 1;
            iter->second.ixMidpointNode[0] = vecFaces[i].ixSplitMidpoints[k];
          }
          vecFaces[i].AddSegment(iter);
        }
      }
    }
/*
    // Verify all faces have a complete set of segments
    //--------------------------------------------------
    for (int i = 0; i < vecFaces.size(); i++) {
      if (vecFaces[i].nSegments != 4) {
        _EXCEPTIONT("Face has an incomplete set of segments.");
      }
    }
*/
/*
    SegmentMapIterator iter;
    for (iter = mapSegments.begin(); iter != mapSegments.end(); iter++) {
      if (iter->second.nFaces != 2) {
        _EXCEPTION2("Segment (%i,%i) has an incomplete set of faces.",
        iter->first.ixNode[0], iter->first.ixNode[1]);
      }
    }
*/

} // ConstructSegmentMap
//=============================================================================

