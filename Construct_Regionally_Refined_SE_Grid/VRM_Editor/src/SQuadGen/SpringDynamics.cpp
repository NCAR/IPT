///////////////////////////////////////////////////////////////////////////////
///
///	\file    SpringDynamics.cpp
///	\author  Paul Ullrich
///	\version March 22, 2013
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

#include "SpringDynamics.h"
#include "CubedSphereTrans.h"

#include <cmath>
#include <iostream>


//=============================================================================
void SpringDynamics( NodeVector & vecNodes         ,
                     FaceVector & vecFaces         ,
                     int          nTransSmoothDist ,
                     int          nSmoothIterations) 
{
    // Spring dynamics parameters
    //-----------------------------
    const double ParamTolerance      = 1.0e-7;
    const double ParamSpringConstant = 1.0e-1;
    
    // Mobile node map
    //-----------------------------
    typedef std::map<int, Node>       MobileNodeMap;
    typedef MobileNodeMap::iterator   MobileNodeIterator;
    typedef MobileNodeMap::value_type MobileNodePair;
    
    // Announce
    //----------
    printf("..Performing spring dynamics smoothing\n");
    
    // Construct the segment map
    //-----------------------------
    SegmentMap mapSegments;
    ConstructSegmentMap(vecFaces, mapSegments, -1);
    
    // Set of mobile nodes
    //-----------------------------
    MobileNodeMap mapMobileNodes;

    // Smooth entire mesh
    //--------------------
    if(nTransSmoothDist == (-1)) {
      for (unsigned int n = 0; n < vecFaces.size(); n++) {
        mapMobileNodes.insert(MobileNodePair(vecFaces[n][0], Node()));
        mapMobileNodes.insert(MobileNodePair(vecFaces[n][1], Node()));
        mapMobileNodes.insert(MobileNodePair(vecFaces[n][2], Node()));
        mapMobileNodes.insert(MobileNodePair(vecFaces[n][3], Node()));
      }
    } else {
      // Add nodes which are connected to the refinement region
      //--------------------------------------------------------
      for (unsigned int n = 0; n < vecFaces.size(); n++) {
        if (vecFaces[n].nRefineLevel == (-1)) {
          mapMobileNodes.insert(MobileNodePair(vecFaces[n][0], Node()));
          mapMobileNodes.insert(MobileNodePair(vecFaces[n][1], Node()));
          mapMobileNodes.insert(MobileNodePair(vecFaces[n][2], Node()));
          mapMobileNodes.insert(MobileNodePair(vecFaces[n][3], Node()));
        }
      }
    
      // Remove nodes on the boundary of the transition region
      //--------------------------------------------------------
      if(nTransSmoothDist == 0) {
        for (unsigned int n = 0; n < vecFaces.size(); n++) {
          if(vecFaces[n].nRefineLevel != (-1)) {
            for (int k = 0; k < 4; k++) {
              MobileNodeIterator iterNode = mapMobileNodes.find(vecFaces[n][k]);
              if(iterNode != mapMobileNodes.end()) { mapMobileNodes.erase(iterNode); }
            }
          }
        }

        // Build up nodes away from boundary of transition region
      } else if(nTransSmoothDist > 1) {
        for (int z = 0; z < nTransSmoothDist-1; z++) {
          MobileNodeMap mapOldMobileNodes = mapMobileNodes;
          SegmentMapIterator iterSeg = mapSegments.begin();
          for (; iterSeg != mapSegments.end(); iterSeg++) {
            MobileNodeIterator iterNode;
            iterNode = mapOldMobileNodes.find(iterSeg->first[0]);
            if(iterNode != mapOldMobileNodes.end()) {
              mapMobileNodes.insert( MobileNodePair(iterSeg->first[1], Node()));
            }
            iterNode = mapOldMobileNodes.find(iterSeg->first[1]);
            if(iterNode != mapOldMobileNodes.end()) {
              mapMobileNodes.insert( MobileNodePair(iterSeg->first[0], Node()));
            }
          }
        }
      }
    }

    // Remove segments which contain no mobile nodes
    //-----------------------------------------------
    SegmentMapIterator iterSeg = mapSegments.begin();
    while (iterSeg != mapSegments.end()) {
      if((mapMobileNodes.find(iterSeg->first[0]) == mapMobileNodes.end()) &&
         (mapMobileNodes.find(iterSeg->first[1]) == mapMobileNodes.end())    ) {
        mapSegments.erase(iterSeg++);
      } else {
        ++iterSeg;
      }
    }

    printf("....Number of active nodes: %i\n", mapMobileNodes.size());
    printf("....Number of active segments: %i\n", mapSegments.size());

    if(mapMobileNodes.size() == 0) {
      printf("....No active nodes to smooth; returning\n");
    }

    // Perform spring dynamics
    //-------------------------
    for (int i = 0; i < nSmoothIterations; i++) {

      MobileNodeIterator iterNode;
    
      // Initialize the vector of displacements
      //----------------------------------------
      iterNode = mapMobileNodes.begin();
      for (; iterNode != mapMobileNodes.end(); iterNode++) {
        iterNode->second.x = 0.0;
        iterNode->second.y = 0.0;
        iterNode->second.z = 0.0;
      }
    
      // Compute displacement
      //-----------------------
      iterSeg = mapSegments.begin();
      for (; iterSeg != mapSegments.end(); iterSeg++) {
        Node disp;
        disp.x = + vecNodes[iterSeg->first[1]].x
                 - vecNodes[iterSeg->first[0]].x;
        disp.y = + vecNodes[iterSeg->first[1]].y
                 - vecNodes[iterSeg->first[0]].y;
        disp.z = + vecNodes[iterSeg->first[1]].z
                 - vecNodes[iterSeg->first[0]].z;
    
        MobileNodeIterator iterNode0 = mapMobileNodes.find(iterSeg->first[0]);
        MobileNodeIterator iterNode1 = mapMobileNodes.find(iterSeg->first[1]);
    
        double dWeight = 1.0;
/*
        Node avg;
        avg.x = 0.5 * ( + vecNodes[iterSeg->first[1]].x
                        + vecNodes[iterSeg->first[0]].x);
        avg.y = 0.5 * ( + vecNodes[iterSeg->first[1]].y
                        - vecNodes[iterSeg->first[0]].y);
        avg.z = 0.5 * ( + vecNodes[iterSeg->first[1]].z
                        - vecNodes[iterSeg->first[0]].z);

        double dR = sqrt(avg.x * avg.x + avg.y * avg.y + avg.z * avg.z);
        avg.x /= dR;
        avg.y /= dR;
        avg.z /= dR;
        
        double dLat = asin(avg.z);
        double dLon = atan2(avg.y, avg.x);
        
        double dA;
        double dB;
        int iP;
        CubedSphereTrans::ABPFromRLL(dLon, dLat, dA, dB, iP);
        
        double dX = tan(dA);
        double dY = tan(dB);

        dWeight = pow(1.0 + dX * dX * dY * dY, 2.0);
*/
        if(iterNode0 != mapMobileNodes.end()) {
          iterNode0->second.x += ParamSpringConstant * dWeight * disp.x;
          iterNode0->second.y += ParamSpringConstant * dWeight * disp.y;
          iterNode0->second.z += ParamSpringConstant * dWeight * disp.z;
        }
        if(iterNode1 != mapMobileNodes.end()) {
          iterNode1->second.x -= ParamSpringConstant * dWeight * disp.x;
          iterNode1->second.y -= ParamSpringConstant * dWeight * disp.y;
          iterNode1->second.z -= ParamSpringConstant * dWeight * disp.z;
        }
      }

      // Apply displacement
      //--------------------
      double dTotalDisplacement = 0.0;

      iterNode = mapMobileNodes.begin();
      for (; iterNode != mapMobileNodes.end(); iterNode++) {
        dTotalDisplacement += sqrt( + iterNode->second.x * iterNode->second.x
                                    + iterNode->second.y * iterNode->second.y
                                    + iterNode->second.z * iterNode->second.z);
      
        vecNodes[iterNode->first].x += iterNode->second.x;
        vecNodes[iterNode->first].y += iterNode->second.y;
        vecNodes[iterNode->first].z += iterNode->second.z;
      
        double dX = vecNodes[iterNode->first].x;
        double dY = vecNodes[iterNode->first].y;
        double dZ = vecNodes[iterNode->first].z;
        double dR = sqrt(dX * dX + dY * dY + dZ * dZ);
      
        vecNodes[iterNode->first].x /= dR;
        vecNodes[iterNode->first].y /= dR;
        vecNodes[iterNode->first].z /= dR;
      }

      printf("....Smoothing iteration %i residual: %1.5f\n", i, dTotalDisplacement);
    }

} // SpringDynamics
//=============================================================================


//=============================================================================
void PressureDynamics( NodeVector & vecNodes         ,
                       FaceVector & vecFaces         ,
                       int          nTransSmoothDist ,
                       int          nSmoothIterations) 
{
    // Spring dynamics parameters
    //-----------------------------
    const double ParamTolerance      = 1.0e-7;
    const double ParamSpringConstant = 0.02;

    // Announce
    //--------------
    printf("..Performing pressure dynamics smoothing\n");

    // Construct the segment map
    //---------------------------
    SegmentMap mapSegments;
    ConstructSegmentMap(vecFaces, mapSegments, -1);

    // Perform spring dynamics
    //---------------------------
    for (int i = 0; i < nSmoothIterations; i++) {

      // Copy the node vector
      //---------------------
      NodeVector vecNodesOld = vecNodes;

      // Find the midpoint of all faces
      //----------------------------------
      for (int j = 0; j < vecFaces.size(); j++) {
        Node nodeCenter;

        nodeCenter.x = 0.25*( + vecNodesOld[vecFaces[j][0]].x
                              + vecNodesOld[vecFaces[j][1]].x
                              + vecNodesOld[vecFaces[j][2]].x
                              + vecNodesOld[vecFaces[j][3]].x);

        nodeCenter.y = 0.25*( + vecNodesOld[vecFaces[j][0]].y
                              + vecNodesOld[vecFaces[j][1]].y
                              + vecNodesOld[vecFaces[j][2]].y
                              + vecNodesOld[vecFaces[j][3]].y);

        nodeCenter.z = 0.25*( + vecNodesOld[vecFaces[j][0]].z
                              + vecNodesOld[vecFaces[j][1]].z
                              + vecNodesOld[vecFaces[j][2]].z
                              + vecNodesOld[vecFaces[j][3]].z);

        for (int k = 0; k < 4; k++) {
          Node disp;
          disp.x = vecNodesOld[vecFaces[j][k]].x - nodeCenter.x;
          disp.y = vecNodesOld[vecFaces[j][k]].y - nodeCenter.y;
          disp.z = vecNodesOld[vecFaces[j][k]].z - nodeCenter.z;

          disp.x *= ParamSpringConstant;
          disp.y *= ParamSpringConstant;
          disp.z *= ParamSpringConstant;

          vecNodes[vecFaces[j][k]].x += disp.x;
          vecNodes[vecFaces[j][k]].y += disp.y;
          vecNodes[vecFaces[j][k]].z += disp.z;
        }
      }

      // Normalize all nodes
      //---------------------
      for (int k = 0; k < vecNodes.size(); k++) {
        double dX = vecNodes[k].x;
        double dY = vecNodes[k].y;
        double dZ = vecNodes[k].z;
        double dR = sqrt(dX * dX + dY * dY + dZ * dZ);

        vecNodes[k].x /= dR;
        vecNodes[k].y /= dR;
        vecNodes[k].z /= dR;
      }
    }
} // PressureDynamics
///////////////////////////////////////////////////////////////////////////////


