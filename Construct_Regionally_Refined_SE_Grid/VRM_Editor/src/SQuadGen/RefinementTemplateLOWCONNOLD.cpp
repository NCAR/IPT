///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefinementTemplateLOWCONNOLD.cpp
///	\author  Paul Ullrich
///	\version March 17, 2013
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

#include "RefinementTemplateLOWCONNOLD.h"

#include <iostream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

//#define NO_CONNECTING_FACES


//=============================================================================
RefinementTemplateLOWCONNOLD::RefinementTemplateLOWCONNOLD() 
{
    dTrisectSpacing[0] = 0.4;
    dTrisectSpacing[1] = 0.8;
}
//=============================================================================


//=============================================================================
void RefinementTemplateLOWCONNOLD::ApplyStandardTemplate(FaceVector & vecFaces    ,
                                                         NodeVector & vecNodes    ,
                                                         int          ixActiveFace,
                                                         int          ixActiveSeg ) 
{
    Face & face = vecFaces[ixActiveFace];

    if(ixActiveSeg == (-1)) { _EXCEPTIONT("Logic error"); }

    SegmentMapIterator iterSeg      = face.iterSegment[ ixActiveSeg     ];
    SegmentMapIterator iterSegLeft  = face.iterSegment[(ixActiveSeg+1)%4];
    SegmentMapIterator iterSegOpp   = face.iterSegment[(ixActiveSeg+2)%4];
    SegmentMapIterator iterSegRight = face.iterSegment[(ixActiveSeg+3)%4];
    
    // Final processed segment
    //--------------------------
    iterSegFinal = iterSegLeft;
    
    // Trisect left and right segments
    //----------------------------------
    if(iterSeg->first[0] != face[ixActiveSeg]) { iterSeg->first.Flip(); }
    
    if(iterSegLeft->first[0] != iterSeg->first[1]) {
      if(iterSegLeft->first[1] != iterSeg->first[1]) { _EXCEPTIONT("Logic error"); }
      iterSegLeft->first.Flip();
    }
    
    if(iterSegRight->first[0] != iterSeg->first[0]) {
      if (iterSegRight->first[1] != iterSeg->first[0]) { _EXCEPTIONT("Logic error"); }
      iterSegRight->first.Flip();
    }
    
    DissectSegment(iterSeg     , vecNodes, 3);
    DissectSegment(iterSegLeft , vecNodes, 2, dTrisectSpacing);
    DissectSegment(iterSegRight, vecNodes, 2, dTrisectSpacing);
    BisectSegment (iterSegOpp  , vecNodes);
    
    // Insert new nodes
    //--------------------
    int ixCentNode = InsertQuadNodeCenter(vecNodes,
                                          vecNodes[face[ixActiveSeg]],
                                          vecNodes[face[(ixActiveSeg+1)%4]],
                                          vecNodes[iterSegRight->second.ixMidpointNode[1]],
                                          vecNodes[iterSegLeft->second.ixMidpointNode[1]]);
    
    int ixNode0 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[iterSeg->second.ixMidpointNode[1]],
                                       vecNodes[face[ixActiveSeg]],
                                       vecNodes[iterSegRight->second.ixMidpointNode[0]],
                                       vecNodes[ixCentNode]);
    
    int ixNode1 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[(ixActiveSeg+1)%4]],
                                       vecNodes[iterSeg->second.ixMidpointNode[1]],
                                       vecNodes[ixCentNode],
                                       vecNodes[iterSegLeft->second.ixMidpointNode[0]]);
    
    int ixNode2 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[iterSegLeft->second.ixMidpointNode[1]],
                                       vecNodes[iterSegOpp->second.ixMidpointNode[0]],
                                       vecNodes[ixCentNode],
                                       vecNodes[iterSegLeft->second.ixMidpointNode[0]]);
    
    int ixNode3 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[iterSegRight->second.ixMidpointNode[0]],
                                       vecNodes[iterSegRight->second.ixMidpointNode[1]],
                                       vecNodes[iterSegOpp->second.ixMidpointNode[0]],
                                       vecNodes[ixCentNode]);
    
#ifndef NO_CONNECTING_FACES
    // Construct connecting face
    //---------------------------
    if(iFirstConnectingNode == InvalidNode) { iFirstConnectingNode = ixNode3; }

    if(iLastConnectingNode != InvalidNode) {
      vecFaces.push_back(Face(iterSegRight->second.ixMidpointNode[1],
                              iLastConnectingNode                   ,
                              iterSegRight->second.ixMidpointNode[0],
                              ixNode3                               ));
    }
    
    // Store connecting node
    //---------------------------
    iLastConnectingNode = ixNode2;
#endif

    // Construct faces
    //----------------------------
    vecFaces.push_back(Face(face[(ixActiveSeg+1)%4]              ,
                            iterSegLeft->second.ixMidpointNode[0],
                            ixNode1                              ,
                            iterSeg->second.ixMidpointNode[2]    ));

    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[2],
                            ixNode1                          ,
                            ixCentNode                       ,
                            iterSeg->second.ixMidpointNode[1]));

    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[1],
                            ixCentNode                       ,
                            ixNode0                          ,
                            iterSeg->second.ixMidpointNode[0]));

    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[0]     ,
                            ixNode0                               ,
                            iterSegRight->second.ixMidpointNode[0],
                            face[ixActiveSeg]                     ));

    vecFaces.push_back(Face(ixCentNode                           ,
                            ixNode1                              ,
                            iterSegLeft->second.ixMidpointNode[0],
                            ixNode2                              ));

    vecFaces.push_back(Face(ixCentNode                          ,
                            ixNode2                             ,
                            iterSegOpp->second.ixMidpointNode[0],
                            ixNode3                             ));

    vecFaces.push_back(Face(ixCentNode                            ,
                            ixNode3                               ,
                            iterSegRight->second.ixMidpointNode[0],
                            ixNode0                               ));

    vecFaces.push_back(Face(ixNode2                              ,
                            iterSegLeft->second.ixMidpointNode[1],
                            face[(ixActiveSeg+2)%4]              ,
                            iterSegOpp->second.ixMidpointNode[0] ));

    vecFaces.push_back(Face(ixNode3                               ,
                            iterSegOpp->second.ixMidpointNode[0]  ,
                            face[(ixActiveSeg+3)%4]               ,
                            iterSegRight->second.ixMidpointNode[1]));

    // Tag this face for deletion
    //----------------------------
    vecFaces[ixActiveFace].nTag = (-1);

} // RefinementTemplateLOWCONNOLD::ApplyStandardTemplate
//=============================================================================


//=============================================================================
void RefinementTemplateLOWCONNOLD::ApplyInnerCornerTemplate(FaceVector & vecFaces    ,
                                                            NodeVector & vecNodes    ,
                                                            int          ixActiveFace,
                                                            int          ixActiveSeg ) 
{
    Face & face = vecFaces[ixActiveFace];

    if(ixActiveSeg == (-1)) { _EXCEPTIONT("Logic error"); }

    SegmentMapIterator iterSeg      = face.iterSegment[ ixActiveSeg     ];
    SegmentMapIterator iterSegLeft  = face.iterSegment[(ixActiveSeg+1)%4];
    SegmentMapIterator iterSegOpp   = face.iterSegment[(ixActiveSeg+2)%4];
    SegmentMapIterator iterSegRight = face.iterSegment[(ixActiveSeg+3)%4];
    
    // Final processed segment
    //--------------------------
    iterSegFinal = iterSegOpp;
    
    // Set up bounding edges
    //--------------------------
    if(iterSeg->first[0] != face[ixActiveSeg]) { iterSeg->first.Flip(); }
    
    if(iterSegLeft->first[0] != iterSeg->first[1]) {
      if(iterSegLeft->first[1] != iterSeg->first[1]) { _EXCEPTIONT("Logic error"); }
      iterSegLeft->first.Flip();
    }
    
    if(iterSegRight->first[0] != iterSeg->first[0]) {
      if(iterSegRight->first[1] != iterSeg->first[0]) { _EXCEPTIONT("Logic error"); }
      iterSegRight->first.Flip();
    }
    
    if(iterSegOpp->first[0] != iterSegLeft->first[1]) {
      if(iterSegOpp->first[1] != iterSegLeft->first[1]) { _EXCEPTIONT("Logic error"); }
      iterSegOpp->first.Flip();
    }
    
    DissectSegment(iterSeg     , vecNodes, 3);
    DissectSegment(iterSegLeft , vecNodes, 3);
    DissectSegment(iterSegRight, vecNodes, 2, dTrisectSpacing);
    DissectSegment(iterSegOpp  , vecNodes, 2, dTrisectSpacing);
    
    // Insert new nodes
    //-------------------
    int ixCentNode = InsertQuadNode(vecNodes,
                                    vecNodes[face[(ixActiveSeg+1)%4]],
                                    vecNodes[face[(ixActiveSeg+2)%4]],
                                    vecNodes[face[(ixActiveSeg+3)%4]],
                                    vecNodes[face[(ixActiveSeg+0)%4]],
                                    0.4, 0.4);
    
    int ixNode0 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[iterSeg->second.ixMidpointNode[1]],
                                       vecNodes[face[(ixActiveSeg+1)%4]],
                                       vecNodes[iterSegLeft->second.ixMidpointNode[1]],
                                       vecNodes[ixCentNode]);
    
    int ixNode1 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[(ixActiveSeg+2)%4]],
                                       vecNodes[iterSegLeft->second.ixMidpointNode[1]],
                                       vecNodes[ixCentNode],
                                       vecNodes[iterSegOpp->second.ixMidpointNode[1]]);
    
    int ixNode2 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[(ixActiveSeg+2)%4]],
                                       vecNodes[face[(ixActiveSeg+3)%4]],
                                       vecNodes[ixCentNode],
                                       vecNodes[face[ixActiveSeg]]);
    
    int ixNode3 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[ixActiveSeg]],
                                       vecNodes[iterSegRight->second.ixMidpointNode[1]],
                                       vecNodes[iterSeg->second.ixMidpointNode[1]],
                                       vecNodes[ixCentNode]);

    int ixNode03 = InsertQuadNodeCenter(vecNodes,
                                        vecNodes[iterSeg->second.ixMidpointNode[0]],
                                        vecNodes[iterSeg->second.ixMidpointNode[2]],
                                        vecNodes[ixNode0],
                                        vecNodes[ixNode3]);

    int ixNode01 = InsertQuadNodeCenter(vecNodes,
                                        vecNodes[iterSegLeft->second.ixMidpointNode[0]],
                                        vecNodes[iterSegLeft->second.ixMidpointNode[2]],
                                        vecNodes[ixNode0],
                                        vecNodes[ixNode1]);

#ifndef NO_CONNECTING_FACES
    // Construct connecting face
    //-----------------------------
    if(iFirstConnectingNode == InvalidNode) {
      iFirstConnectingNode = ixNode2;
    }
    
    if(iLastConnectingNode != InvalidNode) {
      vecFaces.push_back(Face(iterSegRight->second.ixMidpointNode[1],
                              iLastConnectingNode                   ,
                              iterSegRight->second.ixMidpointNode[0],
                              ixNode2                               ));
    }
    
    // Store connecting node
    //-----------------------------
    iLastConnectingNode = ixNode2;
#endif

    // Construct faces along active segment
    //----------------------------------------
    vecFaces.push_back(Face(ixNode3                               ,
                            iterSegRight->second.ixMidpointNode[0],
                            face[ixActiveSeg]                     ,
                            iterSeg->second.ixMidpointNode[0]     ));
    
    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[0],
                            iterSeg->second.ixMidpointNode[1],
                            ixNode03                         ,
                            ixNode3                          ));
    
    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[1],
                            iterSeg->second.ixMidpointNode[2],
                            ixNode0                          ,
                            ixNode03                         ));
    
    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[2]    ,
                            face[(ixActiveSeg+1)%4]              ,
                            iterSegLeft->second.ixMidpointNode[0],
                            ixNode0                              ));
    
    // Construct faces along left segment
    //----------------------------------------
    vecFaces.push_back(Face(iterSegLeft->second.ixMidpointNode[0],
                            iterSegLeft->second.ixMidpointNode[1],
                            ixNode01                             ,
                            ixNode0                              ));
    
    vecFaces.push_back(Face(iterSegLeft->second.ixMidpointNode[1],
                            iterSegLeft->second.ixMidpointNode[2],
                            ixNode1                              ,
                            ixNode01                             ));
    
    vecFaces.push_back(Face(iterSegLeft->second.ixMidpointNode[2],
                            face[(ixActiveSeg+2)%4]              ,
                            iterSegOpp->second.ixMidpointNode[0] ,
                            ixNode1                              ));
    
    // Construct interior faces
    //----------------------------
    vecFaces.push_back(Face(ixNode03  ,
                            ixNode0   ,
                            ixCentNode,
                            ixNode3   ));
    
    vecFaces.push_back(Face(ixNode01  ,
                            ixNode1   ,
                            ixCentNode,
                            ixNode0   ));

    vecFaces.push_back(Face(ixNode1                             ,
                            iterSegOpp->second.ixMidpointNode[0],
                            ixNode2                             ,
                            ixCentNode                          ));
    
    vecFaces.push_back(Face(ixNode3                               ,
                            ixCentNode                            ,
                            ixNode2                               ,
                            iterSegRight->second.ixMidpointNode[0]));
    
    vecFaces.push_back(Face(ixNode2                               ,
                            iterSegOpp->second.ixMidpointNode[1]  ,
                            face[(ixActiveSeg+3)%4]               ,
                            iterSegRight->second.ixMidpointNode[1]));
    
    // Tag this face for deletion
    //----------------------------
    vecFaces[ixActiveFace].nTag = (-1);

} // RefinementTemplateLOWCONNOLD::ApplyInnerCornerTemplate
//=============================================================================


//=============================================================================
void RefinementTemplateLOWCONNOLD::ApplyOuterCornerTemplate(FaceVector & vecFaces    ,
                                                            NodeVector & vecNodes    ,
                                                            int          ixActiveFace,
                                                            int          ixActiveSeg ) 
{
    Face & face = vecFaces[ixActiveFace];

    if(ixActiveSeg == (-1)) { _EXCEPTIONT("Logic error"); }

    SegmentMapIterator iterSeg      = face.iterSegment[ ixActiveSeg     ];
    SegmentMapIterator iterSegLeft  = face.iterSegment[(ixActiveSeg+1)%4];
    SegmentMapIterator iterSegOpp   = face.iterSegment[(ixActiveSeg+2)%4];
    SegmentMapIterator iterSegRight = face.iterSegment[(ixActiveSeg+3)%4];
    
    // Final processed segment
    //-------------------------
    iterSegFinal = iterSegLeft;
    
    // Set up bounding edges
    //-------------------------
    if(iterSeg->first[0] == iterSegLeft->first[0]) {
    } else if(iterSeg->first[1] == iterSegLeft->first[0]) {
      iterSeg->first.Flip();
    } else if(iterSeg->first[0] == iterSegLeft->first[1]) {
      iterSegLeft->first.Flip();
    } else if(iterSeg->first[1] == iterSegLeft->first[1]) {
      iterSeg->first.Flip();
      iterSegLeft->first.Flip();
    }
    
    DissectSegment(iterSegLeft , vecNodes, 2, dTrisectSpacing);
    DissectSegment(iterSeg     , vecNodes, 2, dTrisectSpacing);
    DissectSegment(iterSegOpp  , vecNodes, 1);
    DissectSegment(iterSegRight, vecNodes, 1);
    
    // Insert new nodes
    //------------------
    int ixCentNode = InsertQuadNodeCenter(vecNodes,
                                          vecNodes[face[ixActiveSeg]],
                                          vecNodes[face[(ixActiveSeg+1)%4]],
                                          vecNodes[face[(ixActiveSeg+2)%4]],
                                          vecNodes[face[(ixActiveSeg+3)%4]]);
    
    int ixNode0 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[ixActiveSeg]],
                                       vecNodes[face[(ixActiveSeg+1)%4]],
                                       vecNodes[iterSegRight->second.ixMidpointNode[0]],
                                       vecNodes[ixCentNode]);
    
    int ixNode1 = InsertQuadNodeCenter(vecNodes,
                                       vecNodes[face[(ixActiveSeg+1)%4]],
                                       vecNodes[face[(ixActiveSeg+2)%4]],
                                       vecNodes[iterSegOpp->second.ixMidpointNode[0]],
                                       vecNodes[ixCentNode]);

#ifndef NO_CONNECTING_FACES
    // Construct connecting face
    //---------------------------
    if(iFirstConnectingNode == InvalidNode) { iFirstConnectingNode = ixNode0; }

    if(iLastConnectingNode != InvalidNode) {
      vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[1],
                              iLastConnectingNode              ,
                              iterSeg->second.ixMidpointNode[0],
                              ixNode0                          ));
    }
    iLastConnectingNode = ixNode1;
#endif

    // Construct faces along active segment
    //--------------------------------------
    vecFaces.push_back(Face(ixCentNode                            ,
                            face[(ixActiveSeg+3)%4]               ,
                            iterSegRight->second.ixMidpointNode[0],
                            ixNode0                               ));

    vecFaces.push_back(Face(face[ixActiveSeg]                     ,
                            iterSeg->second.ixMidpointNode[1]     ,
                            ixNode0                               ,
                            iterSegRight->second.ixMidpointNode[0]));

    vecFaces.push_back(Face(iterSeg->second.ixMidpointNode[0],
                            face[(ixActiveSeg+1)%4]          ,
                            ixCentNode                       ,
                            ixNode0                          ));

    vecFaces.push_back(Face(face[(ixActiveSeg+1)%4]              ,
                            iterSegLeft->second.ixMidpointNode[0],
                            ixNode1                              ,
                            ixCentNode                           ));

    vecFaces.push_back(Face(iterSegLeft->second.ixMidpointNode[1],
                            face[(ixActiveSeg+2)%4]              ,
                            iterSegOpp->second.ixMidpointNode[0] ,
                            ixNode1                              ));

    vecFaces.push_back(Face(ixCentNode                          ,
                            ixNode1                             ,
                            iterSegOpp->second.ixMidpointNode[0],
                            face[(ixActiveSeg+3)%4]             ));

    // Tag this face for deletion
    //-----------------------------
    vecFaces[ixActiveFace].nTag = (-1);

} // RefinementTemplateLOWCONNOLD::ApplyOuterCornerTemplate
//=============================================================================


//=============================================================================
void RefinementTemplateLOWCONNOLD::BeginLoop() 
{
    iFirstConnectingNode = InvalidNode;
     iLastConnectingNode = InvalidNode;
}
//=============================================================================


//=============================================================================
void RefinementTemplateLOWCONNOLD::RefineLocalGrid(int                  nRefineLevel  ,
                                                   FaceVector         & vecFaces      ,
                                                   NodeVector         & vecNodes      ,
                                                   bool                 fFirstSegment ,
                                                   bool                 fLastSegment  ,
                                                   int                  ixFaceActive  ,
                                                   int                  ixFaceInactive,
                                                   SegmentMapIterator & iterSegPrev   ,
                                                   SegmentMapIterator & iterSegNext   ,
                                                   SegmentMapIterator & iterSegCurr   )
{
    // Active and inactive faces
    //----------------------------
    Face & faceActive   = vecFaces[ixFaceActive  ];
    Face & faceInactive = vecFaces[ixFaceInactive];
    
    // Segments
    //-----------
    const Segment & segPrev = iterSegPrev->first;
    const Segment & segNext = iterSegNext->first;
    const Segment & segCurr = iterSegCurr->first;
    
    // Check for interior corners
    //-----------------------------
    if(faceActive.ContainsSegment(segNext) != (-1)) {
      ApplyInnerCornerTemplate( vecFaces, vecNodes, ixFaceActive,
      faceActive.ContainsSegment(segCurr));
    } else if (faceActive.ContainsSegment(segPrev) != (-1)) {
      // Ignore interior corner second time it appears
      //------------------------------------------------
    } else {
      // Apply standard refinement pattern
      //------------------------------------
      ApplyStandardTemplate( vecFaces, vecNodes, ixFaceActive,
      faceActive.ContainsSegment(segCurr));
    }
    
    // Entering an exterior corner
    //------------------------------------
    int nCommonNodeNext = segCurr.CommonNode(segNext);
    
    if((faceInactive.ContainsSegment(segNext) != (-1)) && (nCommonNodeNext > 7)) {
    
      // Active segment
      //------------------
      int ixOppSeg = faceActive.ContainsSegment(segCurr);
    
      if(ixOppSeg == (-1)) { _EXCEPTIONT("Logic error."); }
    
      SegmentMapIterator & iterSegActive = faceActive.iterSegment[(ixOppSeg + 1) % 4];
    
      const Segment & segActive = iterSegActive->first;
    
      // Find bounding face
      //---------------------
      int ixBoundingFace;
      if(iterSegActive->second.ixFace[0] == ixFaceActive) {
        ixBoundingFace = iterSegActive->second.ixFace[1];
      } else if(iterSegActive->second.ixFace[1] == ixFaceActive) {
        ixBoundingFace = iterSegActive->second.ixFace[0];
      } else {
        _EXCEPTIONT("Logic error");
      }
    
      Face & faceBounding = vecFaces[ixBoundingFace];
      if(faceBounding.nRefineLevel == nRefineLevel) {
        printf("WARNING (%s, %u) Refined face detected as exterior "
               "corner. Mesh may be corrupted.", __FILE__, __LINE__);
        return;
      }
    
      ApplyOuterCornerTemplate( vecFaces, vecNodes, ixBoundingFace,
      faceBounding.ContainsSegment(segActive));
    }
    
    // Last connecting element
    //--------------------------
    if(fLastSegment) {
      vecFaces.push_back(Face(iterSegFinal->second.ixMidpointNode[1],
                              iLastConnectingNode                   ,
                              iterSegFinal->second.ixMidpointNode[0],
                              iFirstConnectingNode                  ));
    }

} // RefinementTemplateLOWCONNOLD::RefineLocalGrid
//=============================================================================


