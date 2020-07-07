///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefinementTemplateCUBIT.cpp
///	\author  Paul Ullrich
///	\version March 4, 2013
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

#include "RefinementTemplateCUBIT.h"


//=============================================================================
void RefinementTemplateCUBIT::BeginLoop() 
{ 
  iState          = 0;
  fExteriorCorner = false;
  fInteriorCorner = false;
}
//=============================================================================


//=============================================================================
void RefinementTemplateCUBIT::RefineLocalGrid(int                  nRefineLevel  ,
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
    //---------------------------
    Face & faceActive   = vecFaces[ixFaceActive  ];
    Face & faceInactive = vecFaces[ixFaceInactive];

    // Segments
    //----------
    const Segment & segPrev = iterSegPrev->first;
    const Segment & segNext = iterSegNext->first;
    const Segment & segCurr = iterSegCurr->first;

    // Check for an interior corner on the last element
    //--------------------------------------------------
    if ((faceActive.ContainsSegment(segNext) != (-1)) && (fLastSegment)) { return; }

    // Refine the segment
    //--------------------
    BisectSegment(iterSegCurr, vecNodes);

    // Check for interior corners
    //-----------------------------
    if((faceActive.ContainsSegment(segNext) != (-1)) && (iState == 0) 
                                                     && (!fInteriorCorner)) {
      // Tag face for refinement
      //-------------------------
      if(faceActive.nTag != 0) {
        printf("WARNING (%s, %u) Active face already tagged. "
               "Mesh may be corrupted.", __FILE__, __LINE__);
      }
      faceActive.nTag = (+1);

    // Apply standard refinement pattern (full edge previous)
    //-------------------------------------------------------
    } else if((iState == 0) && (!fInteriorCorner)) {
      // Insert centerpoint node
      //-------------------------
      int ixCentroidNode = AddFaceCenterpointNode(faceActive, vecNodes);

      // Bisect opposite face
      //----------------------
      int ixOppSeg = faceActive.ContainsSegment(iterSegCurr->first);
      if(ixOppSeg == (-1)) { _EXCEPTIONT("Logic error."); }

      SegmentMapIterator & iterSegActive = faceActive.iterSegment[(ixOppSeg + 3) % 4];

      int ixBisectNode = BisectSegment(iterSegActive, vecNodes);
      int ixCommonNode = segCurr.CommonNode(iterSegActive->first);

      // Construct new faces
      //---------------------
      vecFaces.push_back(Face( iterSegCurr->second.ixMidpointNode[0], 
                               ixCentroidNode, ixBisectNode, ixCommonNode));

      vecFaces.push_back(Face( faceActive[(ixOppSeg + 2) % 4], 
                               ixCentroidNode, iterSegCurr->second.ixMidpointNode[0], 
                               faceActive[(ixOppSeg + 1) % 4]));

      vecFaces.push_back(Face( faceActive[(ixOppSeg + 3) % 4], ixBisectNode, 
                               ixCentroidNode, faceActive[(ixOppSeg + 2) % 4]));

      // Tag face for deletion
      //-------------------------
      if(faceActive.nTag != 0) { printf("WARNING (%s, %u) Active face already tagged. "
                                        "Mesh may be corrupted.", __FILE__, __LINE__); }
      faceActive.nTag = (-1);

      // Update the state
      //-------------------
      if(faceActive.ContainsSegment(segNext) == (-1)) { iState = 1; }

    // Apply standard refinement pattern (bisected edge previous)
    //------------------------------------------------------------
    } else if(!fInteriorCorner) {

      // Insert centerpoint node
      //-------------------------
      int ixCentroidNode = AddFaceCenterpointNode(faceActive, vecNodes);

      // Bisect opposite face
      //----------------------
      int ixOppSeg = faceActive.ContainsSegment(segCurr);
      if(ixOppSeg == (-1)) { _EXCEPTIONT("Logic error."); }

      SegmentMapIterator & iterSegActive = faceActive.iterSegment[(ixOppSeg + 1) % 4];

      int ixBisectNode = BisectSegment(iterSegActive, vecNodes);
      int ixCommonNode = segCurr.CommonNode(iterSegActive->first);

      // Construct new faces
      //---------------------
      vecFaces.push_back(Face( ixBisectNode, ixCentroidNode, 
                               iterSegCurr->second.ixMidpointNode[0], ixCommonNode));

      vecFaces.push_back(Face( iterSegCurr->second.ixMidpointNode[0], ixCentroidNode, 
                               faceActive[(ixOppSeg + 3) % 4], faceActive[ixOppSeg]));

      vecFaces.push_back(Face( ixCentroidNode, ixBisectNode, 
                               faceActive[(ixOppSeg + 2) % 4], 
                               faceActive[(ixOppSeg + 3) % 4]));

      // Tag face for deletion
      //-------------------------
      if(faceActive.nTag != 0) { printf("WARNING (%s, %u) Active face already tagged. "
                                        "Mesh may be corrupted.", __FILE__, __LINE__); }
      faceActive.nTag = (-1);

      // Update the state
      //--------------------
      if(faceActive.ContainsSegment(segNext) == (-1)) { iState = 0; }
    }

    // Corners of the cubed sphere grid are not exterior corners
    //------------------------------------------------------------
    if((faceInactive.ContainsSegment(segNext) != (-1)) && 
       (           segCurr.CommonNode(segNext) < 8   )   ) {
      fExteriorCorner = false;
      fInteriorCorner = false;
      return;
    }

    // Set interior corner flag
    //---------------------------
    if(faceActive.ContainsSegment(segNext) != (-1)) {
      fInteriorCorner = true;
      return;
    } else {
      fInteriorCorner = false;
    }

    // Check for exterior corners
    //----------------------------
    fExteriorCorner = false;

    // Entering an exterior corner in state 1 (do nothing)
    //-----------------------------------------------------
    if((faceInactive.ContainsSegment(segNext) != (-1)) && (iState == 1)) {
      fExteriorCorner = true;
    }

    // Entering an exterior corner in state 0 (refine)
    //-------------------------------------------------
    if((faceInactive.ContainsSegment(segNext) != (-1)) && (iState == 0)) {

      // Active segment
      //----------------
      int ixOppSeg = faceActive.ContainsSegment(segCurr);
      if(ixOppSeg == (-1)) { _EXCEPTIONT("Logic error."); }

      SegmentMapIterator & iterSegActive = faceActive.iterSegment[(ixOppSeg + 1) % 4];

      const Segment & segActive = iterSegActive->first;

      // Find bounding face
      //--------------------
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

      // Find next edge and bisect
      //----------------------------
      int ixActiveSeg = faceBounding.ContainsSegment(segActive);
      if(ixActiveSeg == (-1)) { _EXCEPTIONT("Logic error"); }

      int   ixNextSeg = (ixActiveSeg + 1) % 4;
      const Segment & segNext = faceBounding.iterSegment[ixNextSeg]->first;

      // Add nodes to active edge and next edge
      //-----------------------------------------
      int ixActiveNode = BisectSegment(iterSegActive                      ,vecNodes);
      int ixNextNode   = BisectSegment(faceBounding.iterSegment[ixNextSeg],vecNodes);

      // Add centroid node
      //-------------------
      int ixCentroidNode = AddFaceCenterpointNode(faceBounding, vecNodes);

      // Common node
      //-------------
      int ixCommonNode = segActive.CommonNode(segNext);

      // Construct new faces
      //-----------------------
      vecFaces.push_back(Face(ixNextNode, ixCentroidNode, ixActiveNode, ixCommonNode));

      vecFaces.push_back(Face(faceBounding[ixActiveSeg], ixActiveNode, ixCentroidNode, 
                              faceBounding[(ixActiveSeg+3)%4]));

      vecFaces.push_back(Face(faceBounding[(ixActiveSeg+2)%4], 
                              faceBounding[(ixActiveSeg+3)%4], 
                              ixCentroidNode, ixNextNode));

      // Tag face for deletion
      //--------------------------
      if(faceBounding.nTag != 0) { 
        printf("WARNING (%s, %u) Bounded face already tagged. "
               "Mesh may be corrupted.", __FILE__, __LINE__);
      }
      faceBounding.nTag = (-1);
    }
}
//=============================================================================

