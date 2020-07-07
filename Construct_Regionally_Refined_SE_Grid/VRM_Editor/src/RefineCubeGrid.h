#ifndef REFINECUBEGRID_H
#define REFINECUBEGRID_H

#include "SQuadGen/GridElements.h"
#include "RefinementCube.h"

class RefinementCube;


//===========================================================================
class RefinementTemplate 
{
    public:
      // Return true if template uses element blocks.
      //---------------------------------------------
      virtual bool UsesElementBlocks() = 0;

      // Initialize the refinement template.
      //-----------------------------------
      virtual void Initialize() { }

      // Initialize a new loop.
      //------------------------
      virtual void BeginLoop() = 0;

      // Perform one grid refinement step.
      //----------------------------------
      virtual void RefineLocalGrid( int                 nRefineLevel  ,
                                    FaceVector         &vecFaces      ,
                                    NodeVector         &vecNodes      ,
                                    bool                fFirstSegment ,
                                    bool                fLastSegment  ,
                                    int                 ixFaceActive  ,
                                    int                 ixFaceInactive,
                                    SegmentMapIterator &segiterPrev   ,
                                    SegmentMapIterator &segiterNext   ,
                                    SegmentMapIterator &segiterCurr   ) = 0;
/*
      // Finalize the refinement template.
      //----------------------------------
      virtual void Finalize( int         nRefineLevel,
                             FaceVector &vecFaces    ,
                             NodeVector &vecNodes    ) { }
*/
      // Finalize the refinement template after all refinement levels.
      //-------------------------------------------------------------
      virtual void FinalizeAll( FaceVector &vecFaces,
                                NodeVector &vecNodes) { }

};
//===========================================================================

//===========================================================================
//  Perform cuts along the given segment at the specified spacing.
void DissectSegment( SegmentMapIterator &iter    ,
                     NodeVector         &vecNodes,
                     int                 nCuts   ,
                     double             *dSpacing);
//===========================================================================


//===========================================================================
//  Perform cuts along the given segment at equiangular spacing.
void DissectSegment( SegmentMapIterator &iter    ,
                     NodeVector         &vecNodes,
                     int                 nCuts   );
//===========================================================================


//===========================================================================
//  Split a segment into two segments at equiangular spacing.
int BisectSegment( SegmentMapIterator &iter    ,
                   NodeVector         &vecNodes);
//===========================================================================


//===========================================================================
//  Add a node at the face centerpoint.
int AddFaceCenterpointNode( const Face &face    ,
                            NodeVector &vecNodes);
//===========================================================================


//===========================================================================
//  Merge two nodes with ixNodeSource > ixNodeTarget.
void MergeNodes( NodeVector &vecNodes    ,
                 FaceVector &vecFaces    ,
                 int         ixNodeTarget,
                 int         ixNodeSource);
//===========================================================================


//===========================================================================
//  Perform 2x2 refinement on the specified face.
void RefineFace( NodeVector &vecNodes,
                 FaceVector &vecFaces,
                 int         ixFace  );
//===========================================================================


//===========================================================================
//  Refine the grid with the given reference map.
void RefineGrid( NodeVector            &vecNodes,
                 FaceVector            &vecFaces,
                 RefinementTemplate    &refTemp ,
                 const RefinementCube  &refCube );
//===========================================================================

#endif // REFINECUBEGRID_H
