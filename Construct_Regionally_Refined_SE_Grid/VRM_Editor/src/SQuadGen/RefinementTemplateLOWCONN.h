///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefinementTemplateLOWCONN.h
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

#ifndef _REFINEMENTTEMPLATELOWCONN_H_
#define _REFINEMENTTEMPLATELOWCONN_H_

#include "GridElements.h"
// #include "RefineGrid.h"
#include "../RefineCubeGrid.h"


//=============================================================================
class RefinementTemplateLOWCONN : public RefinementTemplate 
{
    public: //  Constructor.
            //================
            RefinementTemplateLOWCONN();


    private: //  Apply the standard refinement template along a straight edge.
             //===============================================================
            void ApplyStandardTemplate( FaceVector & vecFaces    ,
                                        NodeVector & vecNodes    ,
                                        int          ixActiveFace,
                                        int          ixActiveSeg );


            //  Apply an inner corner refinement template.
            //=============================================
            void ApplyInnerCornerTemplate( FaceVector & vecFaces    ,
                                           NodeVector & vecNodes    ,
                                           int          ixActiveFace,
                                           int          ixActiveSeg );


            //  Apply an outer corner refinement template.
            //============================================
            void ApplyOuterCornerTemplate( FaceVector & vecFaces    ,
                                           NodeVector & vecNodes    ,
                                           int          ixActiveFace,
                                           int          ixActiveSeg );

    public: //  Template uses 2x2 element blocks.
            //===================================
            virtual bool UsesElementBlocks() 
            {
              return true;
            }


            //  Finalize the grid.
            //=====================
            void Initialize();


            //  Initialize a new loop.
            //========================
            void BeginLoop();


            //  Perform one grid refinement step.
            //=====================================
            void RefineLocalGrid( int                  nRefineLevel  ,
                                  FaceVector         & vecFaces      ,
                                  NodeVector         & vecNodes      ,
                                  bool                 fFirstSegment ,
                                  bool                 fLastSegment  ,
                                  int                  ixFaceActive  ,
                                  int                  ixFaceInactive,
                                  SegmentMapIterator & segiterPrev   ,
                                  SegmentMapIterator & segiterNext   ,
                                  SegmentMapIterator & segiterCurr   );
/*
            //  Finalize the grid.
            //=====================
            void Finalize( int          nRefineLevel,
                           FaceVector & vecFaces    ,
                           NodeVector & vecNodes    );
*/
            //  Finalize the refinement template after all refinement levels.
            //================================================================
            virtual void FinalizeAll( FaceVector & vecFaces, NodeVector & vecNodes);


    private:
            double  dTrisectSpacing[2];      // Spacing of nodes along trisected edges.
            int     iFirstConnectingNode;    // Node index of first connecting node in loop.
            int     iLastConnectingNode;     // Node index of last connecting node.
            SegmentMapIterator iterSegFinal; // Final processed segment.
            typedef std::set<Segment> MergePairSet;          // Pairs of nodes to merge.
            typedef MergePairSet::iterator MergePairSetIterator;
            typedef MergePairSet::value_type MergePairSetElement;
            MergePairSet setMergePairs;
            typedef std::map<int, bool> ActiveMergeMap;               // Map of nodes in merge 
            typedef ActiveMergeMap::iterator ActiveMergeIterator;     // pair map and indicator
            typedef ActiveMergeMap::value_type ActiveMergeMapElement; // of whether it should
            ActiveMergeMap mapMergeActive;                            // be merged or not.
            std::vector<int> vecSpuriousNodes; // Vector of spurious nodes.
};
//=============================================================================

#endif

