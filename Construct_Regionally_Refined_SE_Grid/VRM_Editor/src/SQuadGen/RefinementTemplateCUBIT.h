///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefinementTemplateCUBIT.h
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

#ifndef _REFINEMENTTEMPLATECUBIT_H_
#define _REFINEMENTTEMPLATECUBIT_H_

#include "GridElements.h"
// #include "RefineGrid.h"
#include "../RefineCubeGrid.h"


//=============================================================================
class RefinementTemplateCUBIT : public RefinementTemplate 
{
    public: //  Template does not use element blocks.
            //===================================
            virtual bool UsesElementBlocks() { return false; }


            //  Initialize a new loop.
            //===================================
            virtual void BeginLoop();


            //  Perform one grid refinement step.
            //===================================
            virtual void RefineLocalGrid( int                  nRefineLevel  ,
                                          FaceVector         & vecFaces      ,
                                          NodeVector         & vecNodes      ,
                                          bool                 fFirstSegment ,
                                          bool                 fLastSegment  ,
                                          int                  ixFaceActive  ,
                                          int                  ixFaceInactive,
                                          SegmentMapIterator & segiterPrev   ,
                                          SegmentMapIterator & segiterNext   ,
                                          SegmentMapIterator & segiterCurr   );

    private:
            int  iState;          // Current state.
            bool fExteriorCorner; // Currently working on an exterior corner.
            bool fInteriorCorner; // Currently working on an interior corner.
};
//=============================================================================

#endif

