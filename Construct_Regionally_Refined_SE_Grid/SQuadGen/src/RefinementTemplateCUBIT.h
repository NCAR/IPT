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
#include "RefineGrid.h"

///////////////////////////////////////////////////////////////////////////////

class RefinementTemplateCUBIT : public RefinementTemplate {

public:
	///	<summary>
	///		Template does not use element blocks.
	///	</summary>
	virtual bool UsesElementBlocks() {
		return false;
	}

	///	<summary>
	///		Initialize a new loop.
	///	</summary>
	virtual void BeginLoop();

	///	<summary>
	///		Perform one grid refinement step.
	///	</summary>
	virtual void RefineLocalGrid(
		int nRefineLevel,
		FaceVector & vecFaces,
		NodeVector & vecNodes,
		bool fFirstSegment,
		bool fLastSegment,
		int ixFaceActive,
		int ixFaceInactive,
		SegmentMapIterator & segiterPrev,
		SegmentMapIterator & segiterNext,
		SegmentMapIterator & segiterCurr
	);

private:
	///	<summary>
	///		Current state.
	///	</summary>
	int iState;

	///	<summary>
	///		Currently working on an exterior corner.
	///	</summary>
	bool fExteriorCorner;

	///	<summary>
	///		Currently working on an interior corner.
	///	</summary>
	bool fInteriorCorner;
};

///////////////////////////////////////////////////////////////////////////////

#endif

