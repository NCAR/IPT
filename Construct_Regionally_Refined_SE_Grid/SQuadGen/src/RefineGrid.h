///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefineGrid.h
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

#ifndef _REFINEGRID_H_
#define _REFINEGRID_H_

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

class CSRefinementMap;

///////////////////////////////////////////////////////////////////////////////

class RefinementTemplate {

public:
	///	<summary>
	///		Return true if template uses element blocks.
	///	</summary>
	virtual bool UsesElementBlocks() = 0;

	///	<summary>
	///		Initialize the refinement template.
	///	</summary>
	virtual void Initialize() {
	}

	///	<summary>
	///		Initialize a new loop.
	///	</summary>
	virtual void BeginLoop() = 0;

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
	) = 0;
/*
	///	<summary>
	///		Finalize the refinement template.
	///	</summary>
	virtual void Finalize(
		int nRefineLevel,
		FaceVector & vecFaces,
		NodeVector & vecNodes
	) {
	}
*/
	///	<summary>
	///		Finalize the refinement template after all refinement levels.
	///	</summary>
	virtual void FinalizeAll(
		FaceVector & vecFaces,
		NodeVector & vecNodes
	) {
	}

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Perform cuts along the given segment at the specified spacing.
///	</summary>
void DissectSegment(
	SegmentMapIterator & iter,
	NodeVector & vecNodes,
	int nCuts,
	double * dSpacing
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Perform cuts along the given segment at equiangular spacing.
///	</summary>
void DissectSegment(
	SegmentMapIterator & iter,
	NodeVector & vecNodes,
	int nCuts
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Split a segment into two segments at equiangular spacing.
///	</summary>
int BisectSegment(
	SegmentMapIterator & iter,
	NodeVector & vecNodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Add a node at the face centerpoint.
///	</summary>
int AddFaceCenterpointNode(
	const Face & face,
	NodeVector & vecNodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Merge two nodes with ixNodeSource > ixNodeTarget.
///	</summary>
void MergeNodes(
	NodeVector & vecNodes,
	FaceVector & vecFaces,
	int ixNodeTarget,
	int ixNodeSource
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Perform 2x2 refinement on the specified face.
///	</summary>
void RefineFace(
	NodeVector & vecNodes,
	FaceVector & vecFaces,
	int ixFace
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Refine the grid with the given reference map.
///	</summary>
void RefineGrid(
	NodeVector & vecNodes,
	FaceVector & vecFaces,
	RefinementTemplate & reftemp,
	const CSRefinementMap & refmap
);

///////////////////////////////////////////////////////////////////////////////

#endif

