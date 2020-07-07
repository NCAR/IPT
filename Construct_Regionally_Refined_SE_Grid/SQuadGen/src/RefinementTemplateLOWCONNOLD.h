///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefinementTemplateLOWCONNOLD.h
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

#ifndef _REFINEMENTTEMPLATELOWCONNOLD_H_
#define _REFINEMENTTEMPLATELOWCONNOLD_H_

#include "GridElements.h"
#include "RefineGrid.h"

///////////////////////////////////////////////////////////////////////////////

class RefinementTemplateLOWCONNOLD : public RefinementTemplate {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RefinementTemplateLOWCONNOLD();

private:
	///	<summary>
	///		Apply the standard refinement template along a straight edge.
	///	</summary>
	void ApplyStandardTemplate(
		FaceVector & vecFaces,
		NodeVector & vecNodes,
		int ixActiveFace,
		int ixActiveSeg
	);

	///	<summary>
	///		Apply an inner corner refinement template.
	///	</summary>
	void ApplyInnerCornerTemplate(
		FaceVector & vecFaces,
		NodeVector & vecNodes,
		int ixActiveFace,
		int ixActiveSeg
	);

	///	<summary>
	///		Apply an outer corner refinement template.
	///	</summary>
	void ApplyOuterCornerTemplate(
		FaceVector & vecFaces,
		NodeVector & vecNodes,
		int ixActiveFace,
		int ixActiveSeg
	);

public:
	///	<summary>
	///		Template uses 2x2 element blocks.
	///	</summary>
	virtual bool UsesElementBlocks() {
		return true;
	}

	///	<summary>
	///		Initialize a new loop.
	///	</summary>
	void BeginLoop();

	///	<summary>
	///		Perform one grid refinement step.
	///	</summary>
	void RefineLocalGrid(
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
	///		Spacing of nodes along trisected edges.
	///	</sumary>
	double dTrisectSpacing[2];

	///	<summary>
	///		Node index of first connecting node in loop.
	///	</summary>
	int iFirstConnectingNode;

	///	<summary>
	///		Node index of last connecting node.
	///	</summary>
	int iLastConnectingNode;

	///	<summary>
	///		Final processed segment.
	///	</summary>
	SegmentMapIterator iterSegFinal;
};

///////////////////////////////////////////////////////////////////////////////

#endif

