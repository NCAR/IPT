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
#include "RefineGrid.h"

///////////////////////////////////////////////////////////////////////////////

class RefinementTemplateLOWCONN : public RefinementTemplate {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	RefinementTemplateLOWCONN();

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
	///		Finalize the grid.
	///	</summary>
	void Initialize();

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
/*
	///	<summary>
	///		Finalize the grid.
	///	</summary>
	void Finalize(
		int nRefineLevel,
		FaceVector & vecFaces,
		NodeVector & vecNodes
	);
*/
	///	<summary>
	///		Finalize the refinement template after all refinement levels.
	///	</summary>
	virtual void FinalizeAll(
		FaceVector & vecFaces,
		NodeVector & vecNodes
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

	///	<summary>
	///		Pairs of nodes to merge.
	///	</summary>
	typedef std::set<Segment> MergePairSet;
	typedef MergePairSet::iterator MergePairSetIterator;
	typedef MergePairSet::value_type MergePairSetElement;
	MergePairSet setMergePairs;

	///	<summary>
	///		Map of nodes in the merge pair map and indicator of whether it
	///		should be merged or not.
	///	</summary>
	typedef std::map<int, bool> ActiveMergeMap;
	typedef ActiveMergeMap::iterator ActiveMergeIterator;
	typedef ActiveMergeMap::value_type ActiveMergeMapElement;
	ActiveMergeMap mapMergeActive;

	///	<summary>
	///		Vector of spurious nodes.
	///	</summary>
	std::vector<int> vecSpuriousNodes;
};

///////////////////////////////////////////////////////////////////////////////

#endif

