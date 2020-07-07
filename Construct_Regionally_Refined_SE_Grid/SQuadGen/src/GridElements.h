///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.h
///	\author  Paul Ullrich
///	\version February 21, 2012
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

#ifndef _GRIDELEMENTS_H_
#define _GRIDELEMENTS_H_

///////////////////////////////////////////////////////////////////////////////

#include <climits>
#include <vector>
#include <map>
#include <set>

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

struct Node {
	double x;
	double y;
	double z;

	///	<summary>
	///		Default constructor.
	///	</summary>
	Node() :
		x(0.0),
		y(0.0),
		z(0.0)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	Node(
		double _x,
		double _y,
		double _z
	) :
		x(_x),
		y(_y),
		z(_z)
	{ }
};

typedef std::vector<Node> NodeVector;

static const int InvalidNode = (-1);

///////////////////////////////////////////////////////////////////////////////

struct LonLatNode {
	double lon;
	double lat;

	///	<summary>
	///		Constructor
	///	</summary>
	LonLatNode(
		double _lon,
		double _lat
	) :
		lon(_lon),
		lat(_lat)
	{ }
};

typedef std::vector<LonLatNode> LonLatNodeVector;

///////////////////////////////////////////////////////////////////////////////

class Edge : public std::vector<int> {

public:
	///	<summary>
	///		Flip the edge.
	///	</summary>
	Edge Flip() const {
		Edge edgeFlip;
		for (int i = size()-1; i >= 0; i--) {
			edgeFlip.push_back((*this)[i]);
		}
		return edgeFlip;
	}
};

typedef std::vector<Edge> EdgeVector;

///////////////////////////////////////////////////////////////////////////////

struct Segment {
	int ixNode[2];

	///	<summary>
	///		Constructor.
	///	</summary>
	Segment(
		int ixNode0,
		int ixNode1
	) {
		if (ixNode0 == ixNode1) {
			_EXCEPTIONT("Invalid segment; endpoints must be different nodes.");
		}
		if ((ixNode0 == InvalidNode) && (ixNode1 == InvalidNode)) {
			_EXCEPTIONT("Both nodes invalid");
		}

		ixNode[0] = ixNode0;
		ixNode[1] = ixNode1;
	}

	///	<summary>
	///		Flip the order of the nodes stored in the segment.  Note that this
	///		does not affect the comparator properties of the segment, and so
	///		this method can be treated as const.
	///	</summary>
	void Flip() const {
		int ixTemp = ixNode[0];
		const_cast<int&>(ixNode[0]) = ixNode[1];
		const_cast<int&>(ixNode[1]) = ixTemp;
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return ixNode[i];
	}

	///	<summary>
	///		Get the nodes as an ordered pair.
	///	</summary>
	void GetOrderedNodes(
		int & ixNodeSmall,
		int & ixNodeBig
	) const {
		if (ixNode[0] < ixNode[1]) {
			ixNodeSmall = ixNode[0];
			ixNodeBig   = ixNode[1];
		} else {
			ixNodeSmall = ixNode[1];
			ixNodeBig   = ixNode[0];
		}
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Segment & seg) const {

		// Order the local nodes
		int ixNodeSmall;
		int ixNodeBig;
		GetOrderedNodes(ixNodeSmall, ixNodeBig);

		// Order the nodes in seg
		int ixSegNodeSmall;
		int ixSegNodeBig;
		seg.GetOrderedNodes(ixSegNodeSmall, ixSegNodeBig);

		// Compare
		if (ixNodeSmall < ixSegNodeSmall) {
			return true;
		} else if (ixNodeSmall > ixSegNodeSmall) {
			return false;
		} else if (ixNodeBig < ixSegNodeBig) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Segment & seg) const {
		if ((ixNode[0] == seg.ixNode[0]) &&
			(ixNode[1] == seg.ixNode[1])
		) {
			return true;

		} else if (
			(ixNode[0] == seg.ixNode[1]) &&
			(ixNode[1] == seg.ixNode[0])
		) {
			return true;
		}

		return false;
	}

	///	<summary>
	///		Return the node that is shared between segments.
	///	</summary>
	int CommonNode(
		const Segment & seg
	) const {
		if (seg[0] == ixNode[0]) {
			return ixNode[0];
		} else if (seg[0] == ixNode[1]) {
			return ixNode[1];
		} else if (seg[1] == ixNode[0]) {
			return ixNode[0];
		} else if (seg[1] == ixNode[1]) {
			return ixNode[1];
		} else {
			return InvalidNode;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

struct FacePair {

	int nFaces;
	int ixFace[2];

	int nMidpoints;
	int ixMidpointNode[3];

	///	<summary>
	///		Constructor.
	///	</summary>
	FacePair() :
		nFaces(0),
		nMidpoints(0)
	{
		ixMidpointNode[0] = InvalidNode;
		ixMidpointNode[1] = InvalidNode;
		ixMidpointNode[2] = InvalidNode;
	}

	///	<summary>
	///		Add a face to this segment.
	///	</summary>
	void AddFace(int ixF) {
		if (nFaces == 2) {
			_EXCEPTIONT("Segment already has a full set of faces.");
		}

		ixFace[nFaces] = ixF;
		nFaces++;
	}

	///	<summary>
	///		Does this FacePair have a complete set of faces?
	///	</summary>
	bool IsComplete() const {
		return (nFaces == 2);
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return ixFace[i];
	}
};

///////////////////////////////////////////////////////////////////////////////

typedef std::map<Segment, FacePair> SegmentMap;

typedef SegmentMap::value_type SegmentMapPair;

typedef SegmentMap::iterator SegmentMapIterator;

typedef std::set<Segment> SegmentSet;

typedef std::pair<Segment, FacePair> SegmentPair;

typedef std::vector<SegmentPair> SegmentMapVector;

///////////////////////////////////////////////////////////////////////////////

struct Face {
	int ixNode[4];

	int nSegments;
	SegmentMapIterator iterSegment[4];

	int ixSplitMidpoints[4];

	int nRefineLevel;

	int nColor;

	int nTag;

	int iPanel;
	int iA;
	int iB;

	///	<summary>
	///		Constructor.
	///	</summary>
	Face(
		int ixNode0,
		int ixNode1,
		int ixNode2,
		int ixNode3,
		int _nRefineLevel = (-1)
	) :
		nSegments(0),
		nRefineLevel(_nRefineLevel),
		nColor(0),
		nTag(0)
	{
		ixNode[0] = ixNode0;
		ixNode[1] = ixNode1;
		ixNode[2] = ixNode2;
		ixNode[3] = ixNode3;

		ixSplitMidpoints[0] = InvalidNode;
		ixSplitMidpoints[1] = InvalidNode;
		ixSplitMidpoints[2] = InvalidNode;
		ixSplitMidpoints[3] = InvalidNode;
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	inline int operator[](int ix) const {
		return ixNode[ix];
	}

	inline int & operator[](int ix) {
		return ixNode[ix];
	}

	///	<summary>
	///		Add a new segment to this face.  Segments are stored in
	///		anti-clockwise order starting with the face's first node.
	///	</summary>
	void AddSegment(
		SegmentMapIterator iterSeg
	) {
		int k;

		if (nSegments == 4) {
			_EXCEPTIONT("Invalid number of segments per face.");
		}
		for (k = 0; k < 4; k++) {
			Segment seg(ixNode[k], ixNode[(k+1)%4]);
			if (seg == iterSeg->first) {
				iterSegment[k] = iterSeg;
				break;
			}
		}
		if (k == 4) {
			_EXCEPTIONT("Segment incorrectly matched to face.");
		}
		nSegments++;
	}

	///	<summary>
	///		Returns the index of the segment if it is part of this face,
	///		otherwise returns (-1).
	///	</summary>
	int ContainsSegment(
		const Segment & seg
	) const {
		for (int k = 0; k < 4; k++) {
			if (iterSegment[k]->first == seg) {
				return k;
			}
		}
		return (-1);
	}

	///	<summary>
	///		Returns the orientation of the segment (clockwise or
	///		anti-clockwise) in its natural orientation.
	///	</summary>
	///	<returns>
	///		+1 if the orientation is clockwise
	///		0  if the segment is not part of the face
	///		-1 if the orientation is anti-clockwise
	///	</returns>
	int SegmentOriention(
		const Segment & seg
	) const {
		for (int k = 0; k < 4; k++) {
			if (seg.ixNode[0] == ixNode[k]) {
				int kx = (k+1)%4;
				int ky = (k+3)%4;
				if (seg.ixNode[1] == ixNode[kx]) {
					return (-1);
				} else if (seg.ixNode[1] == ixNode[ky]) {
					return (+1);
				} else {
					return (0);
				}
			}
		}
		return (0);
	}
};

typedef std::vector<Face> FaceVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Reverse the orientation of all faces.
///	</summary>
void ReverseFaceOrientation(
	FaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Insert a node at the given distance between two nodes.
///	</summary>
int InsertSubNode(
	int ix0,
	int ix1,
	double alpha,
	NodeVector & vecNodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Insert a node at the centerpoint of a triangular face.
///	</summary>
int InsertTriFaceCentroidNode(
	int ix0,
	int ix1,
	int ix2,
	NodeVector & vecNodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Insert a node at the centerpoint of a quadrilateral.
///	</summary>
int InsertQuadNodeCenter(
	NodeVector & vecNodes,
	const Node & node0,
	const Node & node1,
	const Node & node2,
	const Node & node3
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Insert a node at a given location within the quadrilateral.
///	</summary>
int InsertQuadNode(
	NodeVector & vecNodes,
	const Node & node0,
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dAlpha,
	double dBeta
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate equiangular spacing of nodes along the given edge.
///	</summary>
void GenerateEdgeVertices(
	int nRefineLevel,
	int ix0,
	int ix1,
	NodeVector & vecNodes,
	Edge & edge
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Construct a map containing all segments connecting faces with at least
///		the specified refinement level.
///	</summary>
void ConstructSegmentMap(
	FaceVector & vecFaces,
	SegmentMap & mapSegments,
	int nRefineLevel
);

///////////////////////////////////////////////////////////////////////////////

#endif

