///////////////////////////////////////////////////////////////////////////////
///
///	\file    Tessellate.cpp
///	\author  Paul Ullrich
///	\version December 5, 2012
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

#include "Tessellate.h"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

void Tessellate(
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	int nInitialNodeListSize = vecNodes.size();

	// Create centerpoint nodes
	for (int i = 0; i < vecFaces.size(); i++) {
		InsertQuadNodeCenter(
			vecNodes,
			vecNodes[vecFaces[i][0]],
			vecNodes[vecFaces[i][1]],
			vecNodes[vecFaces[i][2]],
			vecNodes[vecFaces[i][3]]);
	}

	// Construct tesselation
	SegmentMap mapSegment;
	ConstructSegmentMap(vecFaces, mapSegment, -1);

	vecFaces.clear();

	SegmentMapIterator iter = mapSegment.begin();
	for (; iter != mapSegment.end(); iter++) {
		Face faceNew = 
			Face(
				iter->first[1],
				nInitialNodeListSize + iter->second[0],
				iter->first[0],
				nInitialNodeListSize + iter->second[1]);

		vecFaces.push_back(faceNew);
	}
}

///////////////////////////////////////////////////////////////////////////////

void RefineEverything(
	NodeVector & vecNodes,
	FaceVector & vecFaces,
	int nResolution
) {
	// Generate segment map
	SegmentMap mapSegment;
	ConstructSegmentMap(vecFaces, mapSegment, -1);

	FaceVector vecFacesOld = vecFaces;

	// Loop over all faces
	vecFaces.clear();

	// Construct map from segments to edges
	std::map<Segment, Edge> mapEdge;

	SegmentMapIterator iter = mapSegment.begin();
	for (; iter != mapSegment.end(); iter++) {
		Edge edge;

		GenerateEdgeVertices(
			nResolution,
			iter->first[0],
			iter->first[1],
			vecNodes,
			edge);

		mapEdge.insert(std::pair<Segment, Edge>(iter->first, edge));
	}

	// Loop over all faces and refine
	for (int n = 0 ; n < vecFacesOld.size(); n++) {
		const Segment & seg0 = vecFacesOld[n].iterSegment[0]->first;
		const Segment & seg1 = vecFacesOld[n].iterSegment[1]->first;
		const Segment & seg2 = vecFacesOld[n].iterSegment[2]->first;
		const Segment & seg3 = vecFacesOld[n].iterSegment[3]->first;

		Edge edge0 = mapEdge.find(seg0)->second;
		Edge edge1 = mapEdge.find(seg1)->second;
		Edge edge3 = mapEdge.find(seg2)->second;
		Edge edge2 = mapEdge.find(seg3)->second;

		// Align bottom and left edge
		if (edge0[0] == edge1[0]) {
		} else if (edge0[0] == edge1[edge1.size()-1]) {
			edge1 = edge1.Flip();
		} else if (edge0[edge0.size()-1] == edge1[0]) {
			edge0 = edge0.Flip();
		} else if (edge0[edge0.size()-1] == edge1[edge1.size()-1]) {
			edge0 = edge0.Flip();
			edge1 = edge1.Flip();
		} else {
			_EXCEPTIONT("Logic error");
		}

		// Align bottom and right edge
		if (edge0[edge0.size()-1] == edge2[0]) {
		} else if (edge0[edge0.size()-1] == edge2[edge2.size()-1]) {
			edge2 = edge2.Flip();
		} else {
			_EXCEPTIONT("Logic error");
		}

		// Align top and left edge
		if (edge1[edge1.size()-1] == edge3[0]) {
		} else if (edge1[edge1.size()-1] == edge3[edge3.size()-1]) {
			edge3 = edge3.Flip();
		} else {
			_EXCEPTIONT("Logic error");
		}

		Edge edgeTop;
		Edge edgeBot = edge0;

		for (int j = 0; j < nResolution; j++) {

			// Generate top level edge
			if (j != nResolution-1) {
				int ix0 = edge1[j+1];
				int ix1 = edge2[j+1];

				GenerateEdgeVertices(nResolution, ix0, ix1, vecNodes, edgeTop);

			} else {
				edgeTop = edge3;
			}

			// Generate face
			for (int i = 0; i < nResolution; i++) {
				Face face(
					edgeBot[i+1], edgeBot[i],
					edgeTop[i], edgeTop[i+1],
					vecFacesOld[n].nRefineLevel);

				face.nColor = vecFacesOld[n].nColor;
				face.nTag = vecFacesOld[n].nTag;

				vecFaces.push_back(face);
			}

			// Increment row
			edgeBot = edgeTop;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

