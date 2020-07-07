///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereGrid.cpp
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

#include "CubedSphereGrid.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

int InsertCSSubNode(
	int ix0,
	int ix1,
	double alpha,
	NodeVector & vecNodes
) {
	double dX = vecNodes[ix0].x + (vecNodes[ix1].x - vecNodes[ix0].x) * alpha;
	double dY = vecNodes[ix0].y + (vecNodes[ix1].y - vecNodes[ix0].y) * alpha;
	double dZ = vecNodes[ix0].z + (vecNodes[ix1].z - vecNodes[ix0].z) * alpha;

	// Project to sphere
	double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);

	dX /= dRadius;
	dY /= dRadius;
	dZ /= dRadius;

	// Index
	int ix = vecNodes.size();

	// Insert node
	vecNodes.push_back(Node(dX, dY, dZ));

	return ix;
}

///////////////////////////////////////////////////////////////////////////////

void GenerateCSEdgeVertices(
	int nRefineLevel,
	int ix0,
	int ix1,
	NodeVector & vecNodes,
	Edge & edge
) {
	edge.clear();
	edge.push_back(ix0);

	int i;
	for (i = 1; i < nRefineLevel; i++) {

		// Nodes along line in Cartesian geometry
		double alpha =
			static_cast<double>(i) / static_cast<double>(nRefineLevel);

		alpha = 0.5 * (tan(0.25 * M_PI * (2.0 * alpha - 1.0)) + 1.0);

		// Insert node along edge
		int ixNode = InsertCSSubNode(ix0, ix1, alpha, vecNodes);

		// Add node to edge
		edge.push_back(ixNode);
	}

	edge.push_back(ix1);
}

///////////////////////////////////////////////////////////////////////////////

void GenerateFacesFromQuad(
	int nResolution,
	int iPanel,
	const Edge & edge0,
	const Edge & edge1,
	const Edge & edge2,
	const Edge & edge3,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	Edge edgeTop;
	Edge edgeBot = edge0;

	for (int j = 0; j < nResolution; j++) {

		// Generate top level edge
		if (j != nResolution-1) {
			int ix0 = edge1[j+1];
			int ix1 = edge2[j+1];

			GenerateCSEdgeVertices(nResolution, ix0, ix1, vecNodes, edgeTop);

		} else {
			edgeTop = edge3;
		}

		// Generate face
		for (int i = 0; i < nResolution; i++) {
			Face face(edgeBot[i+1], edgeBot[i], edgeTop[i], edgeTop[i+1], 0);
			face.iPanel = iPanel;
			face.iA = i;
			face.iB = j;

			vecFaces.push_back(face);
		}

		// Increment row
		edgeBot = edgeTop;
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateCubedSphere(
	int nResolution,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	// Generate corner points
	double dInvDeltaX = 1.0 / sqrt(3.0);

	vecNodes.push_back(Node(+dInvDeltaX, -dInvDeltaX, -dInvDeltaX));
	vecNodes.push_back(Node(+dInvDeltaX, +dInvDeltaX, -dInvDeltaX));
	vecNodes.push_back(Node(-dInvDeltaX, +dInvDeltaX, -dInvDeltaX));
	vecNodes.push_back(Node(-dInvDeltaX, -dInvDeltaX, -dInvDeltaX));
	vecNodes.push_back(Node(+dInvDeltaX, -dInvDeltaX, +dInvDeltaX));
	vecNodes.push_back(Node(+dInvDeltaX, +dInvDeltaX, +dInvDeltaX));
	vecNodes.push_back(Node(-dInvDeltaX, +dInvDeltaX, +dInvDeltaX));
	vecNodes.push_back(Node(-dInvDeltaX, -dInvDeltaX, +dInvDeltaX));

	// Generate edges
	EdgeVector vecEdges;
	vecEdges.resize(12);

	GenerateCSEdgeVertices(nResolution, 0, 1, vecNodes, vecEdges[0]);
	GenerateCSEdgeVertices(nResolution, 1, 2, vecNodes, vecEdges[1]);
	GenerateCSEdgeVertices(nResolution, 2, 3, vecNodes, vecEdges[2]);
	GenerateCSEdgeVertices(nResolution, 3, 0, vecNodes, vecEdges[3]);

	GenerateCSEdgeVertices(nResolution, 0, 4, vecNodes, vecEdges[4]);
	GenerateCSEdgeVertices(nResolution, 1, 5, vecNodes, vecEdges[5]);
	GenerateCSEdgeVertices(nResolution, 2, 6, vecNodes, vecEdges[6]);
	GenerateCSEdgeVertices(nResolution, 3, 7, vecNodes, vecEdges[7]);

	GenerateCSEdgeVertices(nResolution, 4, 5, vecNodes, vecEdges[8]);
	GenerateCSEdgeVertices(nResolution, 5, 6, vecNodes, vecEdges[9]);
	GenerateCSEdgeVertices(nResolution, 6, 7, vecNodes, vecEdges[10]);
	GenerateCSEdgeVertices(nResolution, 7, 4, vecNodes, vecEdges[11]);

	// Generate equatorial faces
	GenerateFacesFromQuad(
		nResolution,
		0,
		vecEdges[0],
		vecEdges[4],
		vecEdges[5],
		vecEdges[8],
		vecNodes,
		vecFaces);

	GenerateFacesFromQuad(
		nResolution,
		1,
		vecEdges[1],
		vecEdges[5],
		vecEdges[6],
		vecEdges[9],
		vecNodes,
		vecFaces);

	GenerateFacesFromQuad(
		nResolution,
		2,
		vecEdges[2],
		vecEdges[6],
		vecEdges[7],
		vecEdges[10],
		vecNodes,
		vecFaces);

	GenerateFacesFromQuad(
		nResolution,
		3,
		vecEdges[3],
		vecEdges[7],
		vecEdges[4],
		vecEdges[11],
		vecNodes,
		vecFaces);

	// Generate north polar face
	GenerateFacesFromQuad(
		nResolution,
		4,
		vecEdges[8],
		vecEdges[11].Flip(),
		vecEdges[9],
		vecEdges[10].Flip(),
		vecNodes,
		vecFaces);

	// Generate south polar face
	GenerateFacesFromQuad(
		nResolution,
		5,
		vecEdges[2].Flip(),
		vecEdges[3],
		vecEdges[1].Flip(),
		vecEdges[0],
		vecNodes,
		vecFaces);
}

///////////////////////////////////////////////////////////////////////////////

