///////////////////////////////////////////////////////////////////////////////
///
///	\file    HOMMEGridGen.cpp
///	\author  Paul Ullrich
///	\version February 13, 2012
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

#include <png.h>
#include <netcdfcpp.h>

#include <string>
#include <cmath>
#include <cstdlib>

#include "GridElements.h"
#include "CubedSphereGrid.h"
#include "IcosahedralFlagGrid.h"
#include "CSRefinementMap.h"
#include "RefineGrid.h"
#include "RefinementTemplateCUBIT.h"
#include "RefinementTemplateLOWCONN.h"
#include "RefinementTemplateLOWCONNOLD.h"
#include "SpringDynamics.h"
#include "Tessellate.h"

#include "Exception.h"
#include "CommandLine.h"

///////////////////////////////////////////////////////////////////////////////

void ComputeMeshQuality(
	const NodeVector & vecNodes,
	const FaceVector & vecFaces
) {
	// Calculate the quality of the mesh using the equal edge lengths metric
	double dEdgeLengthQuality = 0.0;

	for (int i = 0; i < vecFaces.size(); i++) {
		double dMinArcLength = +10.0;
		double dMaxArcLength = -10.0;

		for (int k = 0; k < 4; k++) {
			int kx = (k+1)%4;

			double dX0 = vecNodes[vecFaces[i][k]].x;
			double dY0 = vecNodes[vecFaces[i][k]].y;
			double dZ0 = vecNodes[vecFaces[i][k]].z;

			double dX1 = vecNodes[vecFaces[i][kx]].x;
			double dY1 = vecNodes[vecFaces[i][kx]].y;
			double dZ1 = vecNodes[vecFaces[i][kx]].z;

			double dDeltaX = dX1 - dX0;
			double dDeltaY = dY1 - dY0;
			double dDeltaZ = dZ1 - dZ0;

			double dCartLength =
				sqrt(dDeltaX*dDeltaX + dDeltaY*dDeltaY + dDeltaZ*dDeltaZ);

			double dArcLength = 2.0 * asin(0.5 * dCartLength);

			if (dArcLength < dMinArcLength) {
				dMinArcLength = dArcLength;
			}
			if (dArcLength > dMaxArcLength) {
				dMaxArcLength = dArcLength;
			}
		}

		dEdgeLengthQuality += (1.0 - dMinArcLength / dMaxArcLength);
	}

	// Calculate the quality of the mesh using equal angles metric and the
	// area ratio between elements.
	double dAngleQuality = 0.0;

	double dMinArea = 4.0 * M_PI;
	double dMaxArea = 0.0;

	double dTotalArea = 0.0;

	for (int i = 0; i < vecFaces.size(); i++) {
		double dMinAngle = +10.0;
		double dMaxAngle = -10.0;

		double dCurrentArea = 0.0;

		for (int k = 0; k < 4; k++) {
			int kx = (k+1)%4;
			int ky = (k+2)%4;

			// Normals to the faces that intersect at node kx
			//n = ((y1*z2-z1*y2), (z1*x2-x1*z2), (x1*y2-y1*x2))
			double dDeltaX1 =
				+ vecNodes[vecFaces[i][k ]].y * vecNodes[vecFaces[i][kx]].z
				- vecNodes[vecFaces[i][k ]].z * vecNodes[vecFaces[i][kx]].y;
			double dDeltaY1 = 
				+ vecNodes[vecFaces[i][k ]].z * vecNodes[vecFaces[i][kx]].x
				- vecNodes[vecFaces[i][k ]].x * vecNodes[vecFaces[i][kx]].z;
			double dDeltaZ1 = 
				+ vecNodes[vecFaces[i][k ]].x * vecNodes[vecFaces[i][kx]].y
				- vecNodes[vecFaces[i][k ]].y * vecNodes[vecFaces[i][kx]].x;

			double dDeltaX2 =
				+ vecNodes[vecFaces[i][ky]].y * vecNodes[vecFaces[i][kx]].z
				- vecNodes[vecFaces[i][ky]].z * vecNodes[vecFaces[i][kx]].y;
			double dDeltaY2 = 
				+ vecNodes[vecFaces[i][ky]].z * vecNodes[vecFaces[i][kx]].x
				- vecNodes[vecFaces[i][ky]].x * vecNodes[vecFaces[i][kx]].z;
			double dDeltaZ2 = 
				+ vecNodes[vecFaces[i][ky]].x * vecNodes[vecFaces[i][kx]].y
				- vecNodes[vecFaces[i][ky]].y * vecNodes[vecFaces[i][kx]].x;

			// Dot product of intersecting lines
			double dDotProd =
				+ dDeltaX1 * dDeltaX2
				+ dDeltaY1 * dDeltaY2
				+ dDeltaZ1 * dDeltaZ2;

			double dNorm1 = sqrt(
				+ dDeltaX1 * dDeltaX1
				+ dDeltaY1 * dDeltaY1
				+ dDeltaZ1 * dDeltaZ1);

			double dNorm2 = sqrt(
				+ dDeltaX2 * dDeltaX2
				+ dDeltaY2 * dDeltaY2
				+ dDeltaZ2 * dDeltaZ2);

			dDotProd /= (dNorm1 * dNorm2);

			// Angle between these faces
			double dCurrentAngle = acos(dDotProd) - 0.5 * M_PI;

			if (dCurrentAngle < dMinAngle) {
				dMinAngle = dCurrentAngle;
			}
			if (dCurrentAngle > dMaxAngle) {
				dMaxAngle = dCurrentAngle;
			}

			// Unit sphere element area
			dCurrentArea += dCurrentAngle;
		}

		// Verify element area is positive
		if (dCurrentArea < 0.0) {
			printf("WARNING: Negative area element (possibly concave?)\n");
			const_cast<Face &>(vecFaces[i]).nColor = 0;
		}

		// Quality of this element's angles
		dAngleQuality += (1.0 - dMinAngle / dMaxAngle);

		// Compute total area
		dTotalArea += dCurrentArea;

		if (dCurrentArea < dMinArea) {
			dMinArea = dCurrentArea;
		}
		if (dCurrentArea > dMaxArea) {
			dMaxArea = dCurrentArea;
		}
	}

	double dAreaRatio = dMaxArea / dMinArea;

	// Determine if all faces are counter-clockwise
	int nNegCount = 0;
	int nPosCount = 0;
	for (int i = 0; i < vecFaces.size(); i++) {
		double dTotal = 0.0;

		bool fRotatedSpherical = false;
		if (fabs(vecNodes[vecFaces[i][0]].z) > 0.5) {
			fRotatedSpherical = true;
		}

		for (int k = 0; k < 4; k++) {
			int kx = (k+1)%4;

			double dLat1;
			double dLon1;
			double dLat2;
			double dLon2;

			if (fRotatedSpherical) {
				dLat1 = -asin(vecNodes[vecFaces[i][k]].x);
				dLon1 = atan2(
						vecNodes[vecFaces[i][k]].y,
						vecNodes[vecFaces[i][k]].z);

				dLat2 = -asin(vecNodes[vecFaces[i][kx]].x);
				dLon2 = atan2(
					vecNodes[vecFaces[i][kx]].y,
					vecNodes[vecFaces[i][kx]].z);

			} else {
				dLat1 = asin(vecNodes[vecFaces[i][k]].z);
				dLon1 = atan2(
						vecNodes[vecFaces[i][k]].y,
						vecNodes[vecFaces[i][k]].x);

				dLat2 = asin(vecNodes[vecFaces[i][kx]].z);
				dLon2 = atan2(
					vecNodes[vecFaces[i][kx]].y,
					vecNodes[vecFaces[i][kx]].x);
			}

			if ((dLon2 < -0.5 * M_PI) && (dLon1 >  0.5 * M_PI)) {
				dLon2 += 2.0 * M_PI;
			}
			if ((dLon2 >  0.5 * M_PI) && (dLon1 < -0.5 * M_PI)) {
				dLon1 += 2.0 * M_PI;
			}

			dTotal += (dLon2 - dLon1) * (dLat2 + dLat1);
		}

		if (dTotal < 0.0) {
			nNegCount++;
			const_cast<Face&>(vecFaces[i]).nColor = 1;
		} else {
			nPosCount++;
		}
	}

	// Output results
	double dMeshSize = 4.0 * static_cast<double>(vecFaces.size());

	printf("Mesh Quality Results\n");
	printf("----------------------------------------\n");
	printf("(Face Orientation):   %i / %i\n", nPosCount, nNegCount);
	printf("(Norm. Total Area):   %1.10f\n", dTotalArea / (4.0*M_PI));
	printf("(Equal Edge Lengths): %1.5f\n", dEdgeLengthQuality / dMeshSize);
	printf("(Equal Angles):       %1.5f\n", dAngleQuality / dMeshSize);
	printf("(Area Ratio):         %1.5f\n", dAreaRatio);
	printf("----------------------------------------\n");
}

///////////////////////////////////////////////////////////////////////////////

void OutputNetCDFFile(
	std::string strFilename,
	const NodeVector & vecNodes,
	const FaceVector & vecFaces
) {
	const int ParamFour = 4;
	const int ParamLenString = 33;

	// Output to a NetCDF Exodus file
	NcFile ncOut(strFilename.c_str(), NcFile::Replace);

	// Random Exodus dimensions
	NcDim * dimLenString = ncOut.add_dim("len_string", ParamLenString);
	NcDim * dimLenLine = ncOut.add_dim("len_line", 81);
	NcDim * dimFour = ncOut.add_dim("four", ParamFour);
	NcDim * dimTime = ncOut.add_dim("time_step");
	NcDim * dimDimension = ncOut.add_dim("num_dim", 3);

	// Number of nodes
	int nNodeCount = vecNodes.size();
	NcDim * dimNodes = ncOut.add_dim("num_nodes", nNodeCount);

	// Number of elements
	int nElementCount = vecFaces.size();
	NcDim * dimElements = ncOut.add_dim("num_elem", nElementCount);
	NcDim * dimNumElementBlocks = ncOut.add_dim("num_el_blk", 1);
	NcDim * dimNumQARec = ncOut.add_dim("num_qa_rec", 1);
	NcDim * dimElementBlock1 = ncOut.add_dim("num_el_in_blk1", nElementCount);
	NcDim * dimNodesPerElement = ncOut.add_dim("num_nod_per_el1", 4);
	NcDim * dimAttBlock1 = ncOut.add_dim("num_att_in_blk1", 1);

	// Global attributes
	ncOut.add_att("api_version", 4.98f);
	ncOut.add_att("version", 4.98f);
	ncOut.add_att("floating_point_word_size", 8);
	ncOut.add_att("file_size", 0);

	char szTitle[128];
	sprintf(szTitle, "cubit(%s) 01/01/2013: 00:00:00", strFilename.c_str());
	ncOut.add_att("title", szTitle);

	// Time_whole (unused)
	ncOut.add_var("time_whole", ncDouble, dimTime);

	// QA records
	char szQARecord[ParamFour][ParamLenString] = {
		"CUBIT", "13.0", "01/01/2013", "00:00:00"};

	NcVar * varQARecords =
		ncOut.add_var("qa_records", ncChar, dimNumQARec, dimFour, dimLenString);
	varQARecords->set_cur(0,0,0);
	varQARecords->put(&(szQARecord[0][0]), 1, 4, ParamLenString);

	// Coordinate names
	char szCoordNames[3][ParamLenString] = {"x", "y", "z"};
	
	NcVar * varCoordNames =
		ncOut.add_var("coor_names", ncChar, dimDimension, dimLenString);
	varCoordNames->set_cur(0,0,0);
	varCoordNames->put(&(szCoordNames[0][0]), 3, ParamLenString);

	// Element block names
	NcVar * varElementBlockNames =
		ncOut.add_var("eb_names", ncChar, dimNumElementBlocks, dimLenString);

	// Element map
	int * nElementMap = new int[nElementCount];
	for (int i = 0; i < nElementCount; i++) {
		nElementMap[i] = i+1;
	}

	NcVar * varElementMap =
		ncOut.add_var("elem_map", ncInt, dimElements);
	varElementMap->put(nElementMap, nElementCount);

	delete[] nElementMap;

	// Element block status
	int nOne = 1;

	NcVar * varElementBlockStatus =
		ncOut.add_var("eb_status", ncInt, dimNumElementBlocks);
	varElementBlockStatus->put(&nOne, 1);

	NcVar * varElementProperty =
		ncOut.add_var("eb_prop1", ncInt, dimNumElementBlocks);
	varElementProperty->put(&nOne, 1);
	varElementProperty->add_att("name", "ID");

	// Attributes
	double * dAttrib1 = new double[nElementCount];
	for (int i = 0; i < nElementCount; i++) {
		dAttrib1[i] = 1.0;
	}

	NcVar * varAttrib1 =
		ncOut.add_var("attrib1", ncDouble, dimElementBlock1, dimAttBlock1);
	varAttrib1->put(dAttrib1, nElementCount, 1);
	delete[] dAttrib1;

	// Face nodes (1-indexed)
	NcVar * varFaces =
		ncOut.add_var("connect1", ncInt, dimElementBlock1, dimNodesPerElement);

	varFaces->add_att("elem_type", "SHELL4");

	int nNodesPerElement = 4;
	int * nConnect = new int[nNodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		for (int k = 0; k < nNodesPerElement; k++) {
			nConnect[k] = vecFaces[i][k] + 1;
		}
		varFaces->set_cur(i,0);
		varFaces->put(nConnect, 1, nNodesPerElement);
	}

	delete[] nConnect;

	// Node list
	NcVar * varNodes =
		ncOut.add_var("coord", ncDouble, dimDimension, dimNodes);

	double * dCoord = new double[nNodeCount];
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = vecNodes[i].x;
	}
	varNodes->set_cur(0,0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = vecNodes[i].y;
	}
	varNodes->set_cur(1,0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = vecNodes[i].z;
	}
	varNodes->set_cur(2,0);
	varNodes->put(dCoord, 1, nNodeCount);
	delete[] dCoord;

}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	try {

	// Grid type
	std::string strGridType;

	// Type of refinement
	std::string strRefineType;

	// Initial resolution
	int nResolution;

	// Number of levels of refinement
	int nRefinementLevel;

	// Refinement file
	std::string strRefineFile;

	// Output file
	std::string strOutputFile;

	// Refinement map
	bool fLoadCSRefinementMap;

	// Type of smoothing
	std::string strSmoothType;

	// Transition region smoothing distance
	int nTransSmoothDist;

	// Number of iterations for smoothing
	int nSmoothIterations;

	// Reverse orientation
	bool fReverseOrientation = false;

	// Invert the image
	bool fInvertImage;

	// Image longitude shift
	double dImageLonBase;

	// Image latitude shift
	double dImageLatBase;

	// Grid rotation about the X axis
	double dGridXRotate;

	// Grid rotation about the Y axis
	double dGridYRotate;

	// Number of tesselations
	int nTessellations;

	// Sub-cell resolution
	int nSubCellResolution;

	// Block refine
	bool fBlockRefine;

	// Parse the command line
	BeginCommandLine()
		CommandLineStringD(strGridType, "grid_type", "CS",
			"(Options: ICO | CS)");
		CommandLineStringD(strRefineType, "refine_type", "LOWCONN",
			"(Options: LOWCONN | CUBIT | LOWCONNOLD)");
		CommandLineInt(nRefinementLevel, "refine_level", 2);
		CommandLineInt(nResolution, "resolution", 10);
		CommandLineString(strRefineFile, "refine_file", "");
		CommandLineString(strOutputFile, "output", "");
		CommandLineBool(fLoadCSRefinementMap, "loadcsrefinementmap");
		CommandLineStringD(strSmoothType, "smooth_type", "NONE",
			"(Options: NONE | SPRING | PRESSURE)");
		CommandLineIntD(nTransSmoothDist, "smooth_dist", 1,
			"(Smooth distance, -1 = smooth entire mesh)");
		CommandLineInt(nSmoothIterations, "smooth_iter", 10);
		//CommandLineBool(fReverseOrientation, "reverse_orient");
		CommandLineDouble(dImageLonBase, "lon_base", -180.0);
		CommandLineDouble(dImageLatBase, "lat_base", 0.0);
		CommandLineDouble(dGridXRotate, "x_rotate", 0.0);
		CommandLineDouble(dGridYRotate, "y_rotate", 0.0);
		CommandLineInt(nTessellations, "tessellate", 0);
		CommandLineInt(nSubCellResolution, "subcellres", 0);
		CommandLineBool(fInvertImage, "invert");
		CommandLineBool(fBlockRefine, "block_refine");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	printf("----------------------------------------\n");

	// Check for output file
	if (strOutputFile == "") {
		std::cout << argv[0] << ": No output file specified" << std::endl;
		return (-1);
	}

	// Check for odd resolution
	if (strRefineFile != "") {
		if ((strRefineType == "LOWCONN") || (strRefineType == "LOWCONNOLD")) {
			if ((nResolution % 2) == 1) {
				_EXCEPTIONT("\nERROR: "
					"Only even resolutions are supported by LOWCONN template.");
			}
			nResolution /= 2;
		}
	}

	// Uppercase the grid type
	for (int i = 0; i < strGridType.length(); i++) {
		strGridType[i] = toupper(strGridType[i]);
	}

	// Uppercase the refine type
	for (int i = 0; i < strRefineType.length(); i++) {
		strRefineType[i] = toupper(strRefineType[i]);
	}

	// Nodes of the grid
	NodeVector vecNodes;
	FaceVector vecFaces;

	// Generate grid
	printf("Generating mesh\n");
	if (strGridType == "CS") {
		GenerateCubedSphere(nResolution, vecNodes, vecFaces);
	} else if (strGridType == "ICO") {
		GenerateIcosahedralQuadGrid(nResolution, vecNodes, vecFaces);
	} else {
		_EXCEPTION1("Unknown grid type: %s\n", strGridType.c_str());
	}
	printf("..Base mesh complete\n");

	// Perform refinement (if refinement file specified)
	if ((nRefinementLevel > 0) &&
		((strRefineFile != "") || (fLoadCSRefinementMap))
	) {
		CSRefinementMap refmap(nResolution, nRefinementLevel);

		if (fLoadCSRefinementMap) {
			printf("..Loading mesh from refinement file\n");
			refmap.FromFile("refine_map.dat");
		} else {
			refmap.InitializeFromPNG(
				strRefineFile, dImageLonBase, dImageLatBase, fInvertImage, fBlockRefine);
			refmap.Normalize();
			refmap.ToFile("refine_map.dat");
		}

		// CUBIT refinement template
		if (strRefineType == "CUBIT") {
			RefinementTemplateCUBIT reftempCUBIT;
			RefineGrid(vecNodes, vecFaces, reftempCUBIT, refmap);

		// LOWCONN refinement template
		} else if (strRefineType == "LOWCONN") {
			RefinementTemplateLOWCONN reftempLOWCONN;
			RefineGrid(vecNodes, vecFaces, reftempLOWCONN, refmap);

		// LOWCONNOLD refinement template
		} else if (strRefineType == "LOWCONNOLD") {
			RefinementTemplateLOWCONNOLD reftempLOWCONNOLD;
			RefineGrid(vecNodes, vecFaces, reftempLOWCONNOLD, refmap);

		} else {
			_EXCEPTIONT("Invalid refinement type");
		}
	}

	// Tessellate
	if ((nTessellations < 0) || (nTessellations > 100)) {
		_EXCEPTIONT("Invalid number of tesselations; expected [0,100]");
	} else if (nTessellations == 0) {
		// Do nothing
	} else {
		for (int n = 0; n < nTessellations; n++) {
			Tessellate(vecNodes, vecFaces);
		}
	}

	// Add sub-cell resolution
	if (nSubCellResolution < 0) {
		_EXCEPTIONT("--subcellres must be >= 0");
	} else if (nSubCellResolution == 0) {
		// Do nothing
	} else {
		RefineEverything(vecNodes, vecFaces, nSubCellResolution+1);
	}

	// Reverse orientation
	if ((fReverseOrientation) && (strGridType != "CS")) {
		ReverseFaceOrientation(vecFaces);
	} else if ((!fReverseOrientation) && (strGridType == "CS")) {
		ReverseFaceOrientation(vecFaces);
	}

	// No smoothing
	if (strSmoothType == "NONE") {

	// Perform smoothing
	} else if (strSmoothType == "SPRING" ) {
		SpringDynamics(
			vecNodes,
			vecFaces,
			nTransSmoothDist,
			nSmoothIterations);

	// Pressure dynamics
	} else if (strSmoothType == "PRESSURE") {
		PressureDynamics(
			vecNodes,
			vecFaces,
			nTransSmoothDist,
			nSmoothIterations);

	} else {
		_EXCEPTIONT("Invalid smoothing type");
	}

	printf("Mesh refinement complete!\n");


	// Rotate around the X axis
	if (dGridXRotate != 0.0) {
	        printf("Rotating grid around X axis.\n");
		double dCosTheta = cos(dGridXRotate * M_PI / 180.0);
		double dSinTheta = sin(dGridXRotate * M_PI / 180.0);

		for (int i = 0; i < vecNodes.size(); i++) {
			double dTempY = vecNodes[i].y;
			double dTempZ = vecNodes[i].z;

			vecNodes[i].y = dCosTheta * dTempY - dSinTheta * dTempZ;
			vecNodes[i].z = dSinTheta * dTempY + dCosTheta * dTempZ;
		}
	}

	// Rotate around the Y axis
	if (dGridYRotate != 0.0) {
	        printf("Rotating grid around Y axis.\n");
		double dCosTheta = cos(dGridYRotate * M_PI / 180.0);
		double dSinTheta = sin(dGridYRotate * M_PI / 180.0);

		for (int i = 0; i < vecNodes.size(); i++) {
			double dTempX = vecNodes[i].x;
			double dTempZ = vecNodes[i].z;

			vecNodes[i].x = dCosTheta * dTempX - dSinTheta * dTempZ;
			vecNodes[i].z = dSinTheta * dTempX + dCosTheta * dTempZ;
		}
	}

	// Output number of nodes and faces
	printf("----------------------------------------\n");
	printf("Node Count:  %lu\n", vecNodes.size());
	printf("Face Count:  %lu\n", vecFaces.size());
	printf("----------------------------------------\n");

	// Compute mesh quality
	ComputeMeshQuality(vecNodes, vecFaces);
/*
	// Output the nodes and connectivity
	FILE * fpNodes = fopen("nodes.dat", "w");
	for (int i = 0; i < vecNodes.size(); i++) {
		fprintf(fpNodes, "%1.5e %1.5e %1.5e\n",
			vecNodes[i].x, vecNodes[i].y, vecNodes[i].z);
	}
	fclose(fpNodes);

	FILE * fpFaces = fopen("faces.dat", "w");
	for (int i = 0; i < vecFaces.size(); i++) {
		fprintf(fpFaces, "%i %i %i %i %i\n",
			vecFaces[i][0],
			vecFaces[i][1],
			vecFaces[i][2],
			vecFaces[i][3],
			vecFaces[i].nRefineLevel);
	}
	fclose(fpFaces);
*/
	// Write Exodus file
	OutputNetCDFFile(strOutputFile.c_str(), vecNodes, vecFaces);

	// Success
	std::cout << "HOMMEGridGen completed successfully." << std::endl;

	} catch (Exception e) {
		std::cout << e.ToString() << std::endl;
		std::cout << "HOMMEGridGen failed" << std::endl;
	}

	return (0);
}

///////////////////////////////////////////////////////////////////////////////

