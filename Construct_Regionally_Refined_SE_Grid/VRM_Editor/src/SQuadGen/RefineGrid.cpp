///////////////////////////////////////////////////////////////////////////////
///
///	\file    RefineGrid.cpp
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

// #include "RefineGrid.h"
// #include "GridElements.h"
// #include "CSRefinementMap.h"
// 
// #include <cmath>
// #include <iostream>
// #include <fstream>
// 
// 
// //=============================================================================
// void DissectSegment( SegmentMapIterator &iter    ,
//                      NodeVector         &vecNodes,
//                      int                 nCuts   ,
//                      double             *dSpacing) 
// {
//     if ((nCuts < 1) || (nCuts > 3)) { _EXCEPTIONT("Invalid number of cuts"); }
//     if ((iter->second.nMidpoints != 0    ) && 
//         (iter->second.nMidpoints != nCuts)   ) {
//       _EXCEPTION2("\n\tSegment has already been dissected differently" 
//                   "\n\t(Request %i : Current %i)", nCuts, iter->second.nMidpoints);
//       return;
//     }
//     if ((nCuts != 0) && (iter->second.nMidpoints == nCuts)) { return; }
//     
//     // Number of midpoints
//     //----------------------
//     iter->second.nMidpoints = nCuts;
//     
//     // Generate nodes
//     //----------------------
//     for (int i = 0; i < nCuts; i++) {
//     
//       // Nodes along line in Cartesian geometry
//       //---------------------------------------
//       double alpha = dSpacing[i];
//     
//       // Insert node along edge
//       //---------------------------------------
//       const Segment &seg = iter->first;
//     
//       int ixNode = InsertSubNode(seg[0], seg[1], alpha, vecNodes);
//     
//       // Add node to edge
//       //-----------------
//       iter->second.ixMidpointNode[i] = ixNode;
//     }
// 
//     // End Function DissectSegment()
//     //------------------------------
//     return;
// }
// //=============================================================================
// 
// 
// //=============================================================================
// void DissectSegment( SegmentMapIterator &iter    ,
//                      NodeVector         &vecNodes,
//                      int                 nCuts   ) 
// {
//     if ((nCuts < 1) || (nCuts > 3)) { _EXCEPTIONT("Invalid number of cuts"); }
//     if ((iter->second.nMidpoints != 0    ) && 
//         (iter->second.nMidpoints != nCuts)   ) {
//       _EXCEPTION2("\n\tSegment has already been dissected differently"
//                   "\n\t(Request %i : Current %i)", nCuts, iter->second.nMidpoints);
//       return;
//     }
//     if ((nCuts != 0) && (iter->second.nMidpoints == nCuts)) { return; }
//     
//     const Segment &seg = iter->first;
//     Edge  edge;
//     int   ixNewNode = vecNodes.size();
// 
//     GenerateEdgeVertices(nCuts+1, seg[0], seg[1], vecNodes, edge);
//     
//     iter->second.nMidpoints = nCuts;
//     for (int i = 0; i < nCuts; i++) {
//       iter->second.ixMidpointNode[i] = edge[i+1];
//     }
//     
//     // End Function DissectSegment()
//     //------------------------------
//     return;
// }
// //=============================================================================
// 
// 
// //=============================================================================
// int BisectSegment( SegmentMapIterator &iter    ,
//                    NodeVector         &vecNodes) 
// {
//     DissectSegment(iter, vecNodes, 1);
//     return iter->second.ixMidpointNode[0];
// }
// //=============================================================================
// 
// 
// //=============================================================================
// int AddFaceCenterpointNode( const Face &face    ,
//                             NodeVector &vecNodes) 
// {
//     double dX = 0.25*( +vecNodes[face[0]].x +vecNodes[face[1]].x 
//                        +vecNodes[face[2]].x +vecNodes[face[3]].x);
//     double dY = 0.25*( +vecNodes[face[0]].y +vecNodes[face[1]].y 
//                        +vecNodes[face[2]].y +vecNodes[face[3]].y);
//     double dZ = 0.25*( +vecNodes[face[0]].z +vecNodes[face[1]].z 
//                        +vecNodes[face[2]].z +vecNodes[face[3]].z);
// 
//     double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);
//     dX /= dRadius;
//     dY /= dRadius;
//     dZ /= dRadius;
// 
//     int ixNode = vecNodes.size();
//     vecNodes.push_back(Node(dX, dY, dZ));
// 
//     return ixNode;
// }
// //=============================================================================
// 
// 
// //=============================================================================
// void MergeNodes( NodeVector &vecNodes    ,
//                  FaceVector &vecFaces    ,
//                  int         ixNodeTarget,
//                  int         ixNodeSource) 
// {
//     if (ixNodeSource == ixNodeTarget) { _EXCEPTIONT("Invalid source and target node"); }
// 
//     for (unsigned int i = 0; i < vecFaces.size(); i++) {
//     for (         int c = 0; c < 4              ; c++) {
//       if (vecFaces[i][c] == ixNodeSource) { vecFaces[i][c] = ixNodeTarget; }
//     }
//     }
// 
//     // Merge coordinates
//     //------------------
//     vecNodes[ixNodeTarget].x = 0.5*( + vecNodes[ixNodeTarget].x + vecNodes[ixNodeSource].x);
//     vecNodes[ixNodeTarget].y = 0.5*( + vecNodes[ixNodeTarget].y + vecNodes[ixNodeSource].y);
//     vecNodes[ixNodeTarget].z = 0.5*( + vecNodes[ixNodeTarget].z + vecNodes[ixNodeSource].z);
// 
//     double dRadius = sqrt( + vecNodes[ixNodeTarget].x*vecNodes[ixNodeTarget].x 
//                            + vecNodes[ixNodeTarget].y*vecNodes[ixNodeTarget].y 
//                            + vecNodes[ixNodeTarget].z*vecNodes[ixNodeTarget].z);
//     
//     vecNodes[ixNodeTarget].x /= dRadius;
//     vecNodes[ixNodeTarget].y /= dRadius;
//     vecNodes[ixNodeTarget].z /= dRadius;
// }
// //=============================================================================
// 
// 
// //=============================================================================
// void RefineFace( NodeVector &vecNodes,
//                  FaceVector &vecFaces,
//                  int         ixFace  ) 
// {
//     // Reference to face
//     //-----------------------
//     Face & faceOld = vecFaces[ixFace];
// 
//     // Add midpoint nodes between all adjacent segments
//     //--------------------------------------------------
//     for (int k = 0; k < 4; k++) {
//       if ((faceOld.iterSegment[k]->second.nMidpoints == 0          ) &&
//           (faceOld.ixSplitMidpoints[k]               == InvalidNode)   ) {
//         BisectSegment(faceOld.iterSegment[k], vecNodes);
//       }
//     }
// 
//     // Add new faces
//     //-------------------
//     int ixMidpointNodes[4];
//     int ixSplitMidpoints[8];
//     for (int k = 0; k < 8; k++) {
//       ixSplitMidpoints[k] = InvalidNode;
//     }
// 
//     for (int k = 0; k < 4; k++) {
//       if (faceOld.ixSplitMidpoints[k] != InvalidNode) {
//         ixMidpointNodes[k] = faceOld.ixSplitMidpoints[k];
//       } else if (faceOld.iterSegment[k]->second.nMidpoints == 1) {
//         ixMidpointNodes[k] = faceOld.iterSegment[k]->second.ixMidpointNode[0];
//       } else if (faceOld.iterSegment[k]->second.nMidpoints == 3) {
//         ixMidpointNodes[k] = faceOld.iterSegment[k]->second.ixMidpointNode[1];
//     
//         // Orientation of segment relative to face determines
//         // order of split midpoints.
//         //--------------------------------------------------------
//         if (faceOld.iterSegment[k]->first[0] == faceOld[k]) {
//           ixSplitMidpoints[2*k  ] = faceOld.iterSegment[k]->second.ixMidpointNode[0];
//           ixSplitMidpoints[2*k+1] = faceOld.iterSegment[k]->second.ixMidpointNode[2];
//         } else if (faceOld.iterSegment[k]->first[1] == faceOld[k]) {
//           ixSplitMidpoints[2*k  ] = faceOld.iterSegment[k]->second.ixMidpointNode[2];
//           ixSplitMidpoints[2*k+1] = faceOld.iterSegment[k]->second.ixMidpointNode[0];
//         } else {
//           _EXCEPTIONT("Logic error");
//         }
//       } else {
//         _EXCEPTIONT("Logic error");
//       }
//     }
// 
//     // Add new node at face centroid
//     //-----------------------------------
//     int ixCentroidNode  = AddFaceCenterpointNode(faceOld, vecNodes);
//     int nNewRefineLevel = vecFaces[ixFace].nRefineLevel;
//     
//     // Insert new faces by bisecting existing face into quarters
//     //---------------------------------------------------------------
//     Face faceSE(faceOld[0]        , ixMidpointNodes[0], ixCentroidNode    , 
//                                     ixMidpointNodes[3], nNewRefineLevel   );
//     Face faceSW(ixMidpointNodes[0], faceOld[1]        , ixMidpointNodes[1], 
//                                     ixCentroidNode    , nNewRefineLevel   );
//     Face faceNW(ixCentroidNode    , ixMidpointNodes[1], faceOld[2]        , 
//                                     ixMidpointNodes[2], nNewRefineLevel   );
//     Face faceNE(ixMidpointNodes[3], ixCentroidNode    , ixMidpointNodes[2], 
//                                     faceOld[3]        , nNewRefineLevel   );
//     
//     // Add split midpoints to face
//     //---------------------------------------
//     faceSE.ixSplitMidpoints[3] = ixSplitMidpoints[7];
//     faceSE.ixSplitMidpoints[0] = ixSplitMidpoints[0];
//     faceSW.ixSplitMidpoints[0] = ixSplitMidpoints[1];
//     faceSW.ixSplitMidpoints[1] = ixSplitMidpoints[2];
//     faceNW.ixSplitMidpoints[1] = ixSplitMidpoints[3];
//     faceNW.ixSplitMidpoints[2] = ixSplitMidpoints[4];
//     faceNE.ixSplitMidpoints[2] = ixSplitMidpoints[5];
//     faceNE.ixSplitMidpoints[3] = ixSplitMidpoints[6];
//     
//     // Add cubed-sphere coordinates to face
//     //---------------------------------------
//     faceSE.iPanel = faceOld.iPanel;
//     faceSW.iPanel = faceOld.iPanel;
//     faceNW.iPanel = faceOld.iPanel;
//     faceNE.iPanel = faceOld.iPanel;
//     
//     faceSE.iA = (faceOld.iA * 2) + 1;
//     faceSW.iA = (faceOld.iA * 2);
//     faceNW.iA = (faceOld.iA * 2);
//     faceNE.iA = (faceOld.iA * 2) + 1;
//     
//     faceSE.iB = (faceOld.iB * 2);
//     faceSW.iB = (faceOld.iB * 2);
//     faceNW.iB = (faceOld.iB * 2) + 1;
//     faceNE.iB = (faceOld.iB * 2) + 1;
//     
//     vecFaces.push_back(faceSE);
//     vecFaces.push_back(faceSW);
//     vecFaces.push_back(faceNW);
//     vecFaces.push_back(faceNE);
//     
//     // Tag face for deletion
//     //-----------------------
//     vecFaces[ixFace].nTag = (-1);
//
//     // End Function RefineFace()
//     //------------------------------
// }
// //=============================================================================
//  
// 
// //=============================================================================
// void DeleteTaggedFaces( FaceVector & vecFaces) 
// {
//     int nDeleteCount = 0;
// 
//     for (unsigned int i = 0; i < vecFaces.size(); i++) {
//       if (vecFaces[i].nTag == (-1)) { nDeleteCount ++; }
//     }
//     
//     // Delete tagged faces
//     //----------------------
//     int nOffset = 0;
// 
//     for (unsigned int i = 0; i + nOffset < vecFaces.size(); i++) {
//       if (vecFaces[i + nOffset].nTag == (-1)) {
//         nOffset++;
//         for (; i + nOffset < vecFaces.size(); nOffset++) {
//           if (vecFaces[i + nOffset].nTag != (-1)) { break; }
//         }
//       }
//       if (nOffset != 0) { vecFaces[i] = vecFaces[i + nOffset]; }
//     }
//     vecFaces.erase(vecFaces.end() - nOffset, vecFaces.end());
// }
// //=============================================================================
// 
// 
// //=============================================================================
// void RefineGridLevel( NodeVector         &vecNodes    ,
//                       FaceVector         &vecFaces    ,
//                       SegmentMap         &mapSegments ,
//                       RefinementTemplate &reftemp     ,
//                       int                 nRefineLevel) 
// {
//     // Construct refinement boundary
//     //--------------------------------
//     std::vector<SegmentMapIterator> vecRefineSegments;
//     
//     SegmentMapIterator iter;
//     for (iter = mapSegments.begin(); iter != mapSegments.end(); iter++) {
//       if (!(iter->second.IsComplete())) { continue; }
//     
//       int nRefineLevel0 = vecFaces[iter->second.ixFace[0]].nRefineLevel;
//       int nRefineLevel1 = vecFaces[iter->second.ixFace[1]].nRefineLevel;
//     
//       if( (nRefineLevel0 != nRefineLevel1)   &&
//          ((nRefineLevel0 == nRefineLevel ) ||
//           (nRefineLevel1 == nRefineLevel )  )  ) { 
//         vecRefineSegments.push_back(iter); 
//       }
//     }
//     printf("Outline size: %i\n", vecRefineSegments.size());
//     
//     // Loop over all refinement regions
//     //----------------------------------
//     for (;;) {
//     
//       // Check that all refined segments have been accounted for
//       //--------------------------------------------------------
//       if (vecRefineSegments.size() == 0) { break; }
//     
//       // Put segments in order around boundary
//       //------------------------------------------
//       int nRefineSegments = vecRefineSegments.size();
//     
//       // Construct a new vector of ordered segments which traverse the
//       // refinement region in counter-clockwise order.
//       //--------------------------------------------------------------
//       std::vector<SegmentMapIterator> vecOrderedSegments;
//     
//       int ixFirstNode = vecRefineSegments[0]->first.ixNode[0];
//       int ixNextNode  = ixFirstNode;
//       for (int i = 0; i < nRefineSegments; i++) {
//     
//         // Search through refine segments list to find next piece of loop
//         //---------------------------------------------------------------
//         bool fFound = false;
//         for (unsigned int j = 0; j < vecRefineSegments.size(); j++) {
//           if ((vecRefineSegments[j]->first.ixNode[0] == ixNextNode) ||
//               (vecRefineSegments[j]->first.ixNode[1] == ixNextNode)) {
//             if (vecRefineSegments[j]->first.ixNode[0] == ixNextNode) {
//               ixNextNode = vecRefineSegments[j]->first.ixNode[1];
//             } else {
//               ixNextNode = vecRefineSegments[j]->first.ixNode[0];
//               vecRefineSegments[j]->first.Flip();
//             }
//     
//             // Insert a midpoint node
//             //-------------------------
//             vecOrderedSegments.push_back(vecRefineSegments[j]);
//     
//             vecRefineSegments.erase(vecRefineSegments.begin() + j);
//     
//             fFound = true;
//             break;
//           }
//         }
//     
//         if (!fFound) { _EXCEPTIONT("Cannot find connecting edge."); }
//     
//         // Found a complete loop
//         //-------------------------
//         if (ixNextNode == ixFirstNode) { break; }
//       }
//     
//       // Size of the loop
//       //------------------
//       int nOrderedSegmentsCount = vecOrderedSegments.size();
//     
//       // Check orientation of loop
//       //-------------------------------
//       int ixFace0 = vecOrderedSegments[0]->second.ixFace[0];
//       int ixFace1 = vecOrderedSegments[0]->second.ixFace[1];
//       int nOrientation;
// 
//       if (vecFaces[ixFace0].nRefineLevel > vecFaces[ixFace1].nRefineLevel) {
//         nOrientation = vecFaces[ixFace0].SegmentOriention( vecOrderedSegments[0]->first);
//       } else {
//         nOrientation = vecFaces[ixFace1].SegmentOriention( vecOrderedSegments[0]->first);
//       }
//     
//       // Reverse orientation of loop if needed
//       //----------------------------------------
//       if (nOrientation == (-1)) {
//         for (int i = 0; i < nOrderedSegmentsCount/2; i++) {
//           SegmentMapIterator iterTemp = vecOrderedSegments[i];
//           vecOrderedSegments[i]       = vecOrderedSegments[nOrderedSegmentsCount-i-1];
//           vecOrderedSegments[nOrderedSegmentsCount-i-1] = iterTemp;
//         }
//         nOrientation = - nOrientation;
//       }
//     
//       printf("..Loop size: %i\n", vecOrderedSegments.size());
//     //printf("..Orientation: %i\n", nOrientation);
//     
//       // March around exterior and refine
//       //-----------------------------------
//       reftemp.BeginLoop();
//     
//       for (int i = 0; i < nOrderedSegmentsCount; i++) {
//     
//         // Next ordered segment
//         //-------------------------
//         int iPrev = (i + nOrderedSegmentsCount - 1)%nOrderedSegmentsCount;
//         int iNext = (i + 1                        )%nOrderedSegmentsCount;
//     
//     //  const Segment & segPrev = vecOrderedSegments[iPrev]->first;
//     //  const Segment & segCurr = vecOrderedSegments[i]->first;
//     //  const Segment & segNext = vecOrderedSegments[iNext]->first;
//     
//         // Ensure that adequate space is available in the face vector
//         //-----------------------------------------------------------
//         vecFaces.reserve(vecFaces.size() + 16);
//     
//         // Active face is immediately outside the refinement region.
//         // Inactive face is immediately inside the refinement region.
//         //-------------------------------------------------------------
//         int ixInactiveFace;
//         int ixActiveFace;
//         if (vecFaces[vecOrderedSegments[i]->second[0]].nRefineLevel <
//             vecFaces[vecOrderedSegments[i]->second[1]].nRefineLevel   ) {
//           ixActiveFace   = vecOrderedSegments[i]->second[0];
//           ixInactiveFace = vecOrderedSegments[i]->second[1];
//         } else if ( vecFaces[vecOrderedSegments[i]->second[0]].nRefineLevel >
//                     vecFaces[vecOrderedSegments[i]->second[1]].nRefineLevel   ) {
//           ixActiveFace   = vecOrderedSegments[i]->second[1];
//           ixInactiveFace = vecOrderedSegments[i]->second[0];
//         } else {
//           _EXCEPTIONT("Faces on both sides of segment are refined.");
//         }
// 
//         // Perform one grid refinement step
//         //-----------------------------------
//         reftemp.RefineLocalGrid( nRefineLevel             ,
//                                  vecFaces                 ,
//                                  vecNodes                 ,
//                                  (i == 0)                 ,
//                                  (i == nOrderedSegmentsCount-1),
//                                  ixActiveFace             ,
//                                  ixInactiveFace           ,
//                                  vecOrderedSegments[iPrev],
//                                  vecOrderedSegments[iNext],
//                                  vecOrderedSegments[i    ]);
// /*
//         if (nRefineLevel == 1) {
//           const Segment & segCurr = vecOrderedSegments[i]->first;
//           printf("%i %i;", segCurr.ixNode[0], segCurr.ixNode[1]);
//     //    printf("%i %i %i %i;\n", segCurr.ixNode[0], segCurr.ixNode[1],
//     //	  vecOrderedSegments[i]->second.ixFace[0],
//     //	  vecOrderedSegments[i]->second.ixFace[1]);
//         }
// */
//       }
//     } //  for (;;) {
//     
//     // Refine all tagged faces
//     //------------------------------
//     printf("..Performing cleanup\n");
//     
//     for (unsigned int i = 0; i < vecFaces.size(); i++) {
//       if ((vecFaces[i].nColor != 0) && (vecFaces[i].nTag == 0)) 
//       {
//         printf("Invalid face %i\n", i);
//       }
//     }
//     
//     int nVecFacesOriginalSize = vecFaces.size();
//     for (int i = 0; i < nVecFacesOriginalSize; i++) {
//       if (vecFaces[i].nRefineLevel == nRefineLevel) {
//         RefineFace(vecNodes, vecFaces, i);
//         vecFaces[i].nTag = (-1);
//       }
//       if (vecFaces[i].nTag == (+1)) {
//         vecFaces[i].nRefineLevel = vecFaces[i].nRefineLevel + 1;
//         RefineFace(vecNodes, vecFaces, i);
//         vecFaces[i].nTag = (-1);
//       }
//     }
// 
//     // Delete tagged faces
//     //----------------------
// /*
//     nVecFacesOriginalSize = vecFaces.size();
//     for (int i = 0; i < nVecFacesOriginalSize; i++) {
//       if (vecFaces[i].nTag == (-1)) {
//         vecFaces.erase(vecFaces.begin() + i);
//         i--;
//       }
//     }
// */
//   //printf("Size: %i\n", vecFaces.size());
// 
//     DeleteTaggedFaces(vecFaces);
// 
//     // End Function void RefineGridLevel()
//     //---------------------------------------
// }
// //=============================================================================
// 
// 
// //=============================================================================
// void RefineGrid( NodeVector            &vecNodes,
//                  FaceVector            &vecFaces,
//                  RefinementTemplate    &reftemp ,
//                  const CSRefinementMap &refmap  ) 
// {
//     // Verify resolution of base mesh
//     //-------------------------------
//     int nBaseResolution = refmap.GetBaseResolution();
//     if (vecFaces.size() != 6 * nBaseResolution * nBaseResolution) {
//     _EXCEPTIONT("Inconsistency between RefinementMap and base mesh");
//     }
//     
//     // Get refinement map
//     //---------------------
//     const DataMatrix3D<int> & datamap = refmap.GetMap();
//     
//     // Loop over all refinement levels
//     //---------------------------------
//     int nActiveFineRatio = datamap.GetColumns() / nBaseResolution;
//     
//     int iRefineLevel = 0;
//     for (; iRefineLevel < refmap.GetMaxRefineLevel(); iRefineLevel++) {
//     
//       // Refinement level
//       //------------------
//       printf("Refining level %i...\n", iRefineLevel);
//     
//       // Initialize
//       //----------------
//       reftemp.Initialize();
//     
//       // Tag base elements
//       //---------------------
//       for (unsigned int i = 0; i < vecFaces.size(); i++) {
//         if (vecFaces[i].nRefineLevel < iRefineLevel) { continue; }
// 
//         int iPanel = vecFaces[i].iPanel;
//         int iA     = vecFaces[i].iA*nActiveFineRatio;
//         int iB     = vecFaces[i].iB*nActiveFineRatio;
// 
//         if (datamap[iPanel][iA][iB] > iRefineLevel) {
//           vecFaces[i].nRefineLevel = iRefineLevel + 1;
//         }
//       }
// 
//       // Generate segment map
//       //------------------------
//       SegmentMap mapSegments;
//       ConstructSegmentMap(vecFaces, mapSegments, iRefineLevel);
// 
//       printf("Segments: %i\n", mapSegments.size());
// 
//       // Refine this grid level
//       //------------------------
//       RefineGridLevel(vecNodes,vecFaces,mapSegments,reftemp,iRefineLevel+1);
//       
//       // If element blocks are used, refine to true resolution
//       //---------------------------------------------------------
//       if (reftemp.UsesElementBlocks()) {
//         int nVecFacesOriginalSize = vecFaces.size();
// 
//         // Refine lower levels
//         //-----------------------
//         for (int i = 0; i < nVecFacesOriginalSize; i++) {
//           if (vecFaces[i].nRefineLevel == iRefineLevel) {
//             RefineFace(vecNodes, vecFaces, i);
//             vecFaces[i].nTag = (-1);
//           }
//         }
//         
//         // If element blocks are used and we're done refine upper level
//         //-------------------------------------------------------------
//         if (iRefineLevel == refmap.GetMaxRefineLevel()-1) {
//           mapSegments.clear();
//           ConstructSegmentMap(vecFaces, mapSegments, iRefineLevel+1);
//           for (int i = 0; i < nVecFacesOriginalSize; i++) {
//             if (vecFaces[i].nRefineLevel == iRefineLevel+1) {
//               RefineFace(vecNodes, vecFaces, i);
//               vecFaces[i].nTag = (-1);
//             }
//           }
//         }
//         
//         // Remove refined blocks
//         //------------------------
//         DeleteTaggedFaces(vecFaces);
// 
//         /*
//         for (int i = 0; i < nVecFacesOriginalSize; i++) {
//           if (vecFaces[i].nTag == (-1)) {
//             vecFaces.erase(vecFaces.begin() + i);
//             i--;
//           }
//         }
//         */
//       }
//         /*
//       // Finalize (node vector can be changed here)
//       //---------------------------------------------
//       reftemp.Finalize(iRefineLevel, vecFaces, vecNodes);
// 
//       // Remove refined blocks
//       //-----------------------
//       DeleteTaggedFaces(vecFaces);
//         */
// 
//       // Reduce the active fine ratio
//       //------------------------------
//       nActiveFineRatio /= 2;
// 
//     } //  for (; iRefineLevel < refmap.GetMaxRefineLevel(); iRefineLevel++)
// 
//     // Finalize (node vector can be changed here)
//     //---------------------------------------------
//     reftemp.FinalizeAll(vecFaces, vecNodes);
// 
//     // End Function void RefineGrid()
//     //---------------------------------------
// }
// //=============================================================================
// 

