///////////////////////////////////////////////////////////////////////////////
///
///	\file    IcosahedralFlagGrid.h
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

#ifndef _ICOSAHEDRALFLAGGRID_H_
#define _ICOSAHEDRALFLAGGRID_H_

#include "GridElements.h"

//=============================================================================
void GenerateFacesFromTriangle( int          nRefineLevel, 
                                const Edge & edge0       , 
                                const Edge & edge1       , 
                                const Edge & edge2       , 
                                NodeVector & vecNodes    , 
                                FaceVector & vecFaces    );
//=============================================================================


//=============================================================================
void ConvertFromLonLatToCartesian( const LonLatNodeVector & vecLonLatNodes,
                                         NodeVector       & vecNodes      );
//=============================================================================


//=============================================================================
void GenerateIcosahedralQuadGrid( int          nRefineLevel,
                                  NodeVector & vecNodes    ,
                                  FaceVector & vecFaces    );
//=============================================================================
#endif

