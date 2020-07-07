///////////////////////////////////////////////////////////////////////////////
///
///	\file    SpringDynamics.h
///	\author  Paul Ullrich
///	\version March 22, 2013
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

#ifndef _SPRINGDYNAMICS_H_
#define _SPRINGDYNAMICS_H_

#include "GridElements.h"


//=============================================================================
void SpringDynamics(NodeVector & vecNodes         ,
                    FaceVector & vecFaces         ,
                    int          nTransSmoothDist ,
                    int          nSmoothIterations);
//=============================================================================

//=============================================================================
void PressureDynamics(NodeVector & vecNodes         ,
                      FaceVector & vecFaces         ,
                      int          nTransSmoothDist ,
                      int          nSmoothIterations);
//=============================================================================

#endif
