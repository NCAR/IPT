///////////////////////////////////////////////////////////////////////////////
///
///	\file    CSRefinementMap.h
///	\author  Paul Ullrich
///	\version March 12, 2012
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

#ifndef _CSREFINEMENTMAP_H_
#define _CSREFINEMENTMAP_H_

#include "DataMatrix3D.h"

//=====================================================================================
class CSRefinementMap 
  //  Stores the cubed-sphere refinement map, which indicates the
  //  refinement level of all elements on the cubed-sphere down to the
  //  finest resolution.
{
    public: CSRefinementMap( int nBaseResolution,         ///  Constructor.
                             int nMaxRefineLevel);
                                                         
    public: void SetRefineLevel( int iPanel,              ///  Adjust the refinement level  
                                 int iA,                  ///  of a particular element.
                                 int iB,                   
                                 int nRefineLevel);        

    private: void SetMinimumRefineLevel( int iPanel,      ///  Set the minimum refinement 
                                         int iA,          ///  level of a mesh block.
                                         int iB,          
                                         int iRefineLevel,
                                         int nActiveFineRatio);

//    public: void InitializeFromPNG( const std::string & strPNGFile, ///  Initialize the refinement 
//                                    double dLonBase,                ///  map using a PNG image.
//                                    bool fInvert,
//                                    bool fBlockRefine);

    public: void Normalize();                              ///  Normalize the mesh.

    public: 
      int   GetBaseResolution()         const { return m_nBaseResolution; } /// Get the base resolution.
      int   GetMaxRefineLevel()         const { return m_nMaxRefineLevel; } /// Get the maximum refinement level.
      const DataMatrix3D<int> &GetMap() const { return m_nMap;            } /// Get the refinement map.
      void  ToFile  (const char *szFile) const;                             /// Output the refinement map to a file.
      void  FromFile(const char *szFile);                                   /// Input the refinement map from a file.

    private:
      int               m_nBaseResolution;    ///  Base resolution of the mesh.
      int               m_nMaxRefineLevel;    ///  Maximum refinement level of the mesh.
      DataMatrix3D<int> m_nMap;               ///  Refinement map on each cubed sphere panel.
};
//=====================================================================================


#endif

