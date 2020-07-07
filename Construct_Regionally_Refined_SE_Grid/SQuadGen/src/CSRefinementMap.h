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

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Stores the cubed-sphere refinement map, which indicates the
///		refinement level of all elements on the cubed-sphere down to the
///		finest resolution.
///	</summary>
class CSRefinementMap {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	CSRefinementMap(
		int nBaseResolution,
		int nMaxRefineLevel
	);

public:
	///	<summary>
	///		Adjust the refinement level of a particular element.
	///	</summary>
	void SetRefineLevel(
		int iPanel,
		int iA,
		int iB,
		int nRefineLevel
	);

private:
	///	<summary>
	///		Set the minimum refinement level of a mesh block.
	///	</summary>
	void SetMinimumRefineLevel(
		int iPanel,
		int iA,
		int iB,
		int iRefineLevel,
		int nActiveFineRatio
	);

public:
	///	<summary>
	///		Initialize the refinement map using a PNG image.
	///	</smmary>
	void InitializeFromPNG(
		const std::string & strPNGFile,
		double dLonBase,
		double dLatBase,
		bool fInvert,
		bool fBlockRefine
	);

public:
	///	<summary>
	///		Normalize the mesh.
	///	</summary>
	void Normalize();

public:
	///	<summary>
	///		Get the base resolution.
	///	</summary>
	int GetBaseResolution() const {
		return m_nBaseResolution;
	}

	///	<summary>
	///		Get the maximum refinement level.
	///	</summary>
	int GetMaxRefineLevel() const {
		return m_nMaxRefineLevel;
	}

	///	<summary>
	///		Get the refinement map.
	///	</summary>
	const DataMatrix3D<int> & GetMap() const {
		return m_nMap;
	}

	///	<summary>
	///		Output the refinement map to a file.
	///	</summary>
	void ToFile(const char * szFile) const;

	///	<summary>
	///		Input the refinement map from a file.
	///	</summary>
	void FromFile(const char * szFile);

private:
	///	<summary>
	///		Base resolution of the mesh.
	///	</summary>
	int m_nBaseResolution;

	///	<summary>
	///		Maximum refinement level of the mesh.
	///	</summary>
	int m_nMaxRefineLevel;

	///	<summary>
	///		Refinement map on each cubed sphere panel.
	///	</summary>
	DataMatrix3D<int> m_nMap;
};

///////////////////////////////////////////////////////////////////////////////

#endif

