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

#include "CSRefinementMap.h"

#include "Exception.h"
#include "MathHelper.h"
#include "CubedSphereTrans.h"

//*********************************
//#include "lodepng.h"
#ifdef __cplusplus
#include <vector>
//#include <string>
#endif /*__cplusplus*/
//*********************************



//=======================================================================
CSRefinementMap::CSRefinementMap( int nBaseResolution,
                                  int nMaxRefineLevel) : m_nBaseResolution(nBaseResolution),
                                                         m_nMaxRefineLevel(nMaxRefineLevel)
{
    int nMaxResolution = nBaseResolution*IntPow(2,nMaxRefineLevel-1);

    m_nMap.Initialize(6,nMaxResolution,nMaxResolution);
}
//=======================================================================


//=======================================================================
void CSRefinementMap::SetRefineLevel( int iPanel, int iA, int iB, int nRefineLevel) 
{
    if (nRefineLevel > m_nMaxRefineLevel) {
      _EXCEPTIONT("RefineLevel exceed maximum refinement level.");
    }

    m_nMap[iPanel][iA][iB] = nRefineLevel;
}
//=======================================================================


//=======================================================================
void CSRefinementMap::SetMinimumRefineLevel( int iPanel, 
                                             int iA, 
                                             int iB, 
                                             int iRefineLevel, 
                                             int nActiveFineRatio) 
{
    int iABegin = iA    *nActiveFineRatio;
    int iBBegin = iB    *nActiveFineRatio;
    int iAEnd   = (iA+1)*nActiveFineRatio;
    int iBEnd   = (iB+1)*nActiveFineRatio;

    for (int i = iABegin; i < iAEnd; i++) {
    for (int j = iBBegin; j < iBEnd; j++) {
      if (m_nMap[iPanel][i][j] < iRefineLevel) { m_nMap[iPanel][i][j] = iRefineLevel; }
    }
    }
}
//=======================================================================


/*-----------------------------------------------
//=======================================================================
void CSRefinementMap::InitializeFromPNG( const std::string & strPNGFile,
                                         double dLonBase,
                                         bool fInvert,
                                         bool fBlockRefine) 
{
    std::vector<unsigned char> image;
    unsigned int               nWidth;
    unsigned int               nHeight;

//      unsigned int iError = lodepng::decode(image, nWidth, nHeight, strPNGFile.c_str());

//    if (iError) {
//      _EXCEPTION2("\n  Decoder error %i: %s", iError, lodepng_error_text(iError));
//    }

    printf("----------------------------------------\n");
    printf("PNG loaded successfully\n");
    printf("... Dimensions: %i x %i\n", nWidth, nHeight);

    // Convert longitude shift to radians
    //------------------------------------
    dLonBase *= M_PI/180.0;

    printf("... Levels: ");
    int iLastLevel = (-1);
    for (int iC = 0; iC < 255; iC++) { 
      double dC            = static_cast<double>(iC) / 256.0;
      int    iCurrentLevel = floor(dC * (m_nMaxRefineLevel + 1));

      if (iCurrentLevel != iLastLevel) {
        printf("%i ", iC);
        iLastLevel = iCurrentLevel;
      }
    }
    printf("\n");

    // Loop throught all cubed-sphere faces
    //--------------------------------------
    for (int iP = 0; iP < 6; iP++) {
    for (int iA = 0; iA < m_nMap.GetColumns(); iA++) {
    for (int iB = 0; iB < m_nMap.GetColumns(); iB++) {

      int iMaxRefineLevel = 0;

      for (int iS = 0; iS < 2; iS++) {
      for (int iT = 0; iT < 2; iT++) {
        double dA = -0.25*M_PI + 0.5*M_PI*((static_cast<double>(iA+iS))/m_nMap.GetColumns());
        double dB = -0.25*M_PI + 0.5*M_PI*((static_cast<double>(iB+iT))/m_nMap.GetColumns());
        double dLon;
        double dLat;

        CubedSphereTrans::RLLFromABP(dA, dB, iP, dLon, dLat);

        dLon += dLonBase;
        dLon = dLon - 2.0*M_PI*floor(dLon/(2.0*M_PI));
        if ((dLon < 0.0) || (dLon > 2.0*M_PI)) { _EXCEPTIONT("Invalid longitude"); }

        int iLon = floor((dLon / (2.0 * M_PI)) * nWidth);
        int iLat = floor((dLat + 0.5 * M_PI) / M_PI * nHeight);

        if (iLon < 0       ) { iLon = 0;         }
        if (iLat < 0       ) { iLat = 0;         }
        if (iLon >= nWidth ) { iLon = nWidth -1; }
        if (iLat >= nHeight) { iLat = nHeight-1; }

        int iImageIx = ((nHeight-1-iLat)*nWidth + iLon)*4;

        // Calculate luminance
        //-----------------------
        double dY = + 0.2126 * static_cast<double>(image[iImageIx  ])/256.0
                    + 0.7152 * static_cast<double>(image[iImageIx+1])/256.0
                    + 0.0722 * static_cast<double>(image[iImageIx+2])/256.0;
        
        // Find active fine ratio
        //-----------------------------
        int iRefineLevel = floor(dY * static_cast<double>(m_nMaxRefineLevel + 1));
        
        // Invert the image
        //--------------------
        if (fInvert) { iRefineLevel = m_nMaxRefineLevel - iRefineLevel; }
        
        // Store the maximum refinement level among all nodes
        //----------------------------------------------------
        if (iRefineLevel > iMaxRefineLevel) {
          iMaxRefineLevel = iRefineLevel;
        }
      } // for (int iT = 0; iT < 2; iT++) 
      } // for (int iS = 0; iS < 2; iS++) 

      if (iMaxRefineLevel != 0) {
        int nActiveFineRatio;
        if (fBlockRefine) {
          nActiveFineRatio = IntPow(2, m_nMaxRefineLevel - 1);
        } else {
          nActiveFineRatio = IntPow(2, m_nMaxRefineLevel - iMaxRefineLevel);
        }

        SetMinimumRefineLevel(iP, iA/nActiveFineRatio, iB/nActiveFineRatio, 
                                      iMaxRefineLevel, nActiveFineRatio   );
      }
    } // for (int iB = 0; iB < m_nMap.GetColumns(); iB++) 
    } // for (int iA = 0; iA < m_nMap.GetColumns(); iA++) 
    } // for (int iP = 0; iP < 6; iP++) 

    // End Function CSRefinementMap::InitializeFromPNG
    //--------------------------------------------------
}
//=======================================================================
-----------------------------------------------*/


//=======================================================================
void CSRefinementMap::Normalize() 
{
    int iActiveLevel      = m_nMaxRefineLevel;
    int nActiveResolution = m_nMap.GetColumns();

    // Loop over all refinement levels
    //----------------------------------
    while (iActiveLevel > 0) {
      int nActiveFineRatio = m_nMap.GetColumns() / nActiveResolution;

      // Smooth the refinement map at this level
      //----------------------------------------
      bool fChanges = true;
      while (fChanges) {
        fChanges = false;

        for (int p = 0; p < 6; p++) {
        for (int i = 0; i < nActiveResolution; i++) {
        for (int j = 0; j < nActiveResolution; j++) {

          int ix = i*nActiveFineRatio;
          int jx = j*nActiveFineRatio;
  
          if (m_nMap[p][ix][jx] >= iActiveLevel) { continue; }
  
          // Find neighbors in the eight cardinal coordinate directions
          //------------------------------------------------------------
          int iE_src = ix + nActiveFineRatio;
          int iW_src = ix - nActiveFineRatio;
          int jN_src = jx + nActiveFineRatio;
          int jS_src = jx - nActiveFineRatio;
  
          int pE, pW, pN, pS;
          int iE, iW, iN, iS;
          int jE, jW, jN, jS;

          bool fNE, fNW, fSE, fSW;
          int pNE, pNW, pSE, pSW;
          int iNE, iNW, iSE, iSW;
          int jNE, jNW, jSE, jSW;

          CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),p,iE_src,jx    ,pE,iE,jE);
          CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),p,iW_src,jx    ,pW,iW,jW);
          CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),p,ix    ,jN_src,pN,iN,jN);
          CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),p,ix    ,jS_src,pS,iS,jS);

          fNE = CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),
                                                p,iE_src,jN_src,pNE,iNE,jNE);
          fNW = CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),
                                                p,iW_src,jN_src,pNW,iNW,jNW);
          fSE = CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),
                                                p,iE_src,jS_src,pSE,iSE,jSE);
          fSW = CubedSphereTrans::RelativeCoord(m_nMap.GetColumns(),
                                                p,iW_src,jS_src,pSW,iSW,jSW);

          // Check for cross-over corners (although these are not
          // necessarily problematic).  A better refinement algorithm
          // could negate the necessity of refining these points.
          //-------------------------------------------------------------
          if ((fNE) && (m_nMap[pE ][iE ][jE ] >= iActiveLevel) &&
                       (m_nMap[pN ][iN ][jN ] >= iActiveLevel) &&
                       (m_nMap[pNE][iNE][jNE] <  iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ( (fSE) && (m_nMap[pE ][iE ][jE ] >= iActiveLevel) &&
                        (m_nMap[pS ][iS ][jS ] >= iActiveLevel) &&
                        (m_nMap[pSE][iSE][jSE] <  iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          // Check for exterior corner intersections
          //------------------------------------------
          if ((fNE && (m_nMap[pNE][iNE][jNE] >= iActiveLevel))    &&
              (      ((m_nMap[pE ][iE ][jE ] <  iActiveLevel) && 
                      (m_nMap[pS ][iS ][jS ] >= iActiveLevel))  ||
                     ((m_nMap[pN ][iN ][jN ] <  iActiveLevel) && 
                      (m_nMap[pW ][iW ][jW ] >= iActiveLevel))   )  ){
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((fNW && (m_nMap[pNW][iNW][jNW] >= iActiveLevel))    &&
              (      ((m_nMap[pW ][iW ][jW ] <  iActiveLevel) &&
                      (m_nMap[pS ][iS ][jS ] >= iActiveLevel))  ||
                     ((m_nMap[pN ][iN ][jN ] <  iActiveLevel) &&
                      (m_nMap[pE ][iE ][jE ] >= iActiveLevel))   )  ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((fSE && (m_nMap[pSE][iSE][jSE] >= iActiveLevel))    &&
              (      ((m_nMap[pE ][iE ][jE ] <  iActiveLevel) &&
                      (m_nMap[pN ][iN ][jN ] >= iActiveLevel))  ||
                     ((m_nMap[pS ][iS ][jS ] <  iActiveLevel) &&
                      (m_nMap[pW ][iW ][jW ] >= iActiveLevel))   )  ){
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if((fSW && (m_nMap[pSW][iSW][jSW] >= iActiveLevel))    &&
             (      ((m_nMap[pW ][iW ][jW ] <  iActiveLevel) &&
                     (m_nMap[pN ][iN ][jN ] >= iActiveLevel))  ||
                    ((m_nMap[pS ][iS ][jS ] <  iActiveLevel) &&
                     (m_nMap[pE ][iE ][jE ] >= iActiveLevel))   )  ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          // Check for faces bounded on opposite sides
          //--------------------------------------------
          if ((m_nMap[pE][iE][jE] >= iActiveLevel) &&
              (m_nMap[pW][iW][jW] >= iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((m_nMap[pN][iN][jN] >= iActiveLevel) &&
              (m_nMap[pS][iS][jS] >= iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if (fNE && fNW && fSE && fSW) {
            if ((m_nMap[pNE][iNE][jNE] >= iActiveLevel) &&
                (m_nMap[pSW][iSW][jSW] >= iActiveLevel) &&
                (m_nMap[pSE][iSE][jSE] <  iActiveLevel) &&
                (m_nMap[pNW][iNW][jNW] <  iActiveLevel)   ) {
              SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
              fChanges = true;
              continue;
            }
            if ((m_nMap[pNE][iNE][jNE] <  iActiveLevel) &&
                (m_nMap[pSW][iSW][jSW] <  iActiveLevel) &&
                (m_nMap[pSE][iSE][jSE] >= iActiveLevel) &&
                (m_nMap[pNW][iNW][jNW] >= iActiveLevel)   ) {
              SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
              fChanges = true;
              continue;
            }
          }
        } // for (int j = 0; j < nActiveResolution; j++) 
        } // for (int i = 0; i < nActiveResolution; i++) 
        } // for (int p = 0; p < 6; p++) 
      } //  while (fChanges) 

      if (iActiveLevel == 1) { break; }

      // Loop over all elements at this refinement level
      //--------------------------------------------------
      for (int p = 0; p < 6; p++) {
      for (int i = 0; i < nActiveResolution; i++) {
      for (int j = 0; j < nActiveResolution; j++) {
        int iref = i*nActiveFineRatio;
        int jref = j*nActiveFineRatio;

        if (m_nMap[p][iref][jref] < iActiveLevel) { continue; }

        // Loop over all neighbors
        //----------------------------
        for (int ix = -1; ix <= 1; ix++) {
        for (int jx = -1; jx <= 1; jx++) {
          if ((ix == 0) && (jx == 0)) { continue; }

          int  i_src = (i+ix+2)/2-1;
          int  j_src = (j+jx+2)/2-1;
          int  p_dest;
          int  i_dest;
          int  j_dest;
          bool fNeighborExists = CubedSphereTrans::RelativeCoord(nActiveResolution/2,
                                                  p,i_src,j_src,p_dest,i_dest,j_dest);
          if (!fNeighborExists) { continue; }

          SetMinimumRefineLevel(p_dest,i_dest,j_dest,iActiveLevel-1,nActiveFineRatio*2);
        }
        }
      }
      }
      }

      iActiveLevel--;
      nActiveResolution /= 2;
    } //  while (iActiveLevel > 0) 

    // End Function CSRefinementMap::Normalize()
    //-----------------------------------------------------
}
//=======================================================================


//=======================================================================
void CSRefinementMap::ToFile(const char * szFile) const 
{

    FILE *fp = fopen(szFile, "w");

    for (int p = 0                    ; p <  6; p++) {
    for (int j = m_nMap.GetColumns()-1; j >= 0; j--) {
      for (unsigned int i = 0; i < m_nMap.GetColumns(); i++) { fprintf(fp, "%i ", m_nMap[p][i][j]); }
      fprintf(fp, "\n");
    }
    }
    fclose(fp);
}
//=======================================================================


//=======================================================================
void CSRefinementMap::FromFile(const char * szFile) 
{

    FILE *fp = fopen(szFile, "r");

    if (fp == NULL) { _EXCEPTION1("CSRefinement file \"%s\" not found", szFile); }
    
    for (int p = 0                    ; p <  6; p++) {
    for (int j = m_nMap.GetColumns()-1; j >= 0; j--) {
      for (unsigned int i = 0; i < m_nMap.GetColumns(); i++) { fscanf(fp, "%i ", &(m_nMap[p][i][j])); }
    }
    }
    fclose(fp);
    
    // Normalize the input refinement map
    //------------------------------------
    Normalize();
}
//=======================================================================
