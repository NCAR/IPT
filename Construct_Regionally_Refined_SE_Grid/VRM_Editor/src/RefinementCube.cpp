#include "RefinementCube.h"

RefinementCube::RefinementCube()
{

}

RefinementCube::RefinementCube(int I_nBaseResolution, int I_nMaxRefineLevel) : nBaseResolution(I_nBaseResolution), nMaxRefineLevel(I_nMaxRefineLevel)
//
// Constructor:
//==============================================
{
    double dA;
    double dB;
    double dLon;
    double dLat;

    // Set cube rotation to 0.0 for now.
    //------------------------------------
    rLonShift = 0.0;
    rXrot     = 0.0;
    rYrot     = 0.0;

    // Initialize arrays for the cube
    //------------------------------
    if(nMaxRefineLevel > 0) {
        nMaxResolution = nBaseResolution*IntPow(2,nMaxRefineLevel-1);
    } else {
        nMaxResolution = nBaseResolution;
    }
    val.Initialize(6,nMaxResolution,nMaxResolution);
    CubeLon.Initialize(6,nMaxResolution+1,nMaxResolution+1);
    CubeLat.Initialize(6,nMaxResolution+1,nMaxResolution+1);

    // Set lat/lon gridpoints for cube
    //---------------------------------
    for( int iP=0; iP < 6             ; iP++) {
    for( int iA=0; iA < (nMaxResolution+1); iA++) {
    for( int iB=0; iB < (nMaxResolution+1); iB++) {
        dA = M_PI*(0.5*((static_cast<double>(iA)/nMaxResolution)) - 0.25);
        dB = M_PI*(0.5*((static_cast<double>(iB)/nMaxResolution)) - 0.25);
        CubedSphereTrans::RLLFromABP( dA, dB, iP, dLon, dLat);
        dLon *= 180./M_PI;
        dLat *= 180./M_PI;
        if(dLon < 0.) { dLon += 360.; }
        CubeLon[iP][iA][iB] = dLon;
        CubeLat[iP][iA][iB] = dLat;
    }
    }
    }
}

RefinementCube::~RefinementCube()
//
// Destructor:
//============================================
{
    nBaseResolution = 0;
    nMaxRefineLevel = 0;
    nMaxResolution  = 0;
    val.Deinitialize();
    CubeLon.Deinitialize();
    CubeLat.Deinitialize();
}

void RefinementCube::resize(int I_nBaseResolution, int I_nMaxRefineLevel)
//
// resize: Resize the cube grid for the given resolution and
//         reinitialize values.
//===============================================================
{
    double dA;
    double dB;
    double dLon;
    double dLat;

    // Get Rid of the existing map values
    //------------------------------------
    val.Deinitialize();
    CubeLon.Deinitialize();
    CubeLat.Deinitialize();

    // Initialize cube grid arrays
    //-----------------------------
    nBaseResolution = I_nBaseResolution;
    nMaxRefineLevel = I_nMaxRefineLevel;
    if(nMaxRefineLevel > 0) {
        nMaxResolution = nBaseResolution*IntPow(2,nMaxRefineLevel-1);
    } else {
        nMaxResolution = nBaseResolution;
    }
    val.Initialize(6,nMaxResolution,nMaxResolution);
    CubeLon.Initialize(6,nMaxResolution+1,nMaxResolution+1);
    CubeLat.Initialize(6,nMaxResolution+1,nMaxResolution+1);

    // Set lat/lon values for cube gridpoints
    //----------------------------------------
    for( int iP=0; iP < 6             ; iP++) {
    for( int iA=0; iA < (nMaxResolution+1); iA++) {
    for( int iB=0; iB < (nMaxResolution+1); iB++) {
        dA = M_PI*(0.5*((static_cast<double>(iA)/nMaxResolution)) - 0.25);
        dB = M_PI*(0.5*((static_cast<double>(iB)/nMaxResolution)) - 0.25);
        CubedSphereTrans::RLLFromABP( dA, dB, iP, dLon, dLat);
        dLon *= 180./M_PI;
        dLat *= 180./M_PI;
        if( dLon < 0.) { dLon += 360.;}
        CubeLon[iP][iA][iB] = dLon;
        CubeLat[iP][iA][iB] = dLat;
    }
    }
    }
}

void RefinementCube::resize(int I_nBaseResolution, int I_nMaxRefineLevel, double LonShift, double RotX, double RotY)
//
// resize: Resize the cube grid for the given resolution and
//         reinitialize values.
//===============================================================
{
    double dA;
    double dB;
    double dLon;
    double dLat;

    // Get Rid of the existing map values
    //------------------------------------
    val.Deinitialize();
    CubeLon.Deinitialize();
    CubeLat.Deinitialize();

    // Initialize cube grid arrays
    //-----------------------------
    nBaseResolution = I_nBaseResolution;
    nMaxRefineLevel = I_nMaxRefineLevel;
    if(nMaxRefineLevel > 0) {
        nMaxResolution = nBaseResolution*IntPow(2,nMaxRefineLevel-1);
    } else {
        nMaxResolution = nBaseResolution;
    }
    val.Initialize(6,nMaxResolution,nMaxResolution);
    CubeLon.Initialize(6,nMaxResolution+1,nMaxResolution+1);
    CubeLat.Initialize(6,nMaxResolution+1,nMaxResolution+1);

    // set the rotation angles
    //--------------------------
    rLonShift = LonShift * M_PI/180.;
    rXrot     = RotX * M_PI/180.;
    rYrot     = RotY * M_PI/180.;

    // Set lat/lon values for cube gridpoints
    //----------------------------------------
    for( int iP=0; iP < 6             ; iP++) {
    for( int iA=0; iA < (nMaxResolution+1); iA++) {
    for( int iB=0; iB < (nMaxResolution+1); iB++) {
        dA = M_PI*(0.5*((static_cast<double>(iA)/nMaxResolution)) - 0.25);
        dB = M_PI*(0.5*((static_cast<double>(iB)/nMaxResolution)) - 0.25);
        CubedSphereTrans::RLLFromABP( dA, dB, iP, dLon, dLat);


        // convert lat/lon to cartesian
        //------------------------------
        double Xval = cos(dLon)*cos(dLat);
        double Yval = sin(dLon)*cos(dLat);
        double Zval = sin(dLat);

        // apply X-axis rotation
        //-------------------------
        double Yrot = Yval*cos(rXrot) - Zval*sin(rXrot);
        double Zrot = Yval*sin(rXrot) + Zval*cos(rXrot);
        Yval = Yrot;
        Zval = Zrot;

        // apply Y-axis rotation
        //-----------------------
        double Xrot = Xval*cos(rYrot) - Zval*sin(rYrot);
               Zrot = Xval*sin(rYrot) + Zval*cos(rYrot);
        Xval = Xrot;
        Zval = Zrot;

        // convert back to lat/lon
        //-------------------------
        dLat = asin(Zval);
        dLon = atan2(Yval,Xval);

        // Apply LonShit
        //---------------
        dLon += rLonShift;

        dLon *= 180./M_PI;
        dLat *= 180./M_PI;
        if( dLon < 0.) { dLon += 360.;}
        CubeLon[iP][iA][iB] = dLon;
        CubeLat[iP][iA][iB] = dLat;
    }
    }
    }
}

void RefinementCube::loadCubeVals(RefinementMap I_refMap)
{
    double rval;
    double rval2, rval3, rval4;

    // Loop over cube gridpoints and get the corresponding refMap value
    //------------------------------------------------------------------
    if( nMaxRefineLevel > 0 ) {
        for( int iP=0; iP < 6             ; iP++) {
        for( int iA=0; iA < nMaxResolution; iA++) {
        for( int iB=0; iB < nMaxResolution; iB++) {
            rval  = I_refMap.getVal(CubeLon[iP][iA  ][iB  ],CubeLat[iP][iA  ][iB  ]);
            rval2 = I_refMap.getVal(CubeLon[iP][iA+1][iB  ],CubeLat[iP][iA+1][iB  ]);
            rval3 = I_refMap.getVal(CubeLon[iP][iA+1][iB+1],CubeLat[iP][iA+1][iB+1]);
            rval4 = I_refMap.getVal(CubeLon[iP][iA  ][iB+1],CubeLat[iP][iA  ][iB+1]);

            //***********************************************
            if(true) {
                rval = (rval + rval2 + rval3 + rval4)/4.0;
            } else {
                if(rval2 > rval) { rval = rval2;}
                if(rval3 > rval) { rval = rval3;}
                if(rval4 > rval) { rval = rval4;}
            }
            //***********************************************

// ???      val[iP][iA][iB] = floor(      rval*static_cast<double>(nMaxRefineLevel));
            val[iP][iA][iB] = floor(0.5 + rval*static_cast<double>(nMaxRefineLevel));
// ???      val[iP][iA][iB] = floor(0.01+ rval*static_cast<double>(nMaxRefineLevel));
// PFC            if(rval != 0.0 ) {
// PFC            std::cout << " Refcube["<< iP << "," << iA << "," << iB <<"] = " << val[iP][iA][iB] <<" " << rval << std::endl;
// PFC            std::cout << " Lon=" << CubeLon[iP][iA][iB] << " Lat=" << CubeLat[iP][iA][iB] << std::endl;
// PFC            }
        }
        }
        }
    } else {
        for( int iP=0; iP < 6             ; iP++) {
        for( int iA=0; iA < nMaxResolution; iA++) {
        for( int iB=0; iB < nMaxResolution; iB++) {
            val[iP][iA][iB] = 0;
        }
        }
        }
    }
}

void RefinementCube::get_CubeIndices(double lon, double lat, int &iP, int &iA, int &iB)
{
    // Given the lan,lon values, determine the associated cube indices
    //-----------------------------------------------------------------

    // ** If the cube has been rotated Z, then Y, then X the lat/lon
    //    values need to be mapped to lat,lon for un-rotated cube.
    //    Then those lat,lon values can be mapped to get the correct
    //    cube indices.
    //-------------------------------------------------------------
       // Transform lat,lon to cartesian x,y,z relative to rotated cube
       // Un-Rotate X, then Y, then Z to get the cartesian coordinates
       //  relative to unrotated cube.
       // Transform X,Y,Z to lat',lon'.

    double rLon = lon *M_PI/180.;
    double rLat = lat *M_PI/180.;

    // UN-apply LonShit
    //---------------
    rLon -= rLonShift;

    // convert lat/lon to cartesian
    //------------------------------
    double Xval = cos(rLon)*cos(rLat);
    double Yval = sin(rLon)*cos(rLat);
    double Zval = sin(rLat);

    // UN-apply Y-axis rotation
    //-----------------------
    double Xrot = Xval*cos(-rYrot) - Zval*sin(-rYrot);
    double Zrot = Xval*sin(-rYrot) + Zval*cos(-rYrot);
    Xval = Xrot;
    Zval = Zrot;

    // UN-apply X-axis rotation
    //-------------------------
    double Yrot = Yval*cos(-rXrot) - Zval*sin(-rXrot);
           Zrot = Yval*sin(-rXrot) + Zval*cos(-rXrot);
    Yval = Yrot;
    Zval = Zrot;

    // convert back to lat/lon
    //-------------------------
    rLat = asin(Zval);
    rLon = atan2(Yval,Xval);

    double dLon = rLon;
    double dLat = rLat;

    // map lat,lon to gnomic coordinates
    //------------------------------------
    double dA, dB;
 //   double dLon = lon*M_PI/180.;
 //   double dLat = lat*M_PI/180.;
    CubedSphereTrans::ABPFromRLL(dLon,dLat, dA, dB, iP);

    // Get the iA, iB indices for the refinement cube
    //------------------------------------------------
    iA = floor(2.*nMaxResolution*((dA/M_PI) + 0.25));
    iB = floor(2.*nMaxResolution*((dB/M_PI) + 0.25));
}

void RefinementCube::Normalize()
{
    int iActiveLevel      = nMaxRefineLevel;
    int nActiveResolution = val.GetColumns();

    // Loop over all refinement levels
    //----------------------------------
    while (iActiveLevel > 0) {
      int nActiveFineRatio = val.GetColumns() / nActiveResolution;

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

          if (val[p][ix][jx] >= iActiveLevel) { continue; }

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
          int  pNE, pNW, pSE, pSW;
          int  iNE, iNW, iSE, iSW;
          int  jNE, jNW, jSE, jSW;

          CubedSphereTrans::RelativeCoord(val.GetColumns(),p,iE_src,jx    ,pE,iE,jE);
          CubedSphereTrans::RelativeCoord(val.GetColumns(),p,iW_src,jx    ,pW,iW,jW);
          CubedSphereTrans::RelativeCoord(val.GetColumns(),p,ix    ,jN_src,pN,iN,jN);
          CubedSphereTrans::RelativeCoord(val.GetColumns(),p,ix    ,jS_src,pS,iS,jS);

          fNE = CubedSphereTrans::RelativeCoord(val.GetColumns(),
                                                p,iE_src,jN_src,pNE,iNE,jNE);
          fNW = CubedSphereTrans::RelativeCoord(val.GetColumns(),
                                                p,iW_src,jN_src,pNW,iNW,jNW);
          fSE = CubedSphereTrans::RelativeCoord(val.GetColumns(),
                                                p,iE_src,jS_src,pSE,iSE,jSE);
          fSW = CubedSphereTrans::RelativeCoord(val.GetColumns(),
                                                p,iW_src,jS_src,pSW,iSW,jSW);

          // Check for cross-over corners (although these are not
          // necessarily problematic).  A better refinement algorithm
          // could negate the necessity of refining these points.
          //-------------------------------------------------------------
          if ((fNE) && (val[pE ][iE ][jE ] >= iActiveLevel) &&
                       (val[pN ][iN ][jN ] >= iActiveLevel) &&
                       (val[pNE][iNE][jNE] <  iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ( (fSE) && (val[pE ][iE ][jE ] >= iActiveLevel) &&
                        (val[pS ][iS ][jS ] >= iActiveLevel) &&
                        (val[pSE][iSE][jSE] <  iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          // Check for exterior corner intersections
          //------------------------------------------
          if ((fNE && (val[pNE][iNE][jNE] >= iActiveLevel))    &&
              (      ((val[pE ][iE ][jE ] <  iActiveLevel) &&
                      (val[pS ][iS ][jS ] >= iActiveLevel))  ||
                     ((val[pN ][iN ][jN ] <  iActiveLevel) &&
                      (val[pW ][iW ][jW ] >= iActiveLevel))   )  ){
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((fNW && (val[pNW][iNW][jNW] >= iActiveLevel))    &&
              (      ((val[pW ][iW ][jW ] <  iActiveLevel) &&
                      (val[pS ][iS ][jS ] >= iActiveLevel))  ||
                     ((val[pN ][iN ][jN ] <  iActiveLevel) &&
                      (val[pE ][iE ][jE ] >= iActiveLevel))   )  ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((fSE && (val[pSE][iSE][jSE] >= iActiveLevel))    &&
              (      ((val[pE ][iE ][jE ] <  iActiveLevel) &&
                      (val[pN ][iN ][jN ] >= iActiveLevel))  ||
                     ((val[pS ][iS ][jS ] <  iActiveLevel) &&
                      (val[pW ][iW ][jW ] >= iActiveLevel))   )  ){
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if((fSW && (val[pSW][iSW][jSW] >= iActiveLevel))    &&
             (      ((val[pW ][iW ][jW ] <  iActiveLevel) &&
                     (val[pN ][iN ][jN ] >= iActiveLevel))  ||
                    ((val[pS ][iS ][jS ] <  iActiveLevel) &&
                     (val[pE ][iE ][jE ] >= iActiveLevel))   )  ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          // Check for faces bounded on opposite sides
          //--------------------------------------------
          if ((val[pE][iE][jE] >= iActiveLevel) &&
              (val[pW][iW][jW] >= iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if ((val[pN][iN][jN] >= iActiveLevel) &&
              (val[pS][iS][jS] >= iActiveLevel)   ) {
            SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
            fChanges = true;
            continue;
          }

          if (fNE && fNW && fSE && fSW) {
            if ((val[pNE][iNE][jNE] >= iActiveLevel) &&
                (val[pSW][iSW][jSW] >= iActiveLevel) &&
                (val[pSE][iSE][jSE] <  iActiveLevel) &&
                (val[pNW][iNW][jNW] <  iActiveLevel)   ) {
              SetMinimumRefineLevel(p,i,j,iActiveLevel,nActiveFineRatio);
              fChanges = true;
              continue;
            }
            if ((val[pNE][iNE][jNE] <  iActiveLevel) &&
                (val[pSW][iSW][jSW] <  iActiveLevel) &&
                (val[pSE][iSE][jSE] >= iActiveLevel) &&
                (val[pNW][iNW][jNW] >= iActiveLevel)   ) {
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

        if (val[p][iref][jref] < iActiveLevel) { continue; }

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

    // End Function RefinementCube::Normalize()
    //-----------------------------------------------------
}

void RefinementCube::read(const char *szFile)
{

}

void RefinementCube::write(const char *szFile) const
{
    FILE *fp = fopen(szFile,"w");

    for( int p=0; p<6; p++) {
    for( int j=val.GetColumns()-1; j >= 0; j--) {
        for( int i=0; i < val.GetColumns(); i++) {
            fprintf(fp, "%i ", val[p][i][j]);
        }
        fprintf(fp, "\n");
    }
    }
    fclose(fp);
}

void RefinementCube::SetRefineLevel(int iPanel,
                                    int iA,
                                    int iB,
                                    int iRefineLevel)
{
    if (iRefineLevel > nMaxRefineLevel) {
      _EXCEPTIONT("RefineLevel exceed maximum refinement level.");
    }

    val[iPanel][iA][iB] = iRefineLevel;
}

void RefinementCube::SetMinimumRefineLevel(int iPanel, 
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
      if (val[iPanel][i][j] < iRefineLevel) { val[iPanel][i][j] = iRefineLevel; }
    }
    }
}
