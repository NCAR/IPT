#include "RefinementMap.h"
//#include <iostream>

RefinementMap::RefinementMap()
{
    val.Initialize(360,180);
    RefLon.Initialize(360,1);
    RefLat.Initialize(180,1);
    Min = 0.;
    Max = 1.;
}

RefinementMap::RefinementMap(int I_nRefLon, int I_nRefLat) : nRefLon(I_nRefLon), nRefLat(I_nRefLat)
//
// Constructor:
//===================================
{
    // Initialize arrays with the given size
    //---------------------------------------
    val.Initialize(I_nRefLon,I_nRefLat);
    RefLon.Initialize(I_nRefLon,1);
    RefLat.Initialize(I_nRefLat,1);
    Min = 0.;
    Max = 1.;
    dltLon = 360./static_cast<double>(nRefLon);
    dltLat = 180./static_cast<double>(nRefLat);

    // Set lat/lon gridpoints for map
    //---------------------------------
    for( int nn=0; nn < nRefLon; nn++) {
        RefLon[nn][0] = dltLon*(static_cast<double>(nn) + 0.5);
    }

    for( int nn=0; nn< nRefLat; nn++) {
        RefLat[nn][0] = -90.0 + dltLat*(static_cast<double>(nn) + 0.5);
    }
}

RefinementMap::~RefinementMap()
//
// Destructor:
//===============================
{
    nRefLon = 0;
    nRefLat = 0;
    val.Deinitialize();
    RefLon.Deinitialize();
    RefLat.Deinitialize();
}

void RefinementMap::initialize(int I_nRefLon, int I_nRefLat)
{
    // Get rid of the existing map values
    //------------------------------------
    val.Deinitialize();
    RefLon.Deinitialize();
    RefLat.Deinitialize();

    // Initailize map arrays.
    //-----------------------
    nRefLon = I_nRefLon;
    nRefLat = I_nRefLat;
    val.Initialize(nRefLon,nRefLat);
    RefLon.Initialize(nRefLon,1);
    RefLat.Initialize(nRefLat,1);
    Min = 0.;
    Max = 1.;
    dltLon = 360./static_cast<double>(nRefLon);
    dltLat = 180./static_cast<double>(nRefLat);

    // Explicitly set all map values to 0.0
    //-----------------------------------
    for( int Lon=0; Lon < nRefLon; Lon++) {
    for( int Lat=0; Lat < nRefLat; Lat++) {
        val[Lon][Lat] = 0.0;
    }
    }

    // Set Lat/Lon values for map gridpoints
    //---------------------------------------
    for( int nn=0; nn < nRefLon; nn++) {
        RefLon[nn][0] = dltLon*(static_cast<double>(nn) + 0.5);
    }

    for( int nn=0; nn< nRefLat; nn++) {
        RefLat[nn][0] = -90.0 + dltLat*(static_cast<double>(nn) + 0.5);
    }
}


void RefinementMap::resize(int I_nRefLon, int I_nRefLat)
//
// resize: Resize the Refinement Map to the given dimension
//         and reinitialize all values
//==========================================================
{
    // Get rid of the existing map values
    //------------------------------------
    val.Deinitialize();
    RefLon.Deinitialize();
    RefLat.Deinitialize();

    // Initailize map arrays.
    //-----------------------
    nRefLon = I_nRefLon;
    nRefLat = I_nRefLat;
    val.Initialize(I_nRefLon,I_nRefLat);
    RefLon.Initialize(I_nRefLon,1);
    RefLat.Initialize(I_nRefLat,1);
    Min = 0.;
    Max = 1.;
    dltLon = 360./static_cast<double>(nRefLon);
    dltLat = 180./static_cast<double>(nRefLat);

    // Set Lat/Lon values for map gridpoints
    //---------------------------------------
    for( int nn=0; nn < nRefLon; nn++) {
        RefLon[nn][0] = dltLon*(static_cast<double>(nn) + 0.5);
//        std::cout << " RefLon: "<< nn << " = "<< RefLon[nn][0] << std::endl;
    }

    for( int nn=0; nn< nRefLat; nn++) {
        RefLat[nn][0] = -90.0 + dltLat*(static_cast<double>(nn) + 0.5);
//        std::cout << " RefLat: "<< nn << " = "<< RefLat[nn][0] << std::endl;
    }
}

double RefinementMap::getVal(double dLon, double dLat)
//
// getVal: Given the lat,lon values in the range [0,360], [-90,+90] calculate
//         the bi-linear interpolation value of val[][] from the neighboring
//         gridpoints.
//===============================================================================
{
    double nLon;
    double nLat;
    int    iLon0;
    int    iLat0;
    int    iLon1;
    int    iLat1;
    double fLon;
    double fLat;
    double val0;
    double val1;
    double rval;

    // Get the ~real~ gridpoint locations for the given point
    //--------------------------------------------------------
    if(dLon <   0.0) { dLon += 360.; }
    if(dLon > 360.0) { dLon -= 360.; }
    nLon = (     dLon /dltLon) -0.5;
    nLat = ((dLat+90.)/dltLat) -0.5;

    // Set the Longitudes the bound the given dLon value
    //---------------------------------------------------
    if ( nLon >= static_cast<double>(nRefLon-1)) {
        iLon0 = nRefLon-1;
        iLon1 = 0;
        fLon  = nLon - iLon0;
    } else if(nLon < 0.) {
        iLon0 = nRefLon-1;
        iLon1 = 0;
        fLon = 1.0 + nLon;
    } else {
        iLon0 = floor(nLon);
        iLon1 = iLon0 + 1;
        fLon  = nLon - iLon0;
    }

    // Calc bi-linear interpolation value
    //-------------------------------------
    if( nLat < 0.0) {
        iLat1 = 0;
        fLat = 2.0*(nLat + 0.5);
        val0 = 0.0;
        for( int nn=0; nn < nRefLon ; nn++) {
            val0 += val[nn][iLat1];
        }
        val0 /= nRefLon;
        val1 = val[iLon0][iLat1] + fLon*(val[iLon1][iLat1] - val[iLon0][iLat1]);
        rval = val0 + fLat*(val1-val0);
    } else if ( nLat > static_cast<double>(nRefLat-1)) {
        iLat0 = nRefLat-1;
        fLat = 2.0*(nLat - nRefLat + 1.0);
        val0 = val[iLon0][iLat0] + fLon*(val[iLon1][iLat0] - val[iLon0][iLat0]);
        val1 = 0.0;
        for( int nn=0; nn < nRefLon ; nn++) {
            val1 += val[nn][iLat0];
        }
        val1 /= nRefLon;
        rval = val0 + fLat*(val1-val0);
    } else {
        iLat0 = floor(nLat);
        iLat1 = iLat0 + 1;
        fLat  = nLat - iLat0;
        val0  = val[iLon0][iLat0] + fLon*(val[iLon1][iLat0] - val[iLon0][iLat0]);
        val1  = val[iLon0][iLat1] + fLon*(val[iLon1][iLat1] - val[iLon0][iLat1]);
        rval  = val0 + fLat*(val1-val0);
    }

    // Return the desired value
    //---------------------------
    return rval;
}

void RefinementMap::getInd(double dLon, double dLat, int &iLon, int &iLat)
{
    double nLon;
    double nLat;
    if(dLon <   0.0) { dLon += 360.; }
    if(dLon > 360.0) { dLon -= 360.; }
    nLon = (     dLon /dltLon) -0.5;
    nLat = ((dLat+90.)/dltLat) -0.5;
    if ( nLon >= static_cast<double>(nRefLon-1)) {
        iLon = nRefLon-1;
    } else if(nLon < 0.) {
        iLon = nRefLon-1;
    } else {
        iLon = floor(nLon);
    }
    if( nLat < 0.0) {
        iLat = 0;
    } else if ( nLat > static_cast<double>(nRefLat-1)) {
        iLat = nRefLat-1;
    } else {
        iLat = floor(nLat);
    }
}

void RefinementMap::read()
{

}

void RefinementMap::write() const
{

}
