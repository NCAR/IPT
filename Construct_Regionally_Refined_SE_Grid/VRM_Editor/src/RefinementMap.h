#ifndef REFINEMENTMAP_H
#define REFINEMENTMAP_H

#include "SQuadGen/DataMatrix.h"
#include <math.h>

class RefinementMap
{
public:
    RefinementMap();
    RefinementMap(int I_nRefLon, int I_nRefLat);

    ~RefinementMap();

    void initialize(int I_nRefLon, int I_nRefLat);
    void resize( int I_nRefLon, int I_nRefLat);

    double getVal( double dLon, double dLat);
    void   getInd(double dLon, double dLat, int &iLon, int &iLat);
    void   read();
    void   write() const;

public:
    DataMatrix<double> val;
    int    nRefLon;
    int    nRefLat;
    double Min;
    double Max;
    double dltLon;
    double dltLat;
    DataMatrix<double> RefLon;
    DataMatrix<double> RefLat;
};

#endif // REFINEMENTMAP_H
