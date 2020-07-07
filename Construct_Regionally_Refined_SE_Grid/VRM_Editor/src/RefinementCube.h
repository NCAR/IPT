#ifndef REFINEMENTCUBE_H
#define REFINEMENTCUBE_H

#include "SQuadGen/DataMatrix3D.h"
#include "SQuadGen/CubedSphereTrans.h"
#include "SQuadGen/MathHelper.h"
#include "RefinementMap.h"

class RefinementCube
{
public:
    RefinementCube();
    RefinementCube(int I_nBaseResolution, int I_nMaxRefineLevel);
    ~RefinementCube();

    void resize( int I_nBaseResolution, int I_nMaxRefineLevel);
    void resize( int I_nBaseResolution, int I_nMaxRefineLevel, double LonShift, double RotX, double RotY);

    void loadCubeVals( RefinementMap I_refMap);
    void get_CubeIndices( double lon, double lat, int &iP, int &iA, int &iB);

    /// Normqalize the mesh.
    void Normalize();

    /// I/O the refinement map for the cube.
    void read (const char *szFile);
    void write(const char *szFile) const;

    int   GetBaseResolution()         const { return nBaseResolution; }
    int   GetMaxRefineLevel()         const { return nMaxRefineLevel; }
    int   GetMaxResolution()          const { return nMaxResolution;  }
    const DataMatrix3D<int>    &GetMap()     const { return val;      }
    const DataMatrix3D<double> &GetCubeLon() const { return CubeLon;  }
    const DataMatrix3D<double> &GetCubeLat() const { return CubeLat;  }

    ///  Adjust the refinement level of a particular element.
    void SetRefineLevel( int iPanel,
                         int iA,
                         int iB,           
                         int nRefineLevel);

private:
    ///  Set the minimum refinement level of a mesh block.
    void SetMinimumRefineLevel(int iPanel,
                               int iA,
                               int iB,
                               int iRefineLevel,
                               int nActiveFineRatio);

private:
    int    nBaseResolution;
    int    nMaxRefineLevel;
    int    nMaxResolution;
    double rLonShift;
    double rXrot;
    double rYrot;

public:
    DataMatrix3D<int>    val;
    DataMatrix3D<double> CubeLon;
    DataMatrix3D<double> CubeLat;
};

#endif // REFINEMENTCUBE_H
