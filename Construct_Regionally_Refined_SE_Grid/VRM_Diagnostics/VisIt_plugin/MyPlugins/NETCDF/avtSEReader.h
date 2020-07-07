/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef AVT_SE_READER_H
#define AVT_SE_READER_H
#include <avtNETCDFReaderBase.h>
#include <avtVariableCache.h>
#include <vectortypes.h>
#include <map>

class NETCDFFileObject;
class avtFileFormatInterface;

// ****************************************************************************
//  Class: avtSEReader
//
//  Purpose:
//      Reads in NETCDF files as a plugin to VisIt.
//
//  Programmer: Brad Whitlock
//  Creation:   Wed Aug 10 15:21:14 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class avtSEReader : public avtNETCDFReaderBase
{
public: 

                   avtSEReader(const char *filename, 
                                        NETCDFFileObject *);
                   avtSEReader(const char *filename);
                  ~avtSEReader();

    vtkDataSet    *GetMesh(int timeState, const char *);
    vtkDataArray  *GetVar(int timeState, const char *);
    vtkDataArray  *GetVectorVar(int timeState, const char *);
    void           PopulateDatabaseMetaData(int timeState, avtDatabaseMetaData *);

    void           CreateGlobalAttributesString(int nGlobalAtts, std::string &gaString);

    void          *GetAuxiliaryData(const char *var , int timeState,
                                    const char *type, void *args, DestructorFunction &);


protected:
    bool ReturnSpatialDimensionIndices(const intVector &dims, int sDims[3], int &nSDims) const;
    void ReturnDimStartsAndCounts(int timeState, const intVector &dims, 
                                  intVector &dimStarts, intVector &dimCounts) const;
private:

    // DATA Structures For the known Grids
    // --------------------------------------
    struct SEgrid_t { 
      int    ncol;
      int    ncenters;
      int   *connect;
      float *lon;
      float *lat;
      bool   has3D;
      int    nlev;
      float *lev;
      bool   has2Dheight;
      float *Zs;
      bool   has3Dheight;
      float *Z;
      bool   has2Dtangent;
      float *gradx_Zs;
      float *grady_Zs;
      bool   has3Dtangent;
      float *gradx_Z;
      float *grady_Z;
    };
      
    struct SEelem_t { 
      int    np;
      int    nelem;
      float *lon;
      float *lat;
      bool   has3D;
      int    nlev;
      float *lev;
      bool   has2Dheight;
      float *Zs;
      bool   has3Dheight;
      float *Z;
      bool   has2Dtangent;
      float *gradx_Zs;
      float *grady_Zs;
      bool   has3Dtangent;
      float *gradx_Z;
      float *grady_Z;
    };

    // DATA Structures For the known Grids
    // --------------------------------------
    typedef enum {SE_GRID_DATA, SE_ELEM_DATA, UNKNOWN       } GridType_t;
    typedef enum {SPACE_2D , TIME_2D  , SPACE_3D , TIME_3D  } GridName_t;
    typedef enum {ZMAP_GRID, ZMAP_HGHT, ZMAP_FTAN, ZMAP_VTAN} GridZmap_t;
    typedef enum {VOLUME_DATA , SURFACE_DATA                } GridData_t;

    struct SEvariable_t {
      GridType_t  GridType;
      GridName_t  GridName;
      GridZmap_t  GridZmap;
      GridData_t  GridData;
      std::string VarName;
      std::string LongName;
      std::string Units;
      std::string VectorComponent;
      std::string VectorName;
      std::string VectorLongName;
      std::string MeshName;
      std::string MyName;
    };

    struct SEvector_t {
      GridType_t  GridType;
      GridName_t  GridName;
      GridZmap_t  GridZmap;
      GridData_t  GridData;
      std::string VectorName;
      std::string LongName;
      std::string Units;
      std::string XvarName;
      std::string YvarName;
      std::string ZvarName;
      std::string MeshName;
      std::string MyName;
    };
   
    typedef std::map<std::string, SEvariable_t*>           varList_t   ;
    typedef std::map<std::string, SEvector_t*  >           vectorList_t;
    typedef std::map<std::string, SEvector_t*  >::iterator vectorIter_t;

    // Global Data
    // ---------------
    bool         hasSEgrid;
    SEgrid_t     SEgrid;
    bool         hasSEelem;
    SEelem_t     SEelem;
    varList_t    varList;
    vectorList_t vecBaseList;
    vectorList_t vectorList;
    vectorIter_t vecIter;

    // DATA MEMBERS
    bool                   meshNamesCreated;
    StringIntVectorMap     meshNameToDimensionsSizes;
    StringIntVectorMap     meshNameToNCDimensions;
    StringIntVectorMap     varToDimensionsSizes;

    int                    procNum;
    int                    procCount;
};

#endif
