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

#include <vector>
#include <snprintf.h>
#include <netcdf.h>

#include <avtSEReader.h>
#include <NETCDFFileObject.h>
#include <avtMaterial.h>
#include <avtDatabase.h>
#include <avtDatabaseMetaData.h>
#include <avtSTSDFileFormatInterface.h>
#include <DebugStream.h>

#include <vtkCellArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>

#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkVisItUtility.h>

#include <InvalidVariableException.h>
#include <InstallationFunctions.h>

#define TIME_DIMENSION      -1

#ifdef PARALLEL
#include <avtParallel.h>
#endif

// ****************************************************************************
// Method: avtSEReader::avtSEReader
//
// Purpose: 
//   Constructor for the avtSEReader class.
//
// Arguments:
//   filename : The name of the file being read.
//   f        : The file object associated with the file being read.
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 18 18:03:59 PST 2005
//
// Modifications:
//   Brad Whitlock, Wed Apr 26 17:40:20 PST 2006
//   Initialized meshNamesCreated.
//
//   Mark C. Miller, Tue Aug 15 15:28:11 PDT 2006
//   Added procNum, procCount to support on-the-fly parallel decomposition
//
// ****************************************************************************

avtSEReader::avtSEReader(const char *filename) :
    avtNETCDFReaderBase(filename), meshNameToDimensionsSizes(), meshNameToNCDimensions(), varToDimensionsSizes()
{
    meshNamesCreated = false;
    procNum = 0;
    procCount = 1;
#ifdef PARALLEL
    procNum = PAR_Rank();
    procCount = PAR_Size();
#endif
}

avtSEReader::avtSEReader(const char *filename, NETCDFFileObject *f) :
    avtNETCDFReaderBase(filename, f), meshNameToDimensionsSizes(), meshNameToNCDimensions(), varToDimensionsSizes()
{
    meshNamesCreated = false;
    procNum = 0;
    procCount = 1;
#ifdef PARALLEL
    procNum = PAR_Rank();
    procCount = PAR_Size();
#endif
}

// ****************************************************************************
// Method: avtSEReader::~avtSEReader
//
// Purpose: 
//   Destructor for the avtSEReader class.
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 18 18:04:36 PST 2005
//
// Modifications:
//   
// ****************************************************************************

avtSEReader::~avtSEReader()
{
}

// ****************************************************************************
// Method: avtSEReader::CreateGlobalAttributesString
//
// Purpose: 
//   Create a string for the global attributes.
//
// Arguments:
//   nGlobalAtts : The number of global attributes
//   gaString    : global attributes string.
//
// Returns:    
//
// Note:       
//
// Programmer: Brad Whitlock
// Creation:   Fri Oct 30 11:01:08 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

void
avtSEReader::CreateGlobalAttributesString(int nGlobalAtts, std::string &gaString)
{
    for(int i = 0; i < nGlobalAtts; ++i)
    {
        int     status;
        char    attname[NC_MAX_NAME+1];
        nc_type atttype;
        size_t  attsize;
        if((status = nc_inq_attname(fileObject->GetFileHandle(), NC_GLOBAL, i, attname))
            == NC_NOERR)
        {
            if((status = nc_inq_att(fileObject->GetFileHandle(), NC_GLOBAL, attname, &atttype,
                                    &attsize)) == NC_NOERR)
            {
                std::string tmpStr("\t");
                tmpStr += attname;

                if(atttype == NC_CHAR)
                {
                    char *value = new char[attsize+1];
                    nc_get_att_text(fileObject->GetFileHandle(), NC_GLOBAL, attname, value);
                    value[attsize] = '\0';
                    char *c2 = value + attsize - 1;
                    while(c2 >= value && *c2 == ' ')
                        *c2-- = '\0';
                    tmpStr += " = \"";
                    tmpStr += value;
                    tmpStr += "\"";

                    delete [] value;
                }

                tmpStr += "\n";
                gaString += tmpStr;
            }
        }
    }
}

// ****************************************************************************
// Method: avtSEReader::PopulateDatabaseMetaData
//
// Purpose: 
//   Populates the metadata from information in the file.
//
// Arguments:
//   md : The metadata object to populate.
//
// Programmer: Patrick Callaghan
// Creation:   Sept 23 2016
//
// Modifications:
//
// ****************************************************************************

void
avtSEReader::PopulateDatabaseMetaData(int timeState, avtDatabaseMetaData *md)
{
    debug4 << "avtSEReader::PopulateDatabaseMetaData" << endl;
    if(DebugStream::Level4())
        fileObject->PrintFileContents(DebugStream::Stream4());

    // Get the Number of dimensions and variables
    // --------------------------------------------
    int status, nDims, nVars, nGlobalAtts, unlimitedDimension;
    status = nc_inq(fileObject->GetFileHandle(), &nDims, &nVars, 
                               &nGlobalAtts, &unlimitedDimension);
    debug4 << "  nDims="              << nDims
           << ", nVars="              << nVars
           << ", nGlobalAtts="        << nGlobalAtts
           << ", unlimitedDimension=" << unlimitedDimension
           << endl;
    if(status != NC_NOERR)
    {
        fileObject->HandleError(status);
        return;
    }

    // Process Global Attributes
    // ---------------------------
    if(md != 0)
    {
        std::string gaString;
        CreateGlobalAttributesString(nGlobalAtts, gaString);
        md->SetDatabaseComment(gaString);
    }

    // =============================================
    // Check for variables defining the SEgrid mesh
    // =============================================
    hasSEgrid = false;
    SEgrid.ncol         = -1;
    SEgrid.ncenters     = -1;
    SEgrid.nlev         = -1;
    SEgrid.has3D        = false;
    SEgrid.has2Dheight  = false;
    SEgrid.has3Dheight  = false;
    SEgrid.has2Dtangent = false;
    SEgrid.has3Dtangent = false;

    if(!hasSEgrid) 
    {
      size_t  ncol;
      size_t  ncenters;
      size_t  nlev;
      int     Varid_lat;
      int     Varid_lon;
      int     Varid_con;
      int     Varid_lev;
      int     Varid;
      bool hasNcol     = fileObject->GetDimensionInfo("ncol"    , &ncol    );
      bool hasNcenters = fileObject->GetDimensionInfo("ncenters", &ncenters);
      bool hasNlev     = fileObject->GetDimensionInfo("lev"     , &nlev    );
      bool hasLat  = fileObject->GetVarId        ("lat"         , &Varid_lat);
      bool hasLon  = fileObject->GetVarId        ("lon"         , &Varid_lon);
      bool hasCon  = fileObject->GetVarId        ("connectivity", &Varid_con);
      bool hasLev  = fileObject->GetVarId        ("lev"         , &Varid_lev);
      if( hasNcol && hasLat && hasLon) 
      {
        // First Set the basic mesh values
        // ---------------------------------
        if(hasNcenters && hasCon) 
        {
          hasSEgrid = true;
          SEgrid.ncol     = ncol;
          SEgrid.ncenters = ncenters;
          SEgrid.lat      = new float[ncol];
          SEgrid.lon      = new float[ncol];
          SEgrid.connect  = new int[4*ncenters];
          fileObject->ReadVariableInto("lat",FLOATARRAY_TYPE,SEgrid.lat);
          fileObject->ReadVariableInto("lon",FLOATARRAY_TYPE,SEgrid.lon);
          fileObject->ReadVariableInto("connectivity",INTEGERARRAY_TYPE,
                                                         SEgrid.connect);
        } else {
          // TODO: Add the ability to read a  connectivity file
          //       here for 'old' history files.
          //       For now we try to find a SEMapping file in 
          //       the home .visit directory with the currnt value of ncol
          // --------------------------------------------------------------
#ifdef WIN32
          char SlashChar = '\\';
#else
          char SlashChar = '/';
#endif
          char tmp[100];
          int status;
          std::string mapFileName;
          std::string homedir = GetUserVisItDirectory();
          if(!homedir.empty())
          {
            if(homedir[homedir.size() - 1] != SlashChar) homedir += VISIT_SLASH_STRING;
            SNPRINTF(tmp, 100, "SEMapping_%08d.nc", ncol);
            mapFileName = homedir + std::string(tmp);
            int mapid;
            int ncentersID;
            if((status = nc_open(mapFileName.c_str(), NC_NOWRITE, &mapid)) == NC_NOERR)
            {
              debug4 << "avtSEReader::PopulateDatabaseMetaData" << mapFileName.c_str() << " was opened." << endl;
              hasNcenters =                (nc_inq_dimid (mapid,"ncenters"    , &ncentersID) == NC_NOERR);
              hasNcenters = hasNcenters && (nc_inq_dimlen(mapid, ncentersID   , &ncenters  ) == NC_NOERR);
//              hasCon      =                (nc_inq_varid (mapid,"connectivity", &Varid_con ) == NC_NOERR);
              hasCon      =                (nc_inq_varid (mapid,"element_corners", &Varid_con ) == NC_NOERR);
              if( hasNcenters && hasCon )
              {
                SEgrid.ncenters = ncenters;
                SEgrid.connect  = new int[4*ncenters];
                if(nc_get_var_int(mapid, Varid_con, SEgrid.connect) == NC_NOERR)
                {
                  SEgrid.ncol     = ncol;
                  SEgrid.lat      = new float[ncol];
                  SEgrid.lon      = new float[ncol];
                  fileObject->ReadVariableInto("lat",FLOATARRAY_TYPE,SEgrid.lat);
                  fileObject->ReadVariableInto("lon",FLOATARRAY_TYPE,SEgrid.lon);
                  hasSEgrid = true;
                } else {
                  debug4 << "avtSEReader::PopulateDatabaseMetaData Could not read connectivity data" << endl;
                  hasSEgrid = false;
                }
              } else {
                debug4 << "avtSEReader::PopulateDatabaseMetaData: File does not have ncenters/connectivity data" << endl;
                hasSEgrid = false;
              }
            } else {
              debug4 << "avtSEReader::PopulateDatabaseMetaData: Could not open " << mapFileName.c_str() << endl;
              hasSEgrid = false;
            }
          } else {
            hasSEgrid = false;
          }
        }

        // Now determine if optional mesh information is available.
        // ---------------------------------------------------------
        if(hasSEgrid)
        {
          if(hasNlev && hasLev) 
          {
            SEgrid.has3D = true;
            SEgrid.nlev  = nlev;
            SEgrid.lev   = new float[nlev];
            fileObject->ReadVariableInto("lev",FLOATARRAY_TYPE,SEgrid.lev);
          } else {
            SEgrid.has3D = false;
          }
          if(fileObject->GetVarId("PHIS", &Varid))
          {
            SEgrid.has2Dheight = true;
            SEgrid.Zs          = new float[ncol];
          } else {
            SEgrid.has2Dheight = false;
          }
          if(fileObject->GetVarId("Z3", &Varid))
          {
            SEgrid.has3Dheight = true;
            SEgrid.Z           = new float[ncol];
          } else {
            SEgrid.has3Dheight = false;
          }
          if((fileObject->GetVarId("gradx_PHIS", &Varid)) ||
             (fileObject->GetVarId("grady_PHIS", &Varid))   )
          {
            SEgrid.has2Dtangent = true;
            SEgrid.gradx_Zs     = new float[ncol];
            SEgrid.grady_Zs     = new float[ncol];
          } else {
            SEgrid.has2Dtangent = false;
          }
          if((fileObject->GetVarId("gradx_Z3", &Varid)) ||
             (fileObject->GetVarId("grady_Z3", &Varid))   )
          {
            SEgrid.has3Dtangent = true;
            SEgrid.gradx_Z      = new float[ncol];
            SEgrid.grady_Z      = new float[ncol];
          } else {
            SEgrid.has3Dtangent = false;
          }
        }
      } else {
        hasSEgrid = false;
      }
    }
    debug4 << "PFC:  hasSEgrid =" << hasSEgrid
           << ", has3D= "         << SEgrid.has3D
           << ", has2Dheight= "   << SEgrid.has2Dheight 
           << ", has3Dheight= "   << SEgrid.has3Dheight
           << ", has2Dtangent= "  << SEgrid.has2Dtangent
           << ", has3Dtangent= "  << SEgrid.has3Dtangent << endl;

    // ---------------------------------------------------------
    // Now create all of the display meshes supported by SEgrid
    // ---------------------------------------------------------
    if(md != 0 && hasSEgrid) 
    {
      int SpacDims = 3;
      int TopoDims = 3;

      // ===================
      // 2D Surfaces
      // ===================
      avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_2D_sfc_grid", 1, 1, 1, 0, 
                                                 SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
      mmd->xLabel = "Longitude";
      fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
      mmd->yLabel = "Latitude";
      fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
      mmd->zLabel = "Surface Level";
      mmd->zUnits = " ";
      md->Add(mmd);

      // ===================
      // 2D Height Surfaces
      // ===================
      if(SEgrid.has2Dheight)
      {
        avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_2D_sfc_hght", 1, 1, 1, 0, 
                                                   SpacDims, TopoDims, AVT_UNSTRUCTURED_MESH);
        mmd->xLabel = "Longitude";
        fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
        mmd->yLabel = "Latitude";
        fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
        mmd->zLabel = "Height";
        mmd->zUnits = "km";
        md->Add(mmd);
      }

      if(SEgrid.has3D)
      {
        // ===========
        // 3D Surfaces
        // ===========
        {
          avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_3D_sfc_grid", 1, 1, 1, 0, 
                                                     SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
          mmd->xLabel = "Longitude";
          fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
          mmd->yLabel = "Latitude";
          fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
          mmd->zLabel = "Model Level";
          mmd->zUnits = " ";
          md->Add(mmd);

          // Add material MetaData for level subsetting
          // -------------------------------------------
          avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
          matmd_lay->name     = "SEgrid_3D_sfc_grid_Levels";
          matmd_lay->meshName = "SEgrid_3D_sfc_grid";
          matmd_lay->numMaterials=SEgrid.nlev;
          for (int i=1; i <= SEgrid.nlev; ++i)
          {
             char buffer[50];
             int n;
             n=sprintf(buffer, "Level %d",(i)); (void) n;   // TODO: Add lev(i) values
             matmd_lay->materialNames.push_back(buffer);
          }
          md->Add(matmd_lay);
        }

        // ===========
        // 3D Volume
        // ===========
        {
          avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_3D_vol_grid", 1, 1, 1, 0, 
                                                     SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
          mmd->xLabel = "Longitude";
          fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
          mmd->yLabel = "Latitude";
          fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
          mmd->zLabel = "Model Level";
          mmd->zUnits = " ";
          md->Add(mmd);

          // Add material MetaData for level subsetting
          // -------------------------------------------
          avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
          matmd_lay->name     = "SEgrid_3D_vol_grid_Layers";
          matmd_lay->meshName = "SEgrid_3D_vol_grid";
          matmd_lay->numMaterials=(SEgrid.nlev-1);
          for (int i=1; i <= (SEgrid.nlev-1); ++i)
          {
             char buffer[50];
             int n;
             n=sprintf(buffer, "Layer %d",(i)); (void) n;   // TODO: Add lev(i) values
             matmd_lay->materialNames.push_back(buffer);
          }
          md->Add(matmd_lay);
        }

        if(SEgrid.has3Dheight)
        {
          // ===================
          // 3D Height Surfaces
          // ===================
          {
            avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_3D_sfc_hght", 1, 1, 1, 0, 
                                                       SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
            mmd->xLabel = "Longitude";
            fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
            mmd->yLabel = "Latitude";
            fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
            mmd->zLabel = "Height";
            mmd->zUnits = "km";
            md->Add(mmd);
  
            // Add material MetaData for level subsetting
            // -------------------------------------------
            avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
            matmd_lay->name     = "SEgrid_3D_sfc_hght_Levels";
            matmd_lay->meshName = "SEgrid_3D_sfc_hght";
            matmd_lay->numMaterials=SEgrid.nlev;
            for (int i=1; i <= SEgrid.nlev; ++i)
            {
               char buffer[50];
               int n;
               n=sprintf(buffer, "Level %d",(i)); (void) n;   // TODO: Add lev(i) values
               matmd_lay->materialNames.push_back(buffer);
            }
            md->Add(matmd_lay);
          }

          // =================
          // 3D Height Volume
          // =================
          {
            avtMeshMetaData *mmd = new avtMeshMetaData("SEgrid_3D_vol_hght", 1, 1, 1, 0, 
                                                       SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
            mmd->xLabel = "Longitude";
            fileObject->ReadStringAttribute("lon", "units", mmd->xUnits);
            mmd->yLabel = "Latitude";
            fileObject->ReadStringAttribute("lat", "units", mmd->yUnits);
            mmd->zLabel = "Height";
            mmd->zUnits = "km";
            md->Add(mmd);
  
            // Add material MetaData for level subsetting
            // -------------------------------------------
            avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
            matmd_lay->name     = "SEgrid_3D_vol_hght_Layers";
            matmd_lay->meshName = "SEgrid_3D_vol_hght";
            matmd_lay->numMaterials=(SEgrid.nlev-1);
            for (int i=1; i <= (SEgrid.nlev-1); ++i)
            {
               char buffer[50];
               int n;
               n=sprintf(buffer, "Layer %d",(i)); (void) n;   // TODO: Add lev(i) values
               matmd_lay->materialNames.push_back(buffer);
            }
            md->Add(matmd_lay);
          }
        }
      }
    }
  debug4 << "PFC:  Done with SEgrid Meshes " << endl;

    // =============================================
    // Check for variables defining the SEelem mesh
    // =============================================
  debug4 << "PFC:  Checking SEelem Grid " << endl;
    hasSEelem = false;
    SEelem.np           = -1;
    SEelem.nelem        = -1;
    SEelem.nlev         = -1;
    SEelem.has3D        = false;
    SEelem.has2Dheight  = false;
    SEelem.has3Dheight  = false;
    SEelem.has2Dtangent = false;
    SEelem.has3Dtangent = false;

    if(!hasSEelem) 
    {
      size_t  np;
      size_t  nelem;
      size_t  nlev;
      int     Varid_lat;
      int     Varid_lon;
      int     Varid_lev;
      int     Varid;
      bool hasNp    = fileObject->GetDimensionInfo("np"   , &np   );
      bool hasNelem = fileObject->GetDimensionInfo("nelem", &nelem);
      bool hasNlev  = fileObject->GetDimensionInfo("lev"  , &nlev );
      bool hasLat  = fileObject->GetVarId         ("E_lat", &Varid_lat);
      bool hasLon  = fileObject->GetVarId         ("E_lon", &Varid_lon);
      bool hasLev  = fileObject->GetVarId         ("lev"  , &Varid_lev);
      if( hasNp && hasNelem && hasLat && hasLon) 
      {
        // First Set the basic mesh values
        // ---------------------------------
        hasSEelem = true;
        SEelem.np    = np;
        SEelem.nelem = nelem;
        SEelem.lat   = new float[np*np*nelem];
        SEelem.lon   = new float[np*np*nelem];
        fileObject->ReadVariableInto("E_lat",FLOATARRAY_TYPE,SEelem.lat);
        fileObject->ReadVariableInto("E_lon",FLOATARRAY_TYPE,SEelem.lon);

        // Now determine if optional mesh information is available.
        // ---------------------------------------------------------
        if(hasSEelem)
        {
          if(hasNlev && hasLev) 
          {
            SEelem.has3D = true;
            SEelem.nlev  = nlev;
            SEelem.lev   = new float[nlev];
            fileObject->ReadVariableInto("lev",FLOATARRAY_TYPE,SEelem.lev);
          } else {
            SEelem.has3D = false;
          }
          if(fileObject->GetVarId("Zs", &Varid))
          {
            SEelem.has2Dheight = true;
            SEelem.Zs          = new float[np*np*nelem];
          } else {
            SEelem.has2Dheight = false;
          }
          if(fileObject->GetVarId("Z", &Varid))
          {
            SEelem.has3Dheight = true;
            SEelem.Z           = new float[np*np*nelem];
          } else {
            SEelem.has3Dheight = false;
          }
          if((fileObject->GetVarId("gradx_Zs", &Varid)) ||
             (fileObject->GetVarId("grady_Zs", &Varid))   )
          {
            SEelem.has2Dtangent = true;
            SEelem.gradx_Zs     = new float[np*np*nelem];
            SEelem.grady_Zs     = new float[np*np*nelem];
          } else {
            SEelem.has2Dtangent = false;
          }
          if((fileObject->GetVarId("gradx_Z", &Varid)) ||
             (fileObject->GetVarId("grady_Z", &Varid))   )
          {
            SEelem.has3Dtangent = true;
            SEelem.gradx_Z      = new float[np*np*nelem];
            SEelem.grady_Z      = new float[np*np*nelem];
          } else {
            SEelem.has3Dtangent = false;
          }
        }
      } else {
        hasSEelem = false;
      }
    }

    debug4 << "PFC:  hasSEelem =" << hasSEelem
           << ", has3D= "         << SEelem.has3D
           << ", has2Dheight= "   << SEelem.has2Dheight 
           << ", has3Dheight= "   << SEelem.has3Dheight
           << ", has2Dtangent= "  << SEelem.has2Dtangent
           << ", has3Dtangent= "  << SEelem.has3Dtangent << endl;

    // ---------------------------------------------------------
    // Now create all of the display meshes supported by SEelem
    // ---------------------------------------------------------
    if(md != 0 && hasSEelem) 
    {
      int SpacDims = 3;
      int TopoDims = 3;

      // ===================
      // 2D Surfaces
      // ===================
      avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_2D_sfc_grid", 1, 1, 1, 0, 
                                                 SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
      mmd->xLabel = "Longitude";
      fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
      mmd->yLabel = "Latitude";
      fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
      mmd->zLabel = "Surface Level";
      mmd->zUnits = " ";
      md->Add(mmd);

      // ===================
      // 2D Height Surfaces
      // ===================
      if(SEelem.has2Dheight)
      {
        avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_2D_sfc_hght", 1, 1, 1, 0, 
                                                   SpacDims, TopoDims, AVT_UNSTRUCTURED_MESH);
        mmd->xLabel = "Longitude";
        fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
        mmd->yLabel = "Latitude";
        fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
        mmd->zLabel = "Height";
        mmd->zUnits = "km";
        md->Add(mmd);
      }

      if(SEelem.has3D)
      {
        // ===========
        // 3D Surfaces
        // ===========
        {
          avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_3D_sfc_grid", 1, 1, 1, 0, 
                                                     SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
          mmd->xLabel = "Longitude";
          fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
          mmd->yLabel = "Latitude";
          fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
          mmd->zLabel = "Model Level";
          mmd->zUnits = " ";
          md->Add(mmd);

          // Add material MetaData for level subsetting
          // -------------------------------------------
          avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
          matmd_lay->name     = "SEelem_3D_sfc_grid_Levels";
          matmd_lay->meshName = "SEelem_3D_sfc_grid";
          matmd_lay->numMaterials=SEelem.nlev;
          for (int i=1; i <= SEelem.nlev; ++i)
          {
             char buffer[50];
             int n;
             n=sprintf(buffer, "Level %d",(i)); (void) n;   // TODO: Add lev(i) values
             matmd_lay->materialNames.push_back(buffer);
          }
          md->Add(matmd_lay);
        }

        // ===========
        // 3D Volume
        // ===========
        {
          avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_3D_vol_grid", 1, 1, 1, 0, 
                                                     SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
          mmd->xLabel = "Longitude";
          fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
          mmd->yLabel = "Latitude";
          fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
          mmd->zLabel = "Model Level";
          mmd->zUnits = " ";
          md->Add(mmd);
  
          // Add material MetaData for level subsetting
          // -------------------------------------------
          avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
          matmd_lay->name     = "SEelem_3D_vol_grid_Layers";
          matmd_lay->meshName = "SEelem_3D_vol_grid";
          matmd_lay->numMaterials=(SEelem.nlev-1);
          for (int i=1; i <= (SEelem.nlev-1); ++i)
          {
             char buffer[50];
             int n;
             n=sprintf(buffer, "Layer %d",(i)); (void) n;   // TODO: Add lev(i) values
             matmd_lay->materialNames.push_back(buffer);
          }
          md->Add(matmd_lay);
        }

        if(SEelem.has3Dheight)
        {
          // ===================
          // 3D Height Surfaces
          // ===================
          {
            avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_3D_sfc_hght", 1, 1, 1, 0, 
                                                       SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
            mmd->xLabel = "Longitude";
            fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
            mmd->yLabel = "Latitude";
            fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
            mmd->zLabel = "Height";
            mmd->zUnits = "km";
            md->Add(mmd);
  
            // Add material MetaData for level subsetting
            // -------------------------------------------
            avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
            matmd_lay->name     = "SEelem_3D_sfc_hght_Levels";
            matmd_lay->meshName = "SEelem_3D_sfc_hght";
            matmd_lay->numMaterials=SEelem.nlev;
            for (int i=1; i <= SEelem.nlev; ++i)
            {
               char buffer[50];
               int n;
               n=sprintf(buffer, "Level %d",(i)); (void) n;   // TODO: Add lev(i) values
               matmd_lay->materialNames.push_back(buffer);
            }
            md->Add(matmd_lay);
          }

          // =================
          // 3D Height Volume
          // =================
          {
            avtMeshMetaData *mmd = new avtMeshMetaData("SEelem_3D_vol_hght", 1, 1, 1, 0, 
                                                       SpacDims, TopoDims,AVT_UNSTRUCTURED_MESH);
            mmd->xLabel = "Longitude";
            fileObject->ReadStringAttribute("E_lon", "units", mmd->xUnits);
            mmd->yLabel = "Latitude";
            fileObject->ReadStringAttribute("E_lat", "units", mmd->yUnits);
            mmd->zLabel = "Height";
            mmd->zUnits = "km";
            md->Add(mmd);
    
            // Add material MetaData for level subsetting
            // -------------------------------------------
            avtMaterialMetaData *matmd_lay = new avtMaterialMetaData;
            matmd_lay->name     = "SEelem_3D_vol_hght_Layers";
            matmd_lay->meshName = "SEelem_3D_vol_hght";
            matmd_lay->numMaterials=(SEelem.nlev-1);
            for (int i=1; i <= (SEelem.nlev-1); ++i)
            {
              char buffer[50];
              int n;
              n=sprintf(buffer, "Layer %d",(i)); (void) n;   // TODO: Add lev(i) values
              matmd_lay->materialNames.push_back(buffer);
            }
            md->Add(matmd_lay);
          }
        }
      }
    }

    debug4 << "PFC:  Done with SEelem Meshes " << endl;

    // Error check - one of the known grids must be present
    // -----------------------------------------------------
    if( !hasSEgrid && !hasSEelem)
    {
      // ERROR MESSAGE!
      //---------------
      debug4 << "PFC:ERROR: This SE file does not have SEgrid or SEelem Mesh data" << endl;
      return;
    }

    // Get the size of all of the dimensions in the file.
    // ---------------------------------------------------
    size_t *dimSizes = new size_t[nDims];
    std::map< int , std::string> dimNames;
    dimNames.clear();
    for(int i = 0; i < nDims; ++i)
    {
      char   dimName[NC_MAX_NAME+1];
      size_t dimSize;
      if((status = nc_inq_dim(fileObject->GetFileHandle(), i, dimName, &dimSize)) == NC_NOERR)
      {
        dimSizes[i] = dimSize;
      } else {
        dimSizes[i] = 1;
        fileObject->HandleError(status);
      }
      dimNames[i] = dimName;
      debug4 << "PFC: File Dimensions: i= "<<i<<"/"<<nDims<<" name="<<dimName<<" Size="<<dimSizes[i]<<endl;

    }

    // Get the indices for the dimensions we will look for
    // ----------------------------------------------------
    int TIMEdim  = -1;
    int NPdim    = -1;
    int NELEMdim = -1;
    int NLEVdim  = -1;
    int NCOLdim  = -1;
    for(int i = 0; i < nDims; ++i)
    {
      if(strncmp(dimNames[i].c_str(),"time" ,4) == 0) TIMEdim  = i;
      if(strncmp(dimNames[i].c_str(),"np"   ,2) == 0) NPdim    = i;
      if(strncmp(dimNames[i].c_str(),"nelem",5) == 0) NELEMdim = i;
      if(strncmp(dimNames[i].c_str(),"lev"  ,4) == 0) NLEVdim  = i;
      if(strncmp(dimNames[i].c_str(),"ncol" ,4) == 0) NCOLdim  = i;
    }
    debug4 << "PFC: Dim Indices: TIMEdim="<<TIMEdim<<" NELEMdim="<<NELEMdim
           <<" NLEVdim="<<NLEVdim<<" NPdim="<<NPdim<<" NCOLdim="<<NCOLdim<<endl;

    // =============================================
    // Now loop over variables in the file and process
    // the ones that match the known SE grids.
    // =============================================
    varList.clear();
    vecBaseList.clear();
    for( int nv=0; nv < nVars; ++nv)
    {
      SEvariable_t SEvar;
      SEvar.GridType = UNKNOWN;

      // Get the variable Name, Type, Dimensions, etc..
      // ----------------------------------------------
      nc_type vartype;
      char    varname[NC_MAX_NAME+1];
      int     varndims;
      int     vardims[NC_MAX_VAR_DIMS];
      int     varnatts;
      if((status = nc_inq_var(fileObject->GetFileHandle(), nv, 
                              varname,&vartype, &varndims, vardims, &varnatts)) == NC_NOERR)
      {
//        debug4 << "PFC: nv="<<nv<<" Var Type="<<vartype<<" Var Name="<<varname<< endl;
//        debug4 << "PFC:     ndims="<<varndims<<" natts="<<varnatts<< endl;
//        for( int nd=0; nd < varndims; ++nd)
//        {
//          debug4 << "PFC:    nd="<<nd<<" vardim[nd]="<<vardims[nd]<<" Name="<<dimNames[vardims[nd]]<<endl;
//        }

        // Process Variables that conform to known SE templates
        // -----------------------------------------------------
        if(varndims == 1) 
        {
          if(vardims[0] == NCOLdim) {                         // TEMPLATE: [NCOL]
            SEvar.GridType = SE_GRID_DATA;
            SEvar.GridName = SPACE_2D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [NCOL] "<<endl;
          }
        } else if(varndims == 2) {
          if((vardims[0] == TIMEdim) && 
             (vardims[1] == NCOLdim)  ) {                     // TEMPLATE: [TIME,NCOL]
            SEvar.GridType = SE_GRID_DATA;
            SEvar.GridName = TIME_2D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [TIME,NCOL] "<<endl;
          }
          else if((vardims[0] == NLEVdim) && 
                  (vardims[1] == NCOLdim)  ) {                // TEMPLATE: [NLEV,NCOL]
            SEvar.GridType = SE_GRID_DATA;
            SEvar.GridName = SPACE_3D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [NLEV,NCOL] "<<endl;
          }
        } else if(varndims == 3) {
          if((vardims[0] == TIMEdim) && 
             (vardims[1] == NLEVdim) && 
             (vardims[2] == NCOLdim)  ) {                     // TEMPLATE: [TIME,NLEV,NCOL]
            SEvar.GridType = SE_GRID_DATA;
            SEvar.GridName = TIME_3D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [TIME,NLEV,NCOL] "<<endl;
          }
          else if((vardims[0] == NELEMdim) && 
                  (vardims[1] == NPdim   ) && 
                  (vardims[2] == NPdim   )  ) {               // TEMPLATE: [NELEM,NP,NP]
            SEvar.GridType = SE_ELEM_DATA;
            SEvar.GridName = SPACE_2D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [NELEM,NP,NP] "<<endl;
          }
        } else if(varndims == 4) {
          if((vardims[0] == NELEMdim) && 
             (vardims[1] == NLEVdim ) && 
             (vardims[2] == NPdim   ) && 
             (vardims[3] == NPdim   )  ) {                    // TEMPLATE: [NELEM,NLEV,NP,NP]
            SEvar.GridType = SE_ELEM_DATA;
            SEvar.GridName = SPACE_3D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [NELEM,NLEV,NP,NP] "<<endl;
          }
          else if((vardims[0] == TIMEdim ) && 
                  (vardims[1] == NELEMdim) && 
                  (vardims[2] == NPdim   ) && 
                  (vardims[3] == NPdim   )  ) {               // TEMPLATE: [TIME,NELEM,NP,NP]
            SEvar.GridType = SE_ELEM_DATA;
            SEvar.GridName = TIME_2D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [TIME,NELEM,NP,NP] "<<endl;
          }
        } else if(varndims == 5) {
          if((vardims[0] == TIMEdim ) && 
             (vardims[1] == NELEMdim) && 
             (vardims[2] == NLEVdim ) && 
             (vardims[3] == NPdim   ) && 
             (vardims[4] == NPdim   )  ) {                    // TEMPLATE: [TIME,NELEM,NLEV,NP,NP]
            SEvar.GridType = SE_ELEM_DATA;
            SEvar.GridName = TIME_3D;
            SEvar.VarName  = varname;
            debug4 << "PFC: Var Template = [TIME,NELEM,NLEV,NP,NP] "<<endl;
          }
        }

        // For Each Known Variable, Read in attributes
        // -------------------------------------------
        if((md != 0) && (SEvar.GridType != UNKNOWN))
        {
          std::string AttVal;
          if(fileObject->ReadStringAttribute(varname,"units", AttVal))
          {
            SEvar.Units = AttVal;
          } else {
            SEvar.Units = " ";
          }
          if(fileObject->ReadStringAttribute(varname,"long_name", AttVal))
          {
            SEvar.LongName = AttVal;
          } else {
            SEvar.LongName = " ";
          }
          if(fileObject->ReadStringAttribute(varname,"vector_component", AttVal))
          {
            SEvar.VectorComponent = AttVal;
          } else {
            SEvar.VectorComponent = "0";
          }
          if(fileObject->ReadStringAttribute(varname,"vector_name", AttVal))
          {
            SEvar.VectorName = AttVal;
          } else {
            SEvar.VectorName = "0";
          }
          if(fileObject->ReadStringAttribute(varname,"vector_long_name", AttVal))
          {
            SEvar.VectorLongName = AttVal;
          } else {
            SEvar.VectorLongName = " ";
          }

          //**************************************************************
          // Impose Vector Attributes on Certain Variables:
          //   For the existing history outputs from CESM, scalar variables
          //   will not have the vector attributes need to construct vector
          //   variables for display. So here we add vector attribues to
          //   certain 'known' variables. 
          //   e.g.  scalar U and V arrays are the x,y components of
          //         horizontal winds. 
          //**************************************************************
          if(SEvar.GridType == SE_GRID_DATA)
          {
            if((strncmp(SEvar.VarName.c_str()        ,"U",1) == 0) && 
               (strncmp(SEvar.VectorComponent.c_str(),"0",1) == 0)  )
            {
              SEvar.VectorComponent = "x";
              SEvar.VectorName      = "V";
              SEvar.VectorLongName  = "Horizontal Wind";
            }
            else if((strncmp(SEvar.VarName.c_str()        ,"V",1) == 0) && 
                    (strncmp(SEvar.VectorComponent.c_str(),"0",1) == 0)  )
            {
              SEvar.VectorComponent = "y";
              SEvar.VectorName      = "V";
              SEvar.VectorLongName  = "Horizontal Wind";
            }
            if((strncmp(SEvar.VarName.c_str()        ,"FU",2) == 0) && 
               (strncmp(SEvar.VectorComponent.c_str(),"0" ,1) == 0)  )
            {
              SEvar.VectorComponent = "x";
              SEvar.VectorName      = "F_phys";
              SEvar.VectorLongName  = "Physics Horizontal Momentum Tendency";
            }
            else if((strncmp(SEvar.VarName.c_str()        ,"FV",2) == 0) && 
                    (strncmp(SEvar.VectorComponent.c_str(),"0" ,1) == 0)  )
            {
              SEvar.VectorComponent = "y";
              SEvar.VectorName      = "F_phys";
              SEvar.VectorLongName  = "Physics Horizontal Momentum Tendency";
            }
          }

          debug4 << "PFC: SEvar is KNOWN " <<endl;
          debug4 << "PFC: SEvar.GridType        ="<<SEvar.GridType       <<endl;
          debug4 << "PFC: SEvar.GridName        ="<<SEvar.GridName       <<endl;
          debug4 << "PFC: SEvar.VarName         ="<<SEvar.VarName        <<endl;
          debug4 << "PFC: SEvar.LongName        ="<<SEvar.LongName       <<endl;
          debug4 << "PFC: SEvar.Units           ="<<SEvar.Units          <<endl;
          debug4 << "PFC: SEvar.VectorComponent ="<<SEvar.VectorComponent<<endl;
          debug4 << "PFC: SEvar.VectorName      ="<<SEvar.VectorName     <<endl;
          debug4 << "PFC: SEvar.VectorLongName  ="<<SEvar.VectorLongName <<endl;
          debug4 << "PFC:    "<<endl;

          // Accumulate a base list of vectors using any Vector* attributes 
          // available in the scalar arrays.
          //   - Check the attributes for a VectorComponent entry. 
          //     if one is found, then either add the entry to an existing 
          //     vector or create a new vector. 
          // --------------------------------------------------------
          if((strncmp(SEvar.VectorComponent.c_str(),"x",1) == 0) ||
             (strncmp(SEvar.VectorComponent.c_str(),"X",1) == 0)  )
          {
            if(vecBaseList.find(SEvar.VectorName) != vecBaseList.end())
            {
              vecBaseList[SEvar.VectorName]->XvarName = SEvar.VarName;
            } else {
              SEvector_t *SEvec_list = new SEvector_t;
              SEvec_list->GridType        = SEvar.GridType;
              SEvec_list->GridName        = SEvar.GridName;
              SEvec_list->GridZmap        = SEvar.GridZmap;
              SEvec_list->GridData        = SEvar.GridData;
              SEvec_list->VectorName      = SEvar.VectorName;
              SEvec_list->LongName        = SEvar.VectorLongName;
              SEvec_list->Units           = SEvar.Units;
              SEvec_list->XvarName        = SEvar.VarName;
              SEvec_list->YvarName        = "NONE";
              SEvec_list->ZvarName        = "NONE";
              SEvec_list->MeshName        = "0";
              SEvec_list->MyName          = "0";
              vecBaseList[SEvar.VectorName] = SEvec_list;
            }
          }
          else if((strncmp(SEvar.VectorComponent.c_str(),"y",1) == 0) ||
                  (strncmp(SEvar.VectorComponent.c_str(),"Y",1) == 0)  )
          {
            if(vecBaseList.find(SEvar.VectorName) != vecBaseList.end())
            {
              vecBaseList[SEvar.VectorName]->YvarName = SEvar.VarName;
            } else {
              SEvector_t *SEvec_list = new SEvector_t;
              SEvec_list->GridType        = SEvar.GridType;
              SEvec_list->GridName        = SEvar.GridName;
              SEvec_list->GridZmap        = SEvar.GridZmap;
              SEvec_list->GridData        = SEvar.GridData;
              SEvec_list->VectorName      = SEvar.VectorName;
              SEvec_list->LongName        = SEvar.VectorLongName;
              SEvec_list->Units           = SEvar.Units;
              SEvec_list->XvarName        = "NONE";
              SEvec_list->YvarName        = SEvar.VarName;
              SEvec_list->ZvarName        = "NONE";
              SEvec_list->MeshName        = "0";
              SEvec_list->MyName          = "0";
              vecBaseList[SEvar.VectorName] = SEvec_list;
            }
          }
          else if((strncmp(SEvar.VectorComponent.c_str(),"z",1) == 0) ||
                  (strncmp(SEvar.VectorComponent.c_str(),"Z",1) == 0)  )
          {
            if(vecBaseList.find(SEvar.VectorName) != vecBaseList.end())
            {
              vecBaseList[SEvar.VectorName]->ZvarName = SEvar.VarName;
            } else {
              SEvector_t *SEvec_list = new SEvector_t;
              SEvec_list->GridType        = SEvar.GridType;
              SEvec_list->GridName        = SEvar.GridName;
              SEvec_list->GridZmap        = SEvar.GridZmap;
              SEvec_list->GridData        = SEvar.GridData;
              SEvec_list->VectorName      = SEvar.VectorName;
              SEvec_list->LongName        = SEvar.VectorLongName;
              SEvec_list->Units           = SEvar.Units;
              SEvec_list->XvarName        = "NONE";
              SEvec_list->YvarName        = "NONE";
              SEvec_list->ZvarName        = SEvar.VarName;
              SEvec_list->MeshName        = "0";
              SEvec_list->MyName          = "0";
              vecBaseList[SEvar.VectorName] = SEvec_list;
            }
          }

          // For the Given Variable Information, Construct a variable
          // and mesh names for each supported display option
          // ---------------------------------------------------------
          if((SEvar.GridName == SPACE_2D) || (SEvar.GridName == TIME_2D))
          {
            if(SEvar.GridType == SE_GRID_DATA)
            {
              {
//              SEgrid.VarName_sfc_grid   SEgrid_2D_sfc_grid ZMAP_GRID SURFACE_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEgrid_2D_sfc_grid";
                SEvar.MyName   = SEvar.VarName+"/sfc_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEgrid.has2Dheight)
              {
//              SEgrid.VarName_sfc_hght   SEgrid_2D_sfc_hght ZMAP_HGHT SURFACE_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEgrid_2D_sfc_hght";
                SEvar.MyName   = SEvar.VarName+"/sfc_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
            }
            else if(SEvar.GridType == SE_ELEM_DATA)
            {
              {
//              SEelem.VarName_sfc_grid   SEelem_2D_sfc_grid ZMAP_GRID SURFACE_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEelem_2D_sfc_grid";
                SEvar.MyName   = SEvar.VarName+"/sfc_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEelem.has2Dheight)
              {
//                SEelem.VarName_sfc_hght   SEelem_2D_sfc_hght ZMAP_HGHT SURFACE_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEelem_2D_sfc_hght";
                SEvar.MyName   = SEvar.VarName+"/sfc_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
            }
          }
          else if((SEvar.GridName == SPACE_3D) || (SEvar.GridName == TIME_3D))
          {
            if(SEvar.GridType == SE_GRID_DATA)
            {
              {
//              SEgrid.VarName_sfc_grid    SEgrid_3D_sfc_grid ZMAP_GRID SURFACE_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEgrid_3D_sfc_grid";
                SEvar.MyName   = SEvar.VarName+"/sfc_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEgrid.has3Dheight)
              {
//                SEgrid.VarName_sfc_hght    SEgrid_3D_sfc_hght ZMAP_HGHT SURFACE_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEgrid_3D_sfc_hght";
                SEvar.MyName   = SEvar.VarName+"/sfc_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEgrid.has3D)
              {
//                SEgrid.VarName_vol_grid    SEgrid_3D_vol_grid ZMAP_GRID VOLUME_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = VOLUME_DATA;
                SEvar.MeshName = "SEgrid_3D_vol_grid";
                SEvar.MyName   = SEvar.VarName+"/vol_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;// smd->centering = AVT_ZONECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if((SEgrid.has3D)&&(SEgrid.has3Dheight))
              {
//                SEgrid.VarName_vol_hght    SEgrid_3D_vol_hght ZMAP_HGHT VOLUME_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = VOLUME_DATA;
                SEvar.MeshName = "SEgrid_3D_vol_hght";
                SEvar.MyName   = SEvar.VarName+"/vol_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;// smd->centering = AVT_ZONECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
            }
            else if(SEvar.GridType == SE_ELEM_DATA)
            {
              {
//              SEelem.VarName_sfc_grid    SEelem_3D_sfc_grid ZMAP_GRID SURFACE_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEelem_3D_sfc_grid";
                SEvar.MyName   = SEvar.VarName+"/sfc_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEelem.has3Dheight)
              {
//                SEelem.VarName_sfc_hght    SEelem_3D_sfc_hght ZMAP_HGHT SURFACE_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = SURFACE_DATA;
                SEvar.MeshName = "SEelem_3D_sfc_hght";
                SEvar.MyName   = SEvar.VarName+"/sfc_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if(SEelem.has3D)
              {
//                SEelem.VarName_vol_grid    SEelem_3D_vol_grid ZMAP_GRID VOLUME_DATA
                SEvar.GridZmap = ZMAP_GRID;
                SEvar.GridData = VOLUME_DATA;
                SEvar.MeshName = "SEelem_3D_vol_grid";
                SEvar.MyName   = SEvar.VarName+"/vol_grid";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;// smd->centering = AVT_ZONECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
              if((SEelem.has3D)&&(SEelem.has3Dheight))
              {
//                SEelem.VarName_vol_hght    SEelem_3D_vol_hght ZMAP_HGHT VOLUME_DATA
                SEvar.GridZmap = ZMAP_HGHT;
                SEvar.GridData = VOLUME_DATA;
                SEvar.MeshName = "SEelem_3D_vol_hght";
                SEvar.MyName   = SEvar.VarName+"/vol_hght";
            
                // Add This entry to the Variable List and
                // then add Meta Data for Visit
                // --------------------------------------
                SEvariable_t *SEvar_list = new SEvariable_t;
                SEvar_list->GridType        = SEvar.GridType;
                SEvar_list->GridName        = SEvar.GridName;
                SEvar_list->GridZmap        = SEvar.GridZmap;
                SEvar_list->GridData        = SEvar.GridData;
                SEvar_list->VarName         = SEvar.VarName;
                SEvar_list->LongName        = SEvar.LongName;
                SEvar_list->Units           = SEvar.Units;
                SEvar_list->VectorComponent = SEvar.VectorComponent;
                SEvar_list->VectorName      = SEvar.VectorName;
                SEvar_list->MeshName        = SEvar.MeshName;
                SEvar_list->MyName          = SEvar.MyName;
                varList[SEvar.MyName] = SEvar_list;

                avtScalarMetaData *smd = new avtScalarMetaData;
                smd->name      = SEvar.MyName;    
                smd->meshName  = SEvar.MeshName;
                smd->centering = AVT_NODECENT;// smd->centering = AVT_ZONECENT;
                if(strncmp(SEvar.Units.c_str()," ",1) != 0)
                {
                  smd->hasUnits = true;
                  smd->units    = SEvar.Units;
                } else {
                  smd->hasUnits = false;
                }
                HandleMissingData(varname, smd);
                md->Add(smd);
              }
            }
          }
        } // if(SEvar.GridType != UNKNOWN)
      }
    } // for( int nv=0; nv < nVars; ++nv)

    // =============================================
    // Now loop over all variables and construct 
    // vectors using the variables that idendify 
    // themselves as vector components.
    // =============================================
    vectorList.clear();
    int vnum = 0;
    for(vecIter = vecBaseList.begin(); vecIter != vecBaseList.end(); ++vecIter)
    {
      vnum = vnum + 1;
      debug4<<" PFC: vecIter number= "<< vnum                       <<endl;
      debug4<<" PFC: VectorName= " << (*vecIter).second->VectorName <<endl;
      debug4<<" PFC: LongName= "   << (*vecIter).second->LongName   <<endl;
      debug4<<" PFC: XvarName= "   << (*vecIter).second->XvarName   <<endl;
      debug4<<" PFC: YvarName= "   << (*vecIter).second->YvarName   <<endl;
      debug4<<" PFC: ZvarName= "   << (*vecIter).second->ZvarName   <<endl;
      debug4<<" PFC: GridType= "   << (*vecIter).second->GridType   <<endl;
      debug4<<" PFC: GridName= "   << (*vecIter).second->GridName   <<endl;
      debug4<<"  "<<endl;
     
      // For the Given Vectir Information, Construct a vector
      // and mesh names for each supported display option
      // ---------------------------------------------------------
      if(((*vecIter).second->GridName == SPACE_2D) || 
         ((*vecIter).second->GridName == TIME_2D)   )
      {
        if((*vecIter).second->GridType == SE_GRID_DATA)
        {
          {
//                SEgrid.VectorName_srf_grid    SEgrid_2D_sfc_grid ZMAP_GRID SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_GRID;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_2D_sfc_grid";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_grid";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if(SEgrid.has2Dheight)
          {
//                SEgrid.VectorName_srf_hght    SEgrid_2D_sfc_hght ZMAP_HGHT SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_HGHT;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_hght";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((SEgrid.has2Dtangent)&&(SEgrid.has2Dheight))
          {
//                SEgrid.VectorName_srf_ftan    SEgrid_2D_sfc_hght ZMAP_FTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_FTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((false)&&(SEgrid.has2Dheight)) // Optionally use Visit Gradient() for tangent vectors
          {
//                SEgrid.VectorName_srf_vtan    SEgrid_2D_sfc_hght ZMAP_VTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_VTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent *est*";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
        }
        else if((*vecIter).second->GridType == SE_ELEM_DATA)
        {
          {
//                SEelem.VectorName_srf_grid    SEelem_2D_sfc_grid ZMAP_GRID SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_GRID;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_2D_sfc_grid";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_grid";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if(SEelem.has2Dheight)
          {
//                SEelem.VectorName_srf_hght    SEelem_2D_sfc_hght ZMAP_HGHT SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_HGHT;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_hght";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((SEelem.has2Dtangent)&&(SEelem.has2Dheight))
          {
//                SEelem.VectorName_srf_ftan    SEelem_2D_sfc_hght ZMAP_FTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_FTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((false)&&(SEelem.has2Dheight)) // Optionally use Visit Gradient() for tangent vectors
          {
//                SEelem.VectorName_srf_vtan    SEelem_2D_sfc_hght ZMAP_VTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_VTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_2D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent *est*";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
        }
      }
      else if(((*vecIter).second->GridName == SPACE_3D) || 
              ((*vecIter).second->GridName == TIME_3D)   )
      {
        if((*vecIter).second->GridType == SE_GRID_DATA)
        {
          {
//                SEgrid.VectorName_srf_grid    SEgrid_3D_sfc_grid ZMAP_GRID SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_GRID;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_3D_sfc_grid";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_grid";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if(SEgrid.has3Dheight)
          {
//                SEgrid.VectorName_srf_hght    SEgrid_3D_sfc_hght ZMAP_HGHT SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_HGHT;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_hght";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((SEgrid.has3Dtangent)&&(SEgrid.has3Dheight))
          {
//                SEgrid.VectorName_srf_ftan    SEgrid_3D_sfc_hght ZMAP_FTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_FTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((false)&&(SEgrid.has3Dheight)) // Optionally use Visit Gradient() for tangent vectors
          {
//                SEgrid.VectorName_srf_vtan    SEgrid_3D_sfc_hght ZMAP_VTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_VTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEgrid_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent *est*";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
        }
        else if((*vecIter).second->GridType == SE_ELEM_DATA)
        {
          {
//                SEelem.VectorName_srf_grid    SEelem_3D_sfc_grid ZMAP_GRID SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_GRID;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_3D_sfc_grid";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_grid";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if(SEelem.has3Dheight)
          {
//                SEelem.VectorName_srf_hght    SEelem_3D_sfc_hght ZMAP_HGHT SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_HGHT;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_hght";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((SEelem.has3Dtangent)&&(SEelem.has3Dheight))
          {
//                SEelem.VectorName_srf_ftan    SEelem_3D_sfc_hght ZMAP_FTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_FTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
          if((false)&&(SEelem.has3Dheight)) // Optionally use Visit Gradient() for tangent vectors
          {
//                SEelem.VectorName_srf_vtan    SEelem_3D_sfc_hght ZMAP_VTAN SURFACE_DATA

            // Add This entry to the Vector List and
            // then add Meta Data for Visit
            // --------------------------------------
            SEvector_t *SEvec_list = new SEvector_t;
            SEvec_list->GridType     = (*vecIter).second->GridType;
            SEvec_list->GridName     = (*vecIter).second->GridName;
            SEvec_list->GridZmap     = ZMAP_VTAN;
            SEvec_list->GridData     = SURFACE_DATA;
            SEvec_list->VectorName   = (*vecIter).second->VectorName;
            SEvec_list->LongName     = (*vecIter).second->LongName;
            SEvec_list->Units        = (*vecIter).second->Units;
            SEvec_list->XvarName     = (*vecIter).second->XvarName;
            SEvec_list->YvarName     = (*vecIter).second->YvarName;
            SEvec_list->ZvarName     = (*vecIter).second->ZvarName;
            SEvec_list->MeshName     = "SEelem_3D_sfc_hght";
            SEvec_list->MyName       = SEvec_list->VectorName+"/sfc_tangent *est*";
            vectorList[SEvec_list->MyName] = SEvec_list;

            avtVectorMetaData *smd = new avtVectorMetaData;
            smd->name      = SEvec_list->MyName;    
            smd->varDim    = 3;
            smd->meshName  = SEvec_list->MeshName;
            smd->centering = AVT_NODECENT;
            if(strncmp(SEvec_list->Units.c_str()," ",1) != 0)
            {
              smd->hasUnits = true;
              smd->units    = SEvec_list->Units;
            } else {
              smd->hasUnits = false;
            }
            md->Add(smd);
          }
        }
      }
    }
}

// ****************************************************************************
// Method: avtSEReader::ReturnSpatialDimensionIndices
//
// Purpose: 
//   Returns the indices of the spatial dimensions.
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 29 15:41:42 PDT 2009
//
// Modifications:
//   
// ****************************************************************************

bool
avtSEReader::ReturnSpatialDimensionIndices(const intVector &dims, int sDims[3], int &nSDims) const
{
    size_t i;
    const char *mName = "avtSEReader::ReturnValidDimensions: ";

    // Look for up to 3 valid spatial dimensions.
    nSDims = 0;
    for(i = 0; i < dims.size() && nSDims < 3; ++i)
    {
        if(dims[i] > 1 && dims[i] != TIME_DIMENSION)
        {
            sDims[nSDims++] = i;
        }
    }

    // Count the number of cells that comprise the spatial dimensions
    int nCells = 1;
    debug5 << mName << "validDims=(";
    for(i = 0; i < (size_t)nSDims; ++i)
    {
        nCells *= dims[sDims[i]];
        debug5 << dims[sDims[i]] << ", ";
    }
    debug5 << ")" << endl;

    // Count the number of cells that comprise all dimensions (except time).
    int nValues = 1;
    debug5 << mName << "actualDims=(";
    for(i = 0; i < dims.size(); ++i)
    {
        if(dims[i] != TIME_DIMENSION)
            nValues *= dims[i];
        debug5 << dims[i] << ", ";
    }
    debug5 << ")" << endl;

    return nCells == nValues;
}

// ****************************************************************************
// Method: avtSEReader::ReturnDimStartsAndCounts
//
// Purpose: 
//   Returns the dimStarts and dimCounts for all dimensions.
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 29 15:41:42 PDT 2009
//
// Modifications:
//   Mark C. Miller, Tue Aug 15 15:28:11 PDT 2006
//   Added code to support on-the-fly domain decomposition 
//
//   Mark C. Miller, Wed Aug 16 14:45:22 PDT 2006
//   Fixed possible exit without initializing all contents of starts/counts
//
//   Mark C. Miller, Tue Dec  5 18:14:58 PST 2006
//   Fixed UMRs
//
// ****************************************************************************

void
avtSEReader::ReturnDimStartsAndCounts(int timeState, const intVector &dims, 
    intVector &dimStarts, intVector &dimCounts) const
{
    //
    // Initialize starts/counts, slicing in time.
    //
    dimStarts.clear();
    dimCounts.clear();
    for (size_t i = 0; i < dims.size(); i++)
    {
        if(dims[i] == TIME_DIMENSION)
        {
            dimStarts.push_back(timeState);
            dimCounts.push_back(1);
        }
        else
        {
            dimStarts.push_back(0);
            dimCounts.push_back(dims[i]);
        }
    }

    // Get the indices of the spatial dimensions
    int spatialDimIndices[3] = {0,0,0}, nSpatialDims = 0;
    ReturnSpatialDimensionIndices(dims, spatialDimIndices, nSpatialDims);

    // Compute how many values that makes.
    int spatialSizes[3] = {0,0,0};
    int nValues = 1;
    for(int i = 0; i < nSpatialDims; ++i)
    {
        spatialSizes[i] = dims[spatialDimIndices[i]];
        nValues *= spatialSizes[i];
    }

    //
    // We won't decompose something that is smaller than some threshold
    //
    if (nValues < 100000)
        return;

    //
    // Here's where we let VisIt decide how to divide up the spatial dimensions.
    //

    //
    // Above, we're counting nodes (e.g. mesh lines). The decomposition
    // stuff operates on zones.
    //
    int validZDims[3] = {0,0,0};
    for (int i = 0; i < nSpatialDims; i++)
        validZDims[i] = spatialSizes[i]-1;

    //
    // Ok, now compute the zone-oriented domain decomposition
    //
    int domCount[3] = {0, 0, 0};
    avtDatabase::ComputeRectilinearDecomposition(
        nSpatialDims, procCount,
        validZDims[0], validZDims[1], validZDims[2],
        &domCount[0], &domCount[1], &domCount[2]);

    debug5 << "Decomposition: " << domCount[0] << ", "
           << domCount[1] << ", " << domCount[2] << endl;

    //
    // Determine this processor's logical domain (e.g. domain ijk) indices
    //
    int domLogicalCoords[3] = {0, 0, 0};
    avtDatabase::ComputeDomainLogicalCoords(nSpatialDims, domCount, procNum,
        domLogicalCoords);

    debug5 << "Processor " << procNum << " domain logical coords: "
           << domLogicalCoords[0] << ", " << domLogicalCoords[1] << ", "
           << domLogicalCoords[2] << endl;

    //
    // compute the bounds, in terms of output zone numbers,
    // of this processor's domain. Store the divided spatial domain's starts
    // and counts, overwriting the previous values.
    //
    debug5 << "Processor " << procNum << " zone-centered bounds..." << endl;
    for (int i = 0; i < nSpatialDims; i++)
    {
        avtDatabase::ComputeDomainBounds(validZDims[i], domCount[i], domLogicalCoords[i],
            &dimStarts[spatialDimIndices[i]], &dimCounts[spatialDimIndices[i]]);
        dimCounts[spatialDimIndices[i]]++; // convert to # of zones to # of nodes  
        debug5 << "   start[" << i << "] = " << dimStarts[spatialDimIndices[i]]
               << ",  count[" << i << "] = " << dimCounts[spatialDimIndices[i]]-1 << endl;
    }
}

// ****************************************************************************
// Method: avtSEReader::GetMesh
//
// Purpose: 
//   Returns the specified mesh.
//
// Arguments:
//   timeState : The time state.
//   var       : The name of the mesh to create.
//
// Returns:    A vtkDataSet containing the mesh or 0.
// 
// Programmer: Patrick Callaghan
// Creation:   Oct 18 2016
//
// Modifications:
//
// ****************************************************************************

vtkDataSet *
avtSEReader::GetMesh(int timeState, const char *var)
{
    const char *mName  = "avtSEReader::GetMesh: ";
    vtkDataSet *retval = 0;
    int         origin = 1;

  debug4 << mName << "var=" << var << endl;
  debug4 << "PFC: Entering GetMesh  var =" << var << endl;

    if(strncmp(var,"SEgrid_2D_sfc_grid",18) == 0)
    {
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = SEgrid.lon;
       float *yc = SEgrid.lat;
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = 0;
       }

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters);
       vtkIdType verts[4];
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] - origin;
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEgrid_3D_vol_grid",18) == 0)
    {
       // make a vector of level indices
       // -----------------------------------
       float *lon3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lat3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lev3 = new float [SEgrid.ncol*SEgrid.nlev];
       for( int i = 0; i < SEgrid.ncol ; ++i)
       for( int k = 0; k < SEgrid.nlev ; ++k)
       {
          lon3[i + SEgrid.ncol*k] = SEgrid.lon[i];
          lat3[i + SEgrid.ncol*k] = SEgrid.lat[i];
          lev3[i + SEgrid.ncol*k] = SEgrid.nlev-k;
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // -----------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol*SEgrid.nlev);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = lev3;
       for( int k = 0; k < SEgrid.nlev ; ++k)
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] lon3;
       delete [] lat3;
       delete [] lev3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters*SEgrid.nlev);
       vtkIdType verts[8];
       for(int k = 1; k < SEgrid.nlev    ; ++k)
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] +(k-1)*SEgrid.ncol - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[4] = SEgrid.connect[i                    ] +(k  )*SEgrid.ncol - origin;
          verts[5] = SEgrid.connect[i +   SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          verts[6] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          verts[7] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          ugrid->InsertNextCell(VTK_HEXAHEDRON,8,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEgrid_3D_sfc_grid",18) == 0)
    {
       // make a vector of level indices
       // --------------------------------
       float *lon3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lat3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lev3 = new float [SEgrid.ncol*SEgrid.nlev];
       for( int i = 0; i < SEgrid.ncol ; ++i)
       for( int k = 0; k < SEgrid.nlev ; ++k)
       {
          lon3[i + SEgrid.ncol*k] = SEgrid.lon[i];
          lat3[i + SEgrid.ncol*k] = SEgrid.lat[i];
          lev3[i + SEgrid.ncol*k] = SEgrid.nlev-k;
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol*SEgrid.nlev);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = lev3;
       for( int k = 0; k < SEgrid.nlev ; ++k)
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] lon3;
       delete [] lat3;
       delete [] lev3;

       // Create a vtkUnstructuredGrid to contian the cells
       // -----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters*SEgrid.nlev);
       vtkIdType verts[4];
       for(int k = 0; k < SEgrid.nlev    ; ++k)
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] +(k*SEgrid.ncol) - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEgrid_2D_sfc_hght",18) == 0)
    {
       // Read in PHIS and convert to height in km.
       // ------------------------------------------
       float        *hght = new float[SEgrid.ncol];
       const int start[2] = {timeState,         0 };
       const int count[2] = {       1 ,SEgrid.ncol};
       fileObject->ReadVariableInto("PHIS", FLOATARRAY_TYPE, start, count, hght);
       for(int i = 0; i < SEgrid.ncol; ++i)
       {
          hght[i] = hght[i]/9800.;
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = SEgrid.lon;
       float *yc = SEgrid.lat;
       float *zc = hght;
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;

       // Create a vtkUnstructuredGrid to contian the cells
       // ---------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters);
       vtkIdType verts[4];
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] - origin;
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEgrid_3D_vol_hght",18) == 0)
    {
       // Read in Z3 and convert to height in km.
       // ------------------------------------------
       float        *hght = new float[SEgrid.nlev*SEgrid.ncol];
       const int start[3] = {timeState,         0 ,         0 };
       const int count[3] = {       1 ,SEgrid.nlev,SEgrid.ncol};
       fileObject->ReadVariableInto("Z3", FLOATARRAY_TYPE, start, count, hght);
       for(int i = 0; i < SEgrid.nlev*SEgrid.ncol; ++i)
       {
          hght[i] = hght[i]/1000.;
       }
     
       // make a vector of lats/lons
       // -----------------------------
       float *lon3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lat3 = new float [SEgrid.ncol*SEgrid.nlev];
       for( int i = 0; i < SEgrid.ncol ; ++i)
       for( int k = 0; k < SEgrid.nlev ; ++k)
       {
          lon3[i + SEgrid.ncol*k] = SEgrid.lon[i];
          lat3[i + SEgrid.ncol*k] = SEgrid.lat[i];
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // -----------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol*SEgrid.nlev);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = hght;
       for( int k = 0; k < SEgrid.nlev ; ++k)
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;
       delete [] lon3;
       delete [] lat3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters*SEgrid.nlev);
       vtkIdType verts[8];
       for(int k = 1; k < SEgrid.nlev    ; ++k)
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] +(k-1)*SEgrid.ncol - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k-1)*SEgrid.ncol - origin;
          verts[4] = SEgrid.connect[i                    ] +(k  )*SEgrid.ncol - origin;
          verts[5] = SEgrid.connect[i +   SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          verts[6] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          verts[7] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k  )*SEgrid.ncol - origin;
          ugrid->InsertNextCell(VTK_HEXAHEDRON,8,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEgrid_3D_sfc_hght",18) == 0)
    {
       // Read in Z3 and convert to height in km.
       // ------------------------------------------
       float        *hght = new float[SEgrid.nlev*SEgrid.ncol];
       const int start[3] = {timeState,         0 ,         0 };
       const int count[3] = {       1 ,SEgrid.nlev,SEgrid.ncol};
       fileObject->ReadVariableInto("Z3", FLOATARRAY_TYPE, start, count, hght);
       for(int i = 0; i < SEgrid.nlev*SEgrid.ncol; ++i)
       {
          hght[i] = hght[i]/1000.;
       }
     
       // make a vector of lats/lons
       // -----------------------------
       float *lon3 = new float [SEgrid.ncol*SEgrid.nlev];
       float *lat3 = new float [SEgrid.ncol*SEgrid.nlev];
       for( int i = 0; i < SEgrid.ncol ; ++i)
       for( int k = 0; k < SEgrid.nlev ; ++k)
       {
          lon3[i + SEgrid.ncol*k] = SEgrid.lon[i];
          lat3[i + SEgrid.ncol*k] = SEgrid.lat[i];
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(SEgrid.ncol*SEgrid.nlev);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = hght;
       for( int k = 0; k < SEgrid.nlev ; ++k)
       for( int i = 0; i < SEgrid.ncol ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;
       delete [] lon3;
       delete [] lat3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ---------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate(SEgrid.ncenters*SEgrid.nlev);
       vtkIdType verts[4];
       for(int k = 0; k < SEgrid.nlev    ; ++k)
       for(int i = 0; i < SEgrid.ncenters; ++i)
       {
          verts[0] = SEgrid.connect[i                    ] +(k*SEgrid.ncol) - origin;
          verts[1] = SEgrid.connect[i +   SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          verts[2] = SEgrid.connect[i + 2*SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          verts[3] = SEgrid.connect[i + 3*SEgrid.ncenters] +(k*SEgrid.ncol) - origin;
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    } 
    else if(strncmp(var,"SEelem_2D_sfc_grid",18) == 0)
    {
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       int   npts = SEelem.np*SEelem.np*SEelem.nelem;
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc  = SEelem.lon;
       float *yc  = SEelem.lat;
       for( int i = 0; i < npts ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = 0;
       }

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem);
       vtkIdType verts[4];
       for(int ie = 0; ie <  SEelem.nelem; ++ie)
       for(int j0 = 0; j0 < (SEelem.np-1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np-1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + ie*SEelem.np);
          verts[1] = i0+1 + SEelem.np*(j0   + ie*SEelem.np);
          verts[2] = i0+1 + SEelem.np*(j0+1 + ie*SEelem.np);
          verts[3] = i0   + SEelem.np*(j0+1 + ie*SEelem.np);
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEelem_3D_vol_grid",18) == 0)
    {
       // make a vector of level indices
       // -----------------------------------
       int   npts2 = SEelem.np*SEelem.np*SEelem.nelem;
       int   npts3 = npts2*SEelem.nlev;
       float *lon3 = new float [npts3];
       float *lat3 = new float [npts3];
       float *lev3 = new float [npts3];
       int   ind2, ind3;
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
       for(int j0 = 0; j0 <  SEelem.np     ; ++j0)
       for(int i0 = 0; i0 <  SEelem.np     ; ++i0)
       {
          ind3 = i0+ SEelem.np*(j0+ SEelem.np*(l0+ ie*SEelem.nlev));
          ind2 = i0+ SEelem.np*(j0+ ie*SEelem.np);
          lon3[ind3] = SEelem.lon[ind2];
          lat3[ind3] = SEelem.lat[ind2];
          lev3[ind3] = SEelem.nlev-l0;
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // -----------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts3);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = lev3;
       for( int k = 0; k < SEelem.nlev; ++k)
       for( int i = 0; i < npts2      ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] lon3;
       delete [] lat3;
       delete [] lev3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem*(SEelem.nlev-1));
       vtkIdType verts[8];
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 < (SEelem.nlev-1); ++l0)
       for(int j0 = 0; j0 < (SEelem.np  -1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np  -1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[1] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[2] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[3] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[4] = i0   + SEelem.np*(j0   + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[5] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[6] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[7] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0+1 + ie*SEelem.nlev));
          ugrid->InsertNextCell(VTK_HEXAHEDRON,8,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEelem_3D_sfc_grid",18) == 0)
    {
       // make a vector of level indices
       // --------------------------------
       int   npts2 = SEelem.np*SEelem.np*SEelem.nelem;
       int   npts3 = npts2*SEelem.nlev;
       float *lon3 = new float [npts3];
       float *lat3 = new float [npts3];
       float *lev3 = new float [npts3];
       int   ind2, ind3;
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
       for(int j0 = 0; j0 <  SEelem.np     ; ++j0)
       for(int i0 = 0; i0 <  SEelem.np     ; ++i0)
       {
          ind3 = i0+ SEelem.np*(j0+ SEelem.np*(l0+ ie*SEelem.nlev));
          ind2 = i0+ SEelem.np*(j0+ ie*SEelem.np);
          lon3[ind3] = SEelem.lon[ind2];
          lat3[ind3] = SEelem.lat[ind2];
          lev3[ind3] = SEelem.nlev-l0;
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts3);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = lev3;
       for( int k = 0; k < SEelem.nlev; ++k)
       for( int i = 0; i < npts2      ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] lon3;
       delete [] lat3;
       delete [] lev3;

       // Create a vtkUnstructuredGrid to contian the cells
       // -----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem*SEelem.nlev);
       vtkIdType verts[4];
       for(int ie = 0; ie <  SEelem.nelem; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev ; ++l0)
       for(int j0 = 0; j0 < (SEelem.np-1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np-1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[1] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[2] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[3] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEelem_2D_sfc_hght",18) == 0)
    {
       // Read in Zs and convert to height in km.
       // ------------------------------------------
       int          npts2 = SEelem.np*SEelem.np*SEelem.nelem;
       float        *hght = new float[npts2];
       const int start[4] = {timeState,          0 ,        0,        0};
       const int count[4] = {       1 ,SEelem.nelem,SEelem.np,SEelem.np};
       fileObject->ReadVariableInto("Zs", FLOATARRAY_TYPE, start, count, hght);
//       for(int i = 0; i < npts2; ++i)
//       {
//          hght[i] = hght[i]/1000.;
//       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts2);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = SEelem.lon;
       float *yc = SEelem.lat;
       float *zc = hght;
       for( int i = 0; i < npts2; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;

       // Create a vtkUnstructuredGrid to contian the cells
       // ---------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem);
       vtkIdType verts[4];
       for(int ie = 0; ie <  SEelem.nelem; ++ie)
       for(int j0 = 0; j0 < (SEelem.np-1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np-1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + ie*SEelem.np);
          verts[1] = i0+1 + SEelem.np*(j0   + ie*SEelem.np);
          verts[2] = i0+1 + SEelem.np*(j0+1 + ie*SEelem.np);
          verts[3] = i0   + SEelem.np*(j0+1 + ie*SEelem.np);
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEelem_3D_vol_hght",18) == 0)
    {
       // Read in Z and convert to height in km.
       // ------------------------------------------
       int          npts2 = SEelem.np*SEelem.np*SEelem.nelem;
       int          npts3 = npts2*SEelem.nlev;
       float        *hght = new float[npts3];
       const int start[5] = {timeState,          0 ,          0,        0,       0 };
       const int count[5] = {       1 ,SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};
       fileObject->ReadVariableInto("Z", FLOATARRAY_TYPE, start, count, hght);
//       for(int i = 0; i < npts3; ++i)
//       {
//          hght[i] = hght[i]/1000.;
//       }
     
       // make a vector of lats/lons
       // -----------------------------
       float *lon3 = new float [npts3];
       float *lat3 = new float [npts3];
       int   ind2, ind3;
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
       for(int j0 = 0; j0 <  SEelem.np     ; ++j0)
       for(int i0 = 0; i0 <  SEelem.np     ; ++i0)
       {
          ind3 = i0+ SEelem.np*(j0+ SEelem.np*(l0+ ie*SEelem.nlev));
          ind2 = i0+ SEelem.np*(j0+ ie*SEelem.np);
          lon3[ind3] = SEelem.lon[ind2];
          lat3[ind3] = SEelem.lat[ind2];
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // -----------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts3);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = hght;
       for( int k = 0; k < SEelem.nlev; ++k)
       for( int i = 0; i < npts2      ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;
       delete [] lon3;
       delete [] lat3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ----------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem*(SEelem.nlev-1));
       vtkIdType verts[8];
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 < (SEelem.nlev-1); ++l0)
       for(int j0 = 0; j0 < (SEelem.np  -1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np  -1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[1] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[2] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[3] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[4] = i0   + SEelem.np*(j0   + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[5] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[6] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0+1 + ie*SEelem.nlev));
          verts[7] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0+1 + ie*SEelem.nlev));
          ugrid->InsertNextCell(VTK_HEXAHEDRON,8,verts);
       }
       return ugrid;
    }
    else if(strncmp(var,"SEelem_3D_sfc_hght",18) == 0)
    {
       // Read in Z and convert to height in km.
       // ------------------------------------------
       int          npts2 = SEelem.np*SEelem.np*SEelem.nelem;
       int          npts3 = npts2*SEelem.nlev;
       float        *hght = new float[npts3];
       const int start[5] = {timeState,          0 ,          0,        0,       0 };
       const int count[5] = {       1 ,SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};
       fileObject->ReadVariableInto("Z", FLOATARRAY_TYPE, start, count, hght);
//       for(int i = 0; i < npts3; ++i)
//       {
//          hght[i] = hght[i]/1000.;
//       }
//DIAG       {
//DIAG         int ind3 = 0;
//DIAG         for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
//DIAG         for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
//DIAG         for(int j0 = 0; j0 <  SEelem.np     ; ++j0)
//DIAG         for(int i0 = 0; i0 <  SEelem.np     ; ++i0)
//DIAG         {
//DIAG           hght[ind3] = SEelem.nlev-l0;
//DIAG           ind3 = ind3 + 1;
//DIAG         }
//DIAG       }
     
       // make a vector of lats/lons
       // -----------------------------
       float *lon3 = new float [npts3];
       float *lat3 = new float [npts3];
       int   ind2, ind3;
       for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
       for(int j0 = 0; j0 <  SEelem.np     ; ++j0)
       for(int i0 = 0; i0 <  SEelem.np     ; ++i0)
       {
          ind3 = i0+ SEelem.np*(j0+ SEelem.np*(l0+ ie*SEelem.nlev));
          ind2 = i0+ SEelem.np*(j0+ ie*SEelem.np);
          lon3[ind3] = SEelem.lon[ind2];
          lat3[ind3] = SEelem.lat[ind2];
       }
     
       // Create the vtkPoints object and copy grid locations into it
       // ------------------------------------------------------------
       vtkPoints *points = vtkPoints::New();
       points->SetNumberOfPoints(npts3);
       float *pts = (float *) points->GetVoidPointer(0);
       float *xc = lon3;
       float *yc = lat3;
       float *zc = hght;
       for( int k = 0; k < SEelem.nlev; ++k)
       for( int i = 0; i < npts2      ; ++i)
       {
          *pts++ = *xc++;
          *pts++ = *yc++;
          *pts++ = *zc++;
       }
       delete [] hght;
       delete [] lon3;
       delete [] lat3;

       // Create a vtkUnstructuredGrid to contian the cells
       // ---------------------------------------------------
       vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::New();
       ugrid->SetPoints(points);
       points->Delete();
       ugrid->Allocate((SEelem.np-1)*(SEelem.np-1)*SEelem.nelem*SEelem.nlev);
       vtkIdType verts[4];
       for(int ie = 0; ie <  SEelem.nelem; ++ie)
       for(int l0 = 0; l0 <  SEelem.nlev ; ++l0)
       for(int j0 = 0; j0 < (SEelem.np-1); ++j0)
       for(int i0 = 0; i0 < (SEelem.np-1); ++i0)
       {
          verts[0] = i0   + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[1] = i0+1 + SEelem.np*(j0   + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[2] = i0+1 + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          verts[3] = i0   + SEelem.np*(j0+1 + SEelem.np*(l0   + ie*SEelem.nlev));
          ugrid->InsertNextCell(VTK_QUAD,4,verts);
       }
       return ugrid;
    } 
}

// ****************************************************************************
// Method: avtSEReader::GetVar
//
// Purpose: 
//   Returns the data for the specified variable.
//
// Arguments:
//   timeState : The time state to read.
//   var       : The name of the variable to read.
//
// Returns:    The data or 0.
//
// Note:       
//
// Programmer: Patrick Callaghan
// Creation:   Oct 18 2016
//
// Modifications:
//
// ****************************************************************************

#define READVAR(VTKTYPE) \
        {\
            VTKTYPE *arr = VTKTYPE::New();\
            arr->SetNumberOfComponents(1);\
            arr->SetNumberOfTuples(nValues);\
            debug4 << "Allocated a " << \
                    #VTKTYPE \
                   << " of " << nValues << " elements" << endl; \
            int *rdimStarts = new int[dimStarts.size()]; \
            int *rdimCounts = new int[dimCounts.size()]; \
            for (size_t kk = 0; kk < dimStarts.size(); kk++)\
            {\
                rdimStarts[ndims-kk-1] = dimStarts[kk];\
                rdimCounts[ndims-kk-1] = dimCounts[kk];\
            }\
            if(fileObject->ReadVariableInto(bare_var, t, rdimStarts, rdimCounts,\
                                            arr->GetVoidPointer(0)))\
                retval = arr;\
            else\
                arr->Delete();\
            delete [] rdimStarts;\
            delete [] rdimCounts;\
        }


vtkDataArray *
avtSEReader::GetVar(int timeState, const char *var)
{
    const char *mName = "avtSEReader::GetVar: ";
    debug4 << mName << "var=" << var << endl;
    vtkDataArray *retval = 0;

    debug4 << "PFC: GetVar varList Size: "<< varList.size() << endl;
    debug4 << "PFC: GetVar Lookup: "<< varList[var]->VarName <<" : "<<varList[var]->MyName << endl;

    if(varList[var]->GridType == SE_GRID_DATA)
    {
      if(varList[var]->GridName == TIME_2D )
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEgrid.ncol); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[2] = {timeState,         0 };
        const int count[2] = {        1,SEgrid.ncol};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      else if(varList[var]->GridName == SPACE_2D)
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEgrid.ncol); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[1] = {         0 };
        const int count[1] = {SEgrid.ncol};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      else if(varList[var]->GridName == TIME_3D )
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEgrid.nlev*SEgrid.ncol); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[3] = {timeState,         0 ,         0 };
        const int count[3] = {        1,SEgrid.nlev,SEgrid.ncol};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      else if(varList[var]->GridName == SPACE_3D)
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEgrid.nlev*SEgrid.ncol); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[2] = {         0 ,         0 };
        const int count[2] = {SEgrid.nlev,SEgrid.ncol};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
    }
    else if(varList[var]->GridType == SE_ELEM_DATA)
    {
      if(varList[var]->GridName == TIME_2D )
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEelem.nelem*SEelem.np*SEelem.np); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[4] = {timeState,          0 ,       0 ,       0 };
        const int count[4] = {        1,SEelem.nelem,SEelem.np,SEelem.np};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      if(varList[var]->GridName == SPACE_2D)
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEelem.nelem*SEelem.np*SEelem.np); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[3] = {          0 ,       0 ,       0 };
        const int count[3] = {SEelem.nelem,SEelem.np,SEelem.np};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      else if(varList[var]->GridName == TIME_3D )
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEelem.nelem*SEelem.nlev*SEelem.np*SEelem.np); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[5] = {timeState,          0 ,         0 ,       0 ,       0 };
        const int count[5] = {        1,SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
      else if(varList[var]->GridName == SPACE_3D)
      {
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfTuples(SEelem.nelem*SEelem.nlev*SEelem.np*SEelem.np); 
        float *Sdata = (float *) arr->GetVoidPointer(0);
        const int start[4] = {          0 ,         0 ,       0 ,       0 };
        const int count[4] = {SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};
        fileObject->ReadVariableInto(varList[var]->VarName.c_str(), 
                                     FLOATARRAY_TYPE, start, count, Sdata);
        return arr;
      }
    }
    return retval;
}

// ****************************************************************************
// Method: avtSEReader::GetVectorVar
//
// Purpose: 
//   Returns the Vector data for the specified variable.
//
// Arguments:
//   timeState : The time state to read.
//   var       : The name of the variable to read.
//
// Returns:    The data or 0.
//
// Note:       
//
// Programmer: Patrick Callaghan
// Creation:   Sept 25 2016
//
// Modifications:
//
// ****************************************************************************

vtkDataArray *
avtSEReader::GetVectorVar(int timeState, const char *var)
{
    const char *mName = "avtSEReader::GetVectorVar: ";
    debug4 << mName << "var=" << var << endl;
    vtkDataArray *retval = 0;

    debug4 << "PFC: GetVectorVar Lookup: "<< vectorList[var]->VectorName <<" : "
           << " MyName="                  << vectorList[var]->MyName     << endl;
    debug4 << "PFC:   GridType="          << vectorList[var]->GridType 
           << "       GridName="          << vectorList[var]->GridName   << endl;
    debug4 << "PFC: ZvarName= "<<  vectorList[var]->ZvarName
           << "   : ZvarName= "<<  vectorList[var]->ZvarName.c_str() << endl;

    if(vectorList[var]->GridType == SE_GRID_DATA)
    {
      if(vectorList[var]->GridName == TIME_2D )
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEgrid.ncol;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[2] = {timeState,         0 };
        const int count[2] = {        1,SEgrid.ncol};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else if(vectorList[var]->GridName == SPACE_2D)
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEgrid.ncol;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[1] = {         0 };
        const int count[1] = {SEgrid.ncol};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else if(vectorList[var]->GridName == TIME_3D )
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEgrid.ncol*SEgrid.nlev;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[3] = {timeState,         0 ,         0 };
        const int count[3] = {        1,SEgrid.nlev,SEgrid.ncol};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
     debug4 << "PFC: Setting Wval==0" << endl;
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        debug4 << "PFC: Done reading 3D U,V and setting W=0" << endl;

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        debug4 << "PFC: Done loading vtk Vectors" << endl;
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else if(vectorList[var]->GridName == SPACE_3D)
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEgrid.ncol*SEgrid.nlev;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[2] = {         0 ,         0 };
        const int count[2] = {SEgrid.nlev,SEgrid.ncol};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else
      {
        return retval;
      }
    }
    else if(vectorList[var]->GridType == SE_ELEM_DATA)
    {
      if(vectorList[var]->GridName == TIME_2D )
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEelem.nelem*SEelem.np*SEelem.np;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[4] = {timeState,          0 ,       0 ,       0 };
        const int count[4] = {        1,SEelem.nelem,SEelem.np,SEelem.np};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      if(vectorList[var]->GridName == SPACE_2D)
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEelem.nelem*SEelem.np*SEelem.np;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[3] = {          0 ,       0 ,       0 };
        const int count[3] = {SEelem.nelem,SEelem.np,SEelem.np};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else if(vectorList[var]->GridName == TIME_3D )
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEelem.nelem*SEelem.nlev*SEelem.np*SEelem.np;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[5] = {timeState,          0 ,         0 ,       0 ,       0 };
        const int count[5] = {        1,SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else if(vectorList[var]->GridName == SPACE_3D)
      {
        // Read in the Vector data
        // ----------------------------
        int   nodes = SEelem.nelem*SEelem.nlev*SEelem.np*SEelem.np;
        float *Uval = new float [nodes];
        float *Vval = new float [nodes];
        float *Wval = new float [nodes];
        const int start[4] = {          0 ,         0 ,       0 ,       0 };
        const int count[4] = {SEelem.nelem,SEelem.nlev,SEelem.np,SEelem.np};

        if(strncmp(vectorList[var]->XvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->XvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Uval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Uval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->YvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->YvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Vval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Vval[i] = 0.;
          }
        }
        if(strncmp(vectorList[var]->ZvarName.c_str(),"NONE",4) != 0) 
        {
          fileObject->ReadVariableInto(vectorList[var]->ZvarName.c_str(), 
                                       FLOATARRAY_TYPE, start, count, Wval);
        } else {
          for(int i=0; i<nodes; ++i)
          {
            Wval[i] = 0.;
          }
        }

        // Now optionally transform the vectors relative to the 
        // local height surface tangent.
        // --------------------------------------------------------
        if(vectorList[var]->GridZmap == ZMAP_FTAN )
        {
          // Read in Height Gradients and use them to transform the vectors
          // ---------------------------------------------------------------
        } 
        else if(vectorList[var]->GridZmap == ZMAP_VTAN )
        {
          // Read in Height Values and use the results from VisIt's gradient 
          // operator to transform the vectors
          // ---------------------------------------------------------------
        } 

        // Uval,Vval,Wval now contain the desired vector data, 
        // Create the vtk Vector array for display.
        // ---------------------------------------------------
        vtkFloatArray *arr = vtkFloatArray::New();
        arr->SetNumberOfComponents(3); 
        arr->SetNumberOfTuples(nodes); 
        float *Wind = (float *) arr->GetVoidPointer(0);
        float *c1 = Uval;
        float *c2 = Vval;
        float *c3 = Wval;
        for(int i=0; i<nodes; ++i)
        {
          *Wind++ = *c1++;
          *Wind++ = *c2++;
          *Wind++ = *c3++;
        }
        delete [] Uval;
        delete [] Vval;
        delete [] Wval;
        return arr;
      }
      else
      {
        return retval;
      }
    }
    else
    {
      return retval;
    }
}

// ****************************************************************************
// Method: avtSEReader::GetAuxiliaryData
//
// Purpose: 
//   Returns the material object for a mesh.
//
// Arguments:
//
// Returns:    The data or 0.
//
// Note:       
//
// Programmer: Patrick Callaghan
// Creation:   Sept 24 2016
//
// Modifications:
//
// ****************************************************************************

void *
avtSEReader::GetAuxiliaryData(const char *var , int ts,
                              const char *type, void *args, DestructorFunction &df)
{
   void      *retval = 0;
   const char *mName = "avtSEReader::GetAuxiliaryData: ";

   if(strcmp(type, AUXILIARY_DATA_MATERIAL) == 0)
   {
   debug4 << mName << " PFC: timestate=" << ts   
                   << " VarName="        << var  
                   << " type="           << type 
                   << " args="           << args << endl;
 
      if((strncmp(var,"SEgrid_3D_vol_grid_Layers",25) == 0) ||
         (strncmp(var,"SEgrid_3D_vol_hght_Layers",25) == 0)  )
      {
         // Create material names and matnos array
         // ---------------------------------------
         int     nmat = SEgrid.nlev-1;
         int  *matnos = new int   [nmat];
         char **names = new char *[nmat];
         for(int i = 0; i < nmat; ++i)
         {
            matnos[i] = i + 1;
            char *buffer=new char[20];
            int n;
            n=sprintf(buffer, "Layer %d",matnos[i]); (void) n;
            names[i] = (char *)buffer;
         }

         // Assign mmaterial numbers to mesh
         // ---------------------------------
         int *matlist =new int[nmat*SEgrid.ncenters];
         for (int i = 0; i< SEgrid.ncenters; ++i)
         for (int k = 0; k< nmat           ; ++k)
         {
            matlist[i + k*SEgrid.ncenters] = matnos[nmat-1-k];
         }

         // Create the Materials
         // ------------------------
         int dims[3] = {1,1,1};
         int   ndims = 1;
         dims[0] = SEgrid.ncenters*nmat;
         avtMaterial *mat= new avtMaterial( nmat, matnos, names, ndims, dims,
                                            0, matlist, 0,   // length of mix arrays
                                                        0,   // mix_mat array
                                                        0,   // mix_next array
                                                        0,   // mix_zone array
                                                        0);  // mix_vf array 
         delete [] matlist;
         delete [] matnos;
         for (int i = 0; i < nmat; ++i) delete [] names[i];
         delete [] names;
         retval =(void *)mat;
         df     = avtMaterial::Destruct;
      }
      else if((strncmp(var,"SEgrid_3D_sfc_grid_Levels",25) == 0) ||
              (strncmp(var,"SEgrid_3D_sfc_hght_Levels",25) == 0)  )
      {
         // Create material names and matnos array
         // ---------------------------------------
         int     nmat = SEgrid.nlev;
         int  *matnos = new int   [nmat];
         char **names = new char *[nmat];
         for(int i = 0; i < nmat; ++i)
         {
            matnos[i]    = i + 1;
            char *buffer = new char[20];
            int n;
            n=sprintf(buffer, "Level %d",matnos[i]); (void) n;
            names[i] = (char *)buffer;
         }

         // Assign mmaterial numbers to mesh
         // ---------------------------------
         int *matlist =new int[nmat*SEgrid.ncenters];
         for (int i = 0; i< SEgrid.ncenters; ++i)
         for (int k = 0; k< nmat           ; ++k)
         {
            matlist[i + k*SEgrid.ncenters] = matnos[nmat-1-k];
         }

         // Create the Materials
         // ---------------------
         int dims[3] = {1,1,1};
         int   ndims = 1;
         dims[0] = SEgrid.ncenters*nmat;
         avtMaterial *mat= new avtMaterial( nmat, matnos, names, ndims, dims,
                                            0, matlist, 0,   // length of mix arrays
                                                        0,   // mix_mat array
                                                        0,   // mix_next array
                                                        0,   // mix_zone array
                                                        0);  // mix_vf array 
         delete [] matlist;
         delete [] matnos;
         for (int i = 0; i < nmat; ++i) 
           delete [] names[i];
         delete [] names;
         retval =(void *)mat;
         df     = avtMaterial::Destruct;
      }
      else if((strncmp(var,"SEelem_3D_sfc_grid_Levels",25) == 0) ||
              (strncmp(var,"SEelem_3D_sfc_hght_Levels",25) == 0)  )
      {
         // Create material names and matnos array
         // ---------------------------------------
         int   ncells = (SEelem.np-1)*(SEelem.np-1)*SEelem.nelem;
         int     nmat =  SEelem.nlev;
         int  *matnos = new int   [nmat];
         char **names = new char *[nmat];
         for(int i = 0; i < nmat; ++i)
         {
            matnos[i]    = i + 1;
            char *buffer = new char[20];
            int n;
            n=sprintf(buffer, "Level %d",matnos[i]); (void) n;
            names[i] = (char *)buffer;
         }

         // Assign mmaterial numbers to mesh
         // ---------------------------------
         int *matlist =new int[nmat*ncells];
         int   ind3 = 0;
         for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
         for(int l0 = 0; l0 <  SEelem.nlev   ; ++l0)
         for(int j0 = 0; j0 < (SEelem.np-1)  ; ++j0)
         for(int i0 = 0; i0 < (SEelem.np-1)  ; ++i0)
         {
            matlist[ind3] = matnos[nmat-1-l0];
            ind3 = ind3 + 1;
         }
     
         // Create the Materials
         // ---------------------
         int dims[3] = {1,1,1};
         int   ndims = 1;
         dims[0] = ncells*nmat;
         avtMaterial *mat= new avtMaterial( nmat, matnos, names, ndims, dims,
                                            0, matlist, 0,   // length of mix arrays
                                                        0,   // mix_mat array
                                                        0,   // mix_next array
                                                        0,   // mix_zone array
                                                        0);  // mix_vf array 
         delete [] matlist;
         delete [] matnos;
         for (int i = 0; i < nmat; ++i) 
           delete [] names[i];
         delete [] names;
         retval =(void *)mat;
         df     = avtMaterial::Destruct;
      }
      else if((strncmp(var,"SEelem_3D_vol_grid_Layers",25) == 0) ||
              (strncmp(var,"SEelem_3D_vol_hght_Layers",25) == 0)  )
      {
         // Create material names and matnos array
         // ---------------------------------------
         int   ncells = (SEelem.np-1)*(SEelem.np-1)*SEelem.nelem;
         int     nmat =  SEelem.nlev-1;
         int  *matnos = new int   [nmat];
         char **names = new char *[nmat];
         for(int i = 0; i < nmat; ++i)
         {
            matnos[i] = i + 1;
            char *buffer=new char[20];
            int n;
            n=sprintf(buffer, "Layer %d",matnos[i]); (void) n;
            names[i] = (char *)buffer;
         }

         // Assign mmaterial numbers to mesh
         // ---------------------------------
         int *matlist =new int[nmat*ncells];
//         for (int i = 0; i< ncells; ++i)
//         for (int k = 0; k< nmat  ; ++k)
//         {
//            matlist[i + k*ncells] = matnos[nmat-1-k];
//         }
         int   ind3 = 0;
         for(int ie = 0; ie <  SEelem.nelem  ; ++ie)
         for(int l0 = 0; l0 < (SEelem.nlev-1); ++l0)
         for(int j0 = 0; j0 < (SEelem.np-1)  ; ++j0)
         for(int i0 = 0; i0 < (SEelem.np-1)  ; ++i0)
         {
            matlist[ind3] = matnos[nmat-1-l0];
            ind3 = ind3 + 1;
         }
     

         // Create the Materials
         // ------------------------
         int dims[3] = {1,1,1};
         int   ndims = 1;
         dims[0] = ncells*nmat;
         avtMaterial *mat= new avtMaterial( nmat, matnos, names, ndims, dims,
                                            0, matlist, 0,   // length of mix arrays
                                                        0,   // mix_mat array
                                                        0,   // mix_next array
                                                        0,   // mix_zone array
                                                        0);  // mix_vf array 
         delete [] matlist;
         delete [] matnos;
         for (int i = 0; i < nmat; ++i) delete [] names[i];
         delete [] names;
         retval =(void *)mat;
         df     = avtMaterial::Destruct;
      }
   }

  return retval;
}
