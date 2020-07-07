///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataMatrix3D.h
///	\author  Paul Ullrich
///	\version July 26, 2010
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

#ifndef _DATAMATRIX3D_H_
#define _DATAMATRIX3D_H_

///////////////////////////////////////////////////////////////////////////////

#include "DataVector.h"
#include "Exception.h"

#include <iostream>
#include <cstdlib>
#include <cstring>



//==================================================================================
template <typename DataType> class DataMatrix3D 
  //  A 3D data matrix is a datatype that stores data in a 3D structure.
  //  Arithmatic operations are not supported for this datatype.
  //=======================================================================
{

    public: // Constructor.
            //========================
            DataMatrix3D() : m_data(NULL) 
            {
              m_sSize[0] = 0;
              m_sSize[1] = 0;
              m_sSize[2] = 0;
            }

            // Constructor.
            //========================
            DataMatrix3D( unsigned int sRows,       
                          unsigned int sColumns,
                          unsigned int sSubColumns) : m_data(NULL)
            {
              m_sSize[0] = 0;
              m_sSize[1] = 0;
              m_sSize[2] = 0;
              Initialize(sRows, sColumns, sSubColumns);
            }

            // Copy constructor.
            //========================
            DataMatrix3D(const DataMatrix3D<DataType> & dm) : m_data(NULL)
            {  
              Assign(dm);
            }

            // Destructor.
            //========================
            virtual ~DataMatrix3D() { 
              if (m_data != NULL) { free(reinterpret_cast<void*>(m_data)); }
            }

    public: // Determine if this DataMatrix3D is initialized.
            //===============================================
            bool IsInitialized() const 
            {      
              if (m_data == NULL) {
                return false;
              } else {
                return true;
              }
            }

    public: // Deallocate data for this object.
            //===============================================
            void Deinitialize() 
            {
              if (m_data != NULL) { free(reinterpret_cast<void*>(m_data)); }
              m_data = NULL;
              m_sSize[0] = 0;
              m_sSize[1] = 0;
              m_sSize[2] = 0;
            }

            // Allocate data for this object.
            //===============================================
            void Initialize( unsigned int sRows,
                             unsigned int sColumns,
                             unsigned int sSubColumns,
                             bool         fAutoZero = true)
            {  
              unsigned int sI;
              unsigned int sJ;

              // Check for zero size
              //-----------------------
              if ((sRows == 0) || (sColumns == 0) || (sSubColumns == 0)) {
                Deinitialize();
                return;
              }
 
              // No need to reallocate memory if this matrix already has
              // the correct dimensions.
              //---------------------------------------------------------
              if ((m_sSize[0] == sRows) && (m_sSize[1] == sColumns   ) 
                                        && (m_sSize[2] == sSubColumns)) {
                // Auto zero
                //-----------
                if (fAutoZero) { Zero(); }
                return;
              }
 
              // Deinitialize existing content
              //--------------------------------
              Deinitialize();
 
              // Calculate the per-row footprint
              //----------------------------------
              unsigned int sRowPtrFootprint    = sRows      *sizeof(DataType **);
              unsigned int sColumnPtrFootprint = sColumns   *sizeof(DataType * );
              unsigned int sColumnFootprint    = sSubColumns*sizeof(DataType   );
 
              // Allocate memory
              //------------------
              char *rawdata = reinterpret_cast<char*>( malloc(sRowPtrFootprint
                                                              +sRows*sColumnPtrFootprint 
                                                              +sRows*sColumns*sColumnFootprint));
              if (rawdata == NULL) { _EXCEPTIONT("Out of memory."); }
 
              // Assign memory pointers
              //--------------------------
              char *pColumnPtrs = rawdata     + sRowPtrFootprint;
              char *pDataStart  = pColumnPtrs + sRows*sColumnPtrFootprint;
               
              m_data = reinterpret_cast<DataType***>(rawdata);
               
              for (sI = 0; sI < sRows; sI++) {
                 m_data[sI] = reinterpret_cast<DataType**>(pColumnPtrs + sI*sColumnPtrFootprint);
               
                for (sJ = 0; sJ < sColumns; sJ++) {
                  m_data[sI][sJ] = reinterpret_cast<DataType*>( pDataStart 
                                                               +(sI*sColumns+sJ)*sColumnFootprint);
                }
              }
               
              // Assign dimensions
              //-----------------------
              m_sSize[0] = sRows;
              m_sSize[1] = sColumns;
              m_sSize[2] = sSubColumns;
               
              // Auto zero
              //-----------------
              if (fAutoZero) { Zero(); }
            }

	public: // Assignment operator.
                //========================
                void Assign(const DataMatrix3D<DataType> & dm) 
                {                    
                  // Check initialization status
                  //-----------------------------
                  if (!dm.IsInitialized()) { 
                    Deinitialize(); 
                    return;
                  }
                  
                  // Allocate memory
                  //------------------
                  Initialize(dm.m_sSize[0], dm.m_sSize[1], dm.m_sSize[2]);
                  
                  // Copy data
                  //--------------
                  unsigned int sRowPtrFootprint = m_sSize[0]*sizeof(DataType **);
                  unsigned int sColPtrFootprint = m_sSize[1]*sizeof(DataType * );
                  
                  memcpy(reinterpret_cast<char*>(m_data   )+sRowPtrFootprint+m_sSize[0]*sColPtrFootprint,
                         reinterpret_cast<char*>(dm.m_data)+sRowPtrFootprint+m_sSize[0]*sColPtrFootprint,
                         m_sSize[0]*m_sSize[1]*m_sSize[2]*sizeof(DataType));
                }

                // Assignment operator.
                //========================
                DataMatrix3D & operator= (const DataMatrix3D<DataType> & dm) 
                {
                  Assign(dm);
                  return (*this);
                }

                // Zero the data content of this object.
                //===============================================
                void Zero() 
                {
                  // Check initialization status
                  //-----------------------------
                  if (!IsInitialized()) { 
                    _EXCEPTIONT( "Attempted operation on uninitialized DataMatrix3D."); 
                  }

                  // Set content to zero
                  //---------------------
                  memset( &(m_data[0][0][0]),0,m_sSize[0]*m_sSize[1]*m_sSize[2]*sizeof(DataType));
                }

    public: // Copy the data of this object into a DataVector.
            //===============================================---
            void Vectorize(DataVector<DataType> & vec) const 
            {
              // Check initialization status
              //------------------------------
              if (!IsInitialized()) 
              {
                _EXCEPTIONT("Attempted operation on uninitialized DataMatrix3D.");
              }
              
              // Allocate memory
              //------------------
              vec.Initialize(m_sSize[0] * m_sSize[1] * m_sSize[2]);
              
              // Copy data
              //------------
              memcpy(reinterpret_cast<char*>(&(vec[0])),
                     reinterpret_cast<char*>(&(m_data[0][0][0])),
                     m_sSize[0]*m_sSize[1]*m_sSize[2]*sizeof(DataType));
            }

    public: // Get the number of rows in this matrix.
            //===============================================
            inline unsigned int GetRows() const 
            { 
              return m_sSize[0]; 
            }

            // Get the number of columns in this matrix.
            //===============================================
            inline unsigned int GetColumns() const 
            { 
              return m_sSize[1]; 
            }

            // Get the number of columns in this matrix.
            //===============================================
            inline unsigned int GetSubColumns() const 
            { 
              return m_sSize[2]; 
            }

            // Get the number of rows in this matrix.
            //===============================================
            inline unsigned int GetSize(int dim) const 
            { 
              return m_sSize[dim]; 
            }

            // Get the total number of elements in this matrix.
            //===============================================
            inline unsigned int GetTotalElements() const 
            { 
              return m_sSize[0]*m_sSize[1]*m_sSize[2]; 
            }

    public: // Cast to an array.
            //=======================
            inline operator DataType***() const { return m_data; }

    private: 
           unsigned int m_sSize[3]; // The number of elements in each dimension for matrix
           DataType*** m_data;      // A pointer to the data associated with this matrix
};
//==================================================================================

#endif

