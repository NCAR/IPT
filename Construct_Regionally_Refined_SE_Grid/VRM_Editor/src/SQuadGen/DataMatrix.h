///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataMatrix.h
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

#ifndef _DATAMATRIX_H_
#define _DATAMATRIX_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstring>


//=============================================================================
template <typename DataType> class DataMatrix 
  //  A DataMatrix is a datatype that stores data in a 2D structure.
  //  Arithmatic operations are not supported for this datatype.
  //==================================================================
{

    public: //  Default constructor.
            //=======================
            DataMatrix() : m_sRows(0),m_sColumns(0),m_data(NULL)
            { }


            //  Constructor.
            //    sRows     Number of rows in this matrix.
            //    sColumns  Number of columns in this matrix.
            //=============================================
            DataMatrix(unsigned int sRows,
                       unsigned int sColumns) : m_sRows   (0),
                                                m_sColumns(0),
                                                m_data (NULL)
            {
              Initialize(sRows, sColumns);
            }


            //  Copy constructor.
            //=======================
            DataMatrix(const DataMatrix<DataType> & dm) : m_sRows   (0),
                                                          m_sColumns(0),
                                                          m_data (NULL)
            {
              Assign(dm);
            }


            //  Destructor.
            //=======================
            virtual ~DataMatrix() 
            { 
              if(m_data != NULL) { free(reinterpret_cast<void*>(m_data)); } 
            }


    public: //  Determine if this DataMatrix is initialized.
            //=============================================-
            bool IsInitialized() const 
            {
              if(m_data == NULL) {
                return false;
              } else {
                return true;
              }
            }

    public: //  Deallocate data for this object.  
            //=============================================
            void Deinitialize() 
            {
              if(m_data != NULL) { free(reinterpret_cast<void*>(m_data)); }
              m_data     = NULL;
              m_sRows    = 0;
              m_sColumns = 0;
            }


            //  Allocate memory for this object.
            //=============================================
            void Initialize( unsigned int sRows    ,
                             unsigned int sColumns ,
                             bool fAutoZero = true ) 
            {
              unsigned int sI;
            
              // Check for zero size
              //---------------------
              if((sRows == 0) || (sColumns == 0)) { Deinitialize(); return; }
            
              // No need to reallocate memory if this matrix already has
              // the correct dimensions.
              //-------------------------------------------------------
              if((m_sRows == sRows) && (m_sColumns == sColumns)) {
                // Auto zero
                //------------
                if(fAutoZero) { Zero(); }
                return;
              }
            
              // Deinitialize existing content
              //-------------------------------
              Deinitialize();
            
              // Calculate the per-row footprint
              //---------------------------------
              unsigned int sRowPtrFootprint = sRows   *sizeof(DataType *);
              unsigned int sRowFootprint    = sColumns*sizeof(DataType  );
            
              // Calculate padding to align on DataType boundaries
              //----------------------------------------------------
              unsigned int sPadding;
              if((sRowPtrFootprint % sizeof(DataType)) == 0) {
                sPadding = 0;
              } else {
                sPadding = sizeof(DataType) - sRowPtrFootprint % sizeof(DataType);
              }
            
              // Allocate memory
              //-----------------
              char *rawdata = reinterpret_cast<char*>(malloc( sRowPtrFootprint 
                                                             +sPadding 
                                                             +sRows*sRowFootprint));
              if(rawdata == NULL) { _EXCEPTIONT("Out of memory."); }
            
              // Assign memory pointers
              //------------------------
              char *rawdataRowStart = rawdata + sPadding + sRowPtrFootprint;
            
              m_data = reinterpret_cast<DataType**>(rawdata);
            
              for(sI = 0; sI < sRows; sI++) {
                m_data[sI] = reinterpret_cast<DataType*>(rawdataRowStart + sI*sRowFootprint);
              }
            
              // Assign dimensions
              //-----------------------
              m_sRows    = sRows;
              m_sColumns = sColumns;
            
              // Auto zero
              //------------
              if(fAutoZero) { Zero(); }
            }


    public: //  Assignment operator.  
            //=======================
            void Assign(const DataMatrix<DataType> & dm) 
            {
              // Check initialization status
              //------------------------------
              if(!dm.IsInitialized()) { Deinitialize(); return; }

              // Allocate memory
              //-------------------
              Initialize(dm.m_sRows, dm.m_sColumns);

              // Copy data
              //------------
              memcpy( &(   m_data[0][0]), 
                      &(dm.m_data[0][0]), m_sRows*m_sColumns*sizeof(DataType));
            }


            //  Assignment operator.
            //=======================
            DataMatrix & operator= (const DataMatrix<DataType> & dm) 
            {
              Assign(dm);
              return(*this);
            }


            //  Zero the data content of this object.
            //=============================================
            void Zero() 
            {
              // Check initialization status
              //-----------------------------
              if(!IsInitialized()) { 
                _EXCEPTIONT("Attempted operation on uninitialized DataMatrix3D.");
              }
            
              // Set content to zero
              //----------------------
              memset( &(m_data[0][0]), 0, m_sRows*m_sColumns*sizeof(DataType));
            }


    public: //  Get the number of rows in this matrix.
            //=============================================
            inline unsigned int GetRows() const 
            {
              return m_sRows;
            }


            //  Get the number of columns in this matrix.
            //=============================================
            inline unsigned int GetColumns() const 
            {
              return m_sColumns;
            }


            //  Get the total number of elements in this matrix.
            //=============================================
            inline unsigned int GetTotalElements() const 
            {
              return m_sRows * m_sColumns;
            }


    public: //  Cast to an array.
            //=======================
            inline operator DataType**() const 
            {
              return m_data;
            }


    public: //  Give a string representation of this object.
            //=============================================--
            std::string ToString() const 
            {
              unsigned int i;
              unsigned int j;
            
              std::stringstream strstr;
            
              strstr << "[";
              for(i = 0; i < GetRows(); i++) {
                strstr << "[";
                for(j = 0; j < GetColumns(); j++) {
                  strstr << m_data[i][j];
                  if(j != GetColumns()-1) { strstr << ", "; }
                }
                strstr << "]";
                if(i != GetRows()-1) { strstr << ", "; }
              }
              strstr << "]";
              return (strstr.str());
            }


    private:
            unsigned int m_sRows;    // The number of rows in this matrix.
            unsigned int m_sColumns; // The number of columns in this matrix.
            DataType**   m_data;     // A pointer to the data associated with this matrix
};
//=============================================================================

#endif

