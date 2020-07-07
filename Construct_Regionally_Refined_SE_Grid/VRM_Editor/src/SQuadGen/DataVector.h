///////////////////////////////////////////////////////////////////////////////
///
///	\file    DataVector.h
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

#ifndef _DATAVECTOR_H_
#define _DATAVECTOR_H_

///////////////////////////////////////////////////////////////////////////////

#include "Exception.h"

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

//=============================================================================
template <typename DataType> class DataVector 
  //  A data vector is a datatype that stores data in a 1D structure.
  //  Arithmatic operations are supported for this datatype.
  //
  //  Memory is allocated using malloc and deallocated using free, so no
  //  calls will be made to the constructor or destructor of DataType.  This
  //  class is primarily designed to efficiently handle primitive types.
  //======================================================================
{

    public: //  Constructor.
            //=======================
            DataVector() : m_sRows(0),m_data(NULL)
            { }
            

            //  Constructor.
            //=======================
            DataVector( unsigned int sRows) : m_data(NULL)
            {
              Initialize(sRows);
            }
            

            //  Copy constructor.
            //=======================
            DataVector( const DataVector<DataType> & dv) : m_data(NULL)
            {
              Assign(dv);
            }
            

            //  Destructor.
            //=======================
            virtual ~DataVector() 
            { 
              if(m_data != NULL) { free(reinterpret_cast<void*>(m_data)); }
            }

    public: //  Determine if this DataVector is initialized.
            //==============================================
            bool IsInitialized() const 
            {
              if(m_data == NULL) {
                return false;
              } else {
                return true;
              }
            }


    public: //  Deallocate data for this object.
            //======================================
            void Deinitialize() 
            {
              if(m_data != NULL) {
                free(reinterpret_cast<void*>(m_data));
              }
              m_data  = NULL;
              m_sRows = 0;
            }


            //  Allocate data for this object.
            //=================================
            void Initialize( unsigned int sRows   ,
                             bool fAutoZero = true) 
            {
              // Check for zero size
              //---------------------
              if (sRows == 0) { Deinitialize(); return; }
            
              // No need to reallocate memory if this vector already has
              // the correct dimensions.
              //--------------------------------------------------------
              if(m_sRows == sRows) {
                // Auto zero
                //-----------
                if(fAutoZero) { Zero(); }
                return;
              }
            
              // Deinitialize existing content
              //------------------------------
              Deinitialize();
            
              // Allocate memory
              //-------------------
              m_data = reinterpret_cast<DataType*>( malloc(sRows*sizeof(DataType)));
              if(m_data == NULL) { _EXCEPTIONT("Out of memory."); }
              
              // Assign dimensions
              //------------------
              m_sRows = sRows;
              
              // Auto zero
              //------------
              if(fAutoZero) { Zero(); }
            }


    public: //  Assignment operator.
            //============================
            void Assign(const DataVector<DataType> & dv) 
            {
              // Check initialization status
              //-----------------------------
              if(!dv.IsInitialized()) { Deinitialize(); return; }

              // Allocate memory
              //-----------------
              Initialize(dv.m_sRows);

              // Copy data
              //--------------
              memcpy(m_data, dv.m_data, m_sRows * sizeof(DataType));
            }


            //  Assignment operator.
            //========================
            DataVector & operator= (const DataVector<DataType> & dv) 
            {
              Assign(dv);
              return(*this);
            }


            //  Zero the data content of this object.
            //========================================
            void Zero() 
            {
              // Check initialization status
              //-----------------------------
              if(!IsInitialized()) { 
                _EXCEPTIONT("Attempted operation on uninitialized DataMatrix3D.");
              }

              // Set content to zero
              //---------------------
              memset(m_data, 0, m_sRows*sizeof(DataType));
            }


    public: //  Get the number of rows in this matrix.
            //========================================
            inline unsigned int GetRows() const 
            {
              return m_sRows;
            }

    public: //  Cast to an array.
            //=======================
            inline operator DataType*() 
            {
              return m_data;
            }

            //  Cast to an array.
            //=======================
            inline operator const DataType*() const 
            {
              return m_data;
            }

    private:
           unsigned int m_sRows; // The number of rows in this matrix
           DataType*    m_data;  // A pointer to the data associated with this matrix
};
//=============================================================================

#endif

