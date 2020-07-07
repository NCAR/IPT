///////////////////////////////////////////////////////////////////////////////
///
///	\file    MathHelper.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		This header file provides access to various mathematical functions
///		that are not normally exposed in the standard C++ libraries.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _MATHHELPER_H_
#define _MATHHELPER_H_

//=============================================================================

template <typename _T> _T Max(_T x1, _T x2) 
  //  Calculate the maximum of two values.
  //=======================================
{
  return (x1>x2)?(x1):(x2);
}


template <typename _T> _T Min(_T x1, _T x2) 
  //  Calculate the minimum of two values.
  //=======================================
{
  return (x1<x2)?(x1):(x2);
}


template <typename _T> _T Sign(_T x1) 
  //  Calculate the sign of a value.
  //==================================
{
  return (x1 < static_cast<_T>(0))?(static_cast<_T>(-1)):(static_cast<_T>(1));
}


template <typename _T> _T Clamp(_T y, _T x1, _T x2) 
  //  Clamp a value to be within a given range.
  //=============================================
{
  return (y>x2)?(x2):((y<x1)?(x1):(y));
}


//  Calculate the integer square root.
//    Source: Crenshaw, Jack.  Integer square roots.
//    http://www.embedded.com/98/9802fe2.htm
//==================================================
inline unsigned int ISqrt(unsigned int a) 
{
  unsigned int irem = 0;
  unsigned int iroot = 0;
  for(int i = 0; i < 16; i++) {
    iroot <<= 1;
    irem = ((irem << 2) + (a >> 30));
    a <<= 2;
    iroot++;
    if(iroot <= irem) {
      irem -= iroot;
      iroot++;
    } else {
      iroot--;
    }
  }
  return (static_cast<unsigned int>(iroot >> 1));
}


//  Calculate the integer power of integers function.
//===================================================
inline int IntPow(int d, unsigned int p) 
{
  if(p == 0) { return 1; }

  unsigned int q;

  int iPow = d;
  for(q = 1; q < p; q++) {
    iPow *= d;
  }
  return iPow;
}


//  Calculate the integer power function.
//=========================================
inline double IPow(double d, unsigned int p) 
{
  unsigned int q;

  double dPow = 1.0;

  for(q = 0; q < p; q++) {
    dPow *= d;
  }
  return dPow;
}


//  Calculate the integer factorial function.
//===========================================
inline unsigned int IFact(unsigned int p) 
{
  unsigned int q;
  unsigned int iFact = 1;

  for (q = 2; q <= p; q++) {
    iFact *= q;
  }
  return iFact;
}


/*
//  Calculate the inverse hyperbolic arcsin of a value.
//------------------------------------------------------
double asinh(double x);

//  Calculate the inverse hyperbolic arccos of a value.
//-------------------------------------------------------
double acosh(double x);

//  Calculate the inverse hyperbolic arctan of a value.
//-----------------------------------------------------
double atanh(double x);
*/

//========================================================================


//  Calculate the linear interpolant between two points (X1, Y1) and (X2, Y2).
//    dXout   - X coordinate of the interpolating point
//    dX1     - First X coordinate
//    dX2		- Second X coordinate
//    dY1     - First Y coordinate
//    dY2     - Second Y coordinate
//===================================================================
double LinearInterpolate( double dXout,
                          double dX1,
                          double dX2,
                          double dY1,
                          double dY2  );


//  Calculate the quadratic interpolant among three equally spaced points.
//    dXout   - Location of the interpolating point relative to X1
//    dDX     - Spacing between points in X
//    dY1     - Value of the function at X1
//    dY2     - Value of the function at X1 + dDX
//    dY3		- Value of the function at X1 + 2 dDX
//===================================================================
double QuadraticInterpolate( double dXout,
                             double dDX,
                             double dY1,
                             double dY2,
                             double dY3  );


//  Calculate the cubic interpolant among four equally spaced points.
//    dXout   - Location of the interpolating point relative to X1
//    dDX     - Spacing between points in X
//    dY1     - Value of the function at X1
//    dY2     - Value of the function at X1 + dDX
//    dY3		- Value of the function at X1 + 2 dDX
//    dY4     - Value of the function at X1 + 3 dDX
//===================================================================
double CubicInterpolate( double dXout,
                         double dDX,
                         double dY1,
                         double dY2,
                         double dY3,
                         double dY4  );

//==================================================================================

#endif
