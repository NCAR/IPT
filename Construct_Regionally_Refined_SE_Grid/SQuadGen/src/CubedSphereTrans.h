///////////////////////////////////////////////////////////////////////////////
///
///	\file    CubedSphereTrans.h
///	\author  Paul Ullrich
///	\version August 11, 2010
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

#ifndef _CUBEDSPHERETRANS_H_
#define _CUBEDSPHERETRANS_H_

#include "MathHelper.h"

#include <cmath>
#include <cfloat>

////////////////////////////////////////////////////////////////////////////////

class CubedSphereTrans {

private:
	///	<summary>
	///		No initialization for objects of this class.
	///	</summary>
	CubedSphereTrans()
	{ }

public:

	///	<summary>
	///		Determine the (lat, lon) coordinates of a point in equiangular ABP 
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		dA        - Equiangular alpha coordinate
	///		dB        - Equiangular beta coordinate
	///		nP        - Panel coordinate (0-5, 4 = north, 5 = south)
	///		lat (OUT) - Latitude
	///		lon (OUT) - Longitude
	///	</parameters>
	static void RLLFromABP(
		double dA,
		double dB,
		int nP,
		double & lon,
		double & lat
	) {
		double dX = tan(dA);
		double dY = tan(dB);

		switch (nP) {
			// Equatorial panel 1
			case 0:
				lon = atan(dX);
				lat = atan(dY / sqrt(1.0 + dX * dX));
				break;

			// Equatorial panel 2
			case 1:
				lon = atan(dX) + 0.5 * M_PI;
				lat = atan(dY / sqrt(1.0 + dX * dX));
				break;

			// Equatorial panel 3
			case 2:
				lon = atan(dX) + M_PI;
				lat = atan(dY / sqrt(1.0 + dX * dX));
				break;

			// Equatorial panel 4
			case 3:
				lon = atan(dX) + 1.5 * M_PI;
				lat = atan(dY / sqrt(1.0 + dX * dX));
				break;

			// North polar panel
			case 4:
				if (fabs(dX) > DBL_EPSILON) {
					lon = atan2(dX, -dY);
				} else if (dY <= 0.0) {
					lon = 0.0;
				} else {
					lon = M_PI;
				}
				lat = 0.5 * M_PI - atan(sqrt(dX * dX + dY * dY));
				break;

			// South polar panel
			case 5:
				if (fabs(dX) > DBL_EPSILON) {
					lon = atan2(dX, dY);
				} else if (dY > 0.0) {
					lon = 0.0;
				} else {
					lon = M_PI;
				}
				lat = -0.5 * M_PI + atan(sqrt(dX * dX + dY * dY));
				break;

			// Invalid panel
			default:
				_EXCEPTION1("Invalid nP coordinate.  "
					"Given: %d, Expected: [0-5].\n", nP);
		}
	}

	///	<summary>
	///		Determine the equiangular ABP coordinates of a point in (lat, lon)
	///		coordinates.
	///	</summary>
	///	<parameters>
	///		lat      - Latitude
	///		lon      - Longitude
	///		dA (OUT) - Equiangular alpha coordinate
	///		dB (OUT) - Equiangular beta coordinate
	///		nP (OUT) - Panel coordinate (0-5, 4 = north, 5 = south)
	///	</parameters>
	static void ABPFromRLL(
		double lon,
		double lat,
		double & dA,
		double & dB,
		int & nP
	) {
		// Default panel to unattainable value
		nP = 6;

		// Translate from RLL coordinates to XYZ space
		double xx, yy, zz, pm;

		xx = cos(lon) * cos(lat);
		yy = sin(lon) * cos(lat);
		zz = sin(lat);

		pm = Max(fabs(xx), Max(fabs(yy), fabs(zz)));

		// Check maxmality of the x coordinate
		if (pm == fabs(xx)) {
			if (xx > 0) {
				nP = 0;
			} else {
				nP = 2;
			}
		}

		// Check maximality of the y coordinate
		if (pm == fabs(yy)) {
			if (yy > 0) {
				nP = 1;
			} else {
				nP = 3;
			}
		}

		// Check maximality of the z coordinate
		if (pm == fabs(zz)) {
			if (zz > 0) {
				nP = 4;
			} else {
				nP = 5;
			}
		}

		// Panel assignments
		double sx, sy, sz;
		if (nP == 0) {
			sx = yy;
			sy = zz;
			sz = xx;

		} else if (nP == 1) {
			sx = -xx;
			sy = zz;
			sz = yy;

		} else if (nP == 2) {
			sx = -yy;
			sy = zz;
			sz = -xx;

		} else if (nP == 3) {
			sx = xx;
			sy = zz;
			sz = -yy;

		} else if (nP == 4) {
			sx = yy;
			sy = -xx;
			sz = zz;

		} else if (nP == 5) {
			sx = yy;
			sy = xx;
			sz = -zz;

		} else {
			_EXCEPTIONT("Logic error.");
		}

		// Convert to gnomonic coordinates
		dA = atan(sx / sz);
		dB = atan(sy / sz);
	}

	///	<summary>
	///		Determine the panel id in the given direction relative to
	///		a given panel id.
	///	</summary>
	///	<parameters>
	///		Nc              - Resolution of the cubed sphere grid
	///		p_src           - Source panel
	///		ix_src          - Index in X direction of source element
	///		jx_src          - Index in Y direction of source element
	///		p_dest (OUT)    - Destination panel
	///		ix_dest (OUT)   - Index in the X direction of destination element
	///		jx_dest (OUT)   - Index in the Y direction of destination element
	///	</parameters>
	static bool RelativeCoord(
		int Nc,
		int p_src,
		int ix_src,
		int jx_src,
		int &p_dest,
		int &ix_dest,
		int &jx_dest
	) {
		// Internal points
		if ((jx_src >= 0) && (ix_src >= 0) && (jx_src < Nc) && (ix_src < Nc)) {
			jx_dest = jx_src;
			ix_dest = ix_src;
			p_dest = p_src;

		// Equatorial panel 0
		} else if (p_src == 0) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = jx_src - Nc;
				ix_dest = ix_src;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = jx_src + Nc;
				ix_dest = ix_src;
			}

		// Equatorial panel 1
		} else if (p_src == 1) {
			if (ix_src >= Nc) {
				p_dest = 2;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = ix_src;
				ix_dest = 2 * Nc - 1 - jx_src;

			} else if (ix_src < 0) {
				p_dest = 0;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = Nc - 1 - ix_src;
				ix_dest = Nc + jx_src;
			}

		// Equatorial panel 2
		} else if (p_src == 2) {
			if (ix_src >= Nc) {
				p_dest = 3;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = 2 * Nc - 1 - jx_src;
				ix_dest = Nc - ix_src - 1;

			} else if (ix_src < 0) {
				p_dest = 1;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = - jx_src - 1;
				ix_dest = Nc - ix_src - 1;
			}

		// Equatorial panel 3
		} else if (p_src == 3) {
			if (ix_src >= Nc) {
				p_dest = 0;
				jx_dest = jx_src;
				ix_dest = ix_src - Nc;

			} else if (jx_src >= Nc) {
				p_dest = 4;
				jx_dest = Nc - ix_src - 1;
				ix_dest = jx_src - Nc;

			} else if (ix_src < 0) {
				p_dest = 2;
				jx_dest = jx_src;
				ix_dest = ix_src + Nc;

			} else if (jx_src < 0) {
				p_dest = 5;
				jx_dest = ix_src;
				ix_dest = - jx_src - 1;
			}

		// North polar panel
		} else if (p_src == 4) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = 2 * Nc - 1 - ix_src;
				ix_dest = jx_src;

			} else if (jx_src >= Nc) {
				p_dest = 2;
				jx_dest = 2 * Nc - 1 - jx_src;
				ix_dest = Nc - 1 - ix_src;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = ix_src + Nc;
				ix_dest = Nc - 1 - jx_src;

			} else if (jx_src < 0) {
				p_dest = 0;
				jx_dest = jx_src + Nc;
				ix_dest = ix_src;
			}

		// South polar panel
		} else if (p_src == 5) {
			if (ix_src >= Nc) {
				p_dest = 1;
				jx_dest = ix_src - Nc;
				ix_dest = Nc - 1 - jx_src;

			} else if (jx_src >= Nc) {
				p_dest = 0;
				jx_dest = jx_src - Nc;
				ix_dest = ix_src;

			} else if (ix_src < 0) {
				p_dest = 3;
				jx_dest = - ix_src - 1;
				ix_dest = jx_src;

			} else if (jx_src < 0) {
				p_dest = 2;
				jx_dest = - jx_src - 1;
				ix_dest = Nc - 1 - ix_src;
			}
		}

		// Check bounds to see if destination element is valid
		if ((ix_dest < 0) || (ix_dest >= Nc) ||
			(jx_dest < 0) || (jx_dest >= Nc)
		) {
			return false;
		}

		return true;
	}
};

////////////////////////////////////////////////////////////////////////////////

#endif

