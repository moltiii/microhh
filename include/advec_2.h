/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ADVEC_2
#define ADVEC_2

#include "advec.h"

// forward declaration
class cmodel;

/**
 * Derived class for 2nd order advection scheme.
 */
class cadvec_2 : public cadvec
{
  public:
    cadvec_2(cmodel *); ///< Constructor of the advection class.
    ~cadvec_2();        ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
    double getcfl(double);                               ///< Get the CFL number.
    int exec();                                          ///< Execute the advection scheme.

  private:
    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.

    inline double interp2(const double, const double); ///< 2nd order interpolation function.
};
#endif
