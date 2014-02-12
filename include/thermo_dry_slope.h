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

#ifndef THERMO_DRY_SLOPE
#define THERMO_DRY_SLOPE

#include "thermo.h"

// forward declarations to speed up build time
class cmaster;
class cgrid;
class cfields;

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
class cthermo_dry_slope : public cthermo
{
  public:
    cthermo_dry_slope(cmodel *); ///< Constructor of the dry thermodynamics class.
    ~cthermo_dry_slope();        ///< Destructor of the dry thermodynamics class.

    int readinifile(cinput *);               ///< Processing data of the input file.
    int exec();                              ///< Add the tendencies belonging to the buoyancy.

    int getbuoyancy(cfield3d *, cfield3d *); ///< Compute the buoyancy for usage in another routine.
    int getbuoyancysurf(cfield3d *);         ///< Compute the near-surface and bottom buoyancy for usage in another routine.
    int getbuoyancyfluxbot(cfield3d *);      ///< Compute the bottom buoyancy flux for usage in another routine.

  private:
    int calcbuoyancy(double *, double *);         ///< Calculation of the buoyancy.
    int calcbuoyancybot(double *, double *,
                        double *, double *);      ///< Calculation of the near-surface and surface buoyancy.
    int calcbuoyancyfluxbot(double *, double *);  ///< Calculation of the buoyancy flux at the bottom.
    int calcbuoyancytendu_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.
    int calcbuoyancytendw_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.
    int calcbuoyancytendb_4th(double *, double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.

    inline double interp4(const double, const double, 
                          const double, const double); ///< 4th order interpolation function.

    double alpha;
    double n2;
};
#endif
