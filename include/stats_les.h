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

#ifndef STATS_LES
#define STATS_LES

#include <netcdfcpp.h>
#include <stats.h>

// forward declarations to reduce compilation time
class cmodel;

class cstats_les : public cstats
{
  public:
    cstats_les(cmodel *);
    ~cstats_les();

    int readinifile(cinput *);
    int init(double);
    int create(int);
    unsigned long gettimelim(unsigned long);
    int exec(int, double, unsigned long);

  private:
    bool allocated;
    bool initialized;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;

    double *umodel, *vmodel;

    profmap profs;
    int addprof(std::string, std::string);

    int calcmean     (double *, double *, double);
    int calcmoment   (double *, double *, double *, double, int);
    int calcdiff     (double *, double *, double *, double *, double *, double *, double);
    int calcgrad     (double *, double *, double *);
    int calcflux     (double *, double *, double *, double *, int, int);
    int addfluxes    (double *, double *, double *);
    int calccount    (double* data, double* prof, double threshold);

    int nstats;
};
#endif
