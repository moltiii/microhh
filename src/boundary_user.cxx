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


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_user.h"
#include "defines.h"
#include "model.h"

using namespace std;

#define NO_VELOCITY 0.
#define NO_OFFSET 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2

cboundary_user::cboundary_user(cmodel *modelin) : cboundary(modelin)
{
}

int cboundary_user::readinifile(cinput *inputin)
{
  int nerror = 0;

  processbcs(inputin);

  // patch type
  nerror += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , "", 2 );
  nerror += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
  nerror += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
  nerror += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
  nerror += inputin->getItem(&patch_facr, "boundary", "patch_facr", "", 1.);
  nerror += inputin->getItem(&patch_facl, "boundary", "patch_facl", "", 0.);
  
  return nerror;
}

int cboundary_user::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->vtrans);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->vtrans);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, NO_OFFSET, fields->s["tmp1"]->data, patch_facl, patch_facr);
    setbc      (it->second->datatop, it->second->datagradtop, it->second->datafluxtop,
                sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, NO_OFFSET);
  }

  return 0;
}

int cboundary_user::setbc_patch(double * restrict a, double * restrict agrad, double * restrict aflux, int sw, double aval, double visc, double offset,
                                double * restrict tmp, double facl, double facr)
{
  double avall, avalr;
  double xmod, ymod;
  double errvalx, errvaly;

  int ij,jj;
  jj = grid->icells;

  avall = facl*aval;
  avalr = facr*aval;

  // save the pattern
  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij = i + j*jj;
      xmod = fmod(grid->x[i], patch_xh);
      ymod = fmod(grid->y[j], patch_xh);

      errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

      if(patch_dim == 2)
        errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
      else
        errvaly = 1.;

      tmp[ij] = errvalx*errvaly;
    }

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        a[ij] = avall + (avalr-avall)*tmp[ij] - offset;
      }
  }
  else if(sw == BC_NEUMANN)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        agrad[ij] = avall + (avalr-avall)*tmp[ij];
        aflux[ij] = -agrad[ij]*visc;
      }
  }
  else if(sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = avall + (avalr-avall)*tmp[ij];
        agrad[ij] = -aflux[ij]/visc;
      }
  }

    double sum = 0;
    double mean_val;
     for(int k=0;k<grid->icells*grid->jcells;k++)
    {
        sum += agrad[k];
    }

    mean_val = sum/(grid->icells*grid->jcells);
    printf("\n new mean value: %f\n",mean_val);

    ofstream ofile;
    ofile.open("boundary_profile.txt",ios::app); 

    for(int k=0;k<grid->icells;k++)
    {
        ofile << k << " " << agrad[k] << endl;
    }

    ofile.close();

  return 0;
}
