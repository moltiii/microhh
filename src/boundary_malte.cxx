
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

#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_malte.h"
#include "defines.h"
#include "model.h"

#define NO_VELOCITY 0.
#define NO_OFFSET 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2

#define PI 3.14159265358979323846

cboundary_malte::cboundary_malte(cmodel *modelin) : cboundary(modelin)
{
}

int cboundary_malte::readinifile(cinput *inputin)
{
  int nerror = 0;

  nerror += processbcs(inputin);

  // patch type
  nerror += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , "", 2 );
  nerror += inputin->getItem(&cut,  "boundary", "cut" , "", 1 );
  nerror += inputin->getItem(&wl,  "boundary", "wl" , "", 0.5 );
 

  return nerror;
}

int cboundary_malte::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->vtrans);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->vtrans);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, NO_OFFSET, fields->s["tmp1"]->data);
    setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, NO_OFFSET);
  }

  return 0;
}

int cboundary_malte::setbc_patch(double * restrict a, double * restrict agrad, double * restrict aflux, int sw, double aval, double visc, double offset,
                                double * restrict tmp)
{

  double cosX, cosY;

  int ij,jj;
  jj = grid->icells;

  // save the pattern
  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij = i + j*jj;
      
      cosX = sin(2*PI*(grid->x[i]/wl))+1;

      if(patch_dim == 2)
      {
        cosY = sin(2*PI*(grid->y[j]/wl))+1;
        tmp[ij] = cosX*cosY;
        // adjust_mean(tmp);
      }
      else
      {
        cosY = 1.;
        tmp[ij] = cosX*cosY;
      }
    }
    adjust_mean(tmp);

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        a[ij] = cut*tmp[ij] - offset;
      }
  }
  else if(sw == BC_NEUMANN)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        agrad[ij] = cut*tmp[ij]*aval;
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
        aflux[ij] = cut*tmp[ij]*aval;
        agrad[ij] = -aflux[ij]/visc;
      }
  }

  return 0;
}

int cboundary_malte::adjust_mean(double * a)

{ 
  int size;
  size = sizeof(a) / sizeof(a[0]);

  double sum = 0.0;
  double mean_f = 0.0; 

  for (int i = 0; i < size; ++i)
  {
     sum += a[i];
  }

  mean_f = sum/size;

  printf("mean value = %f\n", mean_f );

  // for (int i = 0; i < size; ++i)
  // {
  //    a[i] = a[i] * (mean_val/mean_f);
  // }

  return 0;
}
