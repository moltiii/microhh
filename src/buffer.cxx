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
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "defines.h"
#include "model.h"

cbuffer::cbuffer(cmodel *modelin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  allocated = false;
}

cbuffer::~cbuffer()
{
  if(allocated)
    for(std::map<std::string, double *>::const_iterator it=bufferprofs.begin(); it!=bufferprofs.end(); ++it)
      delete[] it->second;
}

int cbuffer::readinifile(cinput *inputin)
{
  int nerror = 0;

  // optional parameters
  nerror += inputin->getItem(&swbuffer, "buffer", "swbuffer", "", "0");

  if(swbuffer == "1")
  {
    nerror += inputin->getItem(&zstart, "buffer", "zstart", "");
    nerror += inputin->getItem(&sigma , "buffer", "sigma" , "", 2.);
    nerror += inputin->getItem(&beta  , "buffer", "beta"  , "", 2.);
  }

  return nerror;
}

int cbuffer::init()
{
  if(swbuffer == "1")
  {
    // allocate the buffer arrays
    for(fieldmap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
      bufferprofs[it->first] = new double[grid->kcells];

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      bufferprofs[it->first] = new double[grid->kcells];

    allocated = true;
  }

  return 0;
}

int cbuffer::create(cinput *inputin)
{
  int nerror = 0;

  if(swbuffer == "1")
  {
    // set the buffers according to the initial profiles of the variables
    nerror += inputin->getProf(&bufferprofs["u"][grid->kstart], "u", grid->kmax);
    nerror += inputin->getProf(&bufferprofs["v"][grid->kstart], "v", grid->kmax);

    // in case of u and v, subtract the grid velocity
    for(int k=grid->kstart; k<grid->kend; ++k)
    {
      bufferprofs["u"][k] -= grid->utrans;
      bufferprofs["v"][k] -= grid->vtrans;
    }

    // allocate the buffer for w on 0
    for(int k=0; k<grid->kcells; ++k)
      bufferprofs["w"][k] = 0.;
 
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      nerror += inputin->getProf(&bufferprofs[it->first][grid->kstart], it->first, grid->kmax);

    // find the starting points
    bufferkstart  = grid->kstart;
    bufferkstarth = grid->kstart;

    for(int k=grid->kstart; k<grid->kend; ++k)
    {
      // check if the cell center is in the buffer zone
      if(grid->z[k] < this->zstart)
        ++bufferkstart;
      // check if the cell face is in the buffer zone
      if(grid->zh[k] < this->zstart)
        ++bufferkstarth;
    }

    // check whether the lowest of the two levels is contained in the buffer layer
    if(bufferkstarth == grid->kend)
    {
      ++nerror;
      if(master->mpiid == 0) std::printf("ERROR buffer is too close to the model top\n");
    }
  }
  return nerror;
}

int cbuffer::exec()
{
  if(swbuffer == "1")
  {
    // calculate the buffer tendencies
    buffer(fields->mt["u"]->data, fields->mp["u"]->data, bufferprofs["u"], grid->z );
    buffer(fields->mt["v"]->data, fields->mp["v"]->data, bufferprofs["v"], grid->z );
    buffer(fields->mt["w"]->data, fields->mp["w"]->data, bufferprofs["w"], grid->zh);
 
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      buffer(fields->st[it->first]->data, it->second->data, bufferprofs[it->first], grid->z);
  }

  return 0;
}

int cbuffer::buffer(double * const restrict at, const double * const restrict a, 
                    const double * const restrict abuf, const double * const restrict z)
{ 
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double sigmaz;
  double zsizebuf;

  zsizebuf = grid->zsize - this->zstart;

  for(int k=bufferkstart; k<grid->kend; ++k)
  {
    sigmaz = this->sigma*std::pow((z[k]-this->zstart)/zsizebuf, this->beta);
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] -= sigmaz*(a[ijk]-abuf[k]);
      }
  }

  return 0;
}

