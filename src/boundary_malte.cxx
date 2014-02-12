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

//#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <time.h>
#include <fftw3.h>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_malte.h"
#include "defines.h"
#include "model.h"

using namespace std;

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
  nerror += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
  nerror += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
  nerror += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
  nerror += inputin->getItem(&patch_facr, "boundary", "patch_facr", "", 1.);
  nerror += inputin->getItem(&patch_facl, "boundary", "patch_facl", "", 0.);
  nerror += inputin->getItem(&cut,  "boundary", "cut" , "", 1 );
  nerror += inputin->getItem(&wl,  "boundary", "wl" , "", 0.5 );
  nerror += inputin->getItem(&lambda,  "boundary", "lambda" , "", 1 );
  nerror += inputin->getItem(&sigma,  "boundary", "sigma" , "", 10 );

 

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
  int ij,jj;
  int go, stop;
  double cosX, cosY;
  int nyq = grid->icells/2 + 1;

  double sum = 0;

  double mean_val, std_val;
  double mean_n, std_n;

  noise = new double[grid->icells];
  step = new double[grid->icells*grid->jcells];

  double mean_step, std_step;
  calc_step();
  calc_stats(step, grid->icells*grid->jcells, &mean_step, &std_step);
  // printf("\nmean value step: %lf\tstd step: %lf\n",mean_step, std_step);

  jj = grid->icells;
  go = grid->jgc*grid->icells;
  stop = (grid->jgc+1)*grid->icells-1;

  // printf("\ngo: %d\tstop: %d\n",go,stop);

  // save the pattern
  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij = i + j*jj;
      
      cosX = -cos(2*PI*(grid->x[i]/wl))+1;

      if(patch_dim == 2)
        cosY = -cos(2*PI*(grid->y[j]/wl))+1;
      else
        cosY = 1.;

      tmp[ij] = cosX*cosY;
    }
    // adjust mean to step_mean
    calc_stats(tmp, grid->icells*grid->jcells, &mean_val, &std_val);
    for(int k=0;k<grid->icells*grid->jcells;k++)
    {
        tmp[k] = tmp[k] * (mean_step/mean_val);
    }
    calc_stats(tmp, grid->icells*grid->jcells, &mean_val, &std_val);
    // for latest stretching of combined signal
    double std_cos = std_val;
    printf("\nmean flux value of uncutted cos: %lf\tstd cos: %lf\n",mean_val*aval,std_val*aval);
  
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

    printf("\nmean flux step: %lf\tstd flux step: %lf\n", mean_step*aval, std_step*aval);
    
    calc_stats(aflux, grid->icells*grid->jcells, &mean_val, &std_val);
    printf("\nmean flux cos: %lf\tstd_val: %lf\n", mean_val,std_val);

    // generate noise
    // int static seed = 0;
    srand(time(NULL));
    for(int i=0; i<grid->icells; ++i)
      noise[i] = (double) rand() / (double) RAND_MAX;
    
    // stats of generated noise
    calc_stats(noise, grid->icells, &mean_n, &std_n);
    printf("\ngenerated mean value of noise: %lf\n",mean_n);

    // shifting the noise signal to zero mean value
    for(int k=0;k<grid->icells;k++)
        noise[k] = noise[k] - mean_n;

    // test wether mean value of noise is zero now 
    calc_stats(noise, grid->icells, &mean_n, &std_n);
    printf("\nis noise mean value equal to zero? %lf\n",mean_n);

    // mean value of noise depending on cut
    mean_n = mean_step*aval - mean_val;
    // std of noise depending on cut
    std_n = sqrt(1/pow(cut,2)-1)*std_val;
    printf("\nmean flux noise should be: %lf\tstd_n: %lf\n",mean_n,std_n);

    // generate gaussian filter 
    double * hg;
    hg = new double[grid->icells]; 
    sum = 0;
    for(int i=0;i<nyq;i++)
    {
      hg[i] = (1/(sqrt(PI)*sigma*std_step))*exp(-(pow((lambda/wl-i),2)/(2*pow(sigma*std_step,2)))); // <-------- must be std_stepF!!!!
      sum += hg[i];
    }
    for(int i=0;i<nyq;i++)
      hg[i] = hg[i]/sum;

    // fft of noise
    fftw_complex out[nyq];
    fftw_plan plan_backward;
    fftw_plan plan_forward;
    plan_forward = fftw_plan_dft_r2c_1d ( grid->icells, noise, out, FFTW_ESTIMATE );
    plan_backward = fftw_plan_dft_c2r_1d ( grid->icells, out, noise, FFTW_ESTIMATE );
    fftw_execute ( plan_forward );
    // printf("real part: %f\timag part: %f\n",out[nyq-1][0],out[nyq-1][1]);

    // apply filter
    for (int i=1;i<nyq;i++) // not touching the mean and the nyquist frequency
    {
      out[i][0] = out[i][0] * hg[i];
      out[i][1] = out[i][1] * hg[i];
    }

    // ifft
    fftw_execute ( plan_backward );

    // shifting and stretching noise to correct mean and std
    double mean_temp, std_temp;
    calc_stats(noise, grid->icells, &mean_temp, &std_temp);
    for(int i=0;i<grid->icells;i++)
      noise[i] = noise[i]*(std_n/std_temp) + mean_n;

    calc_stats(noise, grid->icells, &mean_n, &std_n);
    printf("\nmean flux noise really is: %lf\tstd_n final: %lf\n", mean_n, std_n);

    // add noise to cos
    for(int j=0; j<grid->jcells; ++j)
    {
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = aflux[ij] + noise[i];
      }
    }

    // final stretch of combined signal + test
    calc_stats(aflux, grid->icells*grid->jcells, &mean_val, &std_val);
    for(int j=0; j<grid->jcells; ++j)
    {
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = (aflux[ij]-mean_val) * (std_cos*aval/std_val) + mean_val; // <- to avoid manipulating the adjusted mean_value
        agrad[ij] = -aflux[ij]/visc;
      }
    }
    calc_stats(aflux, grid->icells*grid->jcells, &mean_val, &std_val);
    printf("\nmean flux of combined signal: %lf\tstd_val: %lf\n\n", mean_val,std_val);

     // write to file
    ofstream ofile;
    ofile.open("boundary_profile.txt",ios::app); 
    for(int k=0;k<grid->icells;k++)
    {
        ofile << k << " " << aflux[k] << " " << step[k]*aval << endl;
    }
    ofile.close();

    // fftw_destroy_plan ( plan_forward );
    // fftw_destroy_plan ( plan_backward );
    // fftw_free ( out );

  return 0;
}

int cboundary_malte::calc_stats(double * a, int n, double * mean_val, double * std_val)
{ 
  double sum = 0;

  // mean value
  for(int k=0;k<n;k++)
    sum += a[k];  
  *mean_val = sum/(n);

  sum = 0;

  for(int k=0;k<n;k++)
    sum += pow((*mean_val-a[k]),2);
  *std_val = sqrt(sum/n);
  
  return 0;
}

int cboundary_malte::calc_step()
{
  double xmod, ymod;
  double errvalx, errvaly;
  int ij;

  for(int j=0; j<grid->jcells; ++j)
  {
    for(int i=0; i<grid->icells; ++i)
    {
      ij = i + j*(grid->icells);
      xmod = fmod(grid->x[i], patch_xh);
      ymod = fmod(grid->y[j], patch_xh);

      errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

      if(patch_dim == 2)
        errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
      else
        errvaly = 1.;

      step[ij] = errvalx*errvaly;

    }
  }

  return 0;
}



