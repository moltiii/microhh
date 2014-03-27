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
#include <sstream>
#include <string>
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
#include "master.h"

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
                                double * restrict tmp1)
{
  int ij, jj, iref, jref, igjg;
  double cosX, cosY;
  int tot = grid->itot*grid->jtot;
  double sum = 0;
  double mean_val, var_val;
  double mean_n, var_n;
  double mean_step, var_step;

  noise = new double[tot];

  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " imax: " << grid->imax << " jmax: " << grid->jmax << endl;

  // ---------------------generate step function and cos in physical space----------------------

  // getting the statistics of the step function
  calc_step(tmp1, aval);
  calc_stats(tmp1, tot, &mean_step, &var_step);
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " mean step: " << mean_step << " var step: " << var_step << endl;

  // save the cos pattern for the whole field 
  for(int j=0; j<grid->jtot; j++)
  {
    for(int i=0; i<grid->itot; i++)
    { 
      ij = i + j*grid->itot;
      if(patch_dim == 2)
      {
        cosX = -cos(2*PI*((grid->xsize/grid->itot)*i/wl+(grid->xsize/grid->itot/2)/wl))+1;
        cosY = -cos(2*PI*((grid->ysize/grid->jtot)*j/wl+(grid->ysize/grid->jtot/2)/wl))+1;
      }
      else
      {
        cosX = -cos(2*PI*((grid->xsize/grid->itot)*i/wl));
        cosY = 1.;
      }

      noise[ij] = cosX*cosY;
    }
  }
  
  calc_stats(noise, tot, &mean_val, &var_val);
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " mean of original cos: " << mean_val << " variance of original cos: " << var_val << endl;

  // // adjusting amplitude of cos to mean_step and shifting to zero mean
  // for(int k=0; k<tot; k++)
  // { 
  //   if (patch_dim == 2)
  //     noise[k] = (noise[k]) - mean_val;
  //   else
  //     noise[k] = (noise[k]) - mean_val;
  // }

  // adjusting amplitude to max B of step
  for(int k=0; k<tot; k++)
    noise[k] = (noise[k]) * ((4/3)*(aval-mean_step)/4);
  calc_stats(noise, tot, &mean_val, &var_val);
  for(int k=0; k<tot; k++)
    noise[k] = (noise[k]) -mean_val;

  calc_stats(noise, tot, &mean_val, &var_val);
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " mean of adjusted local cos: " << mean_val << " variance of adjusted local cos: " << var_val << endl;

//  // write the two functions to file
//  if (master->mpiid == 0)
//  {
//  ofstream nfile;
//  nfile.open("surf_cos.txt",ios::app);
//  for (int j=0;j<grid->jtot;j++)
//  {
//    for (int i=0; i<grid->itot;i++)
//    {
//      ij = i + j*grid->itot;
//      nfile << noise[ij] << " ";
//    }
//    nfile << endl;
//  }
//  nfile.close();
//  }
  // write to file
  if (master->mpiid == 0)
  {
  ofstream kfile;
  kfile.open("step.txt",ios::app);
  for (int j=0;j<grid->jtot;j++)
  {
    for (int i=0; i<grid->itot;i++)
    {
      ij = i + j*grid->itot;
      kfile << tmp1[ij]*aval << " ";
    }
    kfile << endl;
  }
  kfile.close();
  }
  // --------------------------------------------------------------------------------------------------------



  // --------------------------------------------------fft for 3D -------------------------------------------
  if (patch_dim == 2)
  {
  // fft of noise
  int nyh = grid->jtot/2 +1;
  double sqrt2 = sqrt(2.0);
  fftw_complex out[grid->itot * nyh];
  fftw_plan plan_forward = fftw_plan_dft_r2c_2d ( grid->itot, grid->jtot, noise, out, FFTW_ESTIMATE );
  fftw_plan plan_backward = fftw_plan_dft_c2r_2d ( grid->itot, grid->jtot, out, noise, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );

  // initialize temp field and scaling coefficents with grid size
  double abs_val;
  fftw_complex temp[grid->itot * nyh];
  sum = 0;
  for (int k=0;k<nyh*grid->itot;k++)
    { 
      temp[k][0] = 0;
      temp[k][1] = 0;
      out[k][0] = out[k][0]/tot;
      out[k][1] = out[k][1]/tot;
    }

  sum = 0;
  for (int k=0;k<nyh*grid->itot;k++)
  { 
    sum += (pow(out[k][0],2) + pow(out[k][1],2));
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum before scaling: " << sum << endl;

  // computing the sum to be cutted
  double cut_sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        cut_sum += 2*(pow(out[ij][0],2) + pow(out[ij][1],2));

      else
        cut_sum += (pow(out[ij][0],2) + pow(out[ij][1],2));  
    }
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " fft coefficents after scaling summed up: " << cut_sum << endl;


  // cutting the spectrum and storing the cutted variance
  cut_sum = (1-cut)*cut_sum*0.5;
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " cut_sum: " << cut_sum << endl;

  for (int k=0;k<nyh*grid->itot;k++)
  {
    out[k][0] = sqrt(cut) * out[k][0];
    out[k][1] = sqrt(cut) * out[k][1];
  }

  // dominant wavenmuber and noise wavenumber
  int dom_wav = (int) grid->xsize/wl;
  int n_wav = lambda * dom_wav;
  if (master->mpiid == 0)
    cout << "dominant wavenumber: " << dom_wav << " noise wavenumber: " << n_wav << endl;

  // around which wavenumber the cutted variance should be located, sure not in wavenumber 0
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    {
      if(j != 0)
      {
      ij = i*nyh + j;
      temp[ij][0] = (1/(2*PI*pow(sigma,2))*exp(-(abs((n_wav*n_wav)-(pow(i,2) + pow(j,2)))/(2*pow(sigma,2)))));
      }
    }
  }

      // control
  sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        sum += 2*(pow(out[ij][0],2) + pow(out[ij][1],2));

      else
        sum += (pow(out[ij][0],2) + pow(out[ij][1],2));  
    }
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum remained in out before 2 ring: " << sum << endl;

  // second ring
   for (int i=1;i<grid->itot/2;i++)
  {
    for (int j=0;j<nyh;j++)
    {
      ij = i*nyh + j;
      temp[(grid->itot-i)*nyh+j][0] = temp[ij][0];
    }
  }
  

    // write to file
  if (master->mpiid == 0)
  {
    ofstream qfile;
    qfile.open("spectrum.txt",ios::app); 
    for (int i=0;i<grid->itot;i++)
    {
      for (int j=0;j<nyh;j++)
      {
        ij = i*nyh + j;
        qfile << pow(temp[ij][0],2) + pow(temp[ij][1],2) << " ";
      }
      qfile << endl;
    }
    qfile.close();
  }   

  // adjusting variance of noise
  sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        sum += 2*(pow(temp[ij][0],2) + pow(temp[ij][1],2));

      else
        sum += (pow(temp[ij][0],2) + pow(temp[ij][1],2));  
    }
  }

  for (int k=0;k<nyh*grid->itot;k++) 
  { 
    temp[k][0] = temp[k][0] * sqrt(2*cut_sum/sum);
  }

  sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        sum += 2*(pow(temp[ij][0],2) + pow(temp[ij][1],2));

      else
        sum += (pow(temp[ij][0],2) + pow(temp[ij][1],2));  
    }
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum before randomizing: " << sum << endl;

  
  // randomizing phases
  srand(2);
  double x,y1,y2;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if ((pow(out[ij][0],2) + pow(out[ij][1],2)) < 0.00001)
      {
      abs_val = (pow(temp[ij][0],2) + pow(temp[ij][1],2));
      x = (double) rand() / (double) RAND_MAX;
      temp[ij][0] = sqrt(x*abs_val);
      temp[ij][1] = sqrt((1-x)*abs_val);

      y1 = (double) rand() / (double) RAND_MAX - 0.5;
      temp[ij][0] = temp[ij][0] * (y1/abs(y1));
      y2 = (double) rand() / (double) RAND_MAX - 0.5;
      temp[ij][1] = temp[ij][1] * (y2/abs(y2));
      }
      else
        if (master->mpiid == 0)
          cout << "i: " << i << " j: " << j << " temp[ij]: " << temp[ij][0] << endl;
    }
  }

  // control
  sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        sum += 2*(pow(out[ij][0],2) + pow(out[ij][1],2));

      else
        sum += (pow(out[ij][0],2) + pow(out[ij][1],2));  
    }
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum remained in out: " << sum << endl;


  // adding the origianl variance of dominant wavenumber again keeping the sign
  for (int i=0;i<grid->itot;i++)
    {
      for (int j=0;j<nyh;j++)
      {
        ij = i*nyh + j;
        if ((pow(out[ij][0],2) + pow(out[ij][1],2)) > 0.00001)
        {
          out[ij][0] = sqrt(pow(out[ij][0],2) + pow(temp[ij][0],2))*(out[ij][0]/abs(out[ij][0]));
          out[ij][1] = sqrt(pow(out[ij][1],2) + pow(temp[ij][1],2))*(out[ij][1]/abs(out[ij][1]));
        }
        else
        {
          out[ij][0] = out[ij][0] + temp[ij][0];
          out[ij][1] = out[ij][1] + temp[ij][1];
        }
      }
    }

  // control
  sum = 0;
  for (int i=0;i<grid->itot;i++)
  {
    for (int j=0;j<nyh;j++)
    { 
      ij = i*nyh + j;
      if (j != 0 && j != nyh-1)
        sum += 2*(pow(out[ij][0],2) + pow(out[ij][1],2));

      else
        sum += (pow(out[ij][0],2) + pow(out[ij][1],2));  
    }
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum after adding: " << sum << endl;
  
  // ifft
  fftw_execute ( plan_backward );
  }
  // -----------------------------------------------------------------------------------------------------



  // ------------------------------fft for 2D computation ------------------------------------------------
  else
  {
  // fft of noise
  double abs_val;
  int nyq = grid->itot/2+1;
  fftw_complex temp[nyq];
  fftw_complex out[nyq];

  fftw_plan plan_forward = fftw_plan_dft_r2c_1d ( tot, noise, out, FFTW_ESTIMATE );
  fftw_plan plan_backward = fftw_plan_dft_c2r_1d (  tot, out, noise, FFTW_ESTIMATE );
  fftw_execute ( plan_forward );
  
  sum = 0;
  for (int k=0;k<nyq;k++)
    { 
      temp[k][0] = 0;
      temp[k][1] = 0;
      out[k][0] = out[k][0]/tot;
      out[k][1] = out[k][1]/tot;
    }

  sum = 0;
  for (int k=0;k<nyq;k++)
  { 
    sum += (pow(out[k][0],2) + pow(out[k][1],2));
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum before scaling: " << sum << endl; 
 
   // computing the sum to be cutted
  double cut_sum = 0;
  for (int k=0;k<nyq;k++)
  {
      if (k != 0 && k != nyq-1)
        cut_sum += 2*(pow(out[k][0],2) + pow(out[k][1],2));

      else
        cut_sum += (pow(out[k][0],2) + pow(out[k][1],2));  
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " fft coefficents after scaling summed up: " << cut_sum << endl;


  // cutting the spectrum and storing the cutted variance
  cut_sum = (1-cut)*cut_sum*0.5;
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " cut_sum: " << cut_sum << endl;

  for (int k=0;k<nyq;k++)
  {
    out[k][0] = sqrt(cut) * out[k][0];
    out[k][1] = sqrt(cut) * out[k][1];
  }

  // dominant wavenmuber and noise wavenumber
  int dom_wav = (int) grid->xsize/wl;
  int n_wav = lambda * dom_wav;
  if (master->mpiid == 0)
    cout << "dominant wavenumber: " << dom_wav << " noise wavenumber: " << n_wav << endl;


  // around which wavenumber the cutted variance should be located, sure not in wavenumber 0
  for (int k=0;k<nyq;k++)
  {
    if(k != 0)
    {
      temp[k][0] = (1/(sqrt(2*PI)*sigma))*exp(-(abs(pow(n_wav-k,2))/(2*pow(sigma,2))));
    }
  }

  // control
  sum = 0;
  for (int k=0;k<nyq;k++)
  {
      if (k != 0 && k != nyq-1)
        sum += 2*(pow(out[k][0],2) + pow(out[k][1],2));

      else
        sum += (pow(out[k][0],2) + pow(out[k][1],2));  
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " fft coefficents remaining in out: " << sum << endl;


    // write to file
  if (master->mpiid == 0)
  {
    ofstream gfile;
    gfile.open("spectrum.txt",ios::app); 
    for (int k=0;k<nyq;k++)
    {
      gfile << pow(temp[k][0],2) + pow(temp[k][1],2) << " ";
    }
    gfile.close();
  }   

  // adjusting variance of noise
  sum = 0;
  for (int k=0;k<nyq;k++)
  {
      if (k != 0 && k != nyq-1)
        sum += 2*(pow(temp[k][0],2) + pow(temp[k][1],2));

      else
        sum += (pow(temp[k][0],2) + pow(temp[k][1],2));  
  }

  for (int k=0;k<nyq;k++) 
  { 
    temp[k][0] = temp[k][0] * sqrt(2*cut_sum/sum);
  }

   // randomizing phases
  srand(2);
  double x,y1,y2;
  for (int k=0;k<nyq;k++)
  {   
      if ((pow(out[k][0],2) + pow(out[k][1],2)) < 0.00001)
      {
      abs_val = (pow(temp[k][0],2) + pow(temp[k][1],2));
      x = (double) rand() / (double) RAND_MAX;
      temp[k][0] = sqrt(x*abs_val);
      temp[k][1] = sqrt((1-x)*abs_val);

      y1 = (double) rand() / (double) RAND_MAX - 0.5;
      temp[k][0] = temp[k][0] * (y1/abs(y1));
      y2 = (double) rand() / (double) RAND_MAX - 0.5;
      temp[k][1] = temp[k][1] * (y2/abs(y2));
      }
      else
        if (master->mpiid == 0)
          cout << "k: " << k <<" temp[k]: " << temp[k][0] << endl;
  }

  sum = 0;
  for (int k=0;k<nyq;k++)
  {
      if (k != 0 && k != nyq-1)
        sum += 2*(pow(temp[k][0],2) + pow(temp[k][1],2));

      else
        sum += (pow(temp[k][0],2) + pow(temp[k][1],2));  
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum after randomizing: " << sum << endl;

  // adding the origianl variance of dominant wavenumber again keeping the sign
  for (int k=0;k<nyq;k++)
    {
        if ((pow(out[k][0],2) + pow(out[k][1],2)) > 0.00001)
        {
          out[k][0] = sqrt(pow(out[k][0],2) + pow(temp[k][0],2))*(out[k][0]/abs(out[k][0]));
          out[k][1] = sqrt(pow(out[k][1],2) + pow(temp[k][1],2))*(out[k][1]/abs(out[k][1]));
        }
        else
        {
          out[k][0] = out[k][0] + temp[k][0];
          out[k][1] = out[k][1] + temp[k][1];
        }
    }

 sum = 0;
  for (int k=0;k<nyq;k++)
  {
      if (k != 0 && k != nyq-1)
        sum += 2*(pow(out[k][0],2) + pow(out[k][1],2));

      else
        sum += (pow(out[k][0],2) + pow(out[k][1],2));  
  }
  if (master->mpiid == 0)
    cout << "id: " << master->mpiid << " sum before ifft: " << sum << endl;

  // ifft
  fftw_execute ( plan_backward );
  }
  // -----------------------------------------------------------------------------------------------------



  // -----------------------------------------------writing the data in array ---------------------------- 

  // adding mean_step
  for(int k=0; k<tot; k++)
    noise[k] = (noise[k] + mean_step);

   calc_stats(noise, tot, &mean_val, &var_val);
  if (master->mpiid == 0)
  cout << "id: " << master->mpiid << " mean after adding step mean: " << mean_val << " var: " << var_val << endl << endl;
  
  // // noise test
  // int z = 1;
  // for (int j=0; j<grid->jtot; j++)
  // {
  //   for (int i=0; i<grid->itot; i++)
  //   {
  //     ij = i + j*grid->itot;
  //     noise[ij] = z;
  //     z = z+1;
  //   }
  // }

      // write to file
  if (master->mpiid == 0)
  {
  ofstream mfile;
  mfile.open("noise.txt",ios::app);
  for (int j=0;j<grid->jtot;j++)
  {
    for (int i=0; i<grid->itot;i++)
    {
      ij = i + j*grid->itot;
      mfile << noise[ij] << " ";
    }
    mfile << endl;
  }
  mfile.close();
  }

  // split whole noise signal in parts of different processors
  jj = grid->icells;
  for (int j=grid->jstart; j<grid->jend; ++j)
  {
    for (int i=grid->istart; i<grid->iend; ++i)
    { 
      iref = i - grid->igc;
      jref = j - grid->jgc;
      ij = i + j*jj;
      igjg = master->mpicoordy*grid->jmax*grid->itot + jref*grid->itot + master->mpicoordx*grid->imax + iref;
      tmp1[ij] = noise[igjg];
      // if (master->mpiid == 0)
      //   cout << "ij " << ij << " " << master->mpicoordx << " " << master->mpicoordy << endl;
    }
  }

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij = i + j*jj;
        a[ij] = tmp1[ij] - offset;
      }
  }
  else if(sw == BC_NEUMANN)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij = i + j*jj;
        agrad[ij] = tmp1[ij];
        aflux[ij] = -agrad[ij]*visc;
      }
  }
  else if(sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ij = i + j*jj;
        aflux[ij] = tmp1[ij];
        agrad[ij] = -aflux[ij]/visc;
      }
  }

  grid->boundary_cyclic2d(a);
  grid->boundary_cyclic2d(aflux);
  grid->boundary_cyclic2d(agrad);

  

  ofstream lfile;
  std::ostringstream fileNameStream("field_p");
  fileNameStream << "field_p" << master->mpiid << ".txt";
  std::string fileName = fileNameStream.str();
  lfile.open(fileName.c_str(),ios::app);

  for (int j=0;j<grid->jcells;++j)
  {
    for (int i=0; i<grid->icells;++i)
    {
      ij = i + j*jj;
      lfile << aflux[ij] << " ";
    }
    lfile << endl;
  }
  lfile.close(); 


  return 0;
}
// --------------------------------------------------------------------------------------------------------

int cboundary_malte::calc_stats(double * a, int n, double * mean_val, double * var_val)
{ 
  double sum = 0;

  // mean value
  for(int k=0;k<n;k++)
    sum += a[k];  
  *mean_val = sum/(n);

  sum = 0;

   for(int k=0;k<n;k++)
    sum += pow((*mean_val-a[k]),2);
  *var_val = sum/(n);
  
  return 0;
}


int cboundary_malte::calc_step(double * step, double aval)
{
  double xmod, ymod;
  double errvalx, errvaly;
  int ij;
  for(int j=0; j<grid->jtot; j++)
  {
    for(int i=0; i<grid->itot; i++)
    {
      ij = i + j*(grid->itot);
      xmod = fmod((grid->xsize/grid->itot)*i, patch_xh);
      ymod = fmod((grid->ysize/grid->jtot)*j, patch_xh);

      errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

      if(patch_dim == 2)
        errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
      else
        errvaly = 1.;

      step[ij] = errvalx*errvaly*aval;

    }
  }

  return 0;
}



