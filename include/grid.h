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

#ifndef GRID
#define GRID

#ifdef PARALLEL
#include <mpi.h>
#endif
#include <fftw3.h>
#include "input.h"

// forward declaration
class cmodel;
class cmaster;

/**
 * Class for the grid settings and operators.
 * This class contains the grid properties, such as dimensions and resolution.
 * The public funcions in this class contain operations that are called from many routines
 * in order to interpolate, transpose and save data. The MPI operations that work over multiple
 * processes on the entire grid are contained in this class.
 */
class cgrid
{
  public:
    cgrid(cmodel *); ///< Constructor of the grid class.
    ~cgrid();        ///< Destructor of the grid class.

    int readinifile(cinput *); ///< Processes data from the input file.
    int init();                ///< Initialization of the grid arrays.
    int create(cinput *);      ///< Creation of the grid data.
    int calculate();           ///< Computation of dimensions, faces and ghost cells.
    int save();                ///< Saves grid data to file.
    int load();                ///< Loads grid data to file.

    int itot; ///< Total number of grid cells in the x-direction.
    int jtot; ///< Total number of grid cells in the y-direction.
    int ktot; ///< Total number of grid cells in the z-direction.

    int imax; ///< Number of grid cells in the x-direction for one process.
    int jmax; ///< Number of grid cells in the y-direction for one process.
    int kmax; ///< Number of grid cells in the z-direction for one process.

    int iblock; ///< Number of grid cells in the x-direction for calculation blocks in transposes.
    int jblock; ///< Number of grid cells in the y-direction for calculation blocks in transposes.
    int kblock; ///< Number of grid cells in the z-direction for calculation blocks in transposes.

    int igc; ///< Number of ghost cells in the x-direction.
    int jgc; ///< Number of ghost cells in the y-direction.
    int kgc; ///< Number of ghost cells in the z-direction.

    int icells;  ///< Number of grid cells in the x-direction including ghost cells for one process.
    int jcells;  ///< Number of grid cells in the y-direction including ghost cells for one process.
    int ijcells; ///< Number of grid cells in the xy-plane including ghost cells for one process.
    int kcells;  ///< Number of grid cells in the z-direction including ghost cells for one process.
    int ncells;  ///< Total number of grid cells for one process including ghost cells.
    int istart;  ///< Index of the first grid point in the x-direction.
    int jstart;  ///< Index of the first grid point in the y-direction.
    int kstart;  ///< Index of the first grid point in the z-direction.
    int iend;    ///< Index of the last gridpoint+1 in the x-direction.
    int jend;    ///< Index of the last gridpoint+1 in the y-direction.
    int kend;    ///< Index of the last gridpoint+1 in the z-direction.

    double xsize; ///< Size of the domain in the x-direction.
    double ysize; ///< Size of the domain in the y-direction.
    double zsize; ///< Size of the domain in the z-direction.

    double dx;     ///< Distance between the center of two grid cell in the x-direction.
    double dy;     ///< Distance between the center of two grid cell in the y-direction.
    double *dz;    ///< Distance between the center of two grid cell in the z-direction.
    double *dzh;   ///< Distance between the two grid cell faces in the z-direction.
    double *dzi;   ///< Reciprocal of dz.
    double *dzhi;  ///< Reciprocal of dzh.
    double *dzi4;  ///< Fourth order gradient of the distance between cell centers to be used in 4th-order schemes.
    double *dzhi4; ///< Fourth order gradient of the distance between cell faces to be used in 4th-order schemes.
    
    double *x;  ///< Grid coordinate of cell center in x-direction.
    double *y;  ///< Grid coordinate of cell center in y-direction.
    double *z;  ///< Grid coordinate of cell center in z-direction.
    double *xh; ///< Grid coordinate of cell faces in x-direction.
    double *yh; ///< Grid coordinate of cell faces in x-direction.
    double *zh; ///< Grid coordinate of cell faces in x-direction.

    double utrans; ///< Galilean transformation velocity in x-direction.
    double vtrans; ///< Galilean transformation velocity in y-direction.

    std::string swspatialorder; ///< Default spatial order of the operators to be used on this grid.

    // MPI functions
    int initmpi(); ///< Creates the MPI data types used in grid operations.
    int exitmpi(); ///< Destructs the MPI data types used in grid operations.
    int boundary_cyclic  (double *); ///< Fills the ghost cells in the periodic directions.
    int boundary_cyclic2d(double *); ///< Fills the ghost cells of one slice in the periodic direction.
    int transposezx(double *, double *); ///< Changes the transpose orientation from z to x.
    int transposexz(double *, double *); ///< Changes the transpose orientation from x to z.
    int transposexy(double *, double *); ///< Changes the transpose orientation from x to y.
    int transposeyx(double *, double *); ///< Changes the transpose orientation from y to x.
    int transposeyz(double *, double *); ///< Changes the transpose orientation from y to z.
    int transposezy(double *, double *); ///< Changes the transpose orientation from z to y.

    int getmax (double *);      ///< Gets the maximum of a number over all processes.
    int getsum (double *);      ///< Gets the sum of a number over all processes.
    int getprof(double *, int); ///< Averages a vertical profile over all processes.
    int calcmean(double *, const double *, int);

    // IO functions
    int savefield3d(double *, double *, double *, char *, double); ///< Saves a full 3d field.
    int loadfield3d(double *, double *, double *, char *, double); ///< Loads a full 3d field.

    int savexzslice(double *, double *, char *, int);           ///< Saves a xz-slice from a 3d field.
    int savexyslice(double *, double *, char *, int kslice=-1); ///< Saves a xy-slice from a 3d field.
    int loadxyslice(double *, double *, char *, int kslice=-1); ///< Loads a xy-slice.

    // Fourier tranforms
    double *fftini, *fftouti; ///< Help arrays for fast-fourier transforms in x-direction.
    double *fftinj, *fftoutj; ///< Help arrays for fast-fourier transforms in y-direction.
    fftw_plan iplanf, iplanb; ///< FFTW3 plans for forward and backward transforms in x-direction.
    fftw_plan jplanf, jplanb; ///< FFTW3 plans for forward and backward transforms in y-direction.

    int fftforward (double *, double *, double *, double *, double *, double *); ///< Forward fast-fourier transform.
    int fftbackward(double *, double *, double *, double *, double *, double *); ///< Backward fast-fourier transform.

    // interpolation functions
    int interpolatex_2nd(double *, double *, int); ///< Second order interpolation in the x-direction.
    int interpolatey_2nd(double *, double *, int); ///< Second order interpolation in the y-direction.
    int interpolatex_4th(double *, double *, int); ///< Fourth order interpolation in the x-direction.
    int interpolatey_4th(double *, double *, int); ///< Fourth order interpolation in the y-direction.

  private:
    cmaster *master; ///< Pointer to master class.
    bool allocated;  ///< Boolean to check whether grid data is allocated.
    bool mpitypes;   ///< Boolean to check whether MPI datatypes are created.
    bool fftwplan;   ///< Boolean to check whether FFTW3 plans are created.

#ifdef PARALLEL
    // MPI Datatypes
    MPI_Datatype eastwestedge;     ///< MPI datatype containing the ghostcells at the east-west sides.
    MPI_Datatype northsouthedge;   ///< MPI datatype containing the ghostcells at the north-south sides.
    MPI_Datatype eastwestedge2d;   ///< MPI datatype containing the ghostcells for one slice at the east-west sides.
    MPI_Datatype northsouthedge2d; ///< MPI datatype containing the ghostcells for one slice at the north-south sides.

    MPI_Datatype transposez;  ///< MPI datatype containing base blocks for z-orientation in zx-transpose.
    MPI_Datatype transposez2; ///< MPI datatype containing base blocks for z-orientation in zy-transpose.
    MPI_Datatype transposex;  ///< MPI datatype containing base blocks for x-orientation in zx-transpose.
    MPI_Datatype transposex2; ///< MPI datatype containing base blocks for x-orientation in xy-transpose.
    MPI_Datatype transposey;  ///< MPI datatype containing base blocks for y-orientation in xy-transpose.
    MPI_Datatype transposey2; ///< MPI datatype containing base blocks for y-orientation in zy-transpose.

    MPI_Datatype subi;       ///< MPI datatype containing a subset of the entire x-axis.
    MPI_Datatype subj;       ///< MPI datatype containing a subset of the entire y-axis.
    MPI_Datatype subarray;   ///< MPI datatype containing the dimensions of the total array that is contained in one process.
    MPI_Datatype subxzslice; ///< MPI datatype containing only one xz-slice.
    MPI_Datatype subxyslice; ///< MPI datatype containing only one xy-slice.

    double *profl; ///< Help array used in profile writing.
#endif
};
#endif
