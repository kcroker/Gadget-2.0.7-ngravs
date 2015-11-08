#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


/*! \file pm_nonperiodic.c
 *  \brief code for non-periodic FFT to compute long-range PM force
 */


#ifdef PMGRID
#if !defined (PERIODIC) || defined (PLACEHIGHRESREGION)

#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

#include "allvars.h"
#include "proto.h"

#define  GRID  (2*PMGRID)
#define  GRID2 (2*(GRID/2 + 1))

// KC 11/6/15
// PUSH: this is for calculation of the r-space vacuum
// green's function from the *periodic* k-space greens function
//
#define KGRID (2*GRID)
#define KGRID2 (2*(KGRID/2 + 1))
//static int kern_slab_to_task[KGRID];
static rfftwnd_mpi_plan kern_ifft_inverse_plan;
static int *kern_slabs_per_task;
static int *kern_first_slab_of_task;
static int kern_slabstart_x, kern_nslab_x, kern_slabstart_y, kern_nslab_y;
static int kern_fftsize, kern_maxfftsize;
static fftw_real *prekernel, *kern_workspace;
static fftw_real *ifft_of_prekernel;

// Original stuff
static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
static int slab_to_task[GRID];
static int *slabs_per_task;
static int *first_slab_of_task;

static int *meshmin_list, *meshmax_list;

static int slabstart_x, nslab_x, slabstart_y, nslab_y;

static int fftsize, maxfftsize;

static fftw_real *kernel[2][N_GRAVS][N_GRAVS], *rhogrid, *forcegrid, *workspace;
static fftw_complex *fft_of_kernel[2][N_GRAVS][N_GRAVS], *fft_of_rhogrid;

/*! This function determines the particle extension of all particles, and for
 *  those types selected with PLACEHIGHRESREGION if this is used, and then
 *  determines the boundaries of the non-periodic FFT-mesh that can be placed
 *  on this region. Note that a sufficient buffer region at the rim of the
 *  occupied part of the mesh needs to be reserved in order to allow a correct
 *  finite differencing using a 4-point formula. In addition, to allow
 *  non-periodic boundaries, the actual FFT mesh used is twice as large in
 *  each dimension compared with PMGRID.
 */
void pm_init_regionsize(void)
{
  double meshinner[2], xmin[2][3], xmax[2][3];
  int i, j, t;

  /* find enclosing rectangle */

  for(j = 0; j < 3; j++)
    {
      xmin[0][j] = xmin[1][j] = 1.0e36;
      xmax[0][j] = xmax[1][j] = -1.0e36;
    }

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	t = 0;
#ifdef PLACEHIGHRESREGION
	if(((1 << P[i].Type) & (PLACEHIGHRESREGION)))
	  t = 1;
#endif
	if(P[i].Pos[j] > xmax[t][j])
	  xmax[t][j] = P[i].Pos[j];
	if(P[i].Pos[j] < xmin[t][j])
	  xmin[t][j] = P[i].Pos[j];
      }

  MPI_Allreduce(xmin, All.Xmintot, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, All.Xmaxtot, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  for(j = 0; j < 2; j++)
    {
      All.TotalMeshSize[j] = All.Xmaxtot[j][0] - All.Xmintot[j][0];
      All.TotalMeshSize[j] = dmax(All.TotalMeshSize[j], All.Xmaxtot[j][1] - All.Xmintot[j][1]);
      All.TotalMeshSize[j] = dmax(All.TotalMeshSize[j], All.Xmaxtot[j][2] - All.Xmintot[j][2]);
#ifdef ENLARGEREGION
      All.TotalMeshSize[j] *= ENLARGEREGION;
#endif

      /* symmetrize the box onto the center */
      for(i = 0; i < 3; i++)
	{
	  All.Xmintot[j][i] = (All.Xmintot[j][i] + All.Xmaxtot[j][i]) / 2 - All.TotalMeshSize[j] / 2;
	  All.Xmaxtot[j][i] = All.Xmintot[j][i] + All.TotalMeshSize[j];
	}
    }

  /* this will produce enough room for zero-padding and buffer region to
     allow finite differencing of the potential  */

  for(j = 0; j < 2; j++)
    {
      meshinner[j] = All.TotalMeshSize[j];
      All.TotalMeshSize[j] *= 2.001 * (GRID) / ((double) (GRID - 2 - 8));
    }

  /* move lower left corner by two cells to allow finite differencing of the potential by a 4-point function */

  for(j = 0; j < 2; j++)
    for(i = 0; i < 3; i++)
      {
	All.Corner[j][i] = All.Xmintot[j][i] - 2.0005 * All.TotalMeshSize[j] / GRID;
	All.UpperCorner[j][i] = All.Corner[j][i] + (GRID / 2 - 1) * (All.TotalMeshSize[j] / GRID);
      }


#ifndef PERIODIC
  All.Asmth[0] = ASMTH * All.TotalMeshSize[0] / GRID;
  All.Rcut[0] = RCUT * All.Asmth[0];
#endif

#ifdef PLACEHIGHRESREGION
  All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
  All.Rcut[1] = RCUT * All.Asmth[1];
#endif

#ifdef PLACEHIGHRESREGION
  if(2 * All.TotalMeshSize[1] / GRID < All.Rcut[0])
    {
      All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double) (GRID - 2));

      for(i = 0; i < 3; i++)
	{
	  All.Corner[1][i] = All.Xmintot[1][i] - 1.0001 * All.Rcut[0];
	  All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
	}

      if(2 * All.TotalMeshSize[1] / GRID > All.Rcut[0])
	{
	  All.TotalMeshSize[1] = 2 * (meshinner[1] + 2 * All.Rcut[0]) * (GRID) / ((double) (GRID - 10));

	  for(i = 0; i < 3; i++)
	    {
	      All.Corner[1][i] = All.Xmintot[1][i] - 1.0001 * (All.Rcut[0] + 2 * All.TotalMeshSize[1] / GRID);
	      All.UpperCorner[1][i] = All.Corner[1][i] + (GRID / 2 - 1) * (All.TotalMeshSize[1] / GRID);
	    }
	}

      All.Asmth[1] = ASMTH * All.TotalMeshSize[1] / GRID;
      All.Rcut[1] = RCUT * All.Asmth[1];
    }
#endif

  if(ThisTask == 0)
    {
#ifndef PERIODIC
      printf("\nAllowed region for isolated PM mesh (coarse):\n");
      printf("(%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  meshsize=%g\n\n",
	     All.Xmintot[0][0], All.Xmintot[0][1], All.Xmintot[0][2],
	     All.Xmaxtot[0][0], All.Xmaxtot[0][1], All.Xmaxtot[0][2], meshinner[0], All.TotalMeshSize[0],
	     All.TotalMeshSize[0] / GRID);
#endif
#ifdef PLACEHIGHRESREGION
      printf("\nAllowed region for isolated PM mesh (high-res):\n");
      printf("(%g|%g|%g)  -> (%g|%g|%g)   ext=%g  totmeshsize=%g  meshsize=%g\n\n",
	     All.Xmintot[1][0], All.Xmintot[1][1], All.Xmintot[1][2],
	     All.Xmaxtot[1][0], All.Xmaxtot[1][1], All.Xmaxtot[1][2],
	     meshinner[1], All.TotalMeshSize[1], All.TotalMeshSize[1] / GRID);
#endif
    }

}

/*! Initialization of the non-periodic PM routines. The plan-files for FFTW
 *  are created. Finally, the routine to set-up the non-periodic Greens
 *  function is called.
 */
void pm_init_nonperiodic(void)
{
  int i, slab_to_task_local[GRID];
  double bytes_tot = 0;
  size_t bytes;
  int n, m;

  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, GRID, GRID, GRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, GRID, GRID, GRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */
  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < GRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, GRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  slabs_per_task = malloc(NTask * sizeof(int));

  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  // KC 11/8/15
  // Use full complex arrays because we need it to work now
  kern_ifft_inverse_plan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, KGRID, KGRID, KGRID, 
						  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  fftwnd_mpi_local_sizes(kern_ifft_inverse_plan, 
			  &kern_nslab_x, 
			  &kern_slabstart_x, 
			  &kern_nslab_y, 
			  &kern_slabstart_y, 
			  &kern_fftsize);
  
  kern_slabs_per_task = malloc(NTask * sizeof(int));
  MPI_Allgather(&kern_nslab_x, 
		1, 
		MPI_INT, 
		kern_slabs_per_task, 
		1, 
		MPI_INT, 
		MPI_COMM_WORLD);

#ifndef PERIODIC
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	printf("Task=%d  FFT-Slabs=%d\n", i, slabs_per_task[i]);
    }
#endif

  first_slab_of_task = malloc(NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  // KC 11/6/15
  kern_first_slab_of_task = malloc(NTask * sizeof(int));
  MPI_Allgather(&kern_slabstart_x, 
		1, 
		MPI_INT, 
		kern_first_slab_of_task, 
		1, 
		MPI_INT, 
		MPI_COMM_WORLD);

  meshmin_list = malloc(3 * NTask * sizeof(int));
  meshmax_list = malloc(3 * NTask * sizeof(int));

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // KC 11/6/15
  MPI_Allreduce(&kern_fftsize, 
		&kern_maxfftsize, 
		1, 
		MPI_INT, 
		MPI_MAX, 
		MPI_COMM_WORLD);

  /* now allocate memory to hold the FFT fields */

#if !defined(PERIODIC)
  // KC 11/6/15
  // Finally, allocate the memory for my slice of the prekernel enormous transform
  if(!(prekernel = (fftw_real *) malloc(bytes = kern_fftsize * sizeof(fftw_real))))
    {
      printf("ngravs: failed to allocate memory for the non-periodic prekernel (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  for(n = 0; n < N_GRAVS; ++n) {
    for(m = 0; m < N_GRAVS; ++m) {
      
      if(!(kernel[0][n][m] = (fftw_real *) malloc(bytes = fftsize * sizeof(fftw_real))))
	{
	  printf("failed to allocate memory for `FFT-kernel[0]' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;
      fft_of_kernel[0][n][m] = (fftw_complex *) kernel[0][n][m];
    }
  }
#endif

#if defined(PLACEHIGHRESREGION)
  for(n = 0; n < N_GRAVS; ++n) {
    for(m = 0; m < N_GRAVS; ++m) {
      
      if(!(kernel[1][n][m] = (fftw_real *) malloc(bytes = fftsize * sizeof(fftw_real))))
	{
	  printf("failed to allocate memory for `FFT-kernel[1]' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;
      fft_of_kernel[1][n][m] = (fftw_complex *) kernel[1][n][m];
    }
  }
#endif

  if(ThisTask == 0)
    printf("\nAllocated %g MByte for FFT kernel(s) and necessary workspace.\n\n", bytes_tot / (1024.0 * 1024.0));

}


/*! This function allocates the workspace needed for the non-periodic FFT
 *  algorithm. Three fields are used, one for the density/potential fields,
 *  one to hold the force field obtained by finite differencing, and finally
 *  an additional workspace which is used both in the parallel FFT itself, and
 *  as a buffer for the communication algorithm.
 */
void pm_init_nonperiodic_allocate(int dimprod)
{
  static int first_alloc = 1;
  int dimprodmax;
  double bytes_tot = 0;
  size_t bytes;

  // KC 11/6/15
  // Note the strange casts are intentional.
  if(!(prekernel = (fftw_real *) malloc(bytes = kern_fftsize * sizeof(fftw_real)))) {
    
    printf("ngravs: failed to allocate memory for the non-periodic prekernel (%g MB).\n", bytes / (1024.0 * 1024.0));
    endrun(1);
  }
  bytes_tot += bytes;
  
  if(!(kern_workspace = (fftw_real *) malloc(bytes = kern_maxfftsize * sizeof(fftw_real))))
    {
      printf("ngravs: failed to allocate memory for the non-periodic prekernel workspace (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
  
  MPI_Allreduce(&dimprod, &dimprodmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  
  if(!(rhogrid = (fftw_real *) malloc(bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  fft_of_rhogrid = (fftw_complex *) rhogrid;

  if(!(forcegrid = (fftw_real *) malloc(bytes = imax(fftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(workspace = (fftw_real *) malloc(bytes = imax(maxfftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-workspace' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(first_alloc == 1)
    {
      first_alloc = 0;
      if(ThisTask == 0)
	printf("\nUsing %g MByte for non-periodic FFT computation.\n\n", bytes_tot / (1024.0 * 1024.0));
    }
}


/*! This function frees the memory allocated for the non-periodic FFT
 *  computation. (With the exception of the Greens function(s), which are kept
 *  statically in memory for the next force computation.)
 */
void pm_init_nonperiodic_free(void)
{
  /* deallocate memory */
  free(workspace);
  free(forcegrid);
  free(rhogrid);

  // KC 11/6/15
  // Deallocate the enormous prekernel
  free(prekernel);
  free(kern_workspace);
}


/*! This function sets-up the Greens function for the non-periodic potential
 *  in real space, and then converts it to Fourier space by means of a FFT.
 */

// KC 10/17/14
// This function will need initialize N_GRAVS^2 more kernels
// both for the coarse mesh and the fine mesh.
// Now we gotta do N_GRAVS^2 of them
// which is kind of wasteful because we should only need to do N_GRAVS(N_GRAVS + 1)/2 of 
//
// This is actually quite a bit of memory for large meshes and large D.  D=5 is 25 vs 15.
// So we *should* make the shift to symmetric checks.
//

// // A = Passive, B = Active
// #define NGL_MACRO(A, B) ( A > N_GRAVS/2 ? (N_GRAVS - A)*(N_GRAVS-1) + B : A*(N_GRAVS-1) + B )

// KC 11/5/15
// This function is carried out more than once, any time the extent of the mesh changes, this 
// must be recomputed.  Gnarly.
//
// XXX
// Right now, this has a fucked tangle of #if !defined PERIODICs...
//
void pm_setup_nonperiodic_kernel(void)
{
  int i, j, k;
  double x, y, z, r, u;
  double kx, ky, kz, k2, fx, fy, fz, ff;
  int ip;

  // ngravs
  int n, m, d;
  double asmth2[2];
  double fac[2];
  double e[2];

  /* now set up kernel and its Fourier transform */

  pm_init_nonperiodic_allocate(0);

#if !defined(PERIODIC)
  for(n = 0; n < N_GRAVS; ++n)
    for(m = 0; m < N_GRAVS; ++m)
      for(i = 0; i < fftsize; i++)	/* clear local density field */
   	kernel[0][n][m][i] = 0;

  // CONTINUE AS DONE BEFORE, BUT SAMPLE THE PREKERNEL

  // KC 11/4/15
  // Note that we must carefully use the non-periodically determined mesh sizes
  // XXX
  // These are probably wrong
  asmth2[0] = (2 * M_PI) * All.Asmth[0] / All.TotalMeshSize[0];
  asmth2[0] *= asmth2[0];
  fac[0] = All.G / pow(All.TotalMeshSize[0], 4) * pow(All.TotalMeshSize[0] / GRID, 3);	/* to get potential */

  printf("ngravs: fac[0] = %f\n", fac[0]);

#if defined PLACEHIGHRESREGION
  asmth2[1] = (2 * M_PI) * All.Asmth[1] / All.TotalMeshSize[1];
  asmth2[1] *= asmth2[1];
  fac[1] = All.G / pow(All.TotalMeshSize[1], 4) * pow(All.TotalMeshSize[1] / GRID, 3);	/* to get potential */
#endif
  
  // UH oh.  What is asmth2[0] and [1]?  Have they been set yet??
  
  for(n = 0; n < N_GRAVS; ++n) {
    for(m = 0; m < N_GRAVS; ++m) {
      
      // KC 11/6/15
      // Do the very large kernels from the periodic k-space functions
      for(d = 0; d < 2; ++d) {	
#if !defined PLACEHIGHRESREGION
	if(d == 1)
	  break;
#endif

	//////////////////////// BEGIN PREKERNEL CREATION
	//
	for(i = 0; i < kern_fftsize; ++i)
	  prekernel[i] = 0.0;
	
	// KC 11/8/15
	// Note z only goes to half-size.  The rest of the array gets wonk
	// ordered.
	// Note that we are transposed in x and y because TRANSPOSED_ORDER
	for(y = slabstart_y; y < slabstart_y + nslab_y; y++) {
	  for(x = 0; x < KGRID; x++) {
	    for(z = 0; z < KGRID / 2 + 1; z++) {
	      if(x > KGRID / 2)
		kx = x - KGRID;
	      else
		kx = x;
	      if(y > KGRID / 2)
		ky = y - KGRID;
	      else
		ky = y;
	      if(z > KGRID / 2)
		kz = z - KGRID;
	      else
		kz = z;

	      k2 = kx * kx + ky * ky + kz * kz;

	      if(k2 > 0)
		{
		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / KGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / KGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / KGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  ff = ff * ff * ff * ff;

		  ip = KGRID * (KGRID / 2 + 1) * (y - slabstart_y) + (KGRID / 2 + 1) * x + z;
	      
		  // KC 11/6/15
		  // Now we first put the antialiasing in place
		  // And then multiply by the appropriately truncated periodic k-space greens function
#if !defined(PERIODIC)
		  if(d == 0) {
		    prekernel[ip] = ff * (-exp(-k2 * ((double) ASMTH/GRID) * ((double)ASMTH/GRID))) *
  		      (*GreensFxns[n][m])(All.MassTable[n], All.MassTable[m], k2, 0.0, 1) * fac[0] / k2;

		    // KC 11/8/15
		    // And we have to get the upper set properly too...
		    if(ip > 0 && ip < kern_fftsize / 2)
		      prekernel[kern_fftsize - ip] = ff * (-exp(-k2 * ((double) ASMTH/GRID) * ((double)ASMTH/GRID))) *
			(*GreensFxns[n][m])(All.MassTable[n], All.MassTable[m], k2, 0.0, 1) * fac[0] / k2;

		    fprintf("%d %d %d  
		  }
#endif
#if defined(PLACEHIGHRESREGION)
		  //    -(/*exp(-k2 * asmth2[0]) - */ exp(-k2 * pow(All.Asmth[0]/All.Asmth[1], 2))/(All.Asmth[0]/All.Asmth[1]))
		  if(d == 1) {
		    prekernel[ip] = ff;
 /* *  */
 /* 		      (exp(-k2 * pow(All.Asmth[0]/All.Asmth[1],2))*(All.Asmth[0]/All.Asmth[1]) - exp(-k2 * asmth2[0])) *  */
 /* 		      (*GreensFxns[n][m])(All.MassTable[n], All.MassTable[m], k2, 0.0, 1) * fac[d] / k2; */

		    prekernel[ip] = 0.0; 

		    /* // Output difference */
		    /* fprintf(stderr, "%.15e %.15e %.15e\n",  */
		    /* 	    sqrt(k2),  */
		    /* 	    ff * (-exp(-k2*asmth2[0])) * */
		    /* 	    (*GreensFxns[n][m])(All.MassTable[n], All.MassTable[m], k2, 0.0, 1) * fac[0] / k2, */
		    /* 	    (exp(-k2 * pow(All.Asmth[0]/All.Asmth[1], 2))*(All.Asmth[0]/All.Asmth[1]) - exp(-k2 * asmth2[0])) *  */
		    /* 	    (*GreensFxns[n][m])(All.MassTable[n], All.MassTable[m], k2, 0.0, 1) * fac[d] / k2); */

		  }
#endif
		}
	    }
	  }
	}
      
	// At this point, x,y,z and i,j,k are usable indexes...
	// d = 0 if doing low-res mesh, d = 1 if doing high-res mesh.

	// KC 11/6/15
	// This code may be wrong, I got it from pm_periodic.c...
	// Seems legit ;)
	if(slabstart_y == 0)
	  prekernel[0] = 0.0;

	// KC 11/6/15
	// Perform the inverse transform!
	rfftwnd_mpi(kern_ifft_inverse_plan, 1, prekernel, kern_workspace, FFTW_TRANSPOSED_ORDER);
	
	// SANITY CHECK.  Output the prekernel, where OUR SLAB overlaps in the lower octant
	// it should look like erf(u) in *all directions.*
	for(i = kern_slabstart_x; i < GRID && i < (kern_slabstart_x + kern_nslab_x); i++) {
	  for(j = 0; j < GRID; j++) {
	    for(k = 0; k < GRID; k++) {
	      
	      x = ((double) i) / GRID;
	      y = ((double) j) / GRID;
	      z = ((double) k) / GRID;
	      
	      if(x >= 0.5)
		x -= 1.0;
	      if(y >= 0.5)
		y -= 1.0;
	      if(z >= 0.5)
		z -= 1.0;

	      r = sqrt(x * x + y * y + z * z);

	      u = 0.5 * r / (((double) ASMTH) / GRID);
	      e[0] = 1-erfc(u);
	      e[1] = erfc(u * All.Asmth[1] / All.Asmth[0]) - erfc(u); 
	      // Alternatively: -erf(u * All.Asmth[1] / All.Asmth[0]) + 1 - erfc(u)
	      
	      if(r > 0) {
		fprintf(stderr, "%.15e %.15e %.15e\n",
			r,
			fac[d] * prekernel[GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k],
			- e[d] / r);
	      }
	    }
	  }
	}
	fprintf(stderr, "\n# Second prekernel begins\n");
      }
      endrun(5656);
    }
  }
#endif

/*       } */
    

/* 		/\* kernel[d][GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k] = -1/r + fac[0] *  *\/ */
/* 		/\*   prekernel[GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k] ;//-fac / r; *\/ */

/* 	// ... PUSH! */
/* 	// Now the corner octant of dim GRID^3 has the correct kernel, assign it */

/* 	// KC 11/6/15 */
/* 	// What if our slab lies partially in, or entirely outside of the corner octant? */
/* 	// Short circuit the assignment if we sit outside of the active octant */
/* 	// */
/* 	// XXX */
/* 	// There is a misplementation here, maybe.  Right now kern_workspace contains only  */
/* 	// the slab that we were to compute? */
/* 	for(i = kern_slabstart_x; i < GRID && i < (kern_slabstart_x + kern_nslab_x); i++) */
/* 	  for(j = 0; j < GRID; j++) */
/* 	    for(k = 0; k < GRID; k++) */
/* 	      { */
/* 	      x = ((double) i) / GRID; */
/* 	      y = ((double) j) / GRID; */
/* 	      z = ((double) k) / GRID; */
	      
/* 	      if(x >= 0.5) */
/* 		x -= 1.0; */
/* 	      if(y >= 0.5) */
/* 		y -= 1.0; */
/* 	      if(z >= 0.5) */
/* 		z -= 1.0; */

/* 	      r = sqrt(x * x + y * y + z * z); */

/* 	      u = 0.5 * r / (((double) ASMTH) / GRID); */

/* 	      //	      fac = 1 - erfc(u); */

/* 	      if(r > 0) */
/* 		kernel[d][GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k] = -1/r + fac[0] *  */
/* 		  prekernel[GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k] ;//-fac / r; */
/* 	      else // XXX fix me for the PLACEHIRESREGION kernels! */
/* 		kernel[d][GRID * GRID2 * (i - kern_slabstart_x) + GRID2 * j + k] = */
/* 		  -1 / (sqrt(M_PI) * (((double) ASMTH) / GRID)); */
/* 	    } */

/*       ////////// END PREKERNEL CREATION */


/* /\* #if defined(PLACEHIGHRESREGION) *\/ */
/* /\*   for(n = 0; n < N_GRAVS; ++n) *\/ */
/* /\*     for(m = 0; m < N_GRAVS; ++m) { *\/ */
/* /\*       for(i = 0; i < fftsize; i++)	/\\* clear local density field *\\/ *\/ */
/* /\* 	kernel[1][n][m][i] = 0; *\/ */

/* /\*       for(i = slabstart_x; i < (slabstart_x + nslab_x); i++) *\/ */
/* /\* 	for(j = 0; j < GRID; j++) *\/ */
/* /\* 	  for(k = 0; k < GRID; k++) *\/ */
/* /\* 	    { *\/ */
/* /\* 	      x = ((double) i) / GRID; *\/ */
/* /\* 	      y = ((double) j) / GRID; *\/ */
/* /\* 	      z = ((double) k) / GRID; *\/ */
	      
/* /\* 	      if(x >= 0.5) *\/ */
/* /\* 		x -= 1.0; *\/ */
/* /\* 	      if(y >= 0.5) *\/ */
/* /\* 		y -= 1.0; *\/ */
/* /\* 	      if(z >= 0.5) *\/ */
/* /\* 		z -= 1.0; *\/ */

/* /\* 	      r = sqrt(x * x + y * y + z * z); *\/ */
	      
/* /\* 	      u = 0.5 * r / (((double) ASMTH) / GRID); *\/ */

/* /\* 	      // KC 11/1/15 *\/ */
/* /\* 	      // XXX!! *\/ */
/* /\* 	      // See above comments.   *\/ */
/* /\* 	      // This one is really bad though as the smoothing scales, which change, *\/ */
/* /\* 	      // would need to enter the tabulation.  No dice. *\/ */
/* /\* 	      fac = erfc(u * All.Asmth[1] / All.Asmth[0]) - erfc(u); *\/ */

/* /\* 	      if(r > 0) { *\/ */
/* /\* 		kernel[1][n][m][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] =  *\/ */
/* /\* 		  -fac * (*PotentialFxns[n][m])(MassTable[nA], MassTable[nB], 0.0, r, 1); *\/ */
/* /\* 	      } *\/ */
/* /\* 	      else { *\/ */
	    
/* /\* 		fac = 1 - All.Asmth[1] / All.Asmth[0]; *\/ */
/* /\* 	    	kernel[1][n][m][GRID * GRID2 * (i - slabstart_x) + GRID2 * j + k] = *\/ */
/* /\* 		  -fac * PotentialZero[n][m]; /// (sqrt(M_PI) * (((double) ASMTH) / GRID)); *\/ */
/* /\* 	      } *\/ */
/* /\* 	    } *\/ */
    

/*       /\* do the forward transform of the kernel *\/ */
/*       rfftwnd_mpi(fft_forward_plan, 1, kernel[1][n][m], workspace, FFTW_TRANSPOSED_ORDER); */
/*     } */
/* /\* #endif *\/ */

/*   // KC 11/23/14 */
/*   // This function must not even be called if we are not doing non-periodic styles */

/*   /\* deconvolve the Greens function twice with the CIC kernel *\/ */

/*   for(y = slabstart_y; y < slabstart_y + nslab_y; y++) */
/*     for(x = 0; x < GRID; x++) */
/*       for(z = 0; z < GRID / 2 + 1; z++) */
/* 	{ */
/* 	  if(x > GRID / 2) */
/* 	    kx = x - GRID; */
/* 	  else */
/* 	    kx = x; */
/* 	  if(y > GRID / 2) */
/* 	    ky = y - GRID; */
/* 	  else */
/* 	    ky = y; */
/* 	  if(z > GRID / 2) */
/* 	    kz = z - GRID; */
/* 	  else */
/* 	    kz = z; */

/* 	  k2 = kx * kx + ky * ky + kz * kz; */

/* 	  if(k2 > 0) */
/* 	    { */
/* 	      fx = fy = fz = 1; */
/* 	      if(kx != 0) */
/* 		{ */
/* 		  fx = (M_PI * kx) / GRID; */
/* 		  fx = sin(fx) / fx; */
/* 		} */
/* 	      if(ky != 0) */
/* 		{ */
/* 		  fy = (M_PI * ky) / GRID; */
/* 		  fy = sin(fy) / fy; */
/* 		} */
/* 	      if(kz != 0) */
/* 		{ */
/* 		  fz = (M_PI * kz) / GRID; */
/* 		  fz = sin(fz) / fz; */
/* 		} */
/* 	      ff = 1 / (fx * fy * fz); */
/* 	      ff = ff * ff * ff * ff; */

/* 	      ip = GRID * (GRID / 2 + 1) * (y - slabstart_y) + (GRID / 2 + 1) * x + z; */
/* 	      for(n = 0; n < N_GRAVS; ++n) */
/* 		for(m = 0; m < N_GRAVS; ++m) { */
		
/* #if !defined(PERIODIC) */
/* 		  fft_of_kernel[0][n][m][ip].re *= ff; */
/* 		  fft_of_kernel[0][n][m][ip].im *= ff; */
/* #endif */
/* #if defined(PLACEHIGHRESREGION) */
/* 		  fft_of_kernel[1][n][m][ip].re *= ff; */
/* 		  fft_of_kernel[1][n][m][ip].im *= ff; */
/* #endif */
/* 		} */
/* 	    } */
/* 	} */
/*   /\* end deconvolution *\/ */

  pm_init_nonperiodic_free();
}



/*! Calculates the long-range non-periodic forces using the PM method.  The
 *  potential is Gaussian filtered with Asmth, given in mesh-cell units. The
 *  potential is finite differenced using a 4-point finite differencing
 *  formula to obtain the force fields, which are then interpolated to the
 *  particle positions. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The CIC kernel is deconvolved.
 */
int pmforce_nonperiodic(int grnr)
{
  double dx, dy, dz;
  double fac, to_slab_fac;
  double re, im, acc_dim;
  int i, j, slab, level, sendTask, recvTask, flag, flagsum;
  int x, y, z, xl, yl, zl, xr, yr, zr, xll, yll, zll, xrr, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;
  MPI_Status status;

  // ngravs extension
  int nA, nB;
  int offsets[N_GRAVS + 1];
  int k;

  if(ThisTask == 0)
    printf("Starting non-periodic PM calculation (grid=%d).\n", grnr);

  fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);	/* to get potential */
  fac *= 1 / (2 * All.TotalMeshSize[grnr] / GRID);	/* for finite differencing */

  to_slab_fac = GRID / All.TotalMeshSize[grnr];


  /* first, establish the extension of the local patch in GRID (for binning) */

  for(j = 0; j < 3; j++)
    {
      meshmin[j] = GRID;
      meshmax[j] = 0;
    }

  // KC 11/7/14
  // Compute the gravitational offsets here
  offsets[0] = 0;
  for(k = 1; k < N_GRAVS; ++k)
    offsets[k] = offsets[k-1] + NgravLocal[k-1];
  offsets[N_GRAVS] = NumPart; 

  for(i = 0, flag = 0; i < NumPart; i++)
    {
#ifdef PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && ((1 << P[i].Type) & (PLACEHIGHRESREGION))))
#endif
	{
	  for(j = 0; j < 3; j++)
	    {
	      if(P[i].Pos[j] < All.Xmintot[grnr][j] || P[i].Pos[j] > All.Xmaxtot[grnr][j])
		{
		  if(flag == 0)
		    {
		      printf
			("Particle Id=%d on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\nStopping\n",
			 (int)P[i].ID, ThisTask, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		  flag++;
		  break;
		}
	    }
	}

      if(flag > 0)
	continue;

      // KC 11/7/14
      // grnr - GridNumberReference ?
      if(P[i].Pos[0] >= All.Corner[grnr][0] && P[i].Pos[0] < All.UpperCorner[grnr][0])
	if(P[i].Pos[1] >= All.Corner[grnr][1] && P[i].Pos[1] < All.UpperCorner[grnr][1])
	  if(P[i].Pos[2] >= All.Corner[grnr][2] && P[i].Pos[2] < All.UpperCorner[grnr][2])
	    {
	      for(j = 0; j < 3; j++)
		{
		  slab = to_slab_fac * (P[i].Pos[j] - All.Corner[grnr][j]);

		  if(slab < meshmin[j])
		    meshmin[j] = slab;

		  if(slab > meshmax[j])
		    meshmax[j] = slab;
		}
	    }
    }


  MPI_Allreduce(&flag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0)
    {
      if(ThisTask == 0)
	{
	  printf("In total %d particles were outside allowed range.\n", flagsum);
	  fflush(stdout);
	}
      return 1;			/* error - need to return because particle were outside allowed range */
    }

  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);

  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  force_treefree();

  pm_init_nonperiodic_allocate((dimx + 4) * (dimy + 4) * (dimz + 4));

  // Sources
  for(nA = 0; nA < N_GRAVS; ++nA) {

    // Receivers
    for(nB = 0; nB < N_GRAVS; ++nB) {

      for(i = 0; i < dimx * dimy * dimz; i++)
	workspace[i] = 0;

      for(i = offsets[nA]; i < offsets[nA+1]; ++i) {

	if(P[i].Pos[0] < All.Corner[grnr][0] || P[i].Pos[0] >= All.UpperCorner[grnr][0])
	  continue;
	if(P[i].Pos[1] < All.Corner[grnr][1] || P[i].Pos[1] >= All.UpperCorner[grnr][1])
	  continue;
	if(P[i].Pos[2] < All.Corner[grnr][2] || P[i].Pos[2] >= All.UpperCorner[grnr][2])
	  continue;

	slab_x = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]);
	dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
	slab_x -= meshmin[0];
	slab_xx = slab_x + 1;

	slab_y = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]);
	dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
	slab_y -= meshmin[1];
	slab_yy = slab_y + 1;

	slab_z = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]);
	dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;
	slab_z -= meshmin[2];
	slab_zz = slab_z + 1;

	workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
	workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * dy * dz;

	workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (dx) * dy * (1.0 - dz);
	workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (dx) * (1.0 - dy) * dz;
	workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (dx) * dy * dz;
      }
  
      for(i = 0; i < fftsize; i++)	/* clear local density field */
	rhogrid[i] = 0;

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;
	  if(recvTask < NTask)
	    {
	      /* check how much we have to send */
	      sendmin = 2 * GRID;
	      sendmax = -1;
	      for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
		if(slab_to_task[slab_x] == recvTask)
		  {
		    if(slab_x < sendmin)
		      sendmin = slab_x;
		    if(slab_x > sendmax)
		      sendmax = slab_x;
		  }
	      if(sendmax == -1)
		sendmin = 0;

	      /* check how much we have to receive */
	      recvmin = 2 * GRID;
	      recvmax = -1;
	      for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
		if(slab_to_task[slab_x] == sendTask)
		  {
		    if(slab_x < recvmin)
		      recvmin = slab_x;
		    if(slab_x > recvmax)
		      recvmax = slab_x;
		  }
	      if(recvmax == -1)
		recvmin = 0;

	      if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
		{
		  recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
		  recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
		  recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

		  if(level > 0)
		    {
		      MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
				   (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
				   TAG_NONPERIOD_A | (nA << 6) | (nB << 9), forcegrid,
				   (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_A | (nA << 6) | (nB << 9), MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			     (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }

		  for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		    {
		      slab_xx = slab_x - first_slab_of_task[ThisTask];

		      if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
			{
			  for(slab_y = meshmin_list[3 * recvTask + 1];
			      slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			    {
			      slab_yy = slab_y;

			      for(slab_z = meshmin_list[3 * recvTask + 2];
				  slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
				{
				  slab_zz = slab_z;

				  rhogrid[GRID * GRID2 * slab_xx + GRID2 * slab_yy + slab_zz] +=
				    forcegrid[((slab_x - recvmin) * recv_dimy +
					       (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					      (slab_z - meshmin_list[3 * recvTask + 2])];
				}
			    }
			}
		    }
		}
	    }
	}
    


      /* Do the FFT of the density field */

      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);


      /* multiply with the Fourier transform of the Green's function (kernel) */

      for(y = 0; y < nslab_y; y++)
	for(x = 0; x < GRID; x++)
	  for(z = 0; z < GRID / 2 + 1; z++)
	    {
	      ip = GRID * (GRID / 2 + 1) * y + (GRID / 2 + 1) * x + z;

	      re =
		fft_of_rhogrid[ip].re * fft_of_kernel[grnr][nA][nB][ip].re -
		fft_of_rhogrid[ip].im * fft_of_kernel[grnr][nA][nB][ip].im;

	      im =
		fft_of_rhogrid[ip].re * fft_of_kernel[grnr][nA][nB][ip].im +
		fft_of_rhogrid[ip].im * fft_of_kernel[grnr][nA][nB][ip].re;

	      fft_of_rhogrid[ip].re = re;
	      fft_of_rhogrid[ip].im = im;
	    }
    
      /* get the potential by inverse FFT */

      rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* Now rhogrid holds the potential */
      /* construct the potential for the local patch */


      /* if we have a high-res mesh, establish the extension of the local patch in GRID (for reading out the
       * forces) 
       */
      // KC 11/7/14
      // This code is baffling.  Because exactly the same meshmin/max and meshlists will
      // be computed (on all the processors?) as were computed the first time this 
      // same code was executed when we checked for grid OOB problems.  Here we 
      // test for grnr == 1 and consider *only* contributions from the high-resolution type.
      // These are the same conditions above (because grnr will never change within any one 
      // run of this method.)  I am commenting this out and we will see.

/* #ifdef PLACEHIGHRESREGION */
/*       if(grnr == 1) */
/* 	{ */
/* 	  for(j = 0; j < 3; j++) */
/* 	    { */
/* 	      meshmin[j] = GRID; */
/* 	      meshmax[j] = 0; */
/* 	    } */

/* 	  for(i = 0; i < NumPart; i++) */
/* 	    { */
/* 	      if(!((1 << P[i].Type) & (PLACEHIGHRESREGION))) */
/* 		continue; */


/* 	      if(P[i].Pos[0] >= All.Corner[grnr][0] && P[i].Pos[0] < All.UpperCorner[grnr][0]) */
/* 		if(P[i].Pos[1] >= All.Corner[grnr][1] && P[i].Pos[1] < All.UpperCorner[grnr][1]) */
/* 		  if(P[i].Pos[2] >= All.Corner[grnr][2] && P[i].Pos[2] < All.UpperCorner[grnr][2]) */
/* 		    { */
/* 		      for(j = 0; j < 3; j++) */
/* 			{ */
/* 			  slab = to_slab_fac * (P[i].Pos[j] - All.Corner[grnr][j]); */

/* 			  if(slab < meshmin[j]) */
/* 			    meshmin[j] = slab; */

/* 			  if(slab > meshmax[j]) */
/* 			    meshmax[j] = slab; */
/* 			} */
/* 		    } */
/* 	    } */

/* 	  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD); */
/* 	  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD); */
/* 	} */
/* #endif */

      // /XXX? Resume porting...
 
      dimx = meshmax[0] - meshmin[0] + 6;
      dimy = meshmax[1] - meshmin[1] + 6;
      dimz = meshmax[2] - meshmin[2] + 6;

      for(j = 0; j < 3; j++)
	{
	  if(meshmin[j] < 2)
	    endrun(131231);
	  if(meshmax[j] > GRID / 2 - 3)
	    endrun(131288);
	}

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      /* check how much we have to send */
	      sendmin = 2 * GRID;
	      sendmax = -GRID;
	      for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
		if(slab_to_task[slab_x] == sendTask)
		  {
		    if(slab_x < sendmin)
		      sendmin = slab_x;
		    if(slab_x > sendmax)
		      sendmax = slab_x;
		  }
	      if(sendmax == -GRID)
		sendmin = sendmax + 1;


	      /* check how much we have to receive */
	      recvmin = 2 * GRID;
	      recvmax = -GRID;
	      for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
		if(slab_to_task[slab_x] == recvTask)
		  {
		    if(slab_x < recvmin)
		      recvmin = slab_x;
		    if(slab_x > recvmax)
		      recvmax = slab_x;
		  }
	      if(recvmax == -GRID)
		recvmin = recvmax + 1;

	      if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
		{
		  recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
		  recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
		  recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = slab_x - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = slab_y;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z < meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = slab_z;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[GRID * GRID2 * slab_xx + GRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_B | (nA << 6) | (nB << 9),
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_B | (nA << 6) | (nB << 9), MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}

      dimx = meshmax[0] - meshmin[0] + 2;
      dimy = meshmax[1] - meshmin[1] + 2;
      dimz = meshmax[2] - meshmin[2] + 2;

      recv_dimx = meshmax[0] - meshmin[0] + 6;
      recv_dimy = meshmax[1] - meshmin[1] + 6;
      recv_dimz = meshmax[2] - meshmin[2] + 6;


      for(dim = 0; dim < 3; dim++)	/* Calculate each component of the force. */
	{
	  /* get the force component by finite differencing the potential */
	  /* note: "workspace" now contains the potential for the local patch, plus a suffiently large buffer region */

	  for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
	    for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
	      for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
		{
		  xrr = xll = xr = xl = x;
		  yrr = yll = yr = yl = y;
		  zrr = zll = zr = zl = z;

		  switch (dim)
		    {
		    case 0:
		      xr = x + 1;
		      xrr = x + 2;
		      xl = x - 1;
		      xll = x - 2;
		      break;
		    case 1:
		      yr = y + 1;
		      yl = y - 1;
		      yrr = y + 2;
		      yll = y - 2;
		      break;
		    case 2:
		      zr = z + 1;
		      zl = z - 1;
		      zrr = z + 2;
		      zll = z - 2;
		      break;
		    }

		  forcegrid[(x * dimy + y) * dimz + z]
		    =
		    fac * ((4.0 / 3) *
			   (workspace[((xl + 2) * recv_dimy + (yl + 2)) * recv_dimz + (zl + 2)]
			    - workspace[((xr + 2) * recv_dimy + (yr + 2)) * recv_dimz + (zr + 2)]) -
			   (1.0 / 6) *
			   (workspace[((xll + 2) * recv_dimy + (yll + 2)) * recv_dimz + (zll + 2)] -
			    workspace[((xrr + 2) * recv_dimy + (yrr + 2)) * recv_dimz + (zrr + 2)]));
		}


	  /* read out the forces */

	  // KC 11/7/14
	  // Only apply forces to the nB type
	  for(i = offsets[nB]; i < offsets[nB+1]; i++)
	    {
#ifdef PLACEHIGHRESREGION
	      if(grnr == 1)
		if(!((1 << P[i].Type) & (PLACEHIGHRESREGION)))
		  continue;
#endif
	      slab_x = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]);
	      dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
	      slab_x -= meshmin[0];
	      slab_xx = slab_x + 1;

	      slab_y = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]);
	      dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
	      slab_y -= meshmin[1];
	      slab_yy = slab_y + 1;

	      slab_z = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]);
	      dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;
	      slab_z -= meshmin[2];
	      slab_zz = slab_z + 1;

	      acc_dim =
		forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	      acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
	      acc_dim += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
	      acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

	      acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
	      acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
	      acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
	      acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;

	      P[i].GravPM[dim] += acc_dim;
	    }
	}
      // Close ngravs loops
    }
  }

  pm_init_nonperiodic_free();
  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);
  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  if(ThisTask == 0)
    printf("done PM.\n");

  return 0;
}



/*! Calculates the long-range non-periodic potential using the PM method.  The
 *  potential is Gaussian filtered with Asmth, given in mesh-cell units.  We
 *  carry out a CIC charge assignment, and compute the potenial by Fourier
 *  transform methods. The CIC kernel is deconvolved.
 */
int pmpotential_nonperiodic(int grnr)
{
  double dx, dy, dz;
  double fac, to_slab_fac;
  double re, im, pot;
  int i, j, slab, level, sendTask, recvTask, flag, flagsum;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;
  MPI_Status status;

  // ngravs extension
  int nA, nB;
  int offsets[N_GRAVS + 1];
  int k;

  if(ThisTask == 0)
    printf("Starting non-periodic PM-potential calculation.\n");

  fac = All.G / pow(All.TotalMeshSize[grnr], 4) * pow(All.TotalMeshSize[grnr] / GRID, 3);	/* to get potential */

  to_slab_fac = GRID / All.TotalMeshSize[grnr];

  /* first, establish the extension of the local patch in GRID (for binning) */

  for(j = 0; j < 3; j++)
    {
      meshmin[j] = GRID;
      meshmax[j] = 0;
    }

  // KC 11/7/14
  // Compute the gravitational offsets here
  offsets[0] = 0;
  for(k = 1; k < N_GRAVS; ++k)
    offsets[k] = offsets[k-1] + NgravLocal[k-1];
  offsets[N_GRAVS] = NumPart; 

  for(i = 0, flag = 0; i < NumPart; i++)
    {
#ifdef PLACEHIGHRESREGION
      if(grnr == 0 || (grnr == 1 && ((1 << P[i].Type) & (PLACEHIGHRESREGION))))
#endif
	{
	  for(j = 0; j < 3; j++)
	    {
	      if(P[i].Pos[j] < All.Xmintot[grnr][j] || P[i].Pos[j] > All.Xmaxtot[grnr][j])
		{
		  if(flag == 0)
		    {
		      printf
			("Particle Id=%d on task=%d with coordinates (%g|%g|%g) lies outside PM mesh.\nStopping\n",
			 (int)P[i].ID, ThisTask, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		  flag++;
		  break;
		}
	    }
	}

      if(flag > 0)
	continue;

      if(P[i].Pos[0] >= All.Corner[grnr][0] && P[i].Pos[0] < All.UpperCorner[grnr][0])
	if(P[i].Pos[1] >= All.Corner[grnr][1] && P[i].Pos[1] < All.UpperCorner[grnr][1])
	  if(P[i].Pos[2] >= All.Corner[grnr][2] && P[i].Pos[2] < All.UpperCorner[grnr][2])
	    {
	      for(j = 0; j < 3; j++)
		{
		  slab = to_slab_fac * (P[i].Pos[j] - All.Corner[grnr][j]);

		  if(slab < meshmin[j])
		    meshmin[j] = slab;

		  if(slab > meshmax[j])
		    meshmax[j] = slab;
		}
	    }
    }


  MPI_Allreduce(&flag, &flagsum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(flagsum > 0) 
    {
      if(ThisTask == 0)
	{
	  printf("In total %d particles were outside allowed range.\n", flagsum);
	  fflush(stdout);
	}
      return 1;			/* error - need to return because particle were outside allowed range */
    }



  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);

  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;


  force_treefree();

  pm_init_nonperiodic_allocate((dimx + 4) * (dimy + 4) * (dimz + 4));

  // Sources
  for(nA = 0; nA < N_GRAVS; ++nA) {

    // Receivers
    for(nB = 0; nB < N_GRAVS; ++nB) {
      
      for(i = 0; i < dimx * dimy * dimz; i++)
	workspace[i] = 0;

      for(i = offsets[nA]; i < offsets[nA+1]; ++i) {
	
	if(P[i].Pos[0] < All.Corner[grnr][0] || P[i].Pos[0] >= All.UpperCorner[grnr][0])
	  continue;
	if(P[i].Pos[1] < All.Corner[grnr][1] || P[i].Pos[1] >= All.UpperCorner[grnr][1])
	  continue;
	if(P[i].Pos[2] < All.Corner[grnr][2] || P[i].Pos[2] >= All.UpperCorner[grnr][2])
	  continue;

	slab_x = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]);
	dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
	slab_x -= meshmin[0];
	slab_xx = slab_x + 1;

	slab_y = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]);
	dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
	slab_y -= meshmin[1];
	slab_yy = slab_y + 1;

	slab_z = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]);
	dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;
	slab_z -= meshmin[2];
	slab_zz = slab_z + 1;

	workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
	workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * dy * dz;

	workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (dx) * dy * (1.0 - dz);
	workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (dx) * (1.0 - dy) * dz;
	workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (dx) * dy * dz;
      }
    

      for(i = 0; i < fftsize; i++)	/* clear local density field */
	rhogrid[i] = 0;

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;
	  if(recvTask < NTask)
	    {
	      /* check how much we have to send */
	      sendmin = 2 * GRID;
	      sendmax = -1;
	      for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
		if(slab_to_task[slab_x] == recvTask)
		  {
		    if(slab_x < sendmin)
		      sendmin = slab_x;
		    if(slab_x > sendmax)
		      sendmax = slab_x;
		  }
	      if(sendmax == -1)
		sendmin = 0;

	      /* check how much we have to receive */
	      recvmin = 2 * GRID;
	      recvmax = -1;
	      for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
		if(slab_to_task[slab_x] == sendTask)
		  {
		    if(slab_x < recvmin)
		      recvmin = slab_x;
		    if(slab_x > recvmax)
		      recvmax = slab_x;
		  }
	      if(recvmax == -1)
		recvmin = 0;

	      if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
		{
		  recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
		  recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
		  recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

		  if(level > 0)
		    {
		      MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
				   (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
				   TAG_NONPERIOD_C | (nA << 6) | (nB << 9), forcegrid,
				   (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_C | (nA << 6) | (nB << 9), MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			     (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }

		  for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		    {
		      slab_xx = slab_x - first_slab_of_task[ThisTask];

		      if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
			{
			  for(slab_y = meshmin_list[3 * recvTask + 1];
			      slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			    {
			      slab_yy = slab_y;

			      for(slab_z = meshmin_list[3 * recvTask + 2];
				  slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
				{
				  slab_zz = slab_z;

				  rhogrid[GRID * GRID2 * slab_xx + GRID2 * slab_yy + slab_zz] +=
				    forcegrid[((slab_x - recvmin) * recv_dimy +
					       (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					      (slab_z - meshmin_list[3 * recvTask + 2])];
				}
			    }
			}
		    }
		}
	    }
	}


      /* Do the FFT of the density field */

      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);


      /* multiply with the Fourier transform of the Green's function (kernel) */

      for(y = 0; y < nslab_y; y++)
	for(x = 0; x < GRID; x++)
	  for(z = 0; z < GRID / 2 + 1; z++)
	    {
	      ip = GRID * (GRID / 2 + 1) * y + (GRID / 2 + 1) * x + z;

	      re =
		fft_of_rhogrid[ip].re * fft_of_kernel[grnr][nA][nB][ip].re -
		fft_of_rhogrid[ip].im * fft_of_kernel[grnr][nA][nB][ip].im;

	      im =
		fft_of_rhogrid[ip].re * fft_of_kernel[grnr][nA][nB][ip].im +
		fft_of_rhogrid[ip].im * fft_of_kernel[grnr][nA][nB][ip].re;

	      fft_of_rhogrid[ip].re = fac * re;
	      fft_of_rhogrid[ip].im = fac * im;
	    }

      /* get the potential by inverse FFT */

      rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

      /* Now rhogrid holds the potential */
      /* construct the potential for the local patch */


      /* if we have a high-res mesh, establish the extension of the local patch in GRID (for reading out the
   * forces) 
   */

/* #ifdef PLACEHIGHRESREGION */
/*   if(grnr == 1) */
/*     { */
/*       for(j = 0; j < 3; j++) */
/* 	{ */
/* 	  meshmin[j] = GRID; */
/* 	  meshmax[j] = 0; */
/* 	} */

/*       for(i = 0; i < NumPart; i++) */
/* 	{ */
/* 	  if(!((1 << P[i].Type) & (PLACEHIGHRESREGION))) */
/* 	    continue; */


/* 	  if(P[i].Pos[0] >= All.Corner[grnr][0] && P[i].Pos[0] < All.UpperCorner[grnr][0]) */
/* 	    if(P[i].Pos[1] >= All.Corner[grnr][1] && P[i].Pos[1] < All.UpperCorner[grnr][1]) */
/* 	      if(P[i].Pos[2] >= All.Corner[grnr][2] && P[i].Pos[2] < All.UpperCorner[grnr][2]) */
/* 		{ */
/* 		  for(j = 0; j < 3; j++) */
/* 		    { */
/* 		      slab = to_slab_fac * (P[i].Pos[j] - All.Corner[grnr][j]); */

/* 		      if(slab < meshmin[j]) */
/* 			meshmin[j] = slab; */

/* 		      if(slab > meshmax[j]) */
/* 			meshmax[j] = slab; */
/* 		    } */
/* 		} */
/* 	} */

/*       MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD); */
/*       MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD); */
/*     } */
/* #endif */

      // XXX? Again, above code looked absolutely redundant.  Commented out
      dimx = meshmax[0] - meshmin[0] + 6;
      dimy = meshmax[1] - meshmin[1] + 6;
      dimz = meshmax[2] - meshmin[2] + 6;

      for(j = 0; j < 3; j++)
	{
	  if(meshmin[j] < 2)
	    endrun(131231);
	  if(meshmax[j] > GRID / 2 - 3)
	    endrun(131288);
	}

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      /* check how much we have to send */
	      sendmin = 2 * GRID;
	      sendmax = -GRID;
	      for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
		if(slab_to_task[slab_x] == sendTask)
		  {
		    if(slab_x < sendmin)
		      sendmin = slab_x;
		    if(slab_x > sendmax)
		      sendmax = slab_x;
		  }
	      if(sendmax == -GRID)
		sendmin = sendmax + 1;


	      /* check how much we have to receive */
	      recvmin = 2 * GRID;
	      recvmax = -GRID;
	      for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
		if(slab_to_task[slab_x] == recvTask)
		  {
		    if(slab_x < recvmin)
		      recvmin = slab_x;
		    if(slab_x > recvmax)
		      recvmax = slab_x;
		  }
	      if(recvmax == -GRID)
		recvmin = recvmax + 1;

	      if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
		{
		  recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
		  recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
		  recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = slab_x - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = slab_y;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z < meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = slab_z;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[GRID * GRID2 * slab_xx + GRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_NONPERIOD_D | (nA << 6) | (nB << 9),
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_NONPERIOD_D | (nA << 6) | (nB << 9), MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}

      dimx = meshmax[0] - meshmin[0] + 2;
      dimy = meshmax[1] - meshmin[1] + 2;
      dimz = meshmax[2] - meshmin[2] + 2;

      recv_dimx = meshmax[0] - meshmin[0] + 6;
      recv_dimy = meshmax[1] - meshmin[1] + 6;
      recv_dimz = meshmax[2] - meshmin[2] + 6;


      for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
	for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
	  for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
	    {
	      forcegrid[(x * dimy + y) * dimz + z]
		= workspace[((x + 2) * recv_dimy + (y + 2)) * recv_dimz + (z + 2)];
	    }


      /* read out the potential */

      for(i = offsets[nB]; i < offsets[nB+1]; ++i)
	{
#ifdef PLACEHIGHRESREGION
	  if(grnr == 1)
	    if(!((1 << P[i].Type) & (PLACEHIGHRESREGION)))
	      continue;
#endif
	  slab_x = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]);
	  dx = to_slab_fac * (P[i].Pos[0] - All.Corner[grnr][0]) - slab_x;
	  slab_x -= meshmin[0];
	  slab_xx = slab_x + 1;

	  slab_y = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]);
	  dy = to_slab_fac * (P[i].Pos[1] - All.Corner[grnr][1]) - slab_y;
	  slab_y -= meshmin[1];
	  slab_yy = slab_y + 1;

	  slab_z = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]);
	  dz = to_slab_fac * (P[i].Pos[2] - All.Corner[grnr][2]) - slab_z;
	  slab_z -= meshmin[2];
	  slab_zz = slab_z + 1;

	  pot = forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  pot += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
	  pot += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
	  pot += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

	  pot += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
	  pot += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
	  pot += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
	  pot += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;

	  P[i].Potential += pot;

	}
      // close out ngravs loops
    }
  }

  pm_init_nonperiodic_free();
  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);
  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  if(ThisTask == 0)
    printf("done PM-potential.\n");

  return 0;
}

#endif
#endif
