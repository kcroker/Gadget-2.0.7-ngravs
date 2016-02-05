#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
//#include <sys/types.h>
//#include <unistd.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_erf.h>
#include "ngravs.h"
#include "allvars.h"
#include "proto.h"

/*! \file ngravs_core.c
 *  \brief Contains init and the high-precision fourier transform code used for generalied TreePM
 *
 * This file contains portions of the ngravs extension that should not need to be altered by the 
 * end-researcher.  It takes care of initialization, sanity checking, reporting on the main thread,
 * and high-precision computation of the Gaussian filtered generalized short range force.
 *
 */

struct ngravsInterpolant *ngravsPeriodicTable;

double normKtoGridK(double normk) {

  // See comments below for gridKtoNormK
  return normk * All.BoxSize/ (4*M_PI*All.Asmth[0]);
}

double gridKtoNormK(double gridk) {

  // Since the correct coefficient for normk exponential suppression
  // is 1/2 (at least in the Newtonian case, but this should not change)
  // but that for gridk is (2 * M_PI) * All.Asmth[0] / All.BoxSize
  // this is the conversion

  return 4*M_PI*All.Asmth[0] * gridk / All.BoxSize;
}

///////////////// BEGIN FOURIER INTEGRATION ROUTINES /////////////////////
//
// These routines can compute the reqired shortrange tabulations of the generic
// force laws from the k-space greens functions to very near machine accuracy.  
//
// Please see the ngravs paper for detailed discussion of the behaviour here
//

FLOAT jTok(int m, double Z, struct ngravsInterpolant *s) {

  return 2.0 * M_PI * m * s->ntab * 6.0 * s->ol/(3.0 * s->ngravs_tpm_n);
}

FLOAT mTox(int j, struct ngravsInterpolant *s) {
  
  return 3.0*j/(6.0 * s->ntab * s->ol);
}

int gadgetToFourier(int j, struct ngravsInterpolant *s) {

  return s->ol * (6*j + 3);
}

int fourierToGadget(int i, struct ngravsInterpolant *s) {

  return (i - 3 * s->ol)/(6 * s->ol);
}

FLOAT fourierIntegrand(FLOAT k, gravity normKGreen, FLOAT Z) {
  
  FLOAT k2 = k*k;
  
  return (*normKGreen)(1, 1, k2, k, 1) * exp(-k2 * Z * Z);
}

int performConvolution(struct ngravsInterpolant *s, gravity normKGreen, FLOAT Z, FLOAT *oRes, FLOAT *oResI) {
  
  fftw_complex *in, *out;
  int m,j;
  double sum, norm;

  in = (fftw_complex *)malloc(sizeof(fftw_complex) * s->ngravs_tpm_n);
  out = (fftw_complex *)malloc(sizeof(fftw_complex) * s->ngravs_tpm_n);
  if(!in || !out)
    return 1;
  
  /* // Debug your mappings?! */
  /* printf("s->ntab: %d\nNGRAVS_TPM_N: %d\nNGRAVS_TPM_N/2: %d\n3/(2s->ntab): %f\n", s->ntab, s->ngravs_tpm_n, s->ngravs_tpm_n/2, 3.0/(2*s->ntab)); */
  /* for(m = 0; m < s->ntab; ++m) { */
  /*   printf("%d %d %d %f\n", m, gadgetToFourier(m, s), fourierToGadget(gadgetToFourier(m, s), s), mTox(gadgetToFourier(m, s), s)); */
  /* } */
  
  /* printf("\n"); */
  /* for(j = 0; j < NGRAVS_TPM_N/2; ++j) { */
  /*   printf("%d %d %d %f %d\n", j, fourierToGadget(j, s), gadgetToFourier(fourierToGadget(j, s)), mTox(j, s), gadgetToFourier(fourierToGadget(j, s), s) - j); */
  /* } */
  /* exit(0); */

  // Zero out all arrays first
  for(j = 0; j < s->ngravs_tpm_n; ++j) {
    in[j].re = 0;
    in[j].im = 0;
  }

  // 1) FFTW needs this loaded in wonk order
  // Note we need zero power in order to compute the 
  // potential term correctly.
  in[0].re = fourierIntegrand(jTok(0, Z, s), normKGreen, Z);

  for(j = 1; j < s->ngravs_tpm_n/2; ++j) {

    in[j].re = fourierIntegrand(jTok(j, Z, s), normKGreen, Z);
    in[s->ngravs_tpm_n - j].re = fourierIntegrand(jTok(j, Z, s), normKGreen, Z);
  }

  // 2) Xform
  fftw_one(s->plan, in, out);
  norm = 2.0*M_PI * s->ntab * 6.0 * s->ol/(3.0 * s->ngravs_tpm_n);

  /* // Debug */
  /* // Make sure the transform is behaving reasonably */
  /* for(m = 0; m < s->ngravs_tpm_n; ++m) */
  /*   printf("%.15f %.15f\n", mTox(m, s), out[m].re*norm); */
  /* exit(0); */

  // ???
  sum = s->ngravs_tpm_n;

  for(m = 0; m < s->ntab; ++m)
    oRes[m] = out[gadgetToFourier(m, s)].re * norm;

  // 3) Integrate so as to constrain the error correctly:
  // Newton-Cotes 4-point rule
  // Run the sum at double precision, though we may assign to lower precision
  // First term of in[] should be zero.
  in[0].re = 0.0;
  sum = 0.0;
  for(m = 0; m < s->ngravs_tpm_n-3; m += 3) {
    sum += (mTox(m+3, s) - mTox(m, s)) * 0.125 * norm * (out[m].re + 3.0*out[m+1].re + 3.0*out[m+2].re + out[m+3].re);
    
    // Put it where you'd expect to find it
    in[m/3+1].re = sum;
  }

  // 3.5) Downsample
  for(m = 0; m < s->ntab; ++m)
    oResI[m] = in[gadgetToFourier(m, s)/3].re;

  /* // Debug */
  /* // Make sure the transform is behaving reasonably */
  /* for(m = 0; m < s->ntab; ++m) */
  /*   printf("%.15f %.15f\n", mTox(gadgetToFourier(m, s)), oResI[m]); */
  /* exit(0); */

  free(in);
  free(out);
  return 0;
}

/*! Create a workspace for a table with ntab entries, which goes len deep
 *  into x-space, and oversamples at a frequency of ol (so ol 2 is 2x as many 
 *  samples as would be needed to critically sample */
struct ngravsInterpolant *ngravsConvolutionInit(int ntab, int len, int ol) {

  struct ngravsInterpolant *s;

  if(! (s = (struct ngravsInterpolant *) malloc(sizeof(struct ngravsInterpolant)))) {

    printf("ngravs: could not allocate table structure... something really bad is happening because I'm tiny!");
    endrun(1045);
  }
  
  s->ntab = ntab;
  s->len = len;
  s->ol = ol;
  s->ngravs_tpm_n = 12*ntab*ol*len-6*ol*len+2;
  s->plan = fftw_create_plan(s->ngravs_tpm_n, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  if(!ThisTask)
    printf("ngravs: max_k = %.15e\n", jTok(s->ngravs_tpm_n/2, 0.5, s));

  return s;
}

void ngravsConvolutionFree(struct ngravsInterpolant *s) {

  fftw_destroy_plan(s->plan);
  free(s);
}

//////////////////////// END FOURIER INTEGRATION ROUTINES //////////

/*! Establish the mapping between particle type and the native
 *  gravitational force between two particles of this type.
 *  Then call the user modified wiring function to initialize the table of 
 *  AccelFxn function pointers to point to the appropriate compute functions.  
 *  After, perform sanity checks on the user.
 */
void init_grav_maps(void) {

  int i, j;
  int counts[N_GRAVS];
  int n;

#ifdef BAMTEST
  double q;

  if(ThisTask == 0) {

    for(q=1.0e-4; q < 10.0; q += 1.0e-4)
      printf("%e %e\n", q, bambam(2*BAM_EPSILON*M_PI, 2*BAM_EPSILON*M_PI, q, 1));
    
    printf("\n");
    
    for(q=1.0e-4; q < 10.0; q += 1.0e-4)
      printf("%e %e\n", q, sourcebaryonbam(4*BAM_EPSILON*M_PI, BAM_EPSILON*M_PI, q, 1));

    printf("\n");
    
    for(q=1.0e-4; q < 10.0; q += 1.0e-4)
      printf("%e %e\n", q, sourcebambaryon(BAM_EPSILON*M_PI, 4*BAM_EPSILON*M_PI, q, 1));
    
    MPI_Barrier(MPI_COMM_WORLD);
    exit(1001);
  }
  else
    MPI_Barrier(MPI_COMM_WORLD);
#endif

  // KC 10/17/14
  // Code consistency checks
  if(ThisTask == 0) {

#ifdef PLACEHIRESREGION
    printf("ngravs: High-resolution mesh requires non-periodic boundary condition FFT, which has been disabled in ngravs (see paper).");
    endrun(1000);
#endif

#if defined PMGRID && !defined PERIODIC
    printf("ngravs: Non-periodic boundary condition FFT has been disabled in ngravs.  Please disable PMGRID.");
#endif

#if !defined PEANOHILBERT && (defined PMGRID || defined PLACEHIRESREGION || defined PERIODIC)
    printf("ngravs: Gravitational type ordering required for ngravs extension of Fourier mesh code.  This ordering has been implemented on top of Peano-Hilbert ordering within the local task.  Please enable PEANOHILBERT and recompile.\n");
    endrun(1000);
#endif

  }

  memset(counts, 0, sizeof(int)*N_GRAVS);
 
  /* KC 8/9/14 */
  TypeToGrav[0] = All.GravityGas;
#ifdef PMGRID
  if(TypeToGrav[0]) {

    printf("ngravs: Fourier mesh code requires that gas particles must interact under gravity type 0.  Please rewire and recompile.\n");
    endrun(1000);
  }
#endif

  TypeToGrav[1] = All.GravityHalo;
  TypeToGrav[2] = All.GravityDisk;
  TypeToGrav[3] = All.GravityBulge;
  TypeToGrav[4] = All.GravityStars;
  TypeToGrav[5] = All.GravityBndry;

  // Sanity check
  for(i = 0; i < 6; ++i) {
    if(TypeToGrav[i] >= N_GRAVS || TypeToGrav[i] < 0) {
      
      printf("ngravs: native interaction %d declared for type %d does not exist.  We better stop\n", TypeToGrav[i], i);
      endrun(1000);
    }
    
    // See how many times we use this interaction
    counts[TypeToGrav[i]]++;
  }
  
  if(ThisTask == 0) {
    printf("This is Gadget2-ngravs.\n");
    printf("ngravs: timestep reduction scaling parameter set to %f\n", NGRAVS_TIMESTEP_SCALE);
  }

  // Provide some output
  for(i = 0; i < N_GRAVS; ++i) {
    
     if(ThisTask == 0) {
      
       if(counts[i] > 0)
	 printf("ngravs: %d particle types natively interact via gravity %d\n", counts[i], i);
       
       if(counts[i] == 6 && N_GRAVS > 1)
	 printf("ngravs: ALL types natively interact via gravity %d.  Recomple with N_GRAVS=1 for best memory performance\n", i);
     }
  }
    
  // First set the AccelFxns to null, so that we can check for missing settings
  memset(AccelFxns, 0, N_GRAVS*N_GRAVS*sizeof(gravity));
  memset(AccelSplines, 0, N_GRAVS*N_GRAVS*sizeof(gravity));

#if defined PERIODIC
  memset(LatticeForce, 0, N_GRAVS*N_GRAVS*sizeof(latforce));
  memset(LatticePotential, 0, N_GRAVS*N_GRAVS*sizeof(latpot));
#endif

#if defined OUTPUTPOTENTIAL || defined PMGRID
  memset(PotentialFxns, 0, N_GRAVS*N_GRAVS*sizeof(gravity));
  memset(PotentialSplines, 0, N_GRAVS*N_GRAVS*sizeof(gravity));
#endif

#ifdef PMGRID
  memset(GreensFxns, 0, N_GRAVS*N_GRAVS*sizeof(gravity));
#endif

  // This is the function you should alter
  wire_grav_maps();

  // Now run the sanity check
  for(i = 0; i < N_GRAVS; ++i) {

    for(j = 0; j < N_GRAVS; ++j) {

      // Check for null
      if(!AccelFxns[i][j]) {
	printf("ngravs: unwired spot for acceleration [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }

      if(!AccelSplines[i][j]) {
	printf("ngravs: unwired spot for softening spline [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }

#if defined PERIODIC && (!defined PMGRID || defined FORCETEST)
      if(!LatticeForce[i][j]) {
	printf("ngravs: uwired spot for lattice correction force [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }

      if(!LatticePotential[i][j]) {
	printf("ngravs: uwired spot for lattice potential [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }
#endif

#if defined OUTPUTPOTENTIAL || defined PMGRID
      if(!PotentialFxns[i][j]) {
	printf("ngravs: uwired spot for potential [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }

      if(!PotentialSplines[i][j]) {
	printf("ngravs: uwired spot for softening potential spline [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }
#endif

#ifdef PMGRID
      if(!GreensFxns[i][j]) {
	printf("ngravs: unwired spot for Green's function [TARGET][SOURCE]: [%d][%d].  Please fix and recompile.\n", i, j);
	endrun(1000);
      }
#endif

#ifndef NGRAVS_L3VIOLATION

      // Check that the forces are symmetric for unit mass
      // Short circuit on the diagonal
      if( i != j && ( (*AccelFxns[i][j])(1,1,0.5,3,1) != (*AccelFxns[j][i])(1,1,0.5,3,1))) {
	printf("ngravs: force between particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }

      if( i != j && ( (*AccelSplines[i][j])(1,1,0.5,3,1) != (*AccelSplines[j][i])(1,1,0.5,3,1))) {
	printf("ngravs: softened (spline) force between particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }

#if defined PERIODIC
      // Check that the lattice correction forces are symmetric for unit mass
      // Short circuit on the diagonal
      if( i != j && (LatticeForce[i][j] != LatticeForce[j][i])) {
	printf("ngravs: lattice correction force between particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }

      if( LatticePotential[i][j] != LatticePotential[j][i]) {
	
	printf("ngravs: lattice correction potential between particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }
#endif

#ifdef PMGRID
      // Check that Green's functions are symmetric
      // Short circuit on the diagonal
      if( i != j && ((*GreensFxns[i][j])(1,1,0.5,3,1) != (*GreensFxns[j][i])(1,1,0.5,3,1)) ) {
	printf("ngravs: Green's functions for particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }
#endif

#if defined OUTPUT_POTENTIAL || defined PMGRID
      // Check that the Potential functions are symmetric
      // Short circuit on the diagonal
      if( i != j && ((*PotentialFxns[i][j])(1,1,0.5,3,1) != (*PotentialFxns[j][i])(1,1,0.5,3,1))) {
	printf("ngravs: potential functions for particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }

      if( i != j && ((*PotentialSplines[i][j])(1,1,0.5,3,1) != (*PotentialSplines[j][i])(1,1,0.5,3,1)) ) {
	printf("ngravs: softened (spline) potential between particles which natively interact through gravities %d and %d is not symmetric.  Newton's 3rd law violated.  Stopping.\n", i, j);
	endrun(1000);
      }
#endif
#else
      if(!ThisTask)
	printf("ngravs: Newton's 3rd law violations *permitted*, no equality check was performed.\n");
#endif

    }
  }
}
