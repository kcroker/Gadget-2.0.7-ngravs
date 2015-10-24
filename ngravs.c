#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

#include <fftw.h>

#include "allvars.h"
#include "proto.h"

/*! \file ngravs.c
 *  \brief defines multiple types of gravitational interaction between particles 
 *
 * This file contains functions that determine
 * the differing force laws between differing particle types.  These functions are 
 * used to compute acceleration, potential, and softened contributions during the tree walk
 * in forcetree.c and the Fourier kernel for free-space boundary conditions.
 *
 * If you wish to add altered interactions between your 6 particle types in simulation,
 * this is where you define those altered interactions.  You must then 'wire' them manually
 * to the particle types, follow the examples below.
 *
 * Beware that GADGET-2 computes acceleration adjustements to a target due to sources
 * so if your force laws violate SEP, you must be careful about the ordering of source
 * and target.
 *
 */


/*! This function must be modified to point to your desired
 *  extensions to gravity.  It determines which force laws
 *  are used to compute interactions.
 */
void wire_grav_maps(void) {

  int i,j;

  // KC 8/11/14 Wiring
  //
  // NOTE: For all interaction functions
  //
  //       InteractionFunctions[TARGET][SOURCE]
  //       (i.e. InteractionFunctions[PASSIVE][ACTIVE])
  //
  // NOTE: NgravsNames[][] is used to index things used by the simulation
  // code, like lattice correction tables.  So please make a unique identifier!
  //

#if !defined NGRAVS_STOCK_TESTING && !defined NGRAVS_ACCUMULATOR_TESTING
  /////////////////////// WIRING FOR RESEARCH RUNS ///////////////////
  //
  // Here is where you wire for your research runs.  If you wish to verify force accuracy,
  // profile, or otherwise, enable the below defines in the Makefile.
  // Please see comments below in NEWTONIAN COMPARISONS for guidance in
  // defining the wiring appropriately.  There is some subtlty with respect to
  // source and target.
  //
  ////////////////////////////////////////////////////////////////////
  

  ///////////////// END WIRING FOR RESEARCH RUNS ///////////////////////
  
#elif defined NGRAVS_STOCK_TESTING && !defined NGRAVS_ACCUMULATOR_TESTING
  ////////////////////// WIRING FOR NEWTONIAN COMPARISONS //////////////
  //
  // This code automatically populates all gravitational types with
  // the Newtonian interaction.  This is for comparison against unmodified
  // Gadget-2 in consistency, force accuracy, and performance
  //
  //////////////////////////////////////////////////////////////////////

  for(i = 0; i < N_GRAVS; ++i) {
    for(j = 0; j < N_GRAVS; ++j) {

      //////////////////////// NAMES ////////////////////////////
      NgravsNames[i][j] = "Newton";
      AccelFxns[i][j] = newtonian; 
      AccelSplines[i][j] = plummer;

#if defined PERIODIC
      // KC 12/6/14
      // Note that if these are computing lattice forces from exotic objects
      // where scale parameters depend on the mass, the computed value of 
      // the correction should be such that:
      //
      //    (summed total mass in node) * (*LatticeForce[][])(...) = (best approximation)
      //
      // For Newtonian things, this is always the case.
      //
      LatticeForce[i][j] = ewald_force;
      LatticePotential[i][j] = ewald_psi;
      LatticeZero[i][j] = 2.8372975;
#endif

      ////////////////////////// MESH AND POTENTIAL //////////////////////////////
      // KC 11/25/14
      // Note that these are the k-space Greens Functions
      // for a source at the origin
#ifdef PMGRID
      GreensFxns[i][j] = pgdelta;
#endif

      // KC 11/25/14
      // Note that these are also the position space Greens functions
      // for a source at the origin
#if defined OUTPUTPOTENTIAL || defined PMGRID
      PotentialFxns[i][j] = newtonian_pot;
      PotentialSplines[i][j] = plummer_pot;

      // KC 1/26/15
      // This is used by the pm_nonperiodic when setting up the free space
      // fft kernel.  Nothing else should ever encounter it as other computations
      // will hit the softening length.
      PotentialZero[i][j] = -1 / (sqrt(M_PI) * (((double) ASMTH) / (2*PMGRID)));

#endif
    }
  }
  ////////////////////// END NEWTONIAN COMPARISONS WIRING ///////////////////
  
#elif !defined NGRAVS_STOCK_TESTING && defined NGRAVS_ACCUMULATOR_TESTING
  ///////////////// ACCUMULATOR TESTING RUNS ////////////////////////////
  //
  // This code specifically tests the N_\perp addition with a force law
  // (BAM) which violates SEP.  Here the N_\perp correction is an exact
  // correction.
  //
  ///////////////////////////////////////////////////////////////////////
  
  NgravsNames[0][0] = "Newton";
  NgravsNames[0][1] = "SourceBAM";
  NgravsNames[1][0] = "TargetBAM";
  NgravsNames[1][1] = "BAMBAM";

  AccelFxns[0][0] = newtonian;
  AccelSplines[0][0] = plummer;
  AccelFxns[0][1] = AccelSplines[0][1] = sourcebambaryon;
  AccelFxns[1][0] = AccelSplines[1][0] = sourcebaryonbam;
  AccelFxns[1][1] = AccelSplines[1][1] = bambam;

#if defined OUTPUTPOTENTIAL || defined PMGRID
  // These won't be used, as the simulation is non-periodic
  GreensFxns[0][0] = none;
  GreensFxns[0][1] = none;
  GreensFxns[1][0] = none;
  GreensFxns[1][1] = none;

  PotentialFxns[0][0] = newtonian_pot;
  PotentialSplines[0][0] = plummer_pot;
  PotentialFxns[0][1] = PotentialSplines[0][1] = sourcebambaryon_pot;
  PotentialFxns[1][0] = PotentialSplines[1][0] = sourcebaryonbam_pot;
  PotentialFxns[1][1] = PotentialSplines[1][1] = bambam_pot;
  
  // BAM has distinct limiting cases here.  These values correspond to 
  // unit masses and are used in computing the Greens function kernel for 
  // non-periodic fft
  PotentialZero[0][0] = -1 / (sqrt(M_PI) * (((double) ASMTH) / (2*PMGRID)));
  PotentialZero[0][1] = -8*BAM_EPSILON;
  PotentialZero[1][0] = -8*BAM_EPSILON;
  PotentialZero[1][1] = -4*BAM_EPSILON;
#endif

  //////////////// END ACCUMULATOR TESTING WIRING ///////////////////////////

#else
  printf("ngravs: unsupported testing options defined in the Makefile.  Cannot do (explicit) accumulator tests and Newtonian comparions at the same time");
  endrun(1000);
#endif

}


///////////////// BEGIN FOURIER INTEGRATION ROUTINES /////////////////////
// 
// The NGRAVS_TPM_N below gives samples in x-space at (1/2Ntab, 2/2Ntab, ..., (6Ntab + 3)/2Ntab)/2 
// The additional factor of 0.5 (2x as many samples in k) is required to critically sample!
// 
// Gadget-2 uses only a subset of these samples in position space (NTAB of them ;) ), but 
// taking many more is a good idea, especially because we need to be integrating, so that's good
//
//#define NTAB 1000  // AMPLITUDE AND WIDTH STABLE
//#define NGRAVS_TPM_N 2*(6*NTAB+3) // AMPLITUDE AND WIDTH STABLE.  Divided out the FFTW N!

double jTok(int j, double Z) {

  return 2*M_PI * j * NTAB / (double)(Z * NGRAVS_TPM_N);
}

double mTox(int m, double Z) {
  
  return m*Z/(double)NTAB;
}

//
// To keep it consistent with the other code
double fourierIntegrand(double k, gravity normKGreen, double Z) {
  
  double k2 = k*k;

  // Will probably be other facs here, not same as pm_periodic though because this is a 1D xform
  // 2\pi from phi integration, 1/\pi from inverse unnormalized Fourier transform.
  return (*normKGreen)(1, 1, k2, k, 1) * exp(-k2 * Z * Z) * 2; 
}

//
// Note that the k-space greens' functions seem to be dimensionless expressions
// with the Boxsize dimension explicitly reintroduced in the usual Gadget-2
// through fac, Asmth[0], etc. in pm_periodic.c
//
double newtonKGreen(double target, double source, double k2, double k, long N) {
  
  return 1.0; // Normalized Yukawa -> k2 / (k2 + 0.5);
}

void performConvolution(fftw_plan plan, gravity normKGreen, double Z, double *oRes, double *oResI) {
  
  fftw_complex in[NGRAVS_TPM_N], out[NGRAVS_TPM_N];
  int m,j;
  double sum;
 
  // Zero out all arrays first 
  for(j = 0; j < NGRAVS_TPM_N; ++j) {
    in[j].re = 0;
    in[j].im = 0;
  }

  // 1) FFTW needs this loaded in wonk order
  in[0].re = fourierIntegrand(jTok(0, Z), normKGreen, Z);
  for(j = 1; j < NGRAVS_TPM_N/2; ++j) {

    in[j].re = fourierIntegrand(jTok(j, Z), normKGreen, Z);
    in[NGRAVS_TPM_N - j].re = fourierIntegrand(jTok(j, Z), normKGreen, Z);
  }

  // 2) Xform
  fftw_one(plan, in, out);

  // 2.5) Do silly mapping (check ths for errors, jeeez)
  // Adjust normalization!
  //
  // Note that the normalization comes from the change of variables in the
  // non-unitary angular Fourier xform (wiki column 4)
  //
  // We also confirm that it is integrating effectively over symmetric limits!
  for(m = 0; m < NGRAVS_TPM_N; ++m)
    oRes[m] = out[m].re * NTAB / (Z*NGRAVS_TPM_N);

  // 3) Integrate the dumb way so we can keep the sampling interval the same
  //    in the integrated function. 
  sum = 0.0;
  oResI[0] = 0.0;
  for(m = 1; m < NGRAVS_TPM_N; ++m) {
    sum += (mTox(m, Z) - mTox(m-1, Z)) * 0.5 * (oRes[m-1] + oRes[m]);
    oResI[m] = sum;
  }
}

//
// Takes a Gadget-2 i and returns the appropriate m to give
// this value in the table:
//
//    m / 4Ntab = 3(i + 1/2)/Ntab
// So:
//    m = 12i + 6
//
int gadgetToTPM(int i) {

  return (12*i + 6);
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

      // Verify that the TreePM smoothing distance is not inside the tree softening distance!
      for(n = 1; n < 6; ++n) {
	
	// NTAB/asmthfac = NTAB / (0.5 / All.Asmth[0] * (NTAB/3))
	if(All.ForceSoftening[n] > 2*3*All.Asmth[0]) {

	  printf("ngravs: TreePM transition scale %f sits within the spline softened force radius %f for particle type %d.   This is not permitted.",
		 2*3*All.Asmth[0], All.ForceSoftening[n], n);
	  endrun(1034);
	}
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

/*! This is no gravity.  It returns 0.0 regardless of input
 */
double none(double target, double source, double h, double r, long N){ 

  return 0.0;
}

// KC 10/18/14
//
// WARNING: Note that, presumably to eliminate a plethora of unnecessary negations,
//          Gadget-2 works with the positive of the acceleration.  
//
// ALL ACCELERATION SIGNS APPEARING HERE ARE TO BE INVERTED FROM WHAT YOU NORMALLY WOULD WRITE
//
// WARNING: These give the normalized accelerations!  They include an additional factor of 1/r
//          so that taking the vector magnitude of this function's return value times
//          the dx_i gives the correct force.
//
// NOTE OPTIMIZATION: since an AccelFxn does not use h as a softening, we pass in
//                    the r^2 (since it is already computed), so that we only have 
//                    to ever perform multiplications.
//  
/*! This is Newtonian gravity, and is the usual baryon-baryon interaction 
 */
double newtonian(double target, double source, double h, double r, long N) {

  // Note newtonian does not violate SEP
  return source / (h * r);
} 

/*! This is **inverted** Newtonian gravity, for use in the Hohmann & Wolfarth scenario
 */
double neg_newtonian(double target, double source, double h, double r, long N) {

  // Note newtonian does not violate SEP
  return -source / (h * r);
} 

/*! This is the usual Newtonian gravitational potential
 *
 */
double newtonian_pot(double target, double source, double h, double r, long N) {

  return source / r;
}

/*! This is the **inverted** Newtonian potential
 */
double neg_newtonian_pot(double target, double source, double h, double r, long N) {

  return -source / r;
}

// KC 10/18/14
// NOTE: Green's functions are not inverted from their usual sign
// NOTE 2: The factor of 4\pi is missing, and G is 1 in internal Gadget-2 units
//

/*! This is the box periodic NORMALIZED Green's function for a point source of unit mass
 */
double pgdelta(double target, double source, double k2, double k, long N) {

  return 1.0;
}

/*! This is the **inverted** box periodic Green's function for a point source of unit mass, 
 *  for use in the Hohmann & Wolfarth scenario
 */
double neg_pgdelta(double target, double source, double k2, double k, long N) {

  return -1.0;
}

/*! This is the Plummer spline used by GADGET-2
 */
double plummer(double target, double source, double h, double r, long N) {

  double h_inv;

  h_inv = 1/h;
  
  r *= h_inv;
  if(r < 0.5)
    return source * h_inv * h_inv * h_inv * 
      (10.666666666667 + r * r * (32.0 * r - 38.4));
  else
    return source * h_inv * h_inv * h_inv * 
      (21.333333333333 - 48.0 * r +
       38.4 * r * r - 10.666666666667 * r * r * r - 0.066666666667 / (r * r * r));
}

/*! This is the "inverted" Plummer spline for use in the Hohmann & Wolfarth scenario
 */
double neg_plummer(double target, double source, double h, double r, long N) {

  double h_inv;

  // KC 10/26/14
  // It remains a question whether calling a fxn with 5 things on the stack is going to be 
  // slower than just multiplying things out every time and calling something with 3 things
  h_inv = 1/h;
  
  r *= h_inv;
  if(r < 0.5)
    return -source * h_inv * h_inv * h_inv * 
      (10.666666666667 + r * r * (32.0 * r - 38.4));
  else
    return -source * h_inv * h_inv * h_inv * 
      (21.333333333333 - 48.0 * r +
       38.4 * r * r - 10.666666666667 * r * r * r - 0.066666666667 / (r * r * r));
}

/*! This is the potential of the Plummer spline used by GADGET-2
 */
double plummer_pot(double target, double source, double h, double r, long N) {
  
  double h_inv;
  
  h_inv = 1/h;
  r *= h_inv;

  if(r < 0.5)
    return source * h_inv * (-2.8 + r * r * (5.333333333333 + r * r * (6.4 * r - 9.6)));
  else
    return source * h_inv * 
      (-3.2 + 0.066666666667 / r + r * r * (10.666666666667 +
					    r * (-16.0 + r * (9.6 - 2.133333333333 * r))));
}

/*! This is the negative potential of the Plummer spline, for use in Hohmann & Wolfarth scenario
 */
double neg_plummer_pot(double target, double source, double h, double r, long N) {
  
  double h_inv;
  
  h_inv = 1/h;
  r *= h_inv;

  if(r < 0.5)
    return -source * h_inv * (-2.8 + r * r * (5.333333333333 + r * r * (6.4 * r - 9.6)));
  else
    return -source * h_inv * 
      (-3.2 + 0.066666666667 / r + r * r * (10.666666666667 +
					    r * (-16.0 + r * (9.6 - 2.133333333333 * r))));
}

/*! This is Newtonian gravity with a Yukawa modification
 * c.f. Shirata, Yoshida, et. al 2005
 */
double newyukawa(double target, double source, double h, double r, long N) {

#define YUKAWA_ALPHA 1
#define YUKAWA_LAMBDA_INV 1e-2

  return source / (r * h) * -YUKAWA_ALPHA * expm1f(-r*YUKAWA_LAMBDA_INV);
}

/*! This is the BAM-BAM interaction
 * c.f. http://arxiv.org/abs/1408.2702
 */
double bambam(double target, double source, double h, double r, long N) {
  
  // Note apparent SEP violation
  // Note naturally softened
  // Note adjustment of the internal scale by N.  Thus the scale is determined by the average mass content of the cell.
  // In the case where the all BAM halos have the same mass parameter, this correction is the *exact* correction.

  double eta, rho;
  double eta3;
  double reta, reta2;

  eta = 4.0*M_PI*BAM_EPSILON/(target+source/N);
  rho = 2*target*source/M_PI;

  reta = r * eta;
  reta2 = reta * reta;
  eta3 = eta * eta * eta;
 
  if(reta < 0.1) {
    
    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    // the differentiate termwise to get the force
    // also multiply by -1: F \def -\grad pot
    // also must divide by an additional 1/r to give the unit vector in code!
    //
    return rho * eta3 * (2.0/3.0 - 4.0*reta2/5.0 + 6.0*reta2*reta2/7.0);
  }
  else
    return rho * eta3 * (atan(reta)/(reta2*reta) - 1.0/(reta2*(1+reta2)));
}

/*! This is the BAM-Baryon interaction sourced by a BAM.
 * Note that here the target is a baryon!!
 * The force laws are necessarily symmetric, but the computation GADGET-2 uses
 * is not the force, but the adjustment to the acceleration of the target.  
 * So you have to be careful here. 
 */
double sourcebambaryon(double target, double source, double h, double r, long N) {

  double eta, rho;
  double eta3;
  double reta, reta2;

  rho = 2*target*source/M_PI;
  eta = 4.0*M_PI*BAM_EPSILON*N/source;
  
  reta = r * eta;
  reta2 = reta * reta;
  eta3 = eta * eta * eta;
 
  if(reta < 0.1) {
    
    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    // the differentiate termwise to get the force
    // also multiply by -1: F \def -\grad pot
    // also must divide by an additional 1/r to give the unit vector in code!
    //
    return rho * eta3 * (2.0/3.0 - 4.0*reta2/5.0 + 6.0*reta2*reta2/7.0);
  }
  else
    return rho * eta3 * (atan(reta)/(reta2*reta) - 1.0/(reta2*(1+reta2)));
}

/*! This is the BAM-Baryon interaction sourced by a baryon.
  Note that here the target is a BAM!!
  The force laws are necessarily symmetric, but the computation GADGET-2 uses
  is not the force, but the adjustment to the acceleration of the target.  
  So you have to be careful here.
 */
double sourcebaryonbam(double target, double source, double h, double r, long N) {

  // Note apparent SEP violation
  // Note naturally softened
  double eta, rho;
  double eta3;
  double reta, reta2;

  eta = 4.0*M_PI*BAM_EPSILON/target;
  rho = 2*target*source/M_PI;

  reta = r * eta;
  reta2 = reta * reta;
  eta3 = eta * eta * eta;
 
  if(reta < 0.1) {
    
    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    // the differentiate termwise to get the force
    // also multiply by -1: F \def -\grad pot
    // also must divide by an additional 1/r to give the unit vector in code!
    //
    return rho * eta3 * (2.0/3.0 - 4.0*reta2/5.0 + 6.0*reta2*reta2/7.0);
  }
  else
    return rho * eta3 * (atan(reta)/(reta2*reta) - 1.0/(reta2*(1+reta2)));
}

/*! This is the BAM-BAM potential (or free-space Greens fxn in position representation)
 * c.f. http://arxiv.org/abs/1408.2702
 */
double bambam_pot(double target, double source, double h, double r, long N) {

  // Use a 7th order Taylor polynomial if tan(x) x < 1/10
  // The error at x = 1/10 for tan(x) is then < 10^-7
  // Its not beyond double precision, but its almost at float.
  //
  // This will be kinda slow...
  //
  double eta;
  double rho;
  double reta, reta2, reta4;

  rho = 2*target*source/M_PI;

  eta = 4.0*M_PI*BAM_EPSILON/(target+source/N);
  reta = r * eta;
  reta2 = reta * reta;
  reta4 = reta2 * reta2;

  if(reta < 0.1) {

    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    //
    return rho * eta*(1 - reta2/3.0 + reta4/5.0 - reta2*reta4/7.0);
  }
  else
    return rho * atan(reta)/r;
}

// Now note that the advantage of the simple taylor series is that it is easy to compute the 
// derivative and turn that into a force.

/*! This is the BAM-Baryon potential.  Here a baryon acts as a source.
 */
double sourcebaryonbam_pot(double target, double source, double h, double r, long N) {

  double eta;
  double rho;
  double reta, reta2, reta4;
  
  rho = 2*target*source/M_PI;
  eta = 4.0*BAM_EPSILON*M_PI*N/target;
  reta = r * eta;
  reta2 = reta * reta;
  reta4 = reta2 * reta2;

  if(reta < 0.1) {

    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    //
    return rho * eta*(1 - reta2/3.0 + reta4/5.0 - reta2*reta4/7.0);
  }
  else
    return rho * atan(reta)/r;
}

/*! This is the BAM-Baryon potential.  Here a BAM acts as a source.
 */
double sourcebambaryon_pot(double target, double source, double h, double r, long N) {

  double eta;
  double rho;
  double reta, reta2, reta4;
  
  rho = 2*target*source/M_PI;
  eta = 4.0*BAM_EPSILON*M_PI*N/source;
  reta = r * eta;
  reta2 = reta * reta;
  reta4 = reta2 * reta2;

  if(reta < 0.1) {

    // orig taylor: reta - (reta)^3/3 + (reta)^5/5 - (reta)^7/7
    // we want to divide out by the r
    // so: rho*(eta - r^2eta^3/3 + r^4eta^5/5 - r^6eta^7/7)
    //
    return rho * eta*(1 - reta2/3.0 + reta4/5.0 - reta2*reta4/7.0);
  }
  else
    return rho * atan(reta)/r;
}
