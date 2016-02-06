#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
//#include <sys/types.h>
//#include <unistd.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_math.h>
#include "ngravs.h"
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
 * Some notes:
 * ------------------
 * Newtonian is specified completely.  
 * Yukawa is specified completely, lattice potential computation
 *  has not been thoroughly tested (because we didn't need it).  The Yukawa Madelung computation is incomplete
 *  and has been commented out, though completion should be straightforward following the referenced work.
 * The supermacho BAM scenario has been implemented non-periodically.
 */

#define YUKAWA_ALPHA 1
// Make sure to use a decimal here
#ifndef YUKAWA_IMASS
#define YUKAWA_IMASS 60.0
#endif
//
// 1/YUKAWA_IMASS sets the suppression *length* scale wrt 1/2 box length
//
// Examples: 0.5 gives 1/e suppression at 4 half-boxlengths out
//           2 gives 1/e suppression at 1 half-boxlength 
//           10 gives (1/e)^5 suppression at 1 half-boxlength
//           24 gives (1/e)^12 at 1 half-boxlength
//
// (These units are required for an interpolation that is invariant to boxlength)
//

/*! This function must be modified to point to your desired
 *  extensions to gravity.  It determines which force laws
 *  are used to compute interactions.
 */
void wire_grav_maps(void) {

  int i,j;
#ifdef NGRAVS_YUKAWA_FORCETEST
  char *fname;
#endif

  // KC 8/11/14 Wiring
  //
  // NOTE: For all interaction functions
  //
  //       InteractionFunctions[TARGET][SOURCE]
  //       (i.e. InteractionFunctions[PASSIVE][ACTIVE])
  //
  // VERY IMPORTANT: NgravsNames[][] is used to index things used by the simulation
  // code, like lattice correction tables.  So please make a unique identifier!
  // (It can also be used eventually to save memory and startup-time by not 
  // computing redundant tables.)
  //

#if !defined NGRAVS_STOCK_TESTING && !defined NGRAVS_ACCUMULATOR_TESTING && !defined NGRAVS_YUKAWA_FORCETEST && !defined NGRAVS_COMBINED_TESTING_UNIFORM
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
  
#elif defined NGRAVS_STOCK_TESTING
  printf("ngravs: wired in stock comparison test mode\n");

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
      // WARNING: The distance scale in the LatticeForce computation is dimensionless in terms
      //          of the a box grid with EN (see forcetree.c) number of points.  The actual force
      //          is interpolated from this force trilinearly.
      //
      LatticeForce[i][j] = ewald_force;
      LatticePotential[i][j] = ewald_psi;
      LatticeZero[i][j] = 2.8372975;
#endif

      ////////////////////////// MESH AND POTENTIAL //////////////////////////////
      // KC 11/25/14
      // Note that these are the periodic k-space Greens Functions
      // for a source at the origin
#ifdef PMGRID
      GreensFxns[i][j] = pgdelta;
      NormedGreensFxns[i][j] = normed_pgdelta;
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

#elif defined NGRAVS_ACCUMULATOR_TESTING
  printf("ngravs: wired in accumulator test mode\n");
  
  ///////////////// ACCUMULATOR TESTING RUNS ////////////////////////////
  //
  // This code specifically tests the N_\perp addition with a force law
  // (BAM) where the N_\perp correction is an exact
  // correction.
  //
  ///////////////////////////////////////////////////////////////////////
  
  NgravsNames[0][0] = "Newton";
  NgravsNames[0][1] = "SourceBAM";
  NgravsNames[1][0] = "TargetBAM";
  NgravsNames[1][1] = "BAMBAM";

  AccelFxns[0][0] = newtonian;
  AccelFxns[0][1] = sourcebambaryon;
  AccelFxns[1][0] = sourcebaryonbam;
  AccelFxns[1][1] = bambam;

  AccelSplines[0][0] = plummer;
  AccelSplines[0][1] = sourcebambaryon_spline; 
  AccelSplines[1][0] = sourcebaryonbam_spline;
  AccelSplines[1][1] = bambam_spline; 

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

#elif defined NGRAVS_YUKAWA_FORCETEST
  printf("ngravs: wired in yukawa force test mode\n");
  //////////////// BEGIN YUKAWA FORCE TESTING WIRING ///////////////////////////
  //
  // This code is used to examine the force accuracy during the TreePM transition
  // for a more general force, in this case Yukawa, which is a pathological edge case
  // due to the r-space decay rate being exponential. 
  //
  ///////////////////////////////////////////////////////////////////////

  for(i = 0; i < N_GRAVS; ++i) {
    for(j = 0; j < N_GRAVS; ++j) {
      
      // Allocate a new one each time, 
      // because that's what we'd have to be doing anyway.
      // This is not a leak because we need these handles throughout the
      // entire program run.
      fname = (char *)malloc(128);
 
      if(i == j) {
	snprintf(fname, 128, "None");
	AccelFxns[i][j] = none;
	AccelSplines[i][j] = none;
      }
      else {
	snprintf(fname, 128, "Yukawa_%e", YUKAWA_IMASS);

	// We set the Yukawa spline to plummer since
	// the force is Newtonian at small r
	// This being incorrect won't matter for force checking the TreePM 
	// stuff, because the force correction uses the spline too.
      	AccelSplines[i][j] = plummer;
	AccelFxns[i][j] = yukawa;
      }
      
      NgravsNames[i][j] = fname;

#if defined PERIODIC
      // Computed from G. Salin and J.M. Caillol
      // J. Chem. Phys., Vol 113, No. 23, 2000
      if(i != j)
	LatticeForce[i][j] = yukawa_lattice_force;
      else
	LatticeForce[i][j] = lattice_force_none;

      LatticePotential[i][j] = yukawa_lattice_psi;
      LatticeZero[i][j] = yukawa_madelung(YUKAWA_IMASS);
      if(!ThisTask)
	printf("ngravs: Yukawa force Madelung constant for [%d][%d] = %f\n", i, j, LatticeZero[i][j]);
#endif

#if defined OUTPUTPOTENTIAL || defined PMGRID
      
      if(i != j) {
	GreensFxns[i][j] = pgyukawa;
	NormedGreensFxns[i][j] = normed_pgyukawa;
      }
      else {
	GreensFxns[i][j] = none;
	NormedGreensFxns[i][j] = none;
      }
      // We don't care about the potentials because we're
      // not doing non-periodic or gastrophysics in this test
      PotentialFxns[i][j] = none;
      PotentialSplines[i][j] = none;
      PotentialFxns[i][j] = none;
#endif
    }
  }
  //////////////////////// END GENERALIZED FORCE TEST WIRING //////////////////

#elif defined NGRAVS_COMBINED_TESTING_UNIFORM
  printf("ngravs: wired in combined force test mode\n");
  //////////////// BEGIN COMBINED FORCE TESTING WIRING ///////////////////////////
  //
  // This code is used to examine the force accuracy during the TreePM transition
  // for a sum of forces.  This seemed prudent because of the Yukawa tweak 
  // required to get correct behaviour. 
  // 
  ///////////////////////////////////////////////////////////////////////

  for(i = 0; i < N_GRAVS; ++i) {
    for(j = 0; j < N_GRAVS; ++j) {
      
      fname = (char *)malloc(128);
      snprintf(fname, 128, "ColoYuk_%e", YUKAWA_IMASS);
 
      NgravsNames[i][j] = fname;
      AccelFxns[i][j] = coloyuk;
      AccelSplines[i][j] = plummer;

#if defined PERIODIC
      LatticeForce[i][j] = coloyuk_lattice_force;
      LatticePotential[i][j] = ewald_psi; // It doesn't matter, not being used in the tests
      LatticeZero[i][j] = yukawa_madelung(YUKAWA_IMASS) + 2.8372975;
      if(!ThisTask)
	printf("ngravs: Yukawa force Madelung constant for [%d][%d] = %f\n", i, j, LatticeZero[i][j]);
#endif

#if defined OUTPUTPOTENTIAL || defined PMGRID
      GreensFxns[i][j] = pgcoloyuk;
      NormedGreensFxns[i][j] = normed_pgcoloyuk;
      PotentialFxns[i][j] = none;
      PotentialSplines[i][j] = none;
      PotentialFxns[i][j] = none;
#endif

    }
  }
#else
  printf("ngravs: unsupported testing options defined in the Makefile.  Cannot do (explicit) accumulator tests and Newtonian comparions at the same time");
  endrun(1000);
#endif
}

///////////////////// BEGIN GENERALIZED FORCE AND GREENS FUNCTIONS /////////
//
// KC 10/18/14
//
// WARNING: Note that, presumably to eliminate a plethora of unnecessary negations,
//          Gadget-2 works with the positive of the acceleration.  
//
// ALL ACCELERATION *SIGNS* APPEARING HERE ARE TO BE INVERTED FROM WHAT YOU NORMALLY WOULD WRITE
//
// NOTE OPTIMIZATION: since an AccelFxn does not use h as a softening, we pass in
//                    the r^2 (since it is already computed), so that we can perform 
//                    fewer multiplications
//  

/*! This is no gravity.  It returns 0.0 regardless of input
 */
double none(double target, double source, double h, double r, long N){ 

  return 0.0;
}

/*! This is Newtonian gravity, and is the usual baryon-baryon interaction 
 */
double newtonian(double target, double source, double h, double r, long N) {

  // Note newtonian does not violate SEP
  return source / h;
} 

/*! This is **inverted** Newtonian gravity, for use in the Hohmann & Wolfarth scenario
 */
double neg_newtonian(double target, double source, double h, double r, long N) {

  // Note newtonian does not violate SEP
  return -source / h;
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

// KC 10/30/15
//
// The k that gadget uses is dimensionless between [-PMGRID/2, PMGRID/2].  The original form 
// plugged into the convolution is also this dimensionless form.  So, your Greens function will need to 
// be dimensionless.  The length scale is All.BoxSize.  
//
// FACTORS ARE SUCH THAT: 4\pi G/k^2 \becomes 1.0

/*! This is the box periodic NORMALIZED Green's function for a point source of unit mass
 */
double pgdelta(double target, double source, double k2, double k, long N) {

  // Return the NAN, since we either never compute it, or we catch it

  // KC 2/5/16 
  // PPP
  // (This might suck for speed, since we will raise floating point exceptions...)
  return 1.0/k2;
}

double normed_pgdelta(double target, double source, double k2, double k, long N) {

  return 1.0;
}

/*! This is the **inverted** box periodic Green's function for a point source of unit mass, 
 *  for use in the Hohmann & Wolfarth scenario
 */
double neg_pgdelta(double target, double source, double k2, double k, long N) {

  return -1.0/k2;
}

/*! This is the Plummer spline used by GADGET-2
 */
// 
// WARNING: Acceleration splines contain an additional factor of 1/r (or 1/h) 
//          as this division is not carried out for splined forces in forcetree.c
//          This is because splines need to divide by the softening scale instead
//          of the radius when computing their forces.
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
    
    
    // KC 11/3/15
    // r put back in because forcetree.c divides it out now
    return rho * eta3 * (2.0*r/3.0 - 4.0*reta2*r/5.0 + 6.0*reta2*reta2*r/7.0);
  }
  else
    // KC 11/3/15 - corrected radial factor
    return rho * eta3 * (atan(reta)/(reta2*eta) - 1.0/(reta*eta*(1+reta2)));
}

double bambam_spline(double target, double source, double h, double r, long N) {
  
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

  // KC 11/3/15 - corrected radial factor
  if(reta < 0.1)
    return rho * eta3 * (2.0/3.0 - 4.0*reta2/5.0 + 6.0*reta2*reta2/7.0);
  else
    return rho * eta3 * (atan(reta)/(reta2*reta) - 1.0/(reta2*(1+reta2)));
}

/*! This is the BAM-Baryon interaction sourced by a BAM.
 * Note that here the target is a baryon!!
 * The force laws are necessarily symmetric, but the computation GADGET-2 uses
 * is not the force, but the adjustment to the acceleration of the target.  
 * So you have to be careful here. 
 */
double sourcebambaryon_spline(double target, double source, double h, double r, long N) {

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

double sourcebambaryon(double target, double source, double h, double r, long N) {

  double eta, rho;
  double eta3;
  double reta, reta2;

  rho = 2*target*source/M_PI;
  eta = 4.0*M_PI*BAM_EPSILON*N/source;
  
  reta = r * eta;
  reta2 = reta * reta;
  eta3 = eta * eta * eta;

  // KC 11/3/15 - corrected radial factor
  if(reta < 0.1)
    return rho * eta3 * (2.0*r/3.0 - 4.0*reta2*r/5.0 + 6.0*reta2*reta2*r/7.0);
  else
    return rho * eta3 * (atan(reta)/(reta2*eta) - 1.0/(reta*eta*(1+reta2)));
}

/*! This is the BAM-Baryon interaction sourced by a baryon.
  Note that here the target is a BAM!!
  The force laws are necessarily symmetric, but the computation GADGET-2 uses
  is not the force, but the adjustment to the acceleration of the target.  
  So you have to be careful here.
 */
double sourcebaryonbam_spline(double target, double source, double h, double r, long N) {

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
 
  // KC 11/3/15
  // Adjusted radial factor
  if(reta < 0.1)
    return rho * eta3 * (2.0*r/3.0 - 4.0*reta2*r/5.0 + 6.0*reta2*reta2*r/7.0);
  else
    return rho * eta3 * (atan(reta)/(reta2*eta) - 1.0/(reta*eta*(1+reta2)));
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

/*! This function computes the potential correction term by means of Ewald
  *  summation.  Newtonian potential!
 */
double ewald_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int i, n[3], h[3], h2;

  // KC 11/16/15
  // We will figure out the mappings between the variables used here
  // and those in J. Chem. Phys., Vol. 113, No. 23, 2000, Eqn. (3.1)
  // when *their* \alpha \to 0
  //

  alpha = 2.0;

  for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];
	  
	  // r (here) = r* (there)
	  // \alpha (here) = \beta * (there)

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  sum1 += erfc(alpha * r) / r;

	  // Residual distinctions:
	  // None!
	}

  for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  
	  // hdotx (here) = n \dot \r* (there)
	  // other mappings the same!
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
	  if(h2 > 0)
	    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);

	  // Residual distinctions:
	  // None!
	}

  r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

  // Note the embedded neutralizing background
  // The mapping holds, with no residual factors! 
  psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;

  return psi;
}

/* Essential notes on units:
 * ----------------------------------
 * Units to k, k2 in Greens' Functions: k \in [-PMGRID/2, PMGRID/2] (mesh cells)
 * Units to r, r2 in Spline and Accel Functions: r \in [0, BoxLength] or unconstrained (given length unit)
 * Units to x in Lattice functions: x \in [0, 0.5] (fractions of the total side length in one octant)
 *
 */

 double coloyuk(double target, double source, double h, double r, long N) {
   return yukawa(target, source, h, r, N) + newtonian(target, source, h, r, N);
 }

#ifdef PMGRID
 double pgcoloyuk(double target, double source, double k2, double k, long N) {
   return pgyukawa(target, source, k2, k, N) + pgdelta(target, source, k2, k, N);
 }

 double normed_pgcoloyuk(double target, double source, double k2, double k, long N) {

   return normed_pgyukawa(target, source, k2, k, N) + normed_pgdelta(target, source, k2, k, N);
 }
#endif

 void coloyuk_lattice_force(int iii, int jjj, int kkk, double x[3], double force[3]) {
   
   double tmp[3];

   yukawa_lattice_force(iii, jjj, kkk, x, tmp);
   ewald_force(iii, jjj, kkk, x, force);

   for(iii = 0; iii < 3; ++iii)
     force[iii] += tmp[iii];

 }
 
 /*! A pure Yukawa force
  *
  */
 double yukawa(double target, double source, double h, double r, long N) {
  
  double ym;
  ym = YUKAWA_IMASS/All.BoxSize;  
  return source * exp(-r*ym) * (ym/r + 1.0/h);
}

/*! A periodic yukawa k-space Greens function, normalized by the Newtonian interaction
 *  NOTE: k is supplied dimensionlessly in terms of PMGRID so k \in [-PMGRID/2, PMGRID/2]
 */
#ifdef PMGRID
// This #ifdef is required because the structures for tracking the smoothing length
// are not allocated unless running in periodic mode.
double pgyukawa(double target, double source, double k2, double k, long N) {
  
  double ym = YUKAWA_IMASS/(2*M_PI);
  double asmth2;

  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  return 1.0 / (k2 + ym*ym) * exp(-ym*ym*asmth2);
}

double normed_pgyukawa(double target, double source, double k2, double k, long N) {

  // This converts from PMGRID units into shortrange interpolation table units
  double ym = gridKtoNormK(YUKAWA_IMASS/(2*M_PI));
  return k2 / (k2 + ym*ym) * exp(-ym*ym*0.25);
}
#endif

/*! This function computes the Madelung constant for the yukawa potential
 * which depends on the box length interestingly...
 * We follow Eqn (2.19) of G. Salin and Caillol (op. cit)
 *
 * Note that we use the values de-dedimensionalized in the same way
 * as the other yukawa lattice functions.
 *
 */
double yukawa_madelung(double ym) {
  
  /* double sum1, sum2, sum3; */
  /* double k2, m; */
  /* double alpha; */
  /* int n[3]; */

  /* // KC 11/16/15 */
  /* // We again adopt the same notation as that used in Gadget-2, so their beta */
  /* // is our alpha, etc (see comments above) */
  /* // */
  /* alpha = 5.64; */

  /* // Going out to the same distance seems like a good idea, no? */
  /* for(n[0] = -5, sum1 = 0, sum2 = 0; n[0] <= 5; n[0]++) { */
  /*   for(n[1] = -5; n[1] <= 5; n[1]++) { */
  /*     for(n[2] = -5; n[2] <= 5; n[2]++) { */

  /* 	// Here we use n for both the k sum and the n sum because they are both dimensionless */
  /* 	k2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2]; */
	
  /* 	if(k2 > 0) { */
  /* 	  m = sqrt(k2); */
  /* 	  k2 *= 4*M_PI*M_PI; */
  /* 	  sum1 += exp(-(k2 + ym*ym)/(4*alpha*alpha))/(k2 + ym*ym); */
  /* 	  sum2 += (erfc(alpha*m + ym/(2*alpha))*exp(ym*m) + erfc(alpha*m-ym/(2*alpha))*exp(-ym*m))/(2*m); */
  /* 	} */
  /* 	else { */
	  
  /* 	  // XXX Need to explicitly take the limit!  erfc(m)/m */
  /* 	  // will have a finite value. */

  /* 	} */
  /*     } */
  /*   } */
  /* } */
  
  /* // The non-summation terms */
  /* sum3 = -2*alpha/sqrt(M_PI)*exp(-ym*ym/(4*alpha*alpha)) + */
  /*   ym*erfc(ym/(2*alpha)); */
  
  /* // Explicitly the zero yukawa-mass case, limit via l'Hopital */
  /* if(ym > 0) */
  /*   sum3 += 4*M_PI/(ym*ym) * expm1(ym*ym/(4*alpha*alpha)); */
  /* else */
  /*   sum3 += 4*M_PI/(4*alpha*alpha); */
  
  /* return (4*M_PI*sum1 + sum2 + sum3); */

  // KC 2/5/16
  // XXX commented out because the above implementation needs to be finished
  // (if you want to use it, which you probably don't!)
  return 0;
}

/*! This function computes the potential correction term by means of Ewald
  *  summation, adjusted for Yukawa! Thanks Salin and Caillol!
 */
double yukawa_lattice_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int i, n[3], h[3], h2;

  // KC 11/16/15
  // We will figure out the mappings between the variables used here
  // and those in J. Chem. Phys., Vol. 113, No. 23, 2000, Eqn. (3.1)
  // when *their* \alpha \to 0
  
  // KC 11/16/15
  // Notice we use Salin's ideal transition of 5.64 with summations out to |n| = 5
  // So we *really* overcompute it!
  alpha = 5.64;

  for(n[0] = -5, sum1 = 0; n[0] <= 5; n[0]++)
    for(n[1] = -5; n[1] <= 5; n[1]++)
      for(n[2] = -5; n[2] <= 5; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];
	  
	  // r (here) = r* (there)
	  // \alpha (here) = \beta * (there)

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  sum1 += (erfc(alpha * r + YUKAWA_IMASS/(2*alpha)) * exp(YUKAWA_IMASS * r)) / (2*r);
	  sum1 += (erfc(alpha * r - YUKAWA_IMASS/(2*alpha)) * exp(-YUKAWA_IMASS * r)) / (2*r);

	  // Residual distinctions:
	  // None!
	}

  for(h[0] = -5, sum2 = 0; h[0] <= 5; h[0]++)
    for(h[1] = -5; h[1] <= 5; h[1]++)
      for(h[2] = -5; h[2] <= 5; h[2]++)
	{
	  
	  // hdotx (here) = n \dot \r* (there)
	  // other mappings the same!
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
	  if(h2 > 0)
	    sum2 += 1 / (M_PI * h2 + YUKAWA_IMASS*YUKAWA_IMASS/(4*M_PI)) * exp(-M_PI * M_PI * h2 / (alpha * alpha) - YUKAWA_IMASS*YUKAWA_IMASS/(4*alpha*alpha)) * cos(2 * M_PI * hdotx);

	  // Residual distinctions:
	  // None!
	}

  r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

  // Note the embedded neutralizing background
  // The mapping holds, with no residual factors! Eqn. (2.18)
  // Note inverted sign!
  psi = M_PI / (alpha * alpha) - sum1 - sum2 + exp(-YUKAWA_IMASS*r)/r;

  return psi;
}

// Now we have to take the derivative of the above function...
/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void yukawa_lattice_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
  double alpha, r2;
  double r, val, hdotx, dx[3];
  int i, h[3], n[3], h2;
  double ym;
  
  // KC 11/16/15
  // Note our use of Salin's optimal 'alpha', and our excessive momentum-space
  alpha = 5.64;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  //
  //
  // This should never have been called r, but instead u.
  // Its takes values in [0, 0.5], and moves within the corner
  // octant of one unit cell.
  //

  // KC 11/16/15
  // Here we add in the original force, with the displacement vector
  // normalization included. 
  //
  // NOTE:
  //   r \in [0, sqrt(3)*0.5] (indexing possible top octant values)
  //   Below is what was required to make the computed force with a call to yukawa()
  //   equal to that computed here, with the distinct units on each.
  //
  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  r = sqrt(r2);
  ym = YUKAWA_IMASS;
  // r is in dimensionless octant coordinates [0, 1/2]
  for(i = 0; i < 3; i++)
    force[i] = exp(-r*ym) * (ym + 1.0/r) * x[i]/r2; 
  
#ifndef NGRAVS_DEBUG_UNITS_CHECK
  // KC 12/4/14
  // Looks like this takes the first four images out in position space in each direction (so
  // bracketing by 8 overall)
  
  for(n[0] = -5; n[0] <= 5; n[0]++)
    for(n[1] = -5; n[1] <= 5; n[1]++)
      for(n[2] = -5; n[2] <= 5; n[2]++)
  	{
  	  for(i = 0; i < 3; i++)
  	    dx[i] = x[i] - n[i];

  	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  
  	  // Note, as YUKAWA_IMASS \to zero, we regenerate the Ewald for Coloumb
  	  //	  val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);
  	  
	  // TYPE I
	  // 0.5*(A + B) eventually /r^3
	  val = 0.5*( exp(ym*r)*erfc(alpha*r + ym/(2*alpha)) +
		      exp(-ym*r)*erfc(alpha*r - ym/(2*alpha)));

	  // KC 1/7/16
	  // For r > rcut, approx = 0.5*exp(ym*r)*erfc(alpha*r + ym/(2*alpha))
	  // Overall, it will carry a factor of 1/r^2.

	  // 0.5*(A + C*r) eventually /r^3
	  /* val = 0.5*gsl_sf_erfc(alpha*r - ym/(2*alpha))*exp(-ym*r)*(1.0/r + ym);  */
	  /* val += 0.5*gsl_sf_erfc(alpha*r + ym/(2*alpha))*exp(ym*r)*(1.0/r - ym); */
	  
	  // 0.5*(B + D*r) eventually /r^3
	  
	  // TYPE I
  	  for(i = 0; i < 3; i++)
  	    force[i] -= dx[i] / (r * r * r) * val;

	  /* for(i = 0; i < 3; i++) */
  	  /*   force[i] -= dx[i] / (r * r) * val; */

	  // TYPE II
	  // Now E
	  /* val += 2*alpha*exp(-alpha*alpha*r*r-ym*ym/(4*alpha*alpha))/sqrt(M_PI); */

	  // 0.5*ym*(C +1 D) + E eventually /r^2
	  // 
	  // The subtraction here could destroy accuracy?
  	  val = 0.5*ym*(-exp(ym*r)*erfc(alpha*r + ym/(2*alpha)) +
			exp(-ym*r)*erfc(alpha*r - ym/(2*alpha))) +
  	    2*alpha*exp(-alpha*alpha*r*r-ym*ym/(4*alpha*alpha))/sqrt(M_PI);

	  // KC 1/7/16
	  // For r > rcut, approx = -0.5*ym*exp(ym*r)*erfc(alpha*r + ym/(2*alpha)) + 2*alpha*exp(-ym*ym/(4*alpha*alpha))/sqrt(M_PI)
	  // eventually /r

  	  // KC 11/16/110
  	  // Note that these terms enter with one less radial power in the denominator
  	  for(i = 0; i < 3; i++)
  	    force[i] -= dx[i] / (r * r) * val;
  	}

  // KC 12/4/14
  // Looks like this takes the first four images in momentum space in each diretion (again
  // bracketing by 8 overall)
  //
  // KC 11/16/110
  // Take careful note of the relative signs!!
  // If you take the negative grad_r, then the signs are the same!
  //
  // KC 12/31/110
  // My ym is already Caillol's alpha*
  // Volker's alpha is Caillol's beta*
  // 

  // Because it only shows up in this combination now
  ym /= 2*M_PI;

  for(h[0] = -5; h[0] <= 5; h[0]++)
    for(h[1] = -5; h[1] <= 5; h[1]++)
      for(h[2] = -5; h[2] <= 5; h[2]++)
  	{
  	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
  	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0) {

	    //	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);
	    val = 2*exp(-M_PI*M_PI*(h2 + ym*ym)/(alpha*alpha))*sin(2*M_PI*hdotx) /
	      (h2 + ym*ym);

	    for(i = 0; i < 3; i++)
	      force[i] -= h[i] * val;
	  }
	}
#endif
}

double lattice_pot_none(double x[3]) {

  return 0.0;
}

void lattice_force_none(int iii, int jjj, int kkk, double x[3], double force[3]) {
  
  int i;

  for(i = 0; i < 3; ++i)
    force[i] = 0;
}

/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 *
 *  Note that the r (x) values used here are dimensionless....
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
  double alpha, r2;
  double r, val, hdotx, dx[3];
  int i, h[3], n[3], h2;

  alpha = 2.0;

  for(i = 0; i < 3; i++)
    force[i] = 0;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

  // KC 11/16/15
  // Here we add in the newtonian force, with the displacement vector
  // normalization included.  This makes me wonder if I should be including this factor in my
  // short-range tabulations.  It would allow me to make the spline functions and acceleration
  // functions interchangable again, and to remove an additional division in non-splined
  // forces, which are both nice features....
  // 
  for(i = 0; i < 3; i++)
    force[i] += x[i] / (r2 * sqrt(r2));

  // KC 12/4/14
  // Looks like this takes the first four images out in position space in each direction (so 
  // bracketing by 8 overall)
  for(n[0] = -4; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	  val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);

	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r * r) * val;
	}

  // KC 12/4/14
  // Looks like this takes the first four images in momentum space in each diretion (again
  // bracketing by 8 overall)
  for(h[0] = -4; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0)
	    {
	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);

	      for(i = 0; i < 3; i++)
		force[i] -= h[i] * val;
	    }
	}
}
