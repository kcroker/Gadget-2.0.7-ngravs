#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

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
 */

#define YUKAWA_ALPHA 1

#ifdef PERIODIC
#define YUKAWA_IMASS 0.9  // Should be given in dimensionless fraction of the boxsize
#else
#define YUKAWA_IMASS (1e2/All.BoxLength) // Otherwise, do it in terms of normal units
#endif


/*! This function must be modified to point to your desired
 *  extensions to gravity.  It determines which force laws
 *  are used to compute interactions.
 */
void wire_grav_maps(void) {

  int i,j;
  char *fname;

  // KC 8/11/14 Wiring
  //
  // NOTE: For all interaction functions
  //
  //       InteractionFunctions[TARGET][SOURCE]
  //       (i.e. InteractionFunctions[PASSIVE][ACTIVE])
  //
  // VERY IMPORTANT: NgravsNames[][] is used to index things used by the simulation
  // code, like lattice correction tables.  So please make a unique identifier!
  // (It will also be used eventually to save memory and startup-time by not 
  // computing redundant tables.)
  //

#if !defined NGRAVS_STOCK_TESTING && !defined NGRAVS_ACCUMULATOR_TESTING && !defined NGRAVS_GEN_TESTING_UNIFORM
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

#elif defined NGRAVS_GEN_TESTING_UNIFORM
  printf("ngravs: wired in uniform generalized force test mode\n");
  //////////////// BEGIN GENERALIZE FORCE TESTING WIRING ///////////////////////////
  //
  // This code is used to examine the force accuracy during the TreePM transition
  // for a more general force, in this case Yukawa.  Two populations of particles
  // interact via Newtonian, except the self-interactions of one type are Yukawa.
  //
  ///////////////////////////////////////////////////////////////////////

  for(i = 0; i < N_GRAVS; ++i) {
    for(j = 0; j < N_GRAVS; ++j) {
      
      // Allocate a new one each time, 
      // because that's what we'd have to be doing anyway.
      // This is not a leak because we need these handles throughout the
      // entire program run.
      fname = (char *)malloc(128);
      snprintf(fname, 128, "Yukawa_%f", YUKAWA_IMASS);
 
      NgravsNames[i][j] = fname;
      AccelFxns[i][j] = yukawa;

      // We set the Yukawa spline to plummer since
      // the force is Newtonian at small r
      // This being incorrect won't matter for force checking the TreePM 
      // stuff, because the force correction uses the spline too.
      AccelSplines[i][j] = plummer;

#if defined PERIODIC
      // Computed from G. Salin and J.M. Caillol
      // J. Chem. Phys., Vol 113, No. 23, 2000
      LatticeForce[i][j] = yukawa_lattice_force;
      LatticePotential[i][j] = yukawa_lattice_psi;
      LatticeZero[i][j] = yukawa_madelung(YUKAWA_IMASS);
      if(!ThisTask)
	printf("ngravs: Yukawa force Madelung constant for [%d][%d] = %f\n", i, j, LatticeZero[i][j]);
#endif

#if defined OUTPUTPOTENTIAL || defined PMGRID
      GreensFxns[i][j] = pgyukawa;
      
      // We don't care about the potentials because we're
      // not doing non-periodic or gastrophysics
      PotentialFxns[i][j] = none;
      PotentialSplines[i][j] = none;
      PotentialFxns[i][j] = none;
#endif
    }
  }

/*   NgravsNames[0][0] = "Newton"; */
/*   NgravsNames[0][1] = "Newton"; */
/*   NgravsNames[1][0] = "Newton"; */
/*   NgravsNames[1][1] = "Yukawa"; */

/*   AccelFxns[0][0] = newtonian; */
/*   AccelFxns[0][1] = newtonian; */
/*   AccelFxns[1][0] = newtonian; */
/*   AccelFxns[1][1] = yukawa; */

/*   // We set the Yukawa spline to plummer since */
/*   // the force is Newtonian at small r */
/*   AccelSplines[0][0] = plummer; */
/*   AccelSplines[0][1] = plummer;  */
/*   AccelSplines[1][0] = plummer; */
/*   AccelSplines[1][1] = plummer;  */

/* #if defined PERIODIC */
/*   LatticeForce[0][0] = ewald_force; */
/*   LatticePotential[0][0] = ewald_psi; */
/*   LatticeZero[0][0] = 2.8372975; */

/*   LatticeForce[0][1] = LatticeForce[1][0] = ewald_force; */
/*   LatticePotential[0][1] = LatticePotential[1][0] = ewald_psi; */
/*   LatticeZero[0][1] = LatticeZero[1][0] = 2.8372975; */

/*   LatticeForce[1][1] = lattice_force_none; */
/*   LatticePotential[1][1] = lattice_pot_none; */
/*   LatticeZero[1][1] = 0.0; */
/* #endif */

/* #if defined OUTPUTPOTENTIAL || defined PMGRID */
/*   GreensFxns[0][0] = pgdelta; */
/*   GreensFxns[0][1] = pgdelta; */
/*   GreensFxns[1][0] = pgdelta; */
/*   GreensFxns[1][1] = pgyukawa; */

/*   // We don't care about the potentials because we're */
/*   // not doing non-periodic or gastrophysics */
/*   PotentialFxns[0][0] = none; */
/*   PotentialSplines[0][0] = none; */

/*   PotentialSplines[0][1] = none; */
/*   PotentialFxns[0][1] = none; */

/*   PotentialSplines[1][0] = none; */
/*   PotentialFxns[1][0] = none; */

/*   PotentialSplines[1][1] = none; */
/*   PotentialFxns[1][1] = none; */
  
/*   PotentialFxns[0][1] = none; */
/*   PotentialFxns[1][0] = none; */
/*   PotentialFxns[1][1] = none; */
/*   PotentialZero[0][0] = 0.0; */
/*   PotentialZero[0][1] = 0.0; */
/*   PotentialZero[1][0] = 0.0; */
/*   PotentialZero[1][1] = 0.0; */
/* #endif */

  //////////////// END GENERALIZE FORCE TESTING WIRING ///////////////////////////

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


/*! This is a (0) lattice correction for Yukawa force
 */

/*! This is Newtonian gravity with a Yukawa modification
 * c.f. Shirata, Yoshida, et. al 2005
 */
double newyukawa(double target, double source, double h, double r, long N) {

  return source / h * -YUKAWA_ALPHA * expm1f(-r*YUKAWA_IMASS);
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

/*! A pure Yukawa force
 *
 * Radii take values in [0, BoxLength] for periodic, unconstrained otherwise.
 * If we are in PERIODIC mode, so that the tables do not become sensitive to
 * the BoxLength as given in the configuration files, YUKAWA_IMASS is in units of EN.
 */
double yukawa(double target, double source, double h, double r, long N) {
  
#if defined PERIODIC
  return source * YUKAWA_ALPHA * exp(-r*YUKAWA_IMASS/All.BoxSize) * (YUKAWA_IMASS/(All.BoxSize*r) + 1.0/h) * pow(YUKAWA_IMASS/All.BoxSize, 2);
#else
  return source * YUKAWA_ALPHA * exp(-r*YUKAWA_IMASS) * (YUKAWA_IMASS / r + 1.0/h) * YUKAWA_IMASS;
#endif
}

/*! A periodic yukawa k-space Greens function, normalized by the Newtonian interaction
 *  NOTE: k is supplied dimensionlessly in terms of PMGRID so k \in [-PMGRID/2, PMGRID/2]
 */
double pgyukawa(double target, double source, double k2, double k, long N) {

  // This function is only called in periodic mode, so it to PMGRID units
  // Note that we have the beta^2 prefactor because Gadget-2 expects the total charge to be
  // unity!
  //
  // KC 19.11.15
  // Or these need to be normalized such that the DFT gives unit charge (this should already be the case though)
  // XXX HACK 0.5 to see if we can restore the correct behaviour in long-range.  I believe this is due to 
  // discrepancy in box normalization....
  return k2 * (YUKAWA_IMASS*YUKAWA_IMASS*PMGRID*PMGRID) / (k2 + (YUKAWA_IMASS*YUKAWA_IMASS)*(PMGRID*PMGRID));
}

/*! This function computes the Madelung constant for the yukawa potential
 * which depends on the box length interestingly...
 * We follow Eqn (2.19) of G. Salin and Caillol (op. cit)
 *
 * Note that we use the values de-dedimensionalized in the same way
 * as the other yukawa lattice functions.
 *
 */
double yukawa_madelung(double ym) {
  
  double sum1, sum2, sum3;
  double k2, m;
  double alpha;
  int n[3];

  // KC 11/16/15
  // We again adopt the same notation as that used in Gadget-2, so their beta
  // is our alpha, etc (see comments above)
  //
  // Note something very interesting:
  // This function should be *INDEPENDENT* of alpha if ym=0
  // Gotta figure out why its not...
  alpha = 5.64;

  // Going out to the same distance seems like a good idea, no?
  // Notice that reducing all the quantities results in an over factor of (1/L)
  // This factor is NOT present in the originally specified number, so Gadget-2
  // must be putting it back in somewhere.  So we elide it.
  for(n[0] = -5, sum1 = 0, sum2 = 0; n[0] <= 5; n[0]++) {
    for(n[1] = -5; n[1] <= 5; n[1]++) {
      for(n[2] = -5; n[2] <= 5; n[2]++) {

	// Here we use n for both the k sum and the n sum because they are both dimensionless
	k2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
	if(k2 > 0) {
	  m = sqrt(k2);
	  k2 *= 4*M_PI*M_PI;
	  sum1 += exp(-(k2 + ym*ym)/(4*alpha*alpha))/(k2 + ym*ym);
	  sum2 += (erfc(alpha*m + ym/(2*alpha))*exp(ym*m) + erfc(alpha*m-ym/(2*alpha))*exp(-ym*m))/(2*m);
	}
      }
    }
  }
  
  // The non-summation terms
  sum3 = -2*alpha/sqrt(M_PI)*exp(-ym*ym/(4*alpha*alpha)) +
    ym*erfc(ym/(2*alpha));
  
  // Explicitly the zero yukawa-mass case, limit via l'Hopital
  if(ym > 0)
    sum3 += 4*M_PI/(ym*ym) * expm1(ym*ym/(4*alpha*alpha));
  else
    sum3 += 4*M_PI/(4*alpha*alpha);
  
  return (4*M_PI*sum1 + sum2 + sum3);
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
  //
  // NOTE: YUKAWA_IMASS is in terms of 1/NGRAVS_EN since r here is in units of NGRAVS_EN
  //
  
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

  // KC 11/16/15
  // Note our use of Salin's optimal 'alpha', and our excessive momentum-space
  alpha = 5.64;

  for(i = 0; i < 3; i++)
    force[i] = 0;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  r = sqrt(r2);

  // KC 11/16/15
  // Here we add in the original force, with the displacement vector
  // normalization included.  This makes me wonder if I should be including this factor in my
  // short-range tabulations.  It would allow me to make the spline functions and acceleration
  // functions interchangable again, and to remove an additional division in non-splined
  // forces, which are both nice features....
  // 
  // Optimize later (this is only used for preliminary tabulations)
  //
  // Note that r and r2 are in units of NGRAVS_EN, so we don't use the usual Yukawa force
  // function which assumes things are in internal units...
  //
  // (Because YUKAWA_IMASS is defined in NGRAVS_EN units so that tabulations do not need to be repeated
  //  for different box lengths)
  //
  
  // KC 11/19/15
  // Note that mass^2 prefactor so that the integrated charge is unity!
  //
  for(i = 0; i < 3; i++)
    force[i] += YUKAWA_ALPHA * exp(-r*YUKAWA_IMASS) * (YUKAWA_IMASS / r2 + 1.0/(r2*r)) * x[i] * YUKAWA_IMASS*YUKAWA_IMASS; 

  //yukawa(1.0, 1.0, r2, sqrt(r2), 1) * (x[i] / sqrt(r2));

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
	  val = 0.5*( exp(YUKAWA_IMASS*r)*erfc(alpha*r + YUKAWA_IMASS/(2*alpha)) + 
		      exp(-YUKAWA_IMASS*r)*erfc(alpha*r - YUKAWA_IMASS/(2*alpha)));
	  
	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r * r) * val;

	  val += 0.5*YUKAWA_IMASS*(-exp(YUKAWA_IMASS*r)*erfc(alpha*r + YUKAWA_IMASS/(2*alpha)) + 
				  exp(-YUKAWA_IMASS*r)*erfc(alpha*r - YUKAWA_IMASS/(2*alpha))) +
	    2*alpha*exp(-alpha*alpha*r*r-YUKAWA_IMASS*YUKAWA_IMASS/(4*alpha*alpha))/sqrt(M_PI);

	  // KC 11/16/15
	  // Note that these terms enter with one less radial power in the denominator
	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r) * val;
	}

  // KC 12/4/14
  // Looks like this takes the first four images in momentum space in each diretion (again
  // bracketing by 8 overall)
  //
  // KC 11/16/15
  // Take careful note of the relative signs!!
  // If you take the negative grad_r, then the signs are the same!
  for(h[0] = -5; h[0] <= 5; h[0]++)
    for(h[1] = -5; h[1] <= 5; h[1]++)
      for(h[2] = -5; h[2] <= 5; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0)
	    {
	      //	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);
	      val = 2*M_PI*exp(-(4*M_PI*M_PI*h2 + YUKAWA_IMASS*YUKAWA_IMASS)/(4*alpha*alpha))*sin(2*M_PI*hdotx) / 
		(M_PI*h2 + YUKAWA_IMASS*YUKAWA_IMASS/(4*M_PI));
	      
	      // KC 11/16/15
	      // XXX?
	      for(i = 0; i < 3; i++)
		force[i] -= h[i] * val;
	    }
	}
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



///////////////// END GENERALIZED FORCE AND GREENS FUNCTIONS ////////////
//
// Below this are things you should not need to modify
//
/////////////////////////////////////////////////////////////////////////

///////////////// BEGIN FOURIER INTEGRATION ROUTINES /////////////////////
//
// These routines can compute the reqired shortrange tabulations of the generic
// force laws from the k-space greens functions to very near machine accuracy.  
//
// Please see the ngravs paper for very detailed discussion of the behaviour here
//

// Some default values here
// How much to oversample by
// WHY: This determines the magnitude of the error term: 
//#define OL 20

// How much further (by this multiplicative factor) to go in x-space 
// than the required 3-\epsilon
// WHY: This (essentially) determines how deep into k-space to probe the function
//#define LEN 4

// This is then the required number of samples.
// #define NGRAVS_TPM_N (12*NTAB*OL*LEN-6*OL*LEN+2)

// Note: the product of mTox * jTok = 2\pi mj/N_G
//       which is what FFTW expects
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
  
  // XXX?
  // But if we fuck with this, we'll break newton.  The exponential factor
  // seems to be correct because it works with newton, where the normgreens = 1.0
  //
  // But if this were the problem, it would only be on the lowend of the force...
  return (*normKGreen)(1, 1, k2, k, 1) * exp(-k2 * Z * Z);
}

//
// Sigh.  Fourier routines need to be reconsidered from the point of view of 
// units being in MESH CELLS, beacuse that is what the normKGreen works in, 
// and if you look at pm_periodic, what k2 takes values in :/
//
// You didn't see this problem with Newton because normKGreen = 1.
//
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

  sum = s->ngravs_tpm_n;

  for(m = 0; m < s->ntab; ++m)
    oRes[m] = out[gadgetToFourier(m, s)].re * norm;
 
  // 3) Integrate so as to constrain the error correctly:
  // Newton-Cotes 4-point rule
  // Run the sum at double precision, though we may assign to lower precision
  sum = 0.0;
  in[0].re = 0.0;
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

  printf("ngravs: max_k = %.15e\n", jTok(s->ngravs_tpm_n/2, 0.5, s));

  return s;
}

void ngravsConvolutionFree(struct ngravsInterpolant *s) {

  fftw_destroy_plan(s->plan);
  free(s);
}

FLOAT hermiteSpline(FLOAT u, struct ngravsInterpolant s, FLOAT *tab, FLOAT *dtab) {

  int i = 0;

  // Note that u should be rescaled to [0,1] a priori
  // Use factorized form to keep addition of terms of roughly the same order  
  return tab[i] * (1 + 2*u) * (1-u) * (1-u) + 
    dtab[i] * u * (1-u) * (1-u) +
    tab[i+1] * u*u*(3-2*u) +
    dtab[i+1] * u*u*(u - 1);
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
    printf("ngravs: High-resolution mesh implementation in progress.  Email kcroker@phys.hawaii.edu for status.");
    endrun(1000);
#endif

#if defined PMGRID && !defined PERIODIC
    printf("ngravs: In this version, non-periodic k-space greens determined by FFT of *tabulated* values.  Full transformation implementation in progress.  For now, verify sanity of outputs!");
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
