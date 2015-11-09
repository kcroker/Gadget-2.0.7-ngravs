#include "allvars.h"

#ifdef NOTYPEPREFIX_FFTW
#include        <fftw.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <dfftw.h>	/* double precision FFTW */
#else
#include     <sfftw.h>
#endif
#endif


struct ngravsInterpolant {

  fftw_plan plan;
  int ntab;
  int len;
  int ngravs_tpm_n;
  int ol;
};
 
// Gravitational interactions (accelerations)
double none(double target, double source, double h, double r, long N);
double newtonian(double target, double source, double h, double r, long N);
double neg_newtonian(double target, double source, double h, double r, long N);
double bambam(double target, double source, double h, double r, long N);
double sourcebaryonbam(double target, double source, double h, double r, long N);
double sourcebambaryon(double target, double source, double h, double r, long N);
double newyukawa(double target, double source, double h, double r, long N);

// (Periodic) Greens functions
double newtonKGreen(double target, double source, double k2, double k, long N);
double pgdelta(double target, double source, double k2, double k, long N);
double neg_pgdelta(double target, double source, double k2, double k, long N); 

// Potential functions
double newtonian_pot(double target, double source, double h, double r, long N);
double neg_newtonian_pot(double target, double source, double h, double r, long N);
double bambam_pot(double target, double source, double h, double r, long N);
double sourcebambaryon_pot(double target, double source, double h, double r, long N);
double sourcebaryonbam_pot(double target, double source, double h, double r, long N);

// Spline functions
double plummer(double target, double source, double h, double r, long N);
double neg_plummer(double target, double source, double h, double r, long N);
double null_spline(double target, double source, double h, double r, long N);

// Potential spline functions
double plummer_pot(double target, double source, double h, double r, long N);
double neg_plummer_pot(double target, double source, double h, double r, long N);

// Functions required for convolution
int gadgetToFourier(int i, struct ngravsInterpolant *s);
int performConvolution(struct ngravsInterpolant *s, gravity normKGreen, FLOAT Z, FLOAT *oRes, FLOAT *oResI);
struct ngravsInterpolant *ngravsConvolutionInit(int ntab, int len, int ol);
void ngravsConvolutionFree(struct ngravsInterpolant *s);
FLOAT mTox(int m, struct ngravsInterpolant *s);

// Initialization functions for ngravs extension
void wire_grav_maps(void);
void init_grav_maps(void);
