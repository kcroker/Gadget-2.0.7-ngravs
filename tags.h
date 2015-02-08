/*! \file tags.h
 *  \brief declares various tags for labelling MPI messages.
 */

#define TAG_N             10      /*!< Various tags used for labelling MPI messages */ 
#define TAG_HEADER        11
#define TAG_PDATA         12
#define TAG_SPHDATA       13
#define TAG_KEY           14
#define TAG_DMOM          15
#define TAG_NODELEN       16
#define TAG_HMAX          17
#define TAG_GRAV_A        18
#define TAG_GRAV_B        19
#define TAG_DIRECT_A      20
#define TAG_DIRECT_B      21
#define TAG_HYDRO_A       22 
#define TAG_HYDRO_B       23
#define TAG_NFORTHISTASK  24

// KC 10/3/14
// These tags are used to keep the pm_periodic fourier computation in sync
// So we will generalize these tags so that I can use bitwise operations against
// them depending on the N_GRAVS, and keep the entire algorithm in sync
#define TAG_PERIODIC_A    25
#define TAG_PERIODIC_B    26
// KC 10/4/14 
// MPI tags are integers, which should be at least 16 bits.
// Note that maxtag is 37.  There are N_GRAVS(N_GRAVS + 1)/2 types of interaction to compute, so at most
// 6(7)/2 = 42/2 = 21 
// So lets start at bit 6: 2^6 = 64
// --> bits 6,7,8: na
// --> bits 9,10,11: nb
// --> Encode: TAG_PERIODIC_A | na << 6 | nb << 9
#define TAG_PERIODIC_C    27
#define TAG_PERIODIC_D    28

#define TAG_NONPERIOD_A   29 
#define TAG_NONPERIOD_B   30
#define TAG_NONPERIOD_C   31
#define TAG_NONPERIOD_D   32
#define TAG_POTENTIAL_A   33
#define TAG_POTENTIAL_B   34
#define TAG_DENS_A        35
#define TAG_DENS_B        36
#define TAG_LOCALN        37

