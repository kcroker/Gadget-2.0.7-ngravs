Gadget-2.0.7-ngravs
---------------------------------------------

Copyright(c) 2015 Kevin Arthur Schiff Croker
GPL v2

This code implements an extension of the N-body code GADGET-2 by Volker Springel so as to
enable computation with D distinct interacting gravitational species.  

Installation
----------------------------------------------
To use the code, one should first be familiar with GADGET-2, can be found here:

http://www.mpa-garching.mpg.de/gadget/

The contents of this repository are sufficient to GNU make:
1) cd to the repo directory
2) adjust Makefile for your system
3) type: <pre>make</pre>

If one wishes the documentation and test files included with GADGET-2, one may download from
the above URL and then replace the contents of the Gadget2/ subdirectory with the contents of this
repository.

(Note that the included Makefile.reference works on modern Ubuntu systems, and with 
the included Galaxy collision initial condition from GADGET-2)

Configuration
---------------------------------------------
Gadget-2.0.7-ngravs is configured identically to GADGET-2.  The following are new and 
mandatory configuration file options:

GravityGas [int] <br>
GravityHalo [int] <br>
GravityDisk [int] <br>
GravityBulge [int] <br>
GravityStars [int] <br>
GravityBndry [int] <br>

Here [int] will be any integer in [0,5] < N_GRAVS (as defined in Makefile) and specifies 
the gravitational interaction to be used by that particular GADGET-2 particle type.

To use ngravs for your particular model, modify the ngravs.c and ngravs.h files.  Comments in the C
file are quite verbose and there are two completely functional examples used to test the code, 
one may follow them as necessary.

(Note that the included Configuration.reference will run with the included Galaxy collision initial 
condition from GADGET-2)

Running
----------------------------------------------
For example, to run the included Galaxy Collision initial condition, type: <br>
<pre>
mpirun -n [int] ./Gadget2 Configuration.reference
</pre>

Note that the included GalaxyCollision.IC file uses little endian encoding (like x86 platform).  If
you require big endian encoding, please use the analogous file from the original GADGET-2 code package.

File list (copied from GADGET-2 Users' Guide, with additional ngravs.c description)
----------------------------------------------
accel.c: Driver routine to carry out force computation <br>
allocate.c: Routines for allocating particle and tree storage <br>
allvars.c: Provides instances of all global variables<br>
begrun.c: Initial set-up of a simulation run<br>
density.c: SPH density computation and smoothing length determination<br>
domain.c: Code for domain decomposition<br>
driftfac.c: Computes look-up tables for prefactors in cosmological integrations<br>
endrun.c: Terminates run upon severe errors<br>
forcetree.c: Gravitational tree code<br>
global.c: Global energy computation<br>
gravtree.c: Main driver routines for gravitational (short-range) force computation<br>
gravtree_forcetest.c: Routines for direct summation forces<br>
hydra.c: Computation of SPH forces and rate of entropy generation<br>
init.c: Code for initialisation of a simulation from initial conditions<br>
io.c: Routines for producing a snapshot file on disk<br>
longrange.c: Driver routines for computation of long-range gravitational PM force<br>
main.c: Start of the program<br>
ngb.c: Neighbour search by means of the tree<br>
ngravs.c: defines the gravitational force law for real space, k, and lattice sums<br>
peano.c: Routines to compute a Peano-Hilbert order<br>
pm_nonperiodic.c: Code for non-periodic FFT to compute long-range PM force<br>
pm_periodic.c: Routines for periodic PM-force computation<br>
potential.c: Computation of the gravitational potential of particles<br>
predict.c: Drift particles by a small time interval<br>
read_ic.c: Read initial conditions in one of the supported file formats<br>
restart.c: Code for reading and writing restart files<br>
run.c: Iterates over timesteps, main loop<br>
system.c: Miscellaneous routines, e.g. elapsed time measurements<br>
timestep.c: Routines for ‘kicking’ particles and assigning new timesteps<br>
allvars.h: Declares global variables<br>
ngravs.h: Contains ngravs specific prototypes
proto.h: This file contains all function prototypes of the code (non-ngravs)<br>
tags.h: Declares various tags for labelling MPI messages (non-ngravs).  Defines how ngravs
tags are determined from the pre-existing ones <br>

Other Notes
------------------------------------------------
To see how the extension works, as well as for some annotation on how the original Gadget-2 algorithms
operate, I have interspersed various dated comments with 'KC MM/DD/YY' into the source files wherever
I have made changes to the original source.  Sometimes, these comments are quite explicit and reference
technical works justifying the origin of a computation.  I leave them in so that someone might easily
check that what I have done is indeed correct (though it certainly seems to be the case...)

Regions in the code marked with XXX may be bugs.  At present, there are only two sitting in pm_nonperiodic.c.
Here, it seemed that a bunch of code was needlessly repeated, and so I commented it out with an explanation, 
but tagged it as possibly suspect.

Regions in the code marked with PPP may be amenable to performance improvements.  If you'll note, Volker was 
quite rigorous in eschewing divisions, making references to the same variables spatially proximal in the code,
unrolling loops, etc.  It turns out that these things matter, see the comments to get an idea of some of the 
performance gains!

I hope ngravs is useful for investigation of your wild gravitation models!  If you use it, please cite
the relevant code paper, submitted to Comput. Phys. Commun. (Computer Physics Communications).

Word!
-k@Manoa 2/7/15
-k@Pitt 6/4/15 

Acknowledgments:

This material is based upon work supported by the National Science
Foundation under Grant Number 1415111.  Any opinions, findings, and
conclusions or recommendations expressed in this material are those of
the author(s) and do not necessarily reflect the views of the National
Science Foundation.