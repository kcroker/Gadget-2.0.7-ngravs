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
2) adjust Makefile for your system (Makefile.reference works on modern Ubuntu systems)
3) type: make

If one wishes the documentation and test files included with GADGET-2, one may download from
the above URL and then replace the contents of the Gadget2/ subdirectory with the contents of this
repository.

Configuration
---------------------------------------------
Gadget-2.0.7-ngravs is configured identically to GADGET-2.  The following are new and 
mandatory configuration file options:

GravityGas <int>
GravityHalo <int>
GravityDisk <int>
GravityBulge <int>
GravityStars <int>
GravityBndry <int>

Here <int> will be any integer in [0,5] < N_GRAVS (as defined in Makefile) and specifies 
the gravitational interaction to be used by that particular GADGET-2 particle type.

To use ngravs for your particular model, modify the ngravs.c and ngravs.h files.  Comments in the C
file are quite verbose and there are two completely functional examples used to test the code, 
one may follow them as necessary.

Running
----------------------------------------------
Identical to GADGET-2.


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