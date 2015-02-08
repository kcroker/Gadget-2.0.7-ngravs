Gadget-2.0.7-ngravs
Copyright(c) 2015 Kevin Arthur Schiff Croker
GPL v2
---------------------------------------------

This code implements an extension of the N-body code GADGET-2 by Volker Springel so as to
enable computation with D distinct interacting gravitational species.  To use the code, 
one should first be familiar with GADGET-2, cas an be found here

http://www.mpa-garching.mpg.de/gadget/  (Users' Guide)

To use ngravs, all the user need modify is the ngravs.c and ngravs.h files.  Comments in the C
file are quite verbose and there are two completely functional examples used to test the code, 
one may follow them as necessary.

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
the relevant code paper, submitted to MNRAS (Monthly Notices of the Royal Astronomical Society).

Word!
-k@Manoa 2/7/15
