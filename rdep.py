#!/usr/bin/python

# Investigate the computed force law as a function of radial separation

import sys
import math
import os
import random
import numpy as np

# So that the bumper code can set things up
proggyName = "rdep"

# Make datum more useful for parsing forcetest
class Datum:
    # and a constructor to take the straight data without having to do anything
    def __init__(self, ptype, crap2, crap3, x, y, z, fx, fy, fz, jx, jy, jz, cfx, cfy, cfz, crap4):
        self.ptype = int(ptype)
        self.pos = np.array([float(u) for u in [x, y, z]])
        self.daccel = np.array([float(u) for u in [fx, fy, fz]])
        self.tpm_accel = np.array([float(u) for u in [cfx, cfy, cfz]])
        self.tree_accel = np.array([float(u) for u in [jx, jy, jz]])

# Check for the directory, make it if its not there
if not os.path.exists(proggyName):
   os.mkdir(proggyName)

if len(sys.argv) < 4:
    print "Usage: %s <# times> <box width> <label> [retake data?]" % sys.argv[0]
    exit(1)

# Convert boxlength to radius (this will miss the edges though...)
L = float(sys.argv[2]) * 0.5
N = int(sys.argv[1])

# Make sure its even
N = N + (N % 2)

label = sys.argv[3]
testm = 2.5e-1
cenm = 1.0e3

random.seed()

# Fixed center for love
# cen = np.array([4967.0, 4967.0, 4967.0])
cen = np.array([random.random()*L for x in range(3)])

# These separations probe the default TPM transition scale
seps_pmxition = [10 + n*(400 - 10)/(N/2) for n in range(N/2)]
seps_longrange = [400 + n*(L - 400)/(N/2) for n in range(N/2)]

#seps = [L/s for s in range(1,N+1)]
seps = [x for x in seps_pmxition + seps_longrange]
print seps

if len(sys.argv) > 4:

    # Make the target directory for the forcetests
    if not os.path.exists("./%s/%s" % (proggyName, label)):
        os.mkdir("./%s/%s" % (proggyName, label))
    else:
        print "%s run already exists.  Refusing to overwrite."
        sys.exit()

    for i,d in enumerate(seps):

        # 0)
        print "Taking data for separation %f..." % d
        
        # 1) Make an initial condition with this setup
        # Have this be random at first unless there are problems that manifest
        tpos = np.array([random.random() for x in range(3)])
        print "Sanity on tpos: %s" % [float(x) for x in tpos]

        norm = np.linalg.norm(tpos)
        tpos = cen + d*tpos/norm

        # This feels dumb, but works
        fmtargs = tuple(tpos) + tuple([testm]) + tuple(cen) + tuple([cenm])

        stub = "# g2munge: 2\n\n# Group 0: 0\n\n# Group 1: 1\n1 %f %f %f 0 0 0 %f -1.0\n\n# Group 2: 1\n2 %f %f %f 0 0 0 %f -1.0" % fmtargs
        os.system("echo \"%s\" | g2munge - ascii > %s_ic" % (stub, proggyName))

        # 2) Run the sim for one timestep, and only on one processor
        os.system("mpirun -n 1 ./Gadget2 Configuration.rdep")
        os.system("mv ./%s/forcetest.txt ./%s/%s/forcetest_%d.txt" % (proggyName, proggyName, label, i))


# Open and record the mass parameters
recv = open("./%s/%s/%s_output" % (proggyName, label, proggyName), 'w', 0)
recv.write("# testm:%.15e\n# cenm:%.15e\n" % (testm, cenm))

for i in range(N):
    # 3) Process the data
    src = enumerate(open("./%s/%s/forcetest_%d.txt" % (proggyName, label, i)))
    
    # Use the splat operator
    e1 = Datum(*next(src)[1].split())
    e2 = Datum(*next(src)[1].split())

    r = np.linalg.norm(e1.pos - e2.pos)
    
    # OKAY is work
    
    # Now output forces, and try to figure out why GravPM does not satisfy Newton's 3rd law?
    # (with stock Gadget!)
    # Remember, forcetest output contains accelerations, not forces!
    #
    # c is for center, t is for test
    # Make shortcuts to the test and center particles
    if e1.ptype == 1:
        t = e1
        c = e2
    else:
        t = e2
        c = e1

    # Scale by the mass
    for o,mass in zip([c, t], [cenm, testm]):
        o.daccel *= mass
        o.tpm_accel *= mass
        o.tree_accel *= mass
    
    # Record
    recv.write("%f %.15e %.15e %.15e %.15e %.15e %.15e\n" % (r, 
                                                             np.linalg.norm(c.daccel), 
                                                             np.linalg.norm(c.tpm_accel),
                                                             np.linalg.norm(c.tree_accel),
                                                             np.linalg.norm(t.daccel),
                                                             np.linalg.norm(t.tpm_accel),
                                                             np.linalg.norm(t.tree_accel)))
                                      
                            
recv.close()

    
