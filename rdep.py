#!/usr/bin/python

# Investigate the computed force law as a function of radial separation

import sys
import math
import os
import random

# So that the bumper code can set things up
proggyName = "rdep"

# Make datum more useful for parsing forcetest
class Datum:
    # and a constructor to take the straight data without having to do anything
    def __init__(self, crap1, crap2, crap3, x, y, z, fx, fy, fz, jx, jy, jz, cfx, cfy, cfz):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.fx = float(fx)
        self.fy = float(fy)
        self.fz = float(fz)
        self.cfx = float(cfx)
        self.cfy = float(cfy)
        self.cfz = float(cfz)

# Check for the directory, make it if its not there
if not os.path.exists(proggyName):
   os.mkdir(proggyName)

if len(sys.argv) < 3:
    print "Usage: %s <# times> <box width> [skip sim?]" % sys.argv[0]
    exit(1)

# Convert boxlength to radius (this will miss the edges though...)
L = float(sys.argv[2]) * 0.5
N = int(sys.argv[1])
random.seed()

# Fixed center for love
cenx = 4967.0
ceny = cenx
cenz = cenx

seps = [L/s for s in range(1,N)]

for i,d in enumerate(seps):

    # 0)
    print "Taking data for separation %f..." % d

    # 1) Make an initial condition with this setup
    # Have this be random at first unless there are problems that manifest
    x = random.random()
    y = random.random()
    z = random.random()
    norm = math.hypot(x, math.hypot(y, z))
    x = cenx + d*x/norm
    y = ceny + d*y/norm
    z = cenz + d*z/norm

    stub = "# g2munge: 2\n\n# Group 0: 0\n\n# Group 1: 1\n1 %f %f %f 0 0 0 1 -1.0\n\n# Group 2: 1\n2 %f %f %f 0 0 0 1e6 -1.0" % (x, y, z, cenx, ceny, cenz)
    os.system("echo \"%s\" | g2munge - ascii > %s_ic" % (stub, proggyName))

    # 2) Run the sim for one timestep
    os.system("mpirun -n 1 ./Gadget2 Configuration.rdep")
    os.system("mv ./%s/forcetest.txt ./%s/forcetest_%d.txt" % (proggyName, proggyName, i))

    # 3) Process the data
    recv = open("./%s/%s_output" % (proggyName, proggyName), 'w', 0)

    src = enumerate(open("./%s/forcetest_%d.txt" % (proggyName, i)))
    
    # Use the splat operator
    e1 = Datum(*next(src)[1].split())
    e2 = Datum(*next(src)[1].split())
    
    # Compute separation distance in non-overflow way
    r = math.hypot(math.hypot(e1.x - e2.x, e1.y - e2.y), e1.z - e2.z)
    
    # OKAY is work
    
    # Now output forces, and try to figure out why GravPM does not satisfy Newton's 3rd law?
    # (with stock Gadget!)
    #
    
