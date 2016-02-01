#!/usr/bin/python
#
# Uses g2munge to generate a bunch of initial conditions
# suitable for forcetests probing behaviour of 
#
# These initial conditions consist:
# ------------------------------------------
# * a single particle of unit mass placed at the origin
# * a single other particle placed elsewhere, randomly
# * TreePM boundary conditions
# * Forcetest = 1.0
# * a Configuration such that only one force computation is performed
# 
# Runs Gadget2 over and over again, collating all the forcetest data
#
# Takes a stub file as input (which has a zero centered g2munge txt dump)
#
import sys
from math import *
import os
import numpy as np
import random

# So that the bumper code can set things up
proggyName = "tpmfp"

# Make datum more useful for parsing forcetest
class Datum:
    # and a constructor to take the straight data without having to do anything
    def __init__(self, ptype, ts, crap3, x, y, z, fx, fy, fz, jx, jy, jz, cfx, cfy, cfz, crap4):
        self.ptype = int(ptype)
        self.ts = ts
        self.pos = np.array([float(u) for u in [x, y, z]])
        self.daccel = np.array([float(u) for u in [fx, fy, fz]])
        self.tpm_accel = np.array([float(u) for u in [cfx, cfy, cfz]])
        self.tree_accel = np.array([float(u) for u in [jx, jy, jz]])

# Check for the directory, make it if its not there
if not os.path.exists(proggyName):
   os.mkdir(proggyName)

if len(sys.argv) < 4:
    print "Usage: %s <# times> <# bins> <label> [skip sim?]" % sys.argv[0]
    exit(1)

# Parse the command line arguments
(times, bins, label) = [t(s) for t,s in zip( (int, int, str), sys.argv[1:4])]

# Make the stash-path
stash = "./%s/%s" % (proggyName, label)

# Make a stash for this analysis run 
if not os.path.exists(stash):
    os.mkdir(stash)
# else:
#     print "%s labelled run already exists, refusing to overwrite" % label
#     sys.exit(1)

# Some general settings
L = 10000
spheres = 1000
particles = 7

# Good values differ depending on what you're 
# trying to test.  (Ngravs can actually do
# exact force accuracies by turning off forces
# amoung the same type)
testm = 1
cenm = 1e6

# So the variable persists
i = 0
stack = []

# Do it 
# If there is a 3rd argument, skip running the sim and just work down the data
if not len(sys.argv) > 4: 
    for i in range(times):

        # Randomly pick a spot
        
        # 1) Choose a center
        center = np.array(np.random.uniform(0.0, float(L), 3))
        body = []

        # 2) Iterate over radii
        for r in [L/2.0 - N*L/(2*spheres) for N in range(1, spheres)]:
            #for k in range(0, int((1 - (r*2.0/L)**3)*particles)):
            for k in range(0, int(particles/(r*0.5/L)**2)):
                radius = np.random.uniform(r, r+L/(2.0*spheres))
                theta = np.random.uniform(0, pi)
                phi = np.random.uniform(0, 2*pi)
                pos = center + radius*np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
            
                body.append(pos)

        # Now we manually create the stub file and write it out
        stub = open("%s/stub" % stash, 'w')
        stub.write("# g2munge: %d\n\n# Group 0: 0\n\n# Group 1: %d\n" % (len(body)+1, len(body)))

        for n,pos in enumerate(body):
            stub.write("%d %.15e %.15e %.15e 0.0 0.0 0.0 %f -1.0\n" % 
                       (n, pos[0], pos[1], pos[2], testm))
        
        # Now write the center
        stub.write("\n# Group 2: 1\n")
        stub.write("%d %.15e %.15e %.15e 0.0 0.0 0.0 %f -1.0\n" % 
                   (len(body), center[0], center[1], center[2], cenm))
        stub.close()

        # 3) Now produce the initial condition
        os.system("g2munge - ascii > %s/tpmfp_IC < %s/stub" % (proggyName, stash))
        
        # 4) Run the simulation for a single timestep
        os.system("mpirun -n 2 ./Gadget2 utilities/Configuration.tpmfp")

        # Stash it
        os.system("mv ./%s/forcetest.txt %s/forcetest_%d.txt" % (proggyName, stash, i))
else:
    # Set i to the desired number of files (minus 1 of course)
    i = times - 1
    print "Skipping sim runs, looking for %d forcetest files..." % (i+1)

# Now make the forceerrors
# Unbuffered
errors = []

def processStack(found, center):
    # We have transitioned, compute the forces from this batch
    if not found:
        print "Never found the heavy guy in run %d..." % i
        exit(1)

    while len(stack) > 0:
        e = stack.pop()

        # Make sure to get separation from nearest periodic neighbor!
        pos = [x if x < L/2 else L - x for x in center - e.pos]
        pos = np.array([x if x > -L/2 else L + x for x in pos])
        
        r = np.linalg.norm(pos)
        truef = np.linalg.norm(e.daccel)
        compf = np.linalg.norm(e.tpm_accel)

        # Append this datapoint
        errors.append([r, compf/truef - 1])

while i > -1:
    
    # make stack, forget center each time
    first = True 
    found = False
    center = 0

    # Each run only has one time
    for n,line in enumerate(open("%s/forcetest_%d.txt" % (stash, i))):

        # Splat is a wonderful thing
        e = Datum(*line.split())

        # Mark the first time we get, so we can process times in chunks
        if first:
            now = e.ts
            first = False

        if not now == e.ts:
            processStack(found, center)
            
            # Set to the new time
            now = e.ts
            found = False
        else:
            # We are in the same timezone
            # Are we the heavy guy?
            if e.ptype == 2:
                # Set center, but do not add to the stack
                center = e.pos
                print "Found run %d massive center at line %d: (%f, %f, %f)" % ((i, n) + tuple(e.pos))
                found = True
            else:    
                # Stash for later
                stack.append(e)
                
    # We still have the last segment on the stack
    processStack(found, center)

    # Process the next forcetest
    i -= 1

# Write it out
recv = open("%s/tpmfp_output" % stash, 'w', 0)
for error in errors:
    recv.write("%.15e %.15e\n" % (error[0], error[1]))
recv.close()

# Sort it by separation
errors = sorted(errors, lambda x,y: cmp(x[0], y[0]))

print "%d Gadget2 forcetest.txt files worked down!" % times

# Now compute the binned RMS, unbuffered
binned = open("%s/tpmfp_rms" % stash, "w", 0)
i = 1

# Gadget-2 paper looks at things out to half the dimensionless boxsize
# sys.argv[2] is now interpreted as a desired number of bins
# x^sys.argv[2] = 10^4 --> x = pow(10^4, 1.0/sys.argv[2])
binbase = pow(10**4, 1.0/float(bins))
binrights = [pow(binbase, i) for i in range(0, bins)]
here = 0

for r,error in errors:
    
    # Does it belong in a later bin?
    if here == len(binrights):
        print 

    if r > binrights[here]:
        bincen = sqrt(binrights[here] * pow(binbase, i-1))
        
        # Go to the next bin
        here += 1

        # How many are we waiting to pop?
        N = len(stack)
        
        if N > 0:
            # Record the RMS
            # Binned radius is a midpoint
            sum = 0.0
            while len(stack) > 0:
                sum += float(stack.pop())**2
            binned.write("%f %f\n" % (bincen, sqrt(sum/N)))
        else:
            # If we don't have things, that means the previous bin is empty
            binned.write("%f %f\n" % (bincen, 0.0))

    # Always push the newest entry
    stack.append(error)

binned.close()
print "Output binned!"
sys.exit()
