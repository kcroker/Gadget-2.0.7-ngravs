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
import math
import os

class Datum:
    def __init__(self, x, y, z, fx, fy, fz, cfx, cfy, cfz):
        self.x = x
        self.y = y
        self.z = z
        self.fx = fx
        self.fy = fy
        self.fz = fz
        self.cfx = cfx
        self.cfy = cfy
        self.cfz = cfz

# Check for the directory, make it if its not there
if not os.path.exists("tpmfp"):
   os.mkdir("tpmfp")

# Unbuffered
recv = open("./tpmfp/tpmfp_output", 'w', 0)

if len(sys.argv) < 2:
    print "Usage: %s <# times>" % sys.argv[0]
    exit(1)

# Do it 
for i in range(int(sys.argv[1])):

    # If there is a 3rd argument, skip running the sim and just work down the data
    if not len(sys.argv) > 2: 
        # Use this as the merge
        os.system("g2munge - rand 9999 10000 2.5e-4 1 2 > ./tpmfp/tpmfp_tests")
    
        # Create an IC containing one very heavy particle randomly
        os.system("g2munge - rand 1 10000 1.0 1 1 | g2munge - merge ./tpmfp/tpmfp_tests > ./tpmfp/tpmfp_IC")

        os.system("mpirun -n 2 ./Gadget2 Configuration.tpmfp")
    else:
        print "Skipping sim runs..."

    # make stack
    stack = []
    first = True 
    found = False

    # Now make the forceerrors
    for n,input in enumerate(open("./tpmfp/forcetest.txt")):

        # Slick scanf like parsing from Chris Dellin @ stackoverflow!
        (ptype, pt, ptdelta, 
         x, y, z, 
         fx, fy, fz, 
         jx, jy, jz, 
         cfx, cfy, cfz) = [t(s) for t,s in 
                           zip((int, float, float,
                                float, float, float,
                                float, float, float,
                                float, float, float,
                                float, float, float), input.split())]

        # Mark the first time we get, so we can process times in chunks
        if first:
            now = pt
            first = False
        
        if not now == pt:
            # We have transitioned, compute the forces from this batch
            if not found:
                print "Never found the heavy guy..."
                exit(0)

            while len(stack) > 0:
                e = stack.pop()
                    
                # Watch for overflow...
                r = math.hypot(math.hypot(cenx - e.x, ceny - e.y), cenz - e.z)
                truef = math.hypot(math.hypot(e.fx, e.fy), e.fz)
                compf = math.hypot(math.hypot(e.cfx, e.cfy), e.cfz)
        
                # Output a datapoint
                # Randomly take the computed force and flip the sign
                if truef < 2e-14:
                    recv.write("%f nan\n" % r)
                else:
                    recv.write("%f %f\n" % (r, (compf - truef)/truef))

            # Set to the new time
            now = pt
            found = False
        else:
            # We are in the same timezone
            # Are we the heavy guy?
            if ptype == 1:
                # Set centers, but do not add to the stack
                cenx = x
                ceny = y
                cenz = z
                print "Found run %d massive center at line %d: (%f, %f, %f)" % (i, n, x, y, z)
                found = True
            else:    
                # Stash for later
                stack.append(Datum(x, y, z, fx, fy, fz, cfx, cfy, cfz));

recv.close()
sys.exit()
