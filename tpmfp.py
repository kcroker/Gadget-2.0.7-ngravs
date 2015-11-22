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

if len(sys.argv) < 3:
    print "Usage: %s <# times> <bin width> [skip sim?]" % sys.argv[0]
    exit(1)

# So the variable persists
i = 0
stack = []

# Do it 
# If there is a 3rd argument, skip running the sim and just work down the data
if not len(sys.argv) > 3: 
    for i in range(int(sys.argv[1])):

        # Use this as the merge
        os.system("g2munge - rand 9999 10000 2.5e-4 1 2 > ./tpmfp/tpmfp_tests")
    
        # Create an IC containing one very heavy particle randomly
        os.system("g2munge - rand 1 10000 1.0 1 1 | g2munge - merge ./tpmfp/tpmfp_tests > ./tpmfp/tpmfp_IC")

        os.system("mpirun -n 2 ./Gadget2 Configuration.tpmfp")

        # Stash it
        os.system("mv ./tpmfp/forcetest.txt ./tpmfp/forcetest_%d.txt" % i)
else:
    print "Skipping sim runs..."


# Now make the forceerrors
while i > 0:
    
    # make stack
    first = True 
    found = False

    for line in open("./tpmfp/forcetest_%d.txt" % i):

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
                                float, float, float), line.split())]

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

    # Process the next forcetest
    i -= 1

recv.close()
print "%d Gadget2 forcetest.txt files worked down!" % int(sys.argv[1])

# Sort the output
os.system("cat ./tpmfp/tpmfp_output | sort -g > ./tpmfp/tpmfp_sorted")
print "Output sorted!"
 
# Now compute the binned RMS
i = 1
binned = open("./tpmfp/tpmfp_rms", "w", 0)

for line in open("./tpmfp/tpmfp_output", "rt"):
    (r, error) = [t(s) for t(s) in zip((float, float), line.split())]
    
    # While we are in the first bin
    if r < i*float(sys.argv[2]):
        stack.append(error)
    else:
        # We hit the next range, compute an RMS
        sum = 0.0
        N = len(stack)
        while len(stack) > 0:
            sum += stack.pop()
        
        # Record the RMS
        write(binned, "%f %f\n" % ( (i - 1)*float(sys.argv[2]) + 0.5*float(sys.argv[2]), math.sqrt(sum*sum/(N*N))))

        # Start processing next homie
        i += 1

        # Stack is empty, put the line we just read onto it
        stack.append(error)

binned.close()
print "Output binned!"
sys.exit()
