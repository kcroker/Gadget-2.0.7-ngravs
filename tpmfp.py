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

if len(sys.argv) < 3:
    print "Usage: %s <# times> <# bins> [skip sim?]" % sys.argv[0]
    exit(1)


# So the variable persists
i = 0
stack = []

# Do it 
# If there is a 3rd argument, skip running the sim and just work down the data
if not len(sys.argv) > 3: 
    for i in range(int(sys.argv[1])):

        # Randomly pick a spot
        
        # Use this as the merge
        os.system("g2munge - rand 4999 10000 2.5e-4 1 2 > ./tpmfp/tpmfp_tests")
        
        
        # Create an IC containing one very heavy particle randomly
        os.system("g2munge - rand 1 10000 1.0 1 1 | g2munge - merge ./tpmfp/tpmfp_tests > ./tpmfp/tpmfp_IC")

        # The above code does not probe sufficiently densly around the source mass.
        # So what we really need to do is write a g2munge function "wrap"
        # that takes an initial condition and wraps every coordinate to the BoxSize
        # We can then use the previous IC construction with xlate and wrap to 
        # get the arbitrary center.
        
        os.system("mpirun -n 2 ./Gadget2 Configuration.tpmfp")

        # Stash it
        os.system("mv ./tpmfp/forcetest.txt ./tpmfp/forcetest_%d.txt" % i)
else:
    # Set i to the desired number of files (minus 1 of course)
    i = int(sys.argv[1]) - 1
    print "Skipping sim runs, looking for %d forcetest files..." % (i+1)

# Now make the forceerrors
# Unbuffered
recv = open("./tpmfp/tpmfp_output", 'w', 0)

def processStack(found):
    # We have transitioned, compute the forces from this batch
    if not found:
        print "Never found the heavy guy in run %d..." % i
        exit(1)

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

while i > -1:
    
    # make stack
    first = True 
    found = False

    for n,line in enumerate(open("./tpmfp/forcetest_%d.txt" % i)):

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
            processStack(found)
            
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
                
    # We still have the last segment on the stack
    processStack(found)

    # Process the next forcetest
    i -= 1

recv.close()
print "%d Gadget2 forcetest.txt files worked down!" % int(sys.argv[1])

# Sort the output
os.system("sort -g < ./tpmfp/tpmfp_output > ./tpmfp/tpmfp_sorted")
print "Output sorted!"
 
# Now compute the binned RMS, unbuffered
binned = open("./tpmfp/tpmfp_rms", "w", 0)
i = 1

# Gadget-2 paper looks at things out to half the dimensionless boxsize
# sys.argv[2] is now interpreted as a desired number of bins
# x^sys.argv[2] = 10^4 --> x = pow(10^4, 1.0/sys.argv[2])
binbase = pow(10^4, 1.0/float(sys.argv[2]))
binrights = [pow(binbase, i) for i in range(0, int(sys.argv[2]))]

# Uhh, just generate a list of the binrights, and then go down the file
# and accumulate them...


for line in open("./tpmfp/tpmfp_sorted", "rt"):
    # Run the anonymous function (float, float) with the split line
    (r, error) = [t(s) for t,s in zip((float, float), line.split())]

    # Does it belong in a later bin?
    if r > binright:
        
        # How many are we waiting to pop?
        N = len(stack)
        
        # If we had things, that means we just advanced into a new bin
        if N > 0:
            # Record the RMS
            # Binned radius is a midpoint
            sum = 0.0
            while len(stack) > 0:
                sum += float(stack.pop())
            
            # Stack is now empty
            bincen = math.sqrt(binright * pow(binbase, i-1))
            binned.write("%f %f\n" % (bincen, sum/N))
        else:
            # If we don't have things, that means the previous bin is empty
            while r > binright:
                bincen = math.sqrt(binright * pow(binbase, i-1))
                binned.write("%f %f\n" % (bincen, 0))
                i += 1
                binright = pow(binbase, i)
    
    stack.append(error)

binned.close()
print "Output binned!"
sys.exit()
