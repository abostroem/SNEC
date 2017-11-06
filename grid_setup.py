###############################################################################
#
# This file generates a pattern used for the gridding of the SNEC models
#
# Total number of grid points is equal to imax.
# Current pattern consists of two geometric progressions, connected at the
# point itran and ensuring finer resolution of the grid in the inner (explosion)
# region and the surface region of the models.
#
# ratio1 gives the ratio between the size of the first and itran-th cell
# ratio2 gives the ratio between the size of the last  and itran-th cell
#
###############################################################################

from __future__ import division
import math

delta = []
grid_pattern = []

imax = 1000
itran = 100

ratio1 = 0.1
ratio2 = 0.001

f1 = ratio1**(1/(itran-2))
f2 = ratio2**(1/(imax-itran-1))

delta_tran = ( (1-f2)*(1-f1)
                    /((1-f2)*(1-f1**(itran-1))+(1-f1)*(1-f2**(imax-itran))) )

for l in range(0,itran-1):
    delta.append(delta_tran*f1**(itran-2-l))

for l in range(0,imax-itran):
    delta.append(delta_tran*f2**(l))

grid_pattern.append(0.0)
for l in range(0,imax-1):
    grid_pattern.append(grid_pattern[l]+delta[l])

print "grid pattern consists of", len(grid_pattern), "points"

outfile = open("GridPattern_new.dat","w")
for l in range(0,imax):
    outfile.write(str(grid_pattern[l]) + '\n')
outfile.close()

print("done!")

