#!/usr/bin/env python
# ot_rand3Dcoord_rect.py -- Create randomly distributed distribution of points in a rectangular box
# Created: 20171214 (Lu Gan)
#
# If you find this script useful for your work, please cite:
# Cai, 2018, MBoC, Natural chromatin is heterogeneous and self-associates in vitro
# https://www.ncbi.nlm.nih.gov/pubmed/29742050
#
import random, os, sys
import numpy as np
#
# Code adapted from this SO thread:
# https://stackoverflow.com/a/19668720
# https://stackoverflow.com/questions/8466014/how-to-convert-a-python-set-to-a-numpy-array
#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   sizeX  = int(sys.argv[1])
   sizeY  = int(sys.argv[2])
   sizeZ  = int(sys.argv[3])
   npart  = int(sys.argv[4])
   radius = int(sys.argv[5])
   cenX   = int(sys.argv[6])
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_rand3Dcoord_rect.py [dimX] [dimY] [dimZ] [npart] [rad] [cenX]"
   print "================================================================"
   print "dimX/Y/Z : box's X/Y/Z axes lengths, in pixels"
   print "npart:     number of particles"
   print "rad:       size of particles, for overlap exclusion"
   print "cenX:      move center of mass along X axes this much"
   print "----------------------------------------------------------------"
   print "Creates a random array of points within a rectangular box"
   print "Note that some combos might never converge:"
   print "Box too big, radius too big, too many points"
   print "----------------------------------------------------------------"
   print "Output:>> out.txt & out.mod"
   sys.exit()
#----------- END User inputs -------------------------------
rangeX = (0, sizeX)
rangeY = (0, sizeY)
rangeZ = (0, sizeZ)

# Generate a set of all points within radius of the origin, to be used as offsets later
deltas = set()
for x in range(-radius, radius+1):
   for y in range(-radius, radius+1):
       for z in range(-radius, radius+1):
          if x*x + y*y + z*z <= radius*radius:
             deltas.add((x,y,z))

randPoints = []
excluded = set()
i = 0

while i<npart:
   x = random.randrange(*rangeX)
   y = random.randrange(*rangeY)
   z = random.randrange(*rangeZ)
   if (x,y,z) in excluded: continue
   randPoints.append((x,y,z))
   i += 1
   excluded.update((x+dx, y+dy, z+dz) for (dx,dy,dz) in deltas)
results = np.array(list(randPoints))
offX = np.array([cenX,0,0])
results += offX

print(results)
#print(deltas)
np.savetxt("out.txt", results, fmt='%d', delimiter=" ")
os.system('point2model -scat -sphere 6 out.txt out.mod')
#os.system('point2model -scat -sphere %s out.txt out.mod' % (radius))
sys.exit()
