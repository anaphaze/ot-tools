#!/usr/bin/env python
# ot_unskew.py -- Compensate .mod coordinates for compression artifacts
#
# If you find this script useful for your work, please cite:
# Cai, 2018, PNAS, Cryo-ET reveals the macromolecular reorganisation of S. pombe mitotic chromosomes in vivo
# https://www.ncbi.nlm.nih.gov/pubmed/########
#
# Dependencies: python 2.7, IMOD
# Created: 20170715 (Lu Gan)
# Revised: 20180621 (LG) Revised comments
#
# To calculate the correction matrix, multiply the following matrices
# 1) Rotate so that knife marks are parallel w/ X axis
# 2) Undistort along X and Z axes
# 3) Reverse rotation
# http://scipython.com/book/chapter-6-numpy/examples/creating-a-rotation-matrix-in-numpy/
# https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations

import os, subprocess, sys, struct, math
import numpy as np

#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   angle =    np.radians(float(sys.argv[1]))
   compress = float(sys.argv[2])
   coords =   sys.argv[3]
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_unskew.py [angle] [compress] [mod.txt]"
   print "================================================================"
   print "Note 1:  Tomo has to be accurately positioned w/ section surface perpendicular to Z axis"
   print "Note 2:  Measure or estimate the following values first:" 
   print "angle:     Slicer clockwise Z rotation (positive, degrees) to align knife mark to X axis"
   print "compress:  Estimated compression ratio (greater than one)"
   print "mod.txt:   Coordinates of your template-matching hits, ascii format"
   print "----------------------------------------------------------------"
   print "This script assumes your tomo was aligned so section surface is parallel to XY plane"
   print "[Unskew] = [Rz^-1] * [Stretch] * [Rz]"
   print "[Unskew]:  overall correction matrix"
   print "[Rz]:      aligns knife marks to X axis"
   print "[Rz^-1]:   rotates back to original orientation"
   print "[Stretch]: undistorts along knife-mark and perpendicular to section surface"
   print "----------------------------------------------------------------"
   print "Output:>> mod_unskewed.txt & mod_unskewed.mod"
   sys.exit()
#----------- END User inputs -------------------------------
outname = os.path.splitext(coords)[0]
c, s = np.cos(angle), np.sin(angle)
Rz = np.matrix([[c,-s,0],[s,c,0],[0,0,1]])
Rzinv = Rz.I
Stretch = np.matrix([[compress,0,0],[0,1,0],[0,0,1/compress]])
Unskew = Rzinv*Stretch*Rz
Unskewarray = np.array(Unskew)
vectors = np.genfromtxt(coords, delimiter="").T
vectorarray = np.array(vectors)
results = np.dot(Unskewarray,vectorarray).T

#----------- Print diagnostics -----------------------------
np.set_printoptions(precision=2, suppress=True)
print "rotate %s radians" % (angle)
print "Rotation:"
print(Rz)
print "Inverse rotation:"
print(Rzinv)
print "Stretch:"
print(Stretch)
print "Unskew:"
print(Unskew)
print "vectors:"
print(vectorarray)
#print "results:"
#print(results)

#np.savetxt("unskewed.txt", results, fmt='%1.2f', delimiter=" ")
np.savetxt("%s_unskewed.txt" % (outname), results, fmt='%1.2f', delimiter=" ")
os.system('point2model %s_unskewed.txt %s_unskewed.mod' % (outname, outname))
sys.exit()
