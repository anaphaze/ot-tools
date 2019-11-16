#!/usr/bin/env python
# ot_nnd.py -- Nearest neighbor distances from a set of 3-D coordinates
#
# If you find this script useful for your work, please cite:
# Cai, Tan, 2019, bioRxiv, Structural and biochemical changes of G0 S. pombe chromatin
# https://www.biorxiv.org/content/10.1101/######
#
# Created 20180101 (Lu Gan)
#
# This script is beased on the following resources:
# https://scikit-learn.org/stable/modules/neighbors.html  (Simple NND example)
# https://stackoverflow.com/a/2828121                     (How to sort)

import os, sys
import numpy as np
from sklearn.neighbors import NearestNeighbors

#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   coords = sys.argv[1]
   K_th   = int(sys.argv[2])
   pixel  = float(sys.argv[3])
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_nnd.py [coords] [K] [pix]"
   print "================================================================"
   print "coords:  File w/ text 3-D coordinates (X Y Z), from model2point"
   print "K:       the Kth nearest neighbor"
   print "pixel:   pixel size, in nanometers"
   print "Example: ot_nnd.py tm_hits.txt 10 0.91  (Gets tenth nearest neighbors)"
   print "----------------------------------------------------------------"
   print "Output:>>  nn_indx.txt: indices of each point and its nearest neighbor"
   print "Output:>>  nn_dist.txt: list of NN distances, matched to indices"
   print "Output:>>  nn_sort.txt: list of NN distances, small to large"
   sys.exit()
#----------- END User inputs -------------------------------
outname = os.path.splitext(coords)[0]
vectors = np.genfromtxt(coords, delimiter="")
vectorarray = np.array(vectors)
K_corr = K_th + 1

nbrs = NearestNeighbors(n_neighbors=K_corr, algorithm='brute').fit(vectorarray)
results, indices = nbrs.kneighbors(vectorarray)

dist = np.delete(results, np.s_[:K_th:], 1)
dist = np.multiply(dist, pixel)

distsort = dist[dist[:,0].argsort()]

#----------- Print diagnostics -----------------------------
#np.set_printoptions(precision=2, suppress=True)
#print "Coordinates:"
#print(vectors)
#print "vectors:"
#print(vectorarray)
#print "indices:"
#print(indices)
#print(distsort)

#np.savetxt("nn_all.txt", results, fmt='%1.2f', delimiter=" ")
np.savetxt("nn_indx.txt", indices, fmt='%d', delimiter=" ")
np.savetxt("nn_dist.txt", dist, fmt='%1.2f', delimiter=" ")
np.savetxt("nn_sort.txt", distsort, fmt='%1.2f', delimiter=" ")
sys.exit()
