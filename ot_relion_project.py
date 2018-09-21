#!/usr/bin/env python
# ot_relion_project.py -- average the central N slices in a large number of tomos
#
# If you find this script useful for your work, please cite:
# Cai, 2018, PNAS, Cryo-ET reveals the macromolecular reorganisation of S. pombe mitotic chromosomes in vivo
# https://www.ncbi.nlm.nih.gov/pubmed/########
#
# Dependencies: IMOD
# Created: 20170713 (Lu Gan)

# Algorithm:
# 1) Use clip to make projections of subtomorams, N slices thick
# 2) Use newstack to create image stack from projections

import os, subprocess, sys, struct, math

#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   name_tomos = sys.argv[2]
   num_slice = int(sys.argv[1])
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_relion_project.py  [N_slices]  [Subtomos]"
   print "Run within /Extract/extract/Tomograms/tomo1"
   print "Example:   ot_relion_project.py 24 *.mrc"
   print "----------------------------------------------------------------"
   print "N_slices:  Number of slices to average"
   print "Subtomos:  Subtomograms to generate projections with"
   print "----------------------------------------------------------------"
   print "Output:>>  stack_N-slices_thick.mrcs"
   sys.exit()
#----------- END User inputs -------------------------------


#----------- BEGIN File and system parameters --------------
for name_tomo in sys.argv[2:]:
   name_root1, name_extension = os.path.splitext(name_tomo)
   name_root = os.path.basename(name_root1)
   print name_root 

   # Get max number of slices
   with open (name_tomo, 'rb') as file:
      content = file.read()
      pixZ_rec = struct.unpack('3i',content[:12])[2]

      file.seek(76)
      content = file.read()
      amin_rec = struct.unpack('3f',content[:12])[0]

   file.close
   print "Total slices: %s" %(pixZ_rec)

   # Calculate N_lower and N_upper (lower and upper slices to average)
   N_middle = int(math.ceil(pixZ_rec/2))
   if (num_slice % 2 == 0):
      N_lower = N_middle - (num_slice/2)
      N_upper = N_middle + (num_slice/2) - 1
   else:
      N_lower = int(N_middle - math.floor(num_slice/2))
      N_upper = int(N_middle + math.floor(num_slice/2))

   print "mid: %s    bot: %s    top: %s" % (N_middle, N_lower, N_upper)
   print "name %s" % (name_root)

   os.system('clip avg -2d -iz %s-%s %s tmp_%s.mrc' % (N_lower, N_upper, name_tomo, name_root))
os.system('newstack tmp*.mrc stack_%s-slices_thick.mrcs' % (num_slice))
os.system('/bin/rm tmp*.mrc *~')

print "---------------------------------------------------------------------------"
print "\n"
print "Please make backup of your original stack, then rename stack.mrcs to the relion stack"
print "Copy the new stack.mrcs file to  ../../../proj3d/Tomograms/tomo001/blah.mrcs"
print "\n"

sys.exit()
