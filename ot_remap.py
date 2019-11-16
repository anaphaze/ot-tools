#!/usr/bin/env python
# ot_remap.py -- remap subtomos into a synthetic tomogram
#
# If you find this script useful for your work, please cite:
# Cai, 2018, MBoC, Natural chromatin is heterogeneous and self-associates in vitro
# https://www.ncbi.nlm.nih.gov/pubmed/29742050
#
# Dependencies: EMAN2, IMOD, Bsoft
# Created: 20161021 (Shujun Cai)
# Revised: 20161227 (Lu Gan) merged scripts, increased user-friendliness
# Revised: 20161231 (SC & LG) e2proc3d uses SPIDER Euler convention; read _data.star & class ID
# Revised: 20180117 (LG) added origin shifts
# Revised: 20190509 (LG) ignore # comments; fixed range(sub1,sub2) out of index bug
#
# Based on Tanmay Bharat's relion_2Dto3D_star.py:
# http://www.sciencedirect.com/science/article/pii/S0969212615002798
# http://www2.mrc-lmb.cam.ac.uk/relion/index.php/Sub-tomogram_averaging
#
# Algorithm:
# 1) Rotate all subtomograms  (Euler = SPIDER "ZYZ" and rotates average into subtomo)
# 2) Enlarge all subtomograms
# 3) Translate particles to correct positions
# 4) Add subsets of transformed subtomograms
# 5) Merge all tomograms with IMOD's clip (not parallelized)

import os, subprocess, csv, sys, multiprocessing, struct, time, math

#---------- BEGIN User inputs ------------------------------
# stackoverflow.com/questions/2194163/python-empty-argument
try:
   name_avg  = sys.argv[1]
   name_tomo = sys.argv[2]
   name_star = sys.argv[3]
   name_cls  = sys.argv[4]
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_remap.py [Average] [Tomogram] [Starfile] [ClassID]"
   print "----------------------------------------------------------------"
   print "Average:   the class average you want to remap"
   print "Tomogram:  original tomogram that has all particles"
   print "Starfile:  _data.star file that contains all classes"
   print "ClassID:   Class number you want to remap"
   print "----------------------------------------------------------------"
   print "Output:>> syn_tomo_cls[ClassID].mrc"
   print "To invert tomogram, run 'bimg -invert positive.mrc negative.mrc'"
   sys.exit()
#----------- END User inputs -------------------------------


#----------- BEGIN File and system parameters --------------
# Get x,y,z dimensions
with open (name_tomo, 'rb') as file:
   content = file.read()
   dimX = struct.unpack('3i',content[:12])[0]
   dimY = struct.unpack('3i',content[:12])[1]
   dimZ = struct.unpack('3i',content[:12])[2]
file.close

# Get system parameters
# stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-using-python
# stackoverflow.com/questions/4271740/how-can-i-use-python-to-get-the-system-hostname
ncpu = multiprocessing.cpu_count()
myhost = os.uname()[1]
# Calculate free space; must be done after rotate, expand, & translate operations
# stackoverflow.com/questions/4260116/find-size-and-free-space-of-the-filesystem-containing-a-given-file
size_free = os.statvfs('.').f_frsize * os.statvfs('.').f_bavail
size_free_gibi = size_free / 1073741824
size_free_tebi = size_free / 1099511627776
# stackoverflow.com/questions/2104080/how-to-check-file-size-in-python
size_tomo = float(os.stat(name_tomo).st_size)
size_tomo_gibi = size_tomo / 1073741824
#---------- END File and system parameters -----------------


#---------- BEGIN RELION header function -------------------
def read_relion_header(filename):

  ifile = open(filename, 'r')

  # first get the header
  i, XextCol, YextCol, ZextCol, RotCol, TiltCol, PsiCol, ClassCol, XoriCol, YoriCol, ZoriCol = -1,0,0,0,0,0,0,0,0,0,0

  for line in ifile:
    emptycheck = line.isspace()     #skip whitespace, usually lines 1 and 3
    if(emptycheck):
      continue

    fields = line.split()
    if fields[0] == 'data_' or fields[0] == 'loop_' or fields[0] == '#':  # skip the data_, loop_, # lines
      continue

    i= i+1

    firstcol = fields[0]

    if firstcol == '_rlnCoordinateX':
      XextCol = i
    if firstcol == '_rlnCoordinateY':
      YextCol = i
    if firstcol == '_rlnCoordinateZ':
      ZextCol = i
    if firstcol == '_rlnAngleRot':
      RotCol = i
    if firstcol == '_rlnAngleTilt':
      TiltCol = i
    if firstcol == '_rlnAnglePsi':
      PsiCol = i
    if firstcol == '_rlnClassNumber':
      ClassCol = i
    if firstcol == '_rlnOriginX':
      XoriCol = i
    if firstcol == '_rlnOriginY':
      YoriCol = i
    if firstcol == '_rlnOriginZ':
      ZoriCol = i

    if firstcol[0] != '_':
      break

  return(XextCol,YextCol,ZextCol,RotCol,TiltCol,PsiCol,ClassCol,XoriCol,YoriCol,ZoriCol)
#---------- END RELION header function ---------------------


#---------- BEGIN Transformation thread function -----------
def clipadd(j):                # j indexes each core
   sub1=j*num_task             # first iter index
   sub2=(j+1)*num_task-1       # last iter index
   for i in range(sub1,sub2):
      if (i == sub1 and i < (num_line-1)):     # Initialize 1st sum, making sure not to overcount
# Rotate i == sub1
         a=float(rot[i])
         b=float(tlt[i])
         c=float(psi[i])
         os.system('e2proc3d.py --rot=spider:phi=%d:theta=%d:psi=%d %s tmp_ema_%05d.mrc' % (a, b, c, name_avg, i))
         os.system('bimg -rescale 0,5.2 tmp_ema_{0:05d}.mrc tmp_sca_{0:05d}.mrc'.format(i))
         os.system('newstack -q -mode 1 tmp_sca_{0:05d}.mrc tmp_mod_{0:05d}.mrc'.format(i))
         os.system('bimg -rescale 0,5.2 tmp_mod_{0:05d}.mrc tmp_rot_{0:05d}.mrc'.format(i))
#        os.system('touch ycount%05d' % i)  # Uncomment for diagnostics

         a=float(rot[i+1])
         b=float(tlt[i+1])
         c=float(psi[i+1])
         os.system('e2proc3d.py --rot=spider:phi=%d:theta=%d:psi=%d %s tmp_ema_%05d.mrc' % (a, b, c, name_avg, i+1))
         os.system('bimg -rescale 0,5.2 tmp_ema_{0:05d}.mrc tmp_sca_{0:05d}.mrc'.format(i+1))
         os.system('newstack -q -mode 1 tmp_sca_{0:05d}.mrc tmp_mod_{0:05d}.mrc'.format(i+1))
         os.system('bimg -rescale 0,5.2 tmp_mod_{0:05d}.mrc tmp_rot_{0:05d}.mrc'.format(i+1))
#        os.system('touch ycount%05d' % (i+1))  # Uncomment for diagnostics

# Resize i == sub1
         os.system('clip resize -m 1 -ox %s -oy %s -oz %s tmp_rot_%05d.mrc tmp_siz_%05d.mrc' % (dimX, dimY, dimZ, i, i))
         os.system('clip resize -m 1 -ox %s -oy %s -oz %s tmp_rot_%05d.mrc tmp_siz_%05d.mrc' % (dimX, dimY, dimZ, i+1, i+1))
         os.system('rm -f tmp_ema_{0:05d}.mrc tmp_sca_{0:05d}.mrc tmp_mod_{0:05d}.mrc tmp_rot_{0:05d}.mrc tmp_ema_{1:05d}.mrc tmp_sca_{1:05d}.mrc tmp_mod_{1:05d}.mrc tmp_rot_{1:05d}.mrc'.format(i, i+1))

# Translate i == sub1
         dx=dimX - float(x[i])
         dy=dimY - float(y[i])
         dz=dimZ - float(z[i])
         os.system('clip resize -m 1 -cx %d -cy %d -cz %d tmp_siz_%05d.mrc tmp_tra_%05d.mrc' % (dx, dy, dz, i, i))
         os.system('rm -f tmp_siz_%05d.mrc' % (i))

         dx=dimX - float(x[i+1])
         dy=dimY - float(y[i+1])
         dz=dimZ - float(z[i+1])
         os.system('clip resize -m 1 -cx %d -cy %d -cz %d tmp_siz_%05d.mrc tmp_tra_%05d.mrc' % (dx, dy, dz, i+1, i+1))
         os.system('rm -f tmp_siz_%05d.mrc' % (i+1))

# Combine i == sub1
         os.system('clip add -m 1 tmp_tra_%05d.mrc tmp_tra_%05d.mrc sum_%05d_%05d.mrc' % (i, i+1, i, i+1))
         os.system('rm -f tmp_tra_%05d.mrc tmp_tra_%05d.mrc' % (i, i+1))

      elif (i > sub1 and i < (num_line-1)):  # Sequentially append to previous sums; don't overcount
# Rotate i > sub1
         a=float(rot[i+1])
         b=float(tlt[i+1])
         c=float(psi[i+1])
         os.system('e2proc3d.py --rot=spider:phi=%d:theta=%d:psi=%d %s tmp_ema_%05d.mrc' % (a, b, c, name_avg, i+1))
         os.system('bimg -rescale 0,5.2 tmp_ema_{0:05d}.mrc tmp_sca_{0:05d}.mrc'.format(i+1))
         os.system('newstack -q -mode 1 tmp_sca_{0:05d}.mrc tmp_mod_{0:05d}.mrc'.format(i+1))
         os.system('bimg -rescale 0,5.2 tmp_mod_{0:05d}.mrc tmp_rot_{0:05d}.mrc'.format(i+1))
#        os.system('touch zcount%05d' % (i+1))  # Uncomment for diagnostics
# Resize i > sub1
         os.system('clip resize -m 1 -ox %s -oy %s -oz %s tmp_rot_%05d.mrc tmp_siz_%05d.mrc' % (dimX, dimY, dimZ, i+1, i+1))
         os.system('rm -f tmp_ema_{0:05d}.mrc tmp_sca_{0:05d}.mrc tmp_mod_{0:05d}.mrc tmp_rot_{0:05d}.mrc'.format(i+1))

# Translate i > sub1
         dx=dimX - float(x[i+1])
         dy=dimY - float(y[i+1])
         dz=dimZ - float(z[i+1])
         os.system('clip resize -m 1 -cx %d -cy %d -cz %d tmp_siz_%05d.mrc tmp_tra_%05d.mrc' % (dx, dy, dz, i+1, i+1))
         os.system('rm -f tmp_siz_%05d.mrc' % (i+1))

# Combine i > sub1
         os.system('clip add -m 1 tmp_tra_%05d.mrc sum_%05d_%05d.mrc sum_%05d_%05d.mrc' % (i+1, sub1, i, sub1, i+1))
         os.system('rm -f sum_%05d_%05d.mrc' % (sub1, i))
         os.system('rm -f tmp_tra_%05d.mrc' % (i+1) )
#---------- END Transformation thread function -------------


#---------- BEGIN Data handling ----------------------------
# Get column numbers
XextCol,YextCol,ZextCol,RotCol,TiltCol,PsiCol,ClassCol,XoriCol,YoriCol,ZoriCol = read_relion_header(name_star)

# Read the 3D starfile, ignoring the header items;  "grep -v" is an inverse selection
grepline = 'grep ' + 'Tomograms ' + name_star + '| grep -v "#" | grep -v "loop" | grep -v "data" | awk "NF"  > temp.txt'
os.system(grepline)

### Make database of only the relevant values in memory
data = csv.reader(open('temp.txt', 'r'), delimiter=" ", skipinitialspace = True)
rot, tlt, psi, x, y, z, cls = [], [], [], [], [], [], []
num_line = 0
for row in data:
   if row[ClassCol] == name_cls:
      x.append(float(row[XextCol])-float(row[XoriCol]))
      y.append(float(row[YextCol])-float(row[YoriCol]))
      z.append(float(row[ZextCol])-float(row[ZoriCol]))
      rot.append(row[RotCol])
      tlt.append(row[TiltCol])
      psi.append(row[PsiCol])
      cls.append(row[ClassCol])  # Unessential, but keep for future diagnostics
      num_line = num_line + 1
#---------- END Data handling ------------------------------


#---------- BEGIN Bookkeeping calculations -----------------
size_free = os.statvfs('.').f_frsize * os.statvfs('.').f_bavail
num_file = int(math.floor(size_free/(3*size_tomo)))
num_jobs = min(ncpu, num_file)
num_task = int(math.ceil(float(num_line)/num_jobs))  # tasks per core
print ""'"%s"'" has %s members in class %s " % (name_avg, num_line, name_cls)
print ""'"%s"'" is %s x %s x %s and uses %.2f GB" % (name_tomo, dimX, dimY, dimZ, size_tomo_gibi)
print "%s has %s cores and this directory has %s GB free" % (myhost, ncpu, size_free_gibi)
print "There will be %s parallel jobs handling %s subtomograms each" % (num_jobs, num_task)
time.sleep(3)
#---------- END Bookkeeping calculations -------------------


#---------- BEGIN Parallel work ----------------------------
sizadd_time1 = time.time()

if __name__ == '__main__':
      pool = multiprocessing.Pool(num_jobs)
      results = pool.map(clipadd, range(num_jobs))
#     results = pool.map(clipadd, range(17,20))  #Uncomment for diagnostics

sizadd_time2 = time.time()
name_cls = int(name_cls)
os.system('clip add -m 1 sum*.mrc syn_pos_cls%02d.mrc' % (name_cls))
os.system('alterheader syn_pos_cls%02d.mrc -org 0,0,0' % (name_cls))
# Uncomment to get synthetic tomo w/ inverted contrast:
os.system('bimg -invert syn_pos_cls%02d.mrc syn_neg_cls%02d.mrc' % (name_cls, name_cls))

os.system('rm -f temp.txt *.mrc~ sum_*_*.mrc')  # Comment out for diagnostics
#---------- END Parallel work ------------------------------


#---------- Sum all to compare with tomo --------------------
os.system('clip add -m 1 syn_pos*mrc syn_sum_pos.mrc')
os.system('clip add -m 1 syn_neg*mrc syn_sum_neg.mrc')

print("--------------------------------------------------------")
time.sleep(1)
print("Processing time: %01d seconds" % (sizadd_time2 - sizadd_time1))
print ""'"%s"'" has %s members in class %s " % (name_avg, num_line, name_cls)
print ""'"%s"'" is %s x %s x %s and uses %.2f GB" % (name_tomo, dimX, dimY, dimZ, size_tomo_gibi)
print "%s has %s cores and this directory has %s GB free" % (myhost, ncpu, size_free_gibi)
print "There were %s parallel jobs handling %s subtomograms each" % (num_jobs, num_task)
print "done"
sys.exit()
