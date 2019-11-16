#!/usr/bin/env python
# ot_rot-ps.py -- Generate a rotational correlation and power spectrum from an image
#
# If you find this script useful for your work, please cite:
# Ng, 2019, JCB, "Electron cryotomography analysis of Dam1C/DASH at the kinetochore–spindle interface in situ"
# https://www.ncbi.nlm.nih.gov/pubmed/30504246
#
# Dependencies: EMAN2
# Created: 20170510 (Lu Gan)
# Revised: 20171211 tightened comments (LG)
# Revised: 20191116 updated reference (LG)
#
import os
from EMAN2 import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft

# http://blake.bcm.edu/emanwiki/FAQ_EMAN_USING_20
# https://matplotlib.org/users/pyplot_tutorial.html
# http://stackoverflow.com/questions/38812611/numpys-fast-fourier-transform-yields-unexpected-results
# http://stackoverflow.com/questions/24943991/matplotlib-change-grid-interval-and-specify-tick-labels

#---------- BEGIN User inputs ------------------------------
try:
   name_img = sys.argv[1]
   max_sym  = int(sys.argv[2])
except IndexError:
   print "================================================================"
   print "Usage:>>   ot_rot-ps.py [Image] [max_sym]"
   print "----------------------------------------------------------------"
   print "Image:     Grayscale 2-D image, PNG or TIFF format"
   print "max_sym:   Maximum rotational symmetry"
   print "----------------------------------------------------------------"
   print "Output:>>  corr_[Image].txt  (rotational correlation, 0 - 360 degrees)"
   print "Output:>>  ps_[Image].txt  (rotational power spectrum, 0 - max_sym)"
   print "Output:>>  A graphical plot of the rotational power spectrum"
   sys.exit()
name_root, name_extension = os.path.splitext(name_img)
#----------- END User inputs -------------------------------

# Calculate CCCs with EMAN2
a=EMData(name_img,0)
a.process_inplace("normalize.edgemean")
out=open("corr_%s.txt" % name_root,"w")

for angle in range(0,360):
    b=a.copy()
    b.rotate(float(angle),0,0)
    out.write("%d\t%f\n"%(angle,a.cmp("ccc",b,{"negative":0})))
out.close()

# Number of samplepoints and spacing
N = 360
T = 1.0 / N

# Input by numpy, using column 1
y = np.genfromtxt("corr_%s.txt" % name_root, dtype = float, delimiter = '\t', usecols = 1)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

# Power spectrum by numpy
yf = np.abs(np.fft.fft(y, axis=0))**2

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# major ticks every 10, minor ticks every 1
major_ticks = np.arange(0, max_sym, 10)
minor_ticks = np.arange(0, max_sym, 1)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)

# Minor ticks have lighter shading
ax.grid(which='minor', alpha=0.5)
ax.grid(which='major', alpha=1)

ax.plot(xf, (2.0/N)*(yf[0:N//2]))
ax.set_xlim(0,max_sym)

# Save FFT as text file
np.savetxt("ps_%s.txt" % name_root, yf)

# Open up graphical plot
plt.title('Rotational power spectrum')
plt.ylabel('Power')
plt.xlabel('Symmetry')
plt.grid()
plt.show()

sys.exit()
