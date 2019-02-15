#SJL 11/15
#Script to plot a given structure and print out its essential features

#import required modules
import numpy as np
import scipy as sp
from scipy import constants as const
import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab
import matplotlib.gridspec as gridspec

#Here you will need to point to HERCULES structures on your machine
sys.path.append('/Users/simonlock/Dropbox/HERCULES/Development/Python_scripts')
from HERCULES_structures import *

###############################################################
###############################################################
###############################################################
#PARAMS

#the path to the structure file you want to plot
fname="/Users/simonlock/Dropbox/HERCULES/Development/HERCULES_v1.0/Output/tutorial1_L0_N50_Nm400_k6_f020_p10_l1_0_1.5_final"

###############################################################
#define blank planet and parameter classes
params=HERCULES_parameters()
p=HERCULES_planet()

#read in HERCULES output
file = open(fname, "rb")
params.read_binary(file)
p.read_binary(file)
file.close()

###############################################################
#initialise the figure
fig = plt.figure(figsize=(3.5,6.0))
gs = gridspec.GridSpec(2, 1,
                       height_ratios=[1.22,1]
                       )
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1], sharex=ax1)


plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rcParams.update({'font.size': 8})

mpl.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',    # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]  

#plot the layers of the HERCULES structure
ax1.hold(True)
for i in np.arange(0, p.Nlayer):
    ax1.plot(p.layers[i].xi*p.layers[i].a*np.sqrt(np.absolute(1-p.layers[i].mu**2))/1E6, p.layers[i].xi*p.layers[i].a*p.layers[i].mu/1E6, '-', color='k', linewidth=1.5)

ax1.set_xlim([0,7])
ax1.set_aspect('equal', adjustable='box-forced')
plt.setp( ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('Distance from midplane [Mm]')
ax1.hold(False)


#plot the pressure profile
radius=np.asarray([])
press=np.asarray([])
for i in np.arange(p.Nlayer):
    radius=np.append(radius, np.asarray(p.layers[i].a/1.0E6))
    press=np.append(press, np.asarray(p.press[i]/1.0E5))
radius=np.append(radius, np.asarray(0.0))
press=np.append(press, np.asarray(p.pcore/1.0E5))

ax2.hold(True)
ax2.plot(radius, press, '-', color='k', linewidth=1.5)

ax1.set_xlim([0,7])
ax2.set_yscale('log')
ax2.set_xlabel('Radius [Mm]')
ax2.set_ylabel('Pressure [bar]')
ax2.hold(False)

fig.tight_layout()
plt.show()

#if you want to save the figure
#plt.savefig('Tutorial1_structure.pdf', dpi=400)
