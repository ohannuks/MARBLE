import numpy as np
import pylab as plt
import matplotlib
font = {'weight' : 'normal',
        'size'   : 24}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)
import sys

r = np.linspace(1.05,100,10000)
plt.figure(figsize=(16,9))
a = 0.998
FFIOVolume = -((np.sqrt(r**4)*(a**2*(a**2 + (-2 + r)*r) - (a**2 + r**2)**2))/(r**2*(a**2 + (-2 + r)*r)))
r         =r[FFIOVolume>0]
FFIOVolume=FFIOVolume[FFIOVolume>0]
plt.loglog(r,FFIOVolume,label=r"Proper FFIO volume ($a=0.998$)",lw=3, color='black')
a = 0
FFIOVolume = -((np.sqrt(r**4)*(a**2*(a**2 + (-2 + r)*r) - (a**2 + r**2)**2))/(r**2*(a**2 + (-2 + r)*r)))
r         =r[FFIOVolume>0]
FFIOVolume=FFIOVolume[FFIOVolume>0]
plt.loglog(r,FFIOVolume,label=r"Proper FFIO volume ($a=0.0  $)",lw=3, color='red')
plt.ylim([0,10**3])
plt.grid(ls='--')
plt.ylabel(r"Proper volume ($r_g^3$)")
plt.xlabel(r"Radius ($r_g$)")
plt.legend(frameon=False)
plt.savefig(sys.argv[1],bbox_tight=True)
plt.show()
