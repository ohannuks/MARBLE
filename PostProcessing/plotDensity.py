import numpy as np
import sys
import pylab as plt

import matplotlib
font = {'weight' : 'normal',
        'size'   : 22}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

r,Density = np.loadtxt(sys.argv[1])
if sys.argv[3] == "log":
  plt.loglog(r,Density,label=r'Density')
else:
  plt.plot(r,Density,label=r"Density")
plt.xlabel(r"r ($r_g$)")
plt.ylabel(r"Relative density")
plt.legend()
plt.tight_layout()
plt.savefig(sys.argv[2],bbox_tight=True)
plt.show()
