import numpy as np
import sys
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

if len(sys.argv) > 4:
  label = sys.argv[4]
else:
  label = r'$E_{\rm COM}$'

r,Density = np.loadtxt(sys.argv[1])
plt.figure(figsize=(16,9))
if sys.argv[3] == "log":
  plt.loglog(r,Density,'x',label=label)
else:
  plt.plot(r,Density,'x',label=label)
plt.xlabel(r"r ($r_g$)")
plt.ylabel(r"$E_{\rm com}$")
#plt.ylim([3.99,5.51])
plt.legend()
plt.tight_layout()
plt.savefig(sys.argv[2],bbox_tight=True)
plt.show()
