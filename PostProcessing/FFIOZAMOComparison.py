import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR
import matplotlib

print "python calculateDensity.py Data/samples.npz Plots/out.pdf"
if sys.argv>1:
  fname = sys.argv[1]
else:
  fname = "Data/samples"
if sys.argv>2:
  fout = sys.argv[2]
else:
  fout = "Plots/out.pdf"

R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt(fname)
import matplotlib

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 22}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

FFIODensity = FFIODensity/np.max(Density)
Density     = Density/np.max(Density)

plt.figure(figsize=(16,9))
plt.errorbar(R/500.,FFIODensity,fmt='.',label=r'FFIO Density', yerr=FFIODensity/np.sqrt(nSamples))
plt.errorbar(R/500.,Density,    fmt='.',label=r"ZAMO Density", yerr=Density/np.sqrt(nSamples))
plt.xlim([0,10])
plt.legend(frameon=False)
plt.xlabel(r"r ($r_g$)")
plt.ylabel(r"Relative Density")
plt.tight_layout()
plt.savefig(fout)
plt.show()

