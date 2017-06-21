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


plt.figure(figsize=(16,9))
plt.axvline(1,ls='--',color='black',label="Kerr LSO")
plt.axvline(4,ls='--',color='red',label="Schwarzschild LSO")
plt.axvline(6,ls='--',color='blue',label="Schwarzschild ISCO")

r,Ecom = np.loadtxt(sys.argv[1])
r = r[Ecom>0.01]
Ecom=Ecom[Ecom>0.01]
label=r"Max $E_{\rm COM}$ for a=0"
color='red'
if sys.argv[4] == "log":
  plt.loglog(r,Ecom,'.',markersize=10,label=label,color=color)
else:
  plt.plot(r,Ecom,'.',markersize=10,label=label,color=color)
r,Ecom = np.loadtxt(sys.argv[2])
r = r[Ecom>0.01]
Ecom=Ecom[Ecom>0.01]
label=r"Max $E_{\rm COM}$ for a=0.998"
color='black'
if sys.argv[4] == "log":
  plt.loglog(r,Ecom,'x',markersize=10,label=label,color=color)
else:
  plt.plot(r,Ecom,'x',markersize=10,label=label,color=color)

plt.xlabel(r"r ($r_g$)")
plt.ylabel(r"$E_{\rm com}$")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(sys.argv[3],bbox_tight=True)
plt.show()
