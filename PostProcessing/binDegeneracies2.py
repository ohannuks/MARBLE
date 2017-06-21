import sys
import matplotlib
#font = {'family' : 'sans-serif',
#        'weight' : 'normal',
#        'size'   : 22}        # Set font
#matplotlib.rc('font', **font) # Set font
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
#matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rc('text', usetex=True)
import numpy as np
import paris as pp
from paris import GR
import pylab as plt

print "python binDegeneracies2.py ProductionData/M10_r2_spin0998"

plt.figure(figsize=(16,9))
for fname in np.array(sys.argv)[1:]:
  data =pp.loadpdata(fname)
  g_2    = np.ravel(data["g_2"])
  Ntotal    = np.ravel(data["Ntotal"])
  rMin   = data["rMin_2"]
  rMax   = data["rMax_2"]
  M      = data["M_2"]
  plt.hist(g_2[rMin < 4.*M],color='black',label="Degeneracy after growth r<4M, fname = " + str(fname),alpha=0.2)
  plt.hist(g_2[(rMin < 8.*M)&(rMin>4.*M)],label="Degeneracy after growth 4M<r<8M, fname = " + str(fname),alpha=0.2,color='black')
  #plt.hist(g_2[rMax < 12.*M],label="Degeneracy after growth r<12M, fname = " + str(fname),alpha=0.2,color='black')
  #plt.hist(g_2[rMax < 16.*M],label="Degeneracy after growth r<16M, fname = " + str(fname),alpha=0.2,color='black')
  plt.yscale("log")
  plt.xscale("log")
  plt.ylabel("Counts")
  plt.xlabel("Degeneracy")
  plt.legend(frameon=False)
  plt.tight_layout()
plt.show()

plt.figure(figsize=(16,9))
for fname in np.array(sys.argv)[1:]:
  data =pp.loadpdata(fname)
  g_1    = np.ravel(data["g_1"])
  g_2    = np.ravel(data["g_2"])
  Ntotal    = np.ravel(data["Ntotal"])
  N_1    = np.ones(len(g_1))*1./float(Ntotal)
  f_1    = N_1 / g_1
  f_2    = f_1
  N_2    = f_2 * g_2
  rMin   = data["rMin_2"]
  rMax   = data["rMax_2"]
  M      = data["M_2"]
  plt.hist(N_2[rMin < 4.*M],color='black',label="Degeneracy after growth r<4M, fname = " + str(fname),alpha=0.2)
  plt.hist(N_2[(rMin < 8.*M)&(rMin>4.*M)],label="Degeneracy after growth 4M<r<8M, fname = " + str(fname),alpha=0.2,color='black')
  #plt.hist(g_2[rMax < 12.*M],label="Degeneracy after growth r<12M, fname = " + str(fname),alpha=0.2,color='black')
  #plt.hist(g_2[rMax < 16.*M],label="Degeneracy after growth r<16M, fname = " + str(fname),alpha=0.2,color='black')
  plt.yscale("log")
  plt.xscale("log")
  plt.ylabel("Counts")
  plt.xlabel(r"$\mathcal{N}$")
  plt.legend(frameon=False)
  plt.tight_layout()
plt.show()

