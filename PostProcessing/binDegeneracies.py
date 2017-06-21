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

print "python binDegeneracies.py ProductionData/M10_r2_spin0998"
plt.figure(figsize=(16,9))

print np.array(sys.argv)[1:]
for fname in np.array(sys.argv)[1:]:
  data =pp.loadpdata(fname)
  g_1    = np.ravel(data["g_1"])
  g_2    = np.ravel(data["g_2"])
  Ntotal    = np.ravel(data["Ntotal"])
  plt.hist(g_1,color='red',label="Degeneracy before growth, fname = "  + str(fname),alpha=0.5)
  plt.hist(g_2,color='black',label="Degeneracy after growth, fname = " + str(fname),alpha=0.5)
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
  print "Ntotal: " + str(Ntotal)
  N_1 = np.ones(len(g_1))*1. / Ntotal
  print N_1
  f_1 = N_1 / g_1
  f_2 = f_1
  N_2 = f_2 * g_2
  plt.hist(N_1,color='red',label="Number count before growth, fname = "  + str(fname),alpha=0.5,bins=100)
  plt.hist(N_2,color='black',label="Number count after growth, fname = " + str(fname),alpha=0.5,bins=100)
  plt.yscale("log")
  plt.xscale("log")
  plt.ylabel("Counts")
  plt.xlabel("Number count")
  plt.legend(frameon=False)
  plt.tight_layout()
plt.show()



