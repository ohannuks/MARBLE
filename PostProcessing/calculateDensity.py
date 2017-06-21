import matplotlib
matplotlib.use("Agg")
import paris as pp
import numpy as np
#import pylab as plt
import sys
from paris import GR


print "python calculateDensity.py Data/samples.npz Plots/out.pdf"
if sys.argv>1:
  fname = sys.argv[1]
else:
  fname = "Data/samples"
if sys.argv>2:
  fout = sys.argv[2]
else:
  fout = "Plots/out.pdf"

# Calculate density from the data file
data =pp.loadpdata(fname)
Ntotal = data["Ntotal"]
Eps_1  = data["Eps_1"]
Lz_1   = data["Lz_1"]
K_1    = data["K_1"]
M_1    = data["M_1"]
a_1    = data["a_1"]
rMin_1 = data["rMin_1"]
rMax_1 = data["rMax_1"]
g_1    = np.ravel(data["g_1"])

Eps_2  = data["Eps_2"]
Lz_2   = data["Lz_2"]
K_2    = data["K_2"]
M_2    = data["M_2"]
a_2    = data["a_2"]
rMin_2 = data["rMin_2"]
rMax_2 = data["rMax_2"]
thMin_2 = data["thMin_2"]
thMax_2 = data["thMax_2"]
g_2    = np.ravel(data["g_2"])

#if GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1) == False:
#  raise ValueError("Eps_1,Lz_1,Q_1 not bound; something wrong")
#if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
#  raise ValueError("Eps_2,Lz_2,Q_2 not bound; something wrong")

# Calculate the density using weighted f
r = np.linspace(M_2,100*M_2,400)
N_1 = np.ones(len(Eps_1))*1. / Ntotal
f_1 = N_1 / g_1
f_2 = f_1 # Constancy of phase space
# Calculate from the 4-current J^\mu J_\nu = - \rho^2
#################################################
th   = np.pi/2.
Jmu  = np.array([GR.Nmu(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, ri, th) for ri in r])
rho2 = np.zeros(len(Jmu))
for i in range(len(Jmu)):
  gmunu = GR.gmunu(M_2, a_2, r[i], th)
  for mu in xrange(4):
    for nu in xrange(4):
      rho2[i] = rho2[i] -1.*gmunu[mu][nu] * Jmu[i][mu] * Jmu[i][nu]
Density = np.sqrt(rho2)
#################################################

# Save in r_g units
np.savetxt(fout,np.array([r/M_2,Density*M_2**3]))

