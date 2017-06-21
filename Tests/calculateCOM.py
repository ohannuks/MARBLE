import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR


print "python calculateCOM.py Data/samples.npz Plots/out.pdf <0 or 1 for calculating mean (0) or max (1) energy>"
fname = sys.argv[1]
fout = sys.argv[2]
if len(sys.argv) > 3:
  maxenergy = bool(int(sys.argv[3]))
else:
  maxenergy = False

# Calculate density from the data file
data =pp.loadpdata(fname)
Eps_1  = data["Eps_1"]
Lz_1   = data["Lz_1"]
K_1    = data["K_1"]
M_1    = data["M_1"]
a_1    = data["a_1"]
rMin_1 = data["rMin_1"]
rMax_1 = data["rMax_1"]
thMin_1 = data["thMin_1"]
thMax_1 = data["thMax_1"]
g_1    = data["g_1"]

Eps_2  = data["Eps_2"]
Lz_2   = data["Lz_2"]
K_2    = data["K_2"]
M_2    = data["M_2"]
a_2    = data["a_2"]
rMin_2 = data["rMin_2"]
rMax_2 = data["rMax_2"]
thMin_2 = data["thMin_2"]
thMax_2 = data["thMax_2"]
g_2    = data["g_2"]

if GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1) == False:
  raise ValueError("Eps_1,Lz_1,Q_1 not bound; something wrong")
if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  raise ValueError("Eps_2,Lz_2,Q_2 not bound; something wrong")

# Calculate the density using weighted f
r = np.linspace(1.*M_2,10*M_2,100)
COMs = []
f_1 = 1. / g_1
f_2 = f_1 # Constancy of phase space
th = np.pi/2.
for ri in r:
  COM = []
  for i in range(100):
    trial = GR.MeanEcom( f_2, Eps_2, Lz_2, K_2, rMin_2, rMax_2, thMin_2, thMax_2,f_2, Eps_2, Lz_2, K_2, rMin_2, rMax_2, thMin_2, thMax_2,M_2, a_2,  ri, th, maxenergy )
    COM.append(trial)
  if maxenergy == True: # Record maximum energy
    COMs.append(np.max(COM))
  else:
    COMs.append(np.mean(COM))
# Save in r_g units
np.savetxt(fout,np.array([r/M_2,COMs]))




