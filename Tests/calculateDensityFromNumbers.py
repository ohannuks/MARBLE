import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR


print "python calculateDensityFromNumbers.py datafile"
fname = sys.argv[1]
fout = sys.argv[2]
data =pp.loadpdata(fname)
Eps1 = data["Eps_1"]
Lz1  = data["Lz_1"]
K1   = data["K_1"]
M    = data["M_1"]
a    = data["a_1"]
rMin1= data["rMin_1"]
rMax1= data["rMax_1"]
g1   = data["g_1"]
if GR.tests.isBound(Eps1,Lz1,K1,M,a) == False:
  raise ValueError("Eps1,Lz1,Q1 not bound; something wrong")
if GR.tests.isValidRLimits(Eps1,Lz1,K1,M,a,rMin1,rMax1) == False:
  raise ValueError("Invalid rmin and rmax inputted by Gyoto")
# Calculate the density using weighted f
f = 1. / g1
r1 = np.linspace(2,100,100)
Density1 = np.array([GR.FFIOdensity2(f, Eps1, Lz1, K1, M, a, rMin1, rMax1, ri) for ri in r1])
np.savetxt(fout, np.array([r1,Density1]))
