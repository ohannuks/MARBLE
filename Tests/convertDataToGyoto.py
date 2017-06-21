import matplotlib
import paris as pp
import numpy as np
import sys
from paris import GR

fname = sys.argv[1]
fout  = sys.argv[2]

data  = pp.loadpdata(fname)
Eps_2  = np.array(data["Eps_2"],dtype=np.float64)
Lz_2   = np.array(data["Lz_2"],dtype=np.float64)
K_2    = np.array(data["K_2"],dtype=np.float64)
M_2    = np.array(data["M_2"],dtype=np.float64)
a_2    = np.array(data["a_2"],dtype=np.float64)
g_1    = np.ravel(np.array(data["g_1"],dtype=np.float64)) 
g_2    = np.ravel(np.array(data["g_2"],dtype=np.float64))
Ntotal = np.array(data["Ntotal"],dtype=np.float64)
Nin    = np.array(data["Nin"],dtype=np.float64)
Q_2 = K_2-(Eps_2*a_2-Lz_2)**2 

if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  raise ValueError("Eps,Lz,Q not bound; something wrong")

# Calculate the density using weighted f
N_1 = np.ones(len(Eps_2)) / float(Ntotal)
f_1 = 1. / g_1
f_2 = f_1
N_2 = g_2 * f_2 # Convert back to numbers
print np.sum(np.isnan(N_2))
print "Total nans: " + str(np.sum(np.isnan(N_2)))
print "Warning: Letting N_2 -> nan go to zero; need to fix this later"
N_2[np.isnan(N_2)] = 0
# Save to gyoto:
# Note: In the Gyoto simulation we need at least:
# - Number of particles N(Eps,L,Q) to add weight to each particle
# - Constants of motions
# - Spin and mass
# - grid
# We can easily provide here the constants of motions and number of particles

# We need to convert everything into M=1 units:
Eps_2 = Eps_2
Lz_2  = Lz_2 / M_2
K_2   = K_2  / (M_2)**2
Q_2   = Q_2  / (M_2)**2
a_2   = a_2 / M_2
M_2   = M_2 / M_2
if a_2 > 1:
  raise ValueError("Bad a_2")
# Make sure they are still bound
if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  raise ValueError("Eps,Lz,Q not bound after unit conversion; something wrong")

M_2 = np.ones(len(Eps_2))*M_2
a_2 = np.ones(len(Eps_2))*a_2
Ntotal=np.ones(len(Eps_2))*Ntotal
Nin   =np.ones(len(Eps_2))*Nin
print M_2
np.savetxt(fout, np.transpose(np.array([Eps_2,Lz_2,Q_2,N_2,M_2,a_2,Ntotal,Nin])))

