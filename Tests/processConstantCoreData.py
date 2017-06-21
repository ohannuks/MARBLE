import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR


print "python processConstantCoreData.py grid.dat griddegeneracy.dat"
fname = sys.argv[1]
fout = sys.argv[2]
data =pp.loadpdata(fname)
Eps_1  = data["Eps_1"]
Lz_1   = data["Lz_1"]
K_1    = data["K_1"]
M_1    = data["M_1"]
a_1    = data["a_1"]
rMin_1 = data["rMin_1"]
rMax_1 = data["rMax_1"]
g_1    = data["g_1"]

Eps_2  = data["Eps_2"]
Lz_2   = data["Lz_2"]
K_2    = data["K_2"]
M_2    = data["M_2"]
a_2    = data["a_2"]
rMin_2 = data["rMin_2"]
rMax_2 = data["rMax_2"]
g_2    = data["g_2"]


if GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1) == False:
  raise ValueError("Eps_1,Lz_1,Q_1 not bound; something wrong")
if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  raise ValueError("Eps_2,Lz_2,Q_2 not bound; something wrong")

#if data["a"] != 0:
#  raise ValueError("processConstantCore built assuming a = 0; this is test-specific")
if data["M"] != 1:
  raise ValueError("processConstantCore built assuming M=1; this is test-specific")

data_save = np.transpose([data["Eps_1"], data["Lz_1"], data["Q_1"], data["rMin_1"], data["rMax_1"], data["thMin_1"], data["thMax_1"], np.ravel(data["g_1"])])
data_save = np.append(np.array([data["Ntotal"], data["a_2"], data["M_2"], 0, 0, 0, 0, 0]), data_save)
np.savetxt(fout, data_save, fmt="%.16e")

