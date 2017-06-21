import sys
import matplotlib
matplotlib.use("Agg")
import paris as pp
from paris import GR
from multiprocessing import Pool
import numpy as np
print "Usage: python $< Data/degeneracydistribution.dat Data/degeneracydistribution"
fname = sys.argv[1]
fout  = sys.argv[2]

print "Reading data.."

# Get the first 3 lines of the data
attr = pp.loadattributes(fname)
Ntotal            = attr[0]
Nin               = attr[1]
a                 = attr[2]
M                 = attr[3]
a_2               = attr[4]
M_2               = attr[5]

print "Ntotal: " + str(Ntotal)

# Load data
data = pp.loaddata(fname)
Eps_1  = data[:,0]
Lz_1   = data[:,1]
Q_1    = data[:,2]
Eps_2  = data[:,3]
Lz_2   = data[:,4]
Q_2    = data[:,5]
# Note that M_1 and a_1 are different depending on the sampler
M_1    = data[:,6]
a_1    = data[:,7]


# Cilculate K constant of motion
K_1 = Q_1+(Eps_1*a_1-Lz_1)**2
K_2 = Q_2+(Eps_2*a_2-Lz_2)**2

print "Checking bounds.."

GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1)
GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2)

print "Calculating limits.."

thMin_1, thMax_1 = GR.thLim(Eps_1,Lz_1,K_1,M_1,a_1)
rLimits = GR.rLim2(Eps_1,Lz_1,K_1,M_1,a_1)
rMin_1  = rLimits[:,2]
rMax_1  = rLimits[:,3]

thMin_2, thMax_2 = GR.thLim(Eps_2,Lz_2,K_2,M_2,a_2)
rLimits = GR.rLim(Eps_2,Lz_2,K_2,M_2,a_2)
rMin_2  = rLimits[:,2]
rMax_2  = rLimits[:,3]

print "Calculating degeneracy.."

def deg_pool(args):
  Eps_1,Lz_1,K_1,M_1,a_1,rMin_1,rMax_1,thMin_1,thMax_1 = args
  return GR.degeneracy(np.array([Eps_1]),np.array([Lz_1]),np.array([K_1]),np.array([M_1]),np.array([a_1]),np.array([rMin_1]),np.array([rMax_1]),np.array([thMin_1]),np.array([thMax_1]))

P1 = Pool(processes=4)

g_1 = np.array(P1.map(deg_pool, zip(Eps_1,Lz_1,K_1,M_1,a_1,rMin_1,rMax_1,thMin_1,thMax_1), chunksize=64))#GR.degeneracy(Eps_1,Lz_1,K_1,M_1,a_1,rMin_1,rMax_1,thMin_1,thMax_1)
P2 = Pool(processes=4)
g_2 = np.array(P2.map(deg_pool, zip(Eps_2,Lz_2,K_2,np.ones(len(Eps_2))*M_2,np.ones(len(Eps_2))*a_2,rMin_2,rMax_2,thMin_2,thMax_2), chunksize=64))#GR.degeneracy(Eps_2,Lz_2,K_2,M_2,a_2,rMin_2,rMax_2,thMin_2,thMax_2)
P1.close()
P1.join()
P2.close()
P2.join()

print "Saving data.."

data = {"Ntotal": Ntotal, "Nin": Nin, "a": a, "M": M, "a_2": a_2, "M_2": M_2,
        "Eps_1": Eps_1, "Lz_1": Lz_1, "K_1": K_1, "Q_1": Q_1, "rMin_1": rMin_1, "rMax_1": rMax_1, "thMin_1": thMin_1, "thMax_1": thMax_1, "g_1": g_1, "M_1": M_1, "a_1": a_1, 
        "Eps_2": Eps_2, "Lz_2": Lz_2, "K_2": K_2, "Q_2": Q_2, "rMin_2": rMin_2, "rMax_2": rMax_2, "thMin_2": thMin_2, "thMax_2": thMax_2, "g_2": g_2, "M_2": M_2, "a_2": a_2}

pp.savepdata(fout,data)
data2=pp.loadpdata(fout)

