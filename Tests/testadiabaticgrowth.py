# REWRITE THIS WHOLE THING; brent algorithm simply DOES NOT WORK. the integral Sr seems to be functioning well, but otherwise this whole thing is broken
import numpy as np
import sys
import paris as pp
import pylab as plt
from paris import GR
from cubature import cubature
from scipy.optimize import brentq

print "Usage: python $< Data/degeneracydistribution.dat Data/degeneracydistribution"
fname = sys.argv[1]
fout  = sys.argv[2]

# Get the first 3 lines of the data
attr = pp.loadattributes(fname)
Ntotal            = attr[0]
Nin               = attr[1]
a                 = attr[2]
M                 = attr[3]
a_2               = attr[4]
M_2               = attr[5]

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

if ((a == 0) == False) or ((a_2 == 0) == False):
  raise ValueError("This test only works when a = a_2 = 0")



# Cilculate K constant of motion
K_1 = Q_1+(Eps_1*a_1-Lz_1)**2
K_2 = Q_2+(Eps_2*a_2-Lz_2)**2

GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1)

thMin_1, thMax_1 = GR.thLim(Eps_1,Lz_1,K_1,M_1,a_1)
rLimits = GR.rLim2(Eps_1,Lz_1,K_1,M_1,a_1)
rMin_1  = rLimits[:,2]
rMax_1  = rLimits[:,3]

failedGrowthCondition = np.isnan(Eps_2)

Srd = lambda r,Eps,Lz,K,M,a: np.sqrt(Eps**2*r**4 - (-2.*M*r + r**2)*(K + r**2))/(-2*M*r + r**2)
def integrand_v(x,Eps,Lz,K,M,a):
  r = x[:,0]
  result = Srd(r,Eps,Lz,K,M,a)
  result[r==1.] = 0
  return result
def Sr( Eps,Lz,K,M,a ):
  rMin, rMax = GR.rMinMax(Eps,Lz,K,M,a)
  xmin = np.array([rMin[0]],np.float64)
  xmax = np.array([rMax[0]],np.float64)
  return cubature(integrand_v,1,1,xmin,xmax,vectorized=True,args=(Eps,Lz,K,M,a),relerr=1.e-3)[0]

# Hand-pick only failed growths
Eps_2 = Eps_2[failedGrowthCondition]
Eps_1 = Eps_1[failedGrowthCondition]
Lz_1  = Lz_1[failedGrowthCondition]
K_1   = K_1[failedGrowthCondition]

# Make sure that the growths really failed
def growthFailed(Eps_1, Lz_1, K_1, M_1, a_1, M_2, a_2):
  # Calculate Sr for many different possible Eps_2 values:
  EpsTrial = np.linspace(0.98,1,3000)
  SrValues = np.array([Sr(EpsTrial[i],Lz_1,K_1,M_2,a_2) for i in range(len(EpsTrial))])
  return EpsTrial, SrValues

falseAlarms = 0

for i in range(len(Eps_1)):
  Sr_1 = Sr(Eps_1[i], Lz_1[i], K_1[i], M_1[i], a_1[i])
  EpsTrial, SrValues = growthFailed(Eps_1[i], Lz_1[i], K_1[i], M_1[i],a_1[i],M_2,a_2)
  SrValues = np.ravel(SrValues)
  EpsTrial = np.array(EpsTrial)[np.isnan(SrValues)==False]
  SrValues = np.array(SrValues)[np.isnan(SrValues)==False]
  if len(SrValues) < 2:
    continue
  if (SrValues[0] - Sr_1) * (SrValues[len(SrValues)-1] - Sr_1) < 1:
    falseAlarms = falseAlarms + 1

np.savetxt(fout,np.array([falseAlarms, len(Eps_1), len(Q_2)]))




