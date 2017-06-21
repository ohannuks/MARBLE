import numpy as np
import paris as pp
import pylab as plt
from paris import GR
from cubature import cubature
from line_profiler import LineProfiler

@profile
def integrand_v(x, M, a, r, th):
  f0 = 1.
  Eps = x[:,0]
  Lz  = x[:,1]
  K   = x[:,2]
  result = np.zeros(len(Eps))*1.
  rMin, rMax  = GR.rMinMax(Eps, Lz, K, M, a)
  thMin, thMax= GR.thMinMax(Eps,Lz, K, M, a)
  # Filter:
  cond=(np.isnan(rMin)==False)&(np.imag(rMin)==0)&(np.isnan(thMin)==False)&(r>rMin)&(r<rMax)&(th<thMax)&(th>thMin)
  #print cond
  Eps = Eps[cond]
  Lz  = Lz[cond]
  K   = K[cond]
  ut = GR.u(Eps,Lz,K,M,a,r,th)[0]
  # Get determinant
  detg= GR.detg(M,a,r,th)
  # Jacobian
  J   = np.abs(GR.jacobian(Eps,Lz,K,M,a,r,th))
  # Calculate density
  result[cond] = np.sqrt(-1.*detg) * f0 * (4. / J)
  if np.sum(np.isnan(result)) != 0:
    print Eps, Lz, K, M, a, J, result
    raise ValueError("Bad result")
  return result

M = 1.
a = 0.
r = 10.
th = np.pi/2.
xmin = [0.99,-7, 0]
xmax = [0.999,  7, 7*7]
val,er = cubature(integrand_v, 3, 1, xmin, xmax, vectorized=True,args=(M,a,r,th),relerr=1.e-1, maxEval=100000000)

print val
