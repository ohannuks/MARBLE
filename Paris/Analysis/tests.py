import GR
import numpy as np

def isBound(Eps,Lz,K,M,a):
  rMin, rMax  = GR.rMinMax( Eps, Lz, K, M, a )
  thMin, thMax= GR.thMinMax(Eps, Lz, K, M, a )
  # Check if roots exist
  if (np.sum(np.abs(np.imag(rMin))) != 0 or np.sum(np.abs(np.imag(rMax))) != 0):
    raise ValueError("Imaginary roots exist for Eps, Lz, K, M, a, please check")
  M = M * np.ones(len(Eps))
  a = a * np.ones(len(Eps))
  # Check if R is positive
  for i in np.arange(len(rMax)):
    rmin = rMin[i]
    rmax = rMax[i]
    thmin= thMin[i]
    thmax= thMax[i]
    if np.isnan(rmin) == True or np.isnan(thmin) == True:
      Q= K[i] - (Eps[i]*a[i]-Lz[i])**2
      raise ValueError("Bad rmin or thmin with values: %lf %lf %lf %lf %lf"%(rmin,thmin, Q, GR.Theta(Eps[i],Lz[i],K[i],M[i],a[i],np.pi/2.), GR.Theta(Eps[i],Lz[i],K[i],M[i],a[i],np.pi/10.)))
    dr   = rmax - rmin
    dth  = thmax- thmin
    r0   = rmin + dr * np.random.rand(100)
    th0  = thmin+ dth* np.random.rand(100)
    _R   = GR.R(Eps[i],Lz[i],K[i],M[i],a[i],r0)
    _Theta=GR.Theta(Eps[i],Lz[i],K[i],M[i],a[i],th0)
    if np.sum(_R>=0) != len(_R) or np.sum(np.isnan(_R)) != 0:
      r0 = r0[np.argmin(_R)]
      raise ValueError("Particle E, Lz, K not bound with %f %f %f at r %f with rmin rmax %f %f and R value %f"%(Eps[i], Lz[i], K[i], r0, rmin, rmax, _R[np.argmin(_R)]))
    if np.sum(_Theta>=0) != len(_Theta) or np.sum(np.isnan(_Theta)) != 0:
      th0 = th0[np.argmin(_Theta)]
      raise ValueError("Particle E, Lz, K not bound with %f %f %f at th %f with thmin thmax %f %f and Theta value %f"%(Eps[i], Lz[i], K[i], th0, thmin, thmax, _Theta[np.argmin(_Theta)]))
  # Everything seems ok
  return True


