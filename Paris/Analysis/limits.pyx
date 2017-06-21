cimport cython
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False) # turn off bounds-checking for entire function
def rLim2( np.ndarray[DTYPE_t] Eps, np.ndarray[DTYPE_t] Lz, np.ndarray[DTYPE_t] K, double M, double a ):
  assert (Eps.dtype == DTYPE and Lz.dtype == DTYPE) and K.dtype == DTYPE
  if (len(Eps) != len(Lz)) or (len(Eps) != len(K)):
    raise ValueError("Bad length")
  cdef np.ndarray[DTYPE_t,ndim=2] coefficients = np.transpose([Eps**2-1,
                  2*M*np.ones(len(Eps)),
                  -K-a*(a+2.*Lz*Eps-2.*a*Eps**2),
                  2.*M*K,
                  a**2*(-1.*K+(Lz-a*Eps)**2)])
  cdef np.ndarray[DTYPE_t,ndim=2] roots = np.zeros((len(Eps),4), dtype=DTYPE)
  cdef int N = len(Eps)
  cdef int i
  for i in range(N):
    res = np.roots(coefficients[i])
    if len(res) == 4:
      roots[i] = np.sort(res)
    #roots[i] = np.sort(np.roots(coefficients[i]))
  return roots



