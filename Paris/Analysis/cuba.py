import numpy as np
from cubature import cubature

def calculateAnnihilationPython(rho=1,P=np.array([0.,0.,0.8]),sigma=np.array([0.1,0.5,0.00001])):
  def gaussian(x):
      return 1./(np.sqrt(2*np.pi)**3*sigma[0]*sigma[1]*sigma[2]) * np.exp(-1.*np.sum((x-P)**2/(2.*sigma**2)))
  def gaussian_v(x):
      return 1./(np.sqrt(2*np.pi)**3*sigma[0]*sigma[1]*sigma[2]) * np.exp(-1.*np.sum((x-P)**2/(2.*sigma**2),axis=1))
  def gaussian2(x):
      p1 = x[:3]
      p2 = x[3:6]
      ga1= np.sqrt(1.+np.sum(p1**2))
      ga2= np.sqrt(1.+np.sum(p2**2))
      v1 = p1/ga1
      v2 = p2/ga2
      dotprod=np.sum(v1*v2)
      return gaussian(x[:3])*gaussian(x[3:6])*(1.-dotprod)*np.linalg.norm(v2-v1)
  def gaussian2_v(x):
      p1 = x[:,:3]
      p2 = x[:,3:6]
      ga1= np.sqrt(1.+np.sum(p1**2,axis=1))
      ga2= np.sqrt(1.+np.sum(p2**2,axis=1))
      v1 = np.transpose(np.transpose(p1)/ga1) 
      v2 = np.transpose(np.transpose(p2)/ga2)
      dotprod=np.sum(v1*v2,axis=1)
      return gaussian_v(x[:,:3])*gaussian_v(x[:,3:6])*(1.-dotprod)*np.linalg.norm(v2-v1,axis=1)
  xmin = P-sigma*10
  xmin = np.ravel(np.append(xmin,xmin))
  xmax = P+sigma*10
  xmax = np.ravel(np.append(xmax,xmax))
  val, err = cubature(gaussian2_v, 6, 1, xmin, xmax, abserr=1.e-2,vectorized=True)
  print('Approximated: {0}'.format(val))
  return val*rho*rho, err*rho*rho


