# Everything related to kerr metric etc here
#import pyximport; pyximport.install()

import numpy as np
#import pylab as plt
import tests
Csc  = lambda th: 1./np.sin(th)
from cubature import cubature
#from limits import rLim2

def getDelta(M,a,r):
  return r**2-2.*M*r+a**2

def getSigma(M,a,r,th):
  return r**2+a**2*np.cos(th)**2

# Eq. 8.5.10 from Frolov
def Theta(Eps, Lz, K, M, a, th):
  return K - a**2 * np.cos(th)**2 - (Eps*a*np.sin(th)-Lz/np.sin(th))**2

def R(Eps, Lz, K, M, a, r):
  Delta = getDelta(M,a,r)
  return (Eps*(r**2+a**2)-a*Lz)**2-(r**2+K)*Delta

# Computes the u_down
def u_d(Eps, Lz, K, M, a, r, th, dir_r=1, dir_th=1):
  # From Frolov
  Delta = getDelta(M,a,r)
  Sigma = getSigma(M,a,r,th)
  p_t  = -1.*Eps
  p_phi= Lz
  p_r  = np.sqrt(R(Eps,Lz,K,M,a,r))/Delta
  p_th = np.sqrt(Theta(Eps,Lz,K,M,a,th))
  # Make the direction:
  p_r  = dir_r * p_r
  p_th = dir_th* p_th
  return p_t, p_r, p_th, p_phi

# Computes the u_up:
def u(Eps,Lz,K,M, a,r,th, dir_r=1, dir_th=1):
  # Computes u via Mathematica
  Delta = getDelta(M,a,r)
  Sigma = getSigma(M,a,r,th)
  p_r  = np.sqrt(R(Eps,Lz,K,M,a,r))/Delta
  p_th = np.sqrt(Theta(Eps,Lz,K,M,a,th))
  pt   = (-4.*a*Lz*M*r)/((a**2 + r*(-2.*M + r))*(a**2 + 2.*r**2 + a**2*np.cos(2.*th))) - (4.*Eps*(r**2 + a**2.*np.cos(th)**2)*(-(a**2 + r**2)**2 + a**2*(a**2 + r*(-2.*M + r))*np.sin(th)**2))/((a**2 + r*(-2.*M + r))*(a**2 + 2.*r**2 + a**2*np.cos(2.*th))**2)
  grr  = Sigma/Delta
  pr   = p_r/grr
  gthth= Sigma
  pth  = p_th/gthth
  pphi = (4.*a*Eps*M*r)/((a**2 + r*(-2.*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2.*th))) + (Lz*(a**2 + 2.*r*(-2.*M + r) + a**2*np.cos(2.*th))*Csc(th)**2)/((a**2 + r*(-2.*M + r))*(a**2 + 2.*r**2 + a**2*np.cos(2.*th)))
  return pt, pr, pth, pphi 




def rLim2( Eps, Lz, K, M, a ):
  if len(np.shape(Eps)) == 0:
    Eps = np.array([Eps])
    Lz  = np.array([Lz])
    K   = np.array([K])
  if (len(Eps) != len(Lz)) or (len(Eps) != len(K)):
    raise ValueError("Bad length")
  coefficients = np.transpose([Eps**2-1,
                  2*M*np.ones(len(Eps)),
                  -K-a*(a+2.*Lz*Eps-2.*a*Eps**2),
                  2.*M*K,
                  a**2*(-1.*K+(Lz-a*Eps)**2)])
  roots = []
  for i in np.arange(len(Eps)):
    coefficient = coefficients[i]
    result = np.sort(np.roots(coefficient))
    if len(result) != 4:
      roots.append(np.ones(4)*1j)
    else:
      roots.append(result)
  return np.array(roots)

def rLim(Eps, Lz, K, M, a):
  Eps = np.array(Eps,dtype=np.complex128)
  Lz  = np.array(Lz, dtype=np.complex128)
  K   = np.array(K,  dtype=np.complex128)
  a2  = a**2
  a3  = a2*a
  a4  = a2*a2
  Eps2= Eps**2
  Eps3= Eps2*Eps
  Eps4= Eps2*Eps2
  K2= K**2
  K3= K2*K
  K4= K2*K2
  Lz2= Lz**2
  Lz3= Lz2*Lz
  Lz4= Lz2*Lz2
  M2  = M**2
  r3 = -M/(2.*(-1 + Eps2)) + np.sqrt((-2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + M2/(-1 + Eps2)**2 + (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)))/2. - np.sqrt((-4*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + (2*M2)/(-1 + Eps2)**2 - (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) - (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)) + ((-16*K*M)/(-1 + Eps2) + (8*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M)/(-1 + Eps2)**2 - (8*M**3)/(-1 + Eps2)**3)/ (4.*np.sqrt((-2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + M2/(-1 + Eps2)**2 + (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)))))/2.
  r4 = -M/(2.*(-1 + Eps2)) + np.sqrt((-2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + M2/(-1 + Eps2)**2 + (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)))/2. + np.sqrt((-4*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + (2*M2)/(-1 + Eps2)**2 - (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) - (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)) + ((-16*K*M)/(-1 + Eps2) + (8*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M)/(-1 + Eps2)**2 - (8*M**3)/(-1 + Eps2)**3)/ (4.*np.sqrt((-2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz))/(3.*(-1 + Eps2)) + M2/(-1 + Eps2)**2 + (2**0.3333333333333333*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2))/ (3.*(-1 + Eps2)*(2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333) + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2 + np.sqrt(-4*((-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**2 + 12*(-1 + Eps2)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) - 12*K*M2)**3 + (2*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)**3 - 72*(-1 + Eps2)*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2) + 108*(-1 + Eps2)*K2*M2 - 36*K*(-a2 + 2*a2*Eps2 - K - 2*a*Eps*Lz)*M2 + 108*(a4*Eps2 - a2*K - 2*a3*Eps*Lz + a2*Lz2)*M2)**2))**0.3333333333333333/ (3.*2**0.3333333333333333*(-1 + Eps2)))))/2.
  return np.transpose([r3, r3, r3, r4])

rLim = rLim2

def rMinMax( Eps, Lz, K, M, a ):
  ''' Calculates the rMin and rMax of given bound orbit with constants of motions Eps, Lz, K (see Frolov's book for more info on the constants of motion). Returns nan for the orbits which are not bound

  '''
  M = np.ones(len(Eps)) * M
  a = np.ones(len(Eps)) * a
  rHorizon = M+np.sqrt(M**2-a**2)
  rLimits = rLim( Eps, Lz, K, M, a )
  rMin = rLimits[:,2]
  rMax = rLimits[:,3]
  rMean = rMin + (rMax-rMin)/2.
  ut = np.array( u(Eps,Lz,K,M, a,r=np.real(rMean),th=np.pi/2., dir_r=1, dir_th=1)[0])
  notBound = ( ((np.imag(rMin)==0)==False) | ((np.imag(rMax)==0)==False)) | (ut<0) | (rMin<=rHorizon)
  notBound = ( ((np.abs(np.imag(rMin))<1.e-14)==False) | ((np.abs(np.imag(rMax))<1.e-14)==False)) | (ut<0) | (np.real(rMin)<=rHorizon)
  rMin[notBound] = np.nan
  rMax[notBound] = np.nan
  rMin[notBound==False] = np.real(rMin[notBound==False])
  rMax[notBound==False] = np.real(rMax[notBound==False])
  #print rMin[notBound==False]
  if np.sum( R(Eps[notBound==False], Lz[notBound==False], K[notBound==False], M[notBound==False], a[notBound==False], rMean[notBound==False]) < 0 ) != 0:
    print R(Eps[notBound==False], Lz[notBound==False], K[notBound==False], M[notBound==False], a[notBound==False], rMean[notBound==False])[R(Eps[notBound==False], Lz[notBound==False], K[notBound==False], M[notBound==False], a[notBound==False], rMean[notBound==False]) < 0]
    raise ValueError("Bad R value ; we have bound orbit with the incorrect rmin and rmax")
  return np.array(rMin,dtype=np.float64), np.array(rMax,dtype=np.float64)

#def thLim2(Eps,Lz,K,M,a):
#  ''' Calculates the th limits of given orbits
#
#  '''
#  M = M * np.ones(len(Eps))
#  a = a * np.ones(len(Eps))
#  nerr=1.e-7
#  sinth =         np.sqrt(-(a**2/(-a**2 + a**2*Eps**2)) + K/(-a**2 + a**2*Eps**2) + (2*a*Eps*Lz)/(-a**2 + a**2*Eps**2) - np.sqrt(a**4 - 2*a**2*K + K**2 - 4*a**3*Eps*Lz + 4*a*Eps*K*Lz + 4*a**2*Lz**2)/(-a**2 + a**2*Eps**2))/np.sqrt(2)
#  th = np.arcsin(sinth)
#  thMin = np.arcsin(sinth)
#  thMax = thMin + (np.pi-thMin)*2.
#  thMin = thMin+nerr
#  thMax = thMax-nerr
#  if a.all() == 0:
#    # TODO: The theta limits are nan if a=0; this is a temporary fix but should also make sure that there's nothing fishy going on; technically the theta limits should be set by the equation Theta>0; make sure to check these
#    print "Warning: Setting a=0 th limits!"
#    thMin = np.ones(len(Eps))*nerr
#    thMax = np.ones(len(Eps))*(np.pi-nerr)
#  return thMin, thMax

from scipy.optimize import brentq
def thLim(Eps,Lz,K,M,a):
  ''' Calculates the th limits of given orbits using brent algorithm

  '''
  M = M * np.ones(len(Eps))
  a = a * np.ones(len(Eps))
  Eps=Eps*np.ones(len(Eps))
  Lz =Lz *np.ones(len(Eps))
  K  =K  *np.ones(len(Eps))
  Q  =K-(Eps*a-Lz)**2
  def _root(x,Eps,Lz,K,M,a):
    return Theta(Eps,Lz,K,M,a,x)
  #if (np.min(Theta(Eps,Lz,K,M,a,np.pi/2.)) >= 0) == False:
  #  raise ValueError("Bad Theta (not a bound orbit as currently assumed; note that we can actually fix this by deciding what to return when an orbit is not bound) value with Theta(pi/2): %f"%np.min(Theta(Eps,Lz,K,M,a,np.pi/2.)))
  thMin = []
  thMax = []
  ERR = 1.e-8
  for i in range(len(Eps)):
    ThetaMax = Theta(Eps[i],Lz[i],K[i],M[i],a[i],np.pi/2.)
    if ThetaMax < 0 or np.isnan(ThetaMax) == True:
      thMin.append(np.nan)
      thMax.append(np.nan)
      continue
    _thMin = brentq(_root,ERR,np.pi/2.,args=(Eps[i], Lz[i], K[i], M[i], a[i]),maxiter=1000,xtol=1.e-22)
    #if _thMin == np.pi/2.:
    #  print _root(ERR,Eps[i],Lz[i],K[i],M[i],a[i]), _root(np.pi/2.,Eps[i],Lz[i],K[i],M[i],a[i])
    #  raise ValueError("Something weird going on here; brent algorithm returns exactly pi./2")
    _thMax = np.pi-_thMin
    thMin.append(_thMin)
    thMax.append(_thMax)
    if np.abs(Theta(Eps[i],Lz[i],K[i],M[i],a[i],_thMin)) > 1.e-5:
      raise ValueError("Bad th limit with thmin, thmax, Theta(thMin): %f %f %f \n"%(_thMin, _thMax, Theta(Eps[i],Lz[i],K[i],M[i],a[i],_thMin)))
  thMin = np.array(thMin)
  thMax = np.array(thMax)
  thMin[K-(Lz-a*Eps)**2<0] = np.nan
  thMax[K-(Lz-a*Eps)**2<0] = np.nan
  thMin[Q==0] = np.pi/2.
  thMax[Q==0] = np.pi/2.
  return np.array(thMin), np.array(thMax)

thLim2   = thLim

thMinMax = thLim

# Down indices for kerr in BL according to Frolov
def gmunu(M,a,r,th, inverse=False):
  ''' Computes g_{\mu \nu} in the BL metric given mass, spin and position

  '''
  Sigma = getSigma(M,a,r,th)
  Delta = getDelta(M,a,r)
  A = (r**2+a**2)**2 - Delta * a**2 *np.sin(th)**2
  gtt   = -1.*(1.-2.*M*r/Sigma)
  grr   = Sigma/Delta
  gthth = Sigma
  gpp   = A*np.sin(th)**2/Sigma
  gtp   = -2.*M*r*a*np.sin(th)**2 / Sigma
  g = np.zeros((4,4))
  g[0][0] = gtt
  g[1][1] = grr
  g[2][2] = gthth
  g[3][3] = gpp
  g[0][3] = gtp
  g[3][0] = gtp
  if inverse == True:
    inverse = np.linalg.inv(g)
    g       = np.transpose(inverse)
  return g

def detg(M,a,r,th):
  return -((a**2 + 2*r**2 + a**2*np.cos(2*th))**2*np.sin(th)**2)/4.

def jacobian(Eps, Lz, K, M, a, r, th):
  ''' Computes the Jacobian |d(E,K,Lz,mu)/d(p^0,p^1,p^2,p^3)|

  '''
  m = 1.
  x1= 1.
  x2= 1. # Direction
  return (-2.*np.sqrt(-((K + m**2*r**2)*(a**2 + r*(-2.*M + r))) + (a*Lz - Eps*(a**2 + r**2))**2)*x1*x2*(r**2 + a**2*np.cos(th)**2)*np.sin(th)**2* np.sqrt(K - a**2*m**2*np.cos(th)**2 - (Lz*Csc(th) - a*Eps*np.sin(th))**2))/m

def degeneracy(Eps,Lz,K,M,a,rMin,rMax,thMin,thMax):
  ''' Calculates degeneracy for one given orbit
      TODO: Make sure jacobian is non-negative

  '''
  M = M * np.ones(len(Eps))
  a = a * np.ones(len(Eps))
  def integrand_v(x,Eps,Lz,K,M,a):
    r = x[:,0]
    th= x[:,1]
    ut = u(Eps,Lz,K,M,a,r,th)[0]
    jac = np.abs(jacobian(Eps,Lz,K,M,a,r,th))
    result = np.abs(detg(M,a,r,th))* ut * 4./jac
    result[r==1.] = 0
    result[np.isnan(result)] = 0
    result[result>1e60] = 0
    return result
  result = []
  for i in np.arange(len(Eps)):
    if thMin[i] == np.nan:
      raise ValueError("thMin nan: %d %lf"%(i,thMin[i]))
    xmin=np.array([rMin[i],thMin[i]], np.float64)
    xmax=np.array([rMax[i],thMax[i]], np.float64)
    val, err = cubature(integrand_v, 2, 1, xmin, xmax, vectorized=True,args=(Eps[i],Lz[i],K[i],M[i],a[i]), relerr=1e-2)
    result.append(val[0])
  return np.array(result)

# Calculate density assuming th = pi/2. since we dont have thMin and thMax
def FFIOdensity( Eps, Lz, K, M, a, rMin, rMax, r ):
  condition = ((rMin<r)&(r<rMax)) # Make sure the particle actually goes to r
  # Make sure everything is of the same length; something M is a vector, sometimes a scalar. This is to make up for it.
  if len(np.shape(M)) != 0:
    if np.sum(M[M==M[0]]) != len(M):
      raise ValueError("Trying to input M with multiple values to FFIODensity; not good")
    M = M[0]
    a = a[0]
  Eps = Eps[condition]
  Lz  = Lz[condition]
  K   = K[condition]
  rMin= rMin[condition]
  rMax= rMax[condition]
  # _d denotes down indices
  f0 = 1. # Unnormalized for now
  th = np.pi/2.
  u_down = u_d(Eps,Lz,K,M,a,r,th)
  u_up   = u(Eps,Lz,K,M,a,r,th)
  u_d_r  = u_down[1]
  u_th   = u_up[2]
  Delta = getDelta(M,a,r)
  Sigma = getSigma(M,a,r,th)
  # Check consistency:
  #print "Consistency check: " + str(u_th/((1*np.sqrt(K - a**2*np.cos(th)**2 - (Lz*Csc(th) - a*Eps*np.sin(th))**2))/(r**2 + a**2*np.cos(th)**2)))
  # Straight from Sadeghian's article
  # Note: Sigma is defined as Sigma**2 in Sadeghian's article; this means that this integrand is not incorrect even though in Sadeghian's article they use Sigma**2
  J_d_0  = -2. * Eps * f0 / (Sigma * Delta * np.abs(u_d_r)*np.abs(u_th)*np.sin(th))
  J_d_phi=  2. * Lz  * f0 / (Sigma * Delta * np.abs(u_d_r)*np.abs(u_th)*np.sin(th))
  # TODO:"Warning: Cheating here; I do not have proper proof why this works: "
  J_d_0  = J_d_0  * np.sqrt(K)
  J_d_phi= J_d_phi* np.sqrt(K)
  # Integrate via monte carlo sampling
  J_d_0  = np.sum(J_d_0)
  J_d_phi= np.sum(J_d_phi)
  if np.isnan(J_d_0) == True:
    raise ValueError("J_d_0 is nan!")
  if np.isnan(J_d_phi) == True and a != 0:
    raise ValueError("J_d_phi is nan!")
  Omega = J_d_phi / J_d_0
  if a == 0:
    Omega = 0
  g = gmunu(M,a,r,th)
  gpp = g[3][3]
  gtp = g[0][3]
  gtt = g[0][0]
  rho = -1.*J_d_0 * np.sqrt( (gpp+2.*gtp*Omega+gtt*Omega**2)/Delta )
  return rho


# Calculate density assuming th = pi/2. since we dont have thMin and thMax
def FFIOdensity2( f0, Eps, Lz, K, M, a, rMin, rMax, r ):
  condition = ((rMin<r)&(r<rMax)) # Make sure the particle actually goes to r
  if len(np.shape(M)) != 0:
    if np.sum(M[M==M[0]]) != len(M):
      raise ValueError("Trying to input M with multiple values to FFIODensity; not good")
    M = M[0]
    a = a[0]
  Eps = Eps[condition]
  Lz  = Lz[condition]
  K   = K[condition]
  rMin= rMin[condition]
  rMax= rMax[condition]
  f0  = f0[condition]
  # _d denotes down indices
  th = np.pi/2.
  u_down = u_d(Eps,Lz,K,M,a,r,th)
  u_up   = u(Eps,Lz,K,M,a,r,th)
  u_d_r  = u_down[1]
  u_th   = u_up[2]
  Delta = getDelta(M,a,r)
  Sigma = getSigma(M,a,r,th)
  # Check consistency:
  #print "Consistency check: " + str(u_th/((1*np.sqrt(K - a**2*np.cos(th)**2 - (Lz*Csc(th) - a*Eps*np.sin(th))**2))/(r**2 + a**2*np.cos(th)**2)))
  # Straight from Sadeghian's article
  # Note: Sigma is defined as Sigma**2 in Sadeghian's article; this means that this integrand is not incorrect even though in Sadeghian's article they use Sigma**2
  J_d_0  = -2. * Eps * f0 / (Sigma * Delta * np.abs(u_d_r)*np.abs(u_th)*np.sin(th))
  J_d_phi=  2. * Lz  * f0 / (Sigma * Delta * np.abs(u_d_r)*np.abs(u_th)*np.sin(th))
  #TODO: "Warning: Cheating here; I do not have proper proof why this works: "
  J_d_0  = J_d_0  * np.sqrt(K)
  J_d_phi= J_d_phi* np.sqrt(K)
  J_d_0[np.abs(J_d_0)>1e30] = 0
  # Integrate via monte carlo sampling
  J_d_0  = np.sum(J_d_0)
  J_d_phi= np.sum(J_d_phi)
  if a == 0:
    J_d_phi = 0.
  if J_d_phi > 1e-9 and a == 0:
    print "J_d_Phi should be zero, but is not: %f and 0: %f"%(J_d_phi, J_d_0)
  if np.isnan(J_d_0) == True:
    raise ValueError("J_d_0 is nan!")
  if np.isnan(J_d_phi) == True:
    raise ValueError("J_d_phi is nan!")
  Omega = J_d_phi / J_d_0
  g = gmunu(M,a,r,th)
  gpp = g[3][3]
  gtp = g[0][3]
  gtt = g[0][0]
  rho = -1.*J_d_0 * np.sqrt( (gpp+2.*gtp*Omega+gtt*Omega**2)/Delta )
  return rho

# Calculate density assuming th = pi/2. since we dont have thMin and thMax
def FFIOdensity3( f0, Eps, Lz, K, M, a, rMin, rMax, r ):
  condition = ((rMin<r)&(r<rMax)) # Make sure the particle actually goes to r
  if len(np.shape(M)) != 0:
    if np.sum(M[M==M[0]]) != len(M):
      raise ValueError("Trying to input M with multiple values to FFIODensity; not good")
    M = M[0]
    a = a[0]
  Eps = Eps[condition]
  Lz  = Lz[condition]
  K   = K[condition]
  rMin= rMin[condition]
  rMax= rMax[condition]
  f0  = 1.
  # _d denotes down indices
  th = np.pi/2.
  u_down = u_d(Eps,Lz,K,M,a,r,th)
  u_up   = u(Eps,Lz,K,M,a,r,th)
  u_d_r  = u_down[1]
  u_th   = u_up[2]
  Delta = getDelta(M,a,r)
  Sigma = getSigma(M,a,r,th)
  # Check consistency:
  #print "Consistency check: " + str(u_th/((1*np.sqrt(K - a**2*np.cos(th)**2 - (Lz*Csc(th) - a*Eps*np.sin(th))**2))/(r**2 + a**2*np.cos(th)**2)))
  # Straight from Sadeghian's article
  # Note: Sigma is defined as Sigma**2 in Sadeghian's article; this means that this integrand is not incorrect even though in Sadeghian's article they use Sigma**2
  J_d_0  = -2. * Eps * f0 / (Sigma * Delta * np.abs(u_d_r)*np.abs(u_th)*np.sin(th))
  #print J_d_0
  # TODO: print "Warning: Cheating here; I do not have proper proof why this works: "
  J_d_0  = J_d_0  * np.sqrt(K)
  J_d_0[np.abs(J_d_0)>1e30] = 0
  # Integrate via monte carlo sampling
  J_d_0  = np.sum(J_d_0)
  Omega = 0.
  g = gmunu(M,a,r,th)
  gpp = g[3][3]
  gtp = g[0][3]
  gtt = g[0][0]
  rho = -1.*J_d_0 * np.sqrt( (gpp+2.*gtp*Omega+gtt*Omega**2)/Delta )
  return rho

def filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th):
  if len(np.shape(r)) != 0 or len(np.shape(th)) != 0:
    raise ValueError("Expecting scalar r and th")
  condition = ((rMin<r)&(r<rMax)&(thMin<th)&(th<thMax))&(np.isnan(rMin)==False)&(np.isnan(thMin)==False)
  if np.sum(condition) < 2:
    success = False
    return success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0
  Eps    = Eps[condition]
  Lz     = Lz[condition]
  K      = K[condition]
  rMin   = rMin[condition]
  rMax   = rMax[condition]
  thMin   = thMin[condition]
  thMax   = thMax[condition]
  f0     = f0[condition]
  success= (len(Eps)>0)
  return success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0


# Calculate mass-current density assuming th = pi/2. since we dont have thMin and thMax
def Nmu( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the 4-current J^\mu = \int f * (p^\mu / m) * \sqrt{-g} d^4 p

  '''
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return np.zeros(4)
  Delta  = getDelta(M,a,r)
  Sigma  = getSigma(M,a,r,th)
  pmu    = u(Eps,Lz,K,M,a,r,th)
  Detg   = detg(M,a,r,th)
  J      = np.abs(jacobian(Eps,Lz,K,M,a,r,th))
  if np.sum(pmu[0] >= 0) != len(pmu[0]):
    print pmu[0], M, a
    print rMin, r
    print rMax, r
    raise ValueError("Bad 4-velocity component (see above print)")
  integrand = f0 * pmu * np.sqrt(-1.*Detg) * (4./J)
  nmu    = np.sum(integrand,axis=1)
  # Check for invalid values
  if np.sum(np.isnan(nmu)) != 0:
    cond=np.isnan(pmu[2])
    print Eps[cond], Lz[cond], K[cond], r, th, pmu[2][cond]
    print Theta(Eps,Lz,K,M,a,th)[cond], getSigma(M,a,r,th)
    raise ValueError("Nmu is nan; above is the integrand where Nmu = np.sum(integrand)")
  return nmu

# Calculate the energy-momentum tensor
def Tmunu( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the energy momentum tensor T^{\mu \nu} = \int f * (p^\mu / m) (p^\nu / m) * \sqrt{-g} d^4 p = \int f p^\mu p^\nu d\mathcal{V}_p

      See e.g. Andreasson, Hakan. "The Einstein-Vlasov system/kinetic theory." Living Reviews in Relativity 14.1 (2011): 4.
  '''
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return np.zeros((4,4))
  Delta  = getDelta(M,a,r)
  Sigma  = getSigma(M,a,r,th)
  pmu    = u(Eps,Lz,K,M,a,r,th)
  pnu    = np.copy(pmu)
  Detg   = detg(M,a,r,th)
  J      = np.abs(jacobian(Eps,Lz,K,M,a,r,th))
  integrand = f0 * np.outer(pmu,pnu) * np.sqrt(-1.*Detg) * (4./J)
  tmunu  = np.sum(integrand)
  # Check for invalid values
  if np.sum(np.isnan(nmu)) != 0:
    print integrand
    raise ValueError("tmunu is nan; above is the integrand where Tmunu = np.sum(integrand)")
  return tmunu

def ZAMOTetrad( M, a, r, th, inverse=False ):
  ''' Returns the ZAMO tetrad at a given point e_m^{\ \ \mu} (or e^m_{\ \ \mu} if inverse=True)
  '''
  if len(np.shape(M)) != 0 or len(np.shape(r)) != 0:
    raise ValueError("ZAMOTetrad assumes all variables to be scalar")
  r2 = r*r;
  sinth = np.sin(th);
  costh = np.cos(th);
  sinth2=sinth*sinth;
  costh2=costh*costh;
  a2 = a*a;
  router = 1.*M + np.sqrt(M**2-a**2*costh2)
  Sigma = getSigma(M,a,r,th);
  Delta = getDelta(M,a,r);
  if Delta < 0:
    # No observer here
    return np.ones((4,4))*np.nan
  A = pow(r2+a2,2)-Delta*a2*sinth2;
  gtt=-1.+(2.*M*r)/(a2*costh2+r2);
  gpp=(sinth2*(pow(a2+r2,2)-a2*(a2+r*(-2.*M+r))*sinth2))/(a2*costh2+r2);
  gtp=-1.*((2.*a*M*r*sinth2)/(a2*costh2+r2));
  Omega=-gtp/gpp;
  Omega2=Omega*Omega;
  tetrad = np.zeros((4,4))
  tetrad[0][0] = 1./np.sqrt(abs(gtt-Omega2*gpp));
  tetrad[0][3] = Omega/np.sqrt(abs(gtt-Omega2*gpp));
  if __debug__:
    if ( Delta<0 or np.isnan(Delta) ):
      raise ValueError("Bad Delta %f at M, a, r, th %lf %lf %lf %lf "%(Delta,M,a,r/M,th))
    if ( Sigma<0 or np.isnan(Delta) ):
      raise ValueError("Bad Sigma %f at M, a, r, th %lf %lf %lf %lf "%(Sigma,M,a,r/M,th))
  tetrad[1][1] = np.sqrt(Delta/Sigma);
  tetrad[2][2] = 1./np.sqrt(Sigma);
  tetrad[3][3] = 1./np.sqrt(gpp);
  if inverse == True:
    inv = np.linalg.inv(tetrad)
    tetrad  = np.transpose(inv)
  # Make a check for the ZAMO Tetrad:
  if __debug__:
    g = gmunu(M,a,r,th, inverse)
    nmn = np.dot(tetrad,np.dot(g,np.transpose(tetrad)))
    if np.abs(np.sum(np.abs(nmn)) - 4) > 1.e-8:
      print g
      print tetrad
      raise ValueError("Bad tetrad value")
  return tetrad


def ZAMODispersion( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the ZAMO-frame dispersion defined as \int f (e^m_\mu (J^\mu/N_0 - u^\mu))^2 d\mathcal{V}_p , where $N_0$ is defined as the total phase space integral

      NOT TESTED
  '''
  if len(np.shape(r)) != 0:
    raise ValueError("Bad r input to ZAMODispersion; expecting a scalar!")
  if len(np.shape(th)) != 0:
    raise ValueError("Bad th input to ZAMODispersion; expecting a scalar!")
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return np.zeros(4)
  Delta  = getDelta(M,a,r)
  Sigma  = getSigma(M,a,r,th)
  pmu    = u(Eps,Lz,K,M,a,r,th)
  Detg   = detg(M,a,r,th)
  J      = np.abs(jacobian(Eps,Lz,K,M,a,r,th))
  if np.sum(pmu[0] >= 0) != len(pmu[0]):
    raise ValueError("Bad 4-velocity component")
  N0     = np.sum(f0 * np.sqrt(-1.*Detg) * (4./J))
  nmu    = np.sum(f0 * pmu * np.sqrt(-1.*Detg) * (4./J),axis=1)
  # Tetrad
  tetrad = ZAMOTetrad(M,a,r,th,inverse=True)
  dispersion= np.sum((np.dot(tetrad, np.transpose(nmu/N0 - np.transpose(pmu))))**2 * f0 * np.sqrt(-1.*Detg) * (4./J), axis=1) / N0
  # Check for invalid values
  if np.sum(np.isnan(nmu)) != 0:
    print integrand
    raise ValueError("Nmu is nan; above is the integrand where Nmu = np.sum(integrand)")
  for i in range(4):
    if dispersion[i] < 0:
      raise ValueError("Bad dispersion value")
  return dispersion

def Ecom(Eps0, Lz0, K0, Eps1, Lz1, K1, M, a, r, th):
  ''' Calculates the COM energy and returns the outer product of (0,1)
  '''
  print "Processing term 1.."
  term1 = -1.*((4*a*Eps0*M*r + Lz0*(a**2 + 2*r*(-2*M + r) + a**2*np.cos(2*th))*Csc(th)**2)*Lz1)/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*th)))
  print "Processing term 2.."
  term2 = ((-4*a*Lz0*M*r + Eps0*(a**4 + 2*r**4 + a**2*r*(2*M + 3*r) + a**2*(a**2 + r*(-2*M + r))*np.cos(2*th)))*Eps1)/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*th)))
  print "Processing term 3.."
  term3 = -1.*((np.sqrt((K0 + r**2)*(-a**2 - r*(-2*M + r)) + (a*Lz0 - Eps0*(a**2 + r**2))**2)*np.sqrt((K1 + r**2)*(-a**2 - r*(-2*M + r)) + (a*Lz1 - Eps1*(a**2 + r**2))**2))/  ((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(th)**2)))
  print "Processing term 4.."
  term4 = (np.sqrt(K0 - a**2*np.cos(th)**2 - (Lz0*Csc(th) - a*Eps0*np.sin(th))**2)*np.sqrt(K1 - a**2*np.cos(th)**2 - (Lz1*Csc(th) - a*Eps1*np.sin(th))**2))/(r**2 + a**2*np.cos(th)**2)
  return 2.*(1. + term1 + term2 + term3 + term4 )

def MeanEcom( f0, Eps0, Lz0, K0, rMin0, rMax0, thMin0, thMax0,f1, Eps1, Lz1, K1, rMin1, rMax1, thMin1, thMax1,M, a,  r, th, maxEcom = False ):
  ''' Calculates the mean center of mass collision energy
  '''
  if len(np.shape(r)) != 0:
    raise ValueError("Bad r input to Ecom; expecting a scalar!")
  if len(np.shape(th)) != 0:
    raise ValueError("Bad th input to Ecom; expecting a scalar!")
  condition0 = ((rMin0<r)&(r<rMax0)&(thMin0<th)&(th<thMax0)) # Make sure the particle actually goes to r
  Eps0    = Eps0[condition0]
  Lz0     = Lz0[condition0]
  K0      = K0[condition0]
  rMin0   = rMin0[condition0]
  rMax0   = rMax0[condition0]
  f0      = f0[condition0]
  pmu0    = u(Eps0,Lz0,K0,M,a,r,th)
  J0      = np.abs(jacobian(Eps0,Lz0,K0,M,a,r,th))
  condition1 = ((rMin1<r)&(r<rMax1)&(thMin1<th)&(th<thMax1)) # Make sure the particle actually goes to r
  Eps1    = Eps1[condition1]
  Lz1     = Lz1[condition1]
  K1      = K1[condition1]
  rMin1   = rMin1[condition1]
  rMax1   = rMax1[condition1]
  f1      = f1[condition1]
  pmu1    = u(Eps1,Lz1,K1,M,a,r,th)
  J1      = np.abs(jacobian(Eps1,Lz1,K1,M,a,r,th))
  Delta   = getDelta(M,a,r)
  Sigma   = getSigma(M,a,r,th)
  Detg    = detg(M,a,r,th)
  if len(f0) < 3:
    return 0.
  if len(Eps0) == 0:
    return 0.
  if np.sum(pmu0[0] >= 0) != len(pmu0[0]):
    raise ValueError("Bad 4-velocity 0 component")
  if np.sum(pmu1[0] >= 0) != len(pmu1[0]):
    raise ValueError("Bad 4-velocity 1 component")
  # Monte carlo picking
  Mnorm = 10000000.
  cond0 = np.random.randint(len(f0)-1,size=int(Mnorm))
  cond1 = np.random.randint(len(f1)-1,size=int(Mnorm))
  Eps0    = Eps0[cond0]
  Lz0     = Lz0[cond0]
  K0      = K0[cond0]
  rMin0   = rMin0[cond0]
  rMax0   = rMax0[cond0]
  f0      = f0[cond0]
  pmu0    = u(Eps0,Lz0,K0,M,a,r,th)
  J0      = np.abs(jacobian(Eps0,Lz0,K0,M,a,r,th))
  Eps1    = Eps1[cond1]
  Lz1     = Lz1[cond1]
  K1      = K1[cond1]
  rMin1   = rMin1[cond1]
  rMax1   = rMax1[cond1]
  f1      = f1[cond1]
  pmu1    = u(Eps1,Lz1,K1,M,a,r,th)
  J1      = np.abs(jacobian(Eps1,Lz1,K1,M,a,r,th))
  # Tetrad
  tetrad = ZAMOTetrad(M,a,r,th,inverse=True)
  collisionalEnergies = Ecom(Eps0,Lz0,K0,Eps1,Lz1,K1,M,a,r,th)
  if maxEcom == True:
    return np.max(collisionalEnergies) # Return the max instead of mean
  # Mean COM collisions:
  integrand = f0 * (4./J0) * f1 * (4./J1) * (-1.*Detg) * collisionalEnergies
  # Integrate:
  MEcom = np.sum(integrand) / np.sum(f0 * (4./J0) * f1 * (4./J1) * (-1.*Detg))
  # Check for invalid values
  if np.isnan(MEcom) == True:
    print integrand, N0, N1
    raise ValueError("MEcom is nan; above is the integrand where Nmu = np.sum(integrand)")
  return MEcom



def ZAMOAnnihilation( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the ZAMO Annihilation rate
  '''
  if len(np.shape(r)) != 0:
    raise ValueError("Bad r input to ZAMODispersion; expecting a scalar!")
  if len(np.shape(th)) != 0:
    raise ValueError("Bad th input to ZAMODispersion; expecting a scalar!")
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return 0.
  if len(f0) < 4:
    return 0.
  Delta  = getDelta(M,a,r)
  Sigma  = getSigma(M,a,r,th)
  pmu    = np.array(u(Eps,Lz,K,M,a,r,th))
  Detg   = detg(M,a,r,th)
  J      = np.abs(jacobian(Eps,Lz,K,M,a,r,th))
  N0     = np.sum(f0 * np.sqrt(-1.*Detg) * (4./J))
  nmu    = np.sum(f0 * pmu * np.sqrt(-1.*Detg) * (4./J),axis=1)
  if np.sum(pmu[0] >= 0) != len(pmu[0]):
    raise ValueError("Bad 4-velocity component")
  # Tetrad
  tetrad = ZAMOTetrad(M,a,r,th,inverse=True)
  # Pick random f1 f2:
  Mnorm = 1000000.
  cond1 = np.random.randint(len(f0)-1,size=int(Mnorm))
  cond2 = np.random.randint(len(f0)-1,size=int(Mnorm))
  f1    = f0[cond1]
  J1    = J[cond1]
  pmu1  = pmu[:,cond1]
  N1    = np.sum(f1 * np.sqrt(-1.*Detg) * (4./J1))
  f2    = f0[cond2]
  J2    = J[cond2]
  pmu2  = pmu[:,cond2]
  N2    = np.sum(f2 * np.sqrt(-1.*Detg) * (4./J2))
  pm1   = np.dot(tetrad,pmu1)
  if __debug__:
    pmm1 = np.zeros(4)
    for mu in range(4):
      for m in range(4):
        pmm1[m]= pmm1[m] + tetrad[m][mu] * pmu1[mu,0]
    if np.sum(np.abs(pmm1 - pm1[:,0])) / np.sum(np.abs(pm1[:,0])) > 1.e-5:
      print pm1
      print pmm1
      raise ValueError("Bad value for pm")
  pm2   = np.dot(tetrad,pmu2)
  # Gamma factors: u[0] = c * gamma
  ga1   = pm1[0,:]
  ga2   = pm2[0,:]
  # Three-velocity components
  v1 = pm1[1:,:] / ga1
  v2 = pm2[1:,:] / ga2
  if np.sum(np.linalg.norm( v1, axis=0 ) > 1) != 0:
    print v1
    raise ValueError("3-velocity higher than speed of light for v1")
  if np.sum(np.linalg.norm( v2, axis=0 ) > 1) != 0:
    print v2
    raise ValueError("3-velocity higher than speed of light for v2")
  if np.sum(ga1 < 0 ) != 0:
    print ga1
    raise ValueError("Bad gamma factor")
  vrel  = np.linalg.norm( v2 - v1, axis=0 )
  garel = ga1 * ga2 * (1. - np.sum(v1*v2, axis=0))
  # Calculate the annihilation rate:
  sigma = 1.
  unnormalizedAnnihilationRate = (4./J1)*(4./J2) * (-1.*Detg) * f1 * f2 * sigma * vrel * garel
  unnormalizedAnnihilationRate = np.sum(unnormalizedAnnihilationRate)
  # Normalize it:
  Ntot = float(len(f0) * len(f0))
  annihilationrate = unnormalizedAnnihilationRate * Ntot / Mnorm
  if np.isnan(annihilationrate):
    print (4./J1)*(4./J2) * (-1.*Detg) * f1 * f2 * sigma * vrel * garel
    raise ValueError("Nan anniilationrate; see the printout above for the integrand")
  return annihilationrate


def dtaudt(r,th,a,M):
  factor = np.sqrt(2.)*np.sqrt((((a**2)*(np.cos(th)**2) + (r**2))*((a**2) + r*(-2.*M + r)))/ (2.*(a**4)*(np.cos(th)**2) + (a**2)*(3. + np.cos(2.*th))*(r**2) + 2.*(r**4) + 4.*(a**2)*M*r*(np.sin(th)**2)))
  return factor

def FDensity( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the FFIO Density
  '''
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return 0.
  Pmu = Nmu(f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th)
  g = gmunu(M,a,r,th)
  return np.sqrt(-1.*np.dot(Pmu,np.dot(g,Pmu)))

def ZAMODensity( f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th ):
  ''' Calculates the ZAMO Density
  '''
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = filterBound(f0,Eps,Lz,K,rMin, rMax, thMin, thMax, r, th)
  if success == False:
    return 0.
  Pmu = Nmu(f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th)
  # Calculate the density in the local frame
  tetrad = ZAMOTetrad(M,a,r,th,inverse=True) # e^m_\mu J^\mu
  Nm  = np.dot(tetrad, Pmu)
  return Nm[0]







