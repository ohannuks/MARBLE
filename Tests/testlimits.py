import numpy as np
import pylab as plt
import sys

import matplotlib
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 22}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

print "Usage: python testlimits.py limits.dat biglimits.dat .."

for filename in sys.argv[1:]:
  th = np.pi/2.

  Eps, L, Lz, rmin, rmax = np.transpose( np.loadtxt(filename) )
  
  def integrand(Eps, L, Lz, rmin, rmax, r):
    cond = (rmin+0.005<r)&(r<rmax-0.005)
    Eps = Eps[cond]
    L   = L[cond]
    L   = np.abs(L)
    rmin = rmin[cond]
    rmax = rmax[cond]
    sth2 = 1.
    f=1
    return Eps, L, Lz, rmin, rmax, r, (-4.*np.pi/r**2)*Eps*L*f / np.sqrt(Eps**2-(1-2./r)*(1.+L**2/r**2))
  
  def integrand2(Eps,L,Lz,rmin,rmax,r):
    cond = (rmin+0.005<r)&(r<rmax-0.005)
    Eps = Eps[cond]
    L   = L[cond]
    L   = np.abs(L)
    rmin = rmin[cond]
    rmax = rmax[cond]
    sth2 = 1.
    f=1
    V   = Eps**2-(1.-2./r)*(1.+L**2/r**2)
    return Eps, L, Lz, rmin, rmax, r, (-2./r**2)*Eps*L*f / np.sqrt(V*(L**2*sth2-Lz**2))
  
  R = []
  J0 = []
  samples=[]
  for r in np.linspace(4,100,230):
    Eps2, L2, Lz2, rmin2, rmax2, r2, a = integrand(Eps,L,Lz,rmin,rmax, r)
    R.append(r2)
    J0.append(np.sum(a))
    samples.append(len(a))
  R = np.array(R)
  J0= np.array(J0)
  samples = np.array(samples)
  
  r = np.copy(R)

  Delta = r*(r-2.)
  gpp   = r**2*np.sin(th)
  rho   = -J0 * np.sqrt(gpp/Delta)
  # SADEGHIAN INTEGRATION
  #################################################
  import scipy.integrate as itg
  f0 = 1.0 # when you need to change the height of the whole plot, just change this value
  Gm = 1.0
  def sarho(r):
      rho = []
      emax = 0.999
      for i in r:
          term1 = 4.*np.pi*f0/(1.-2.*Gm/i)**1.5
          if (i >= 6.*Gm):
              emin = (1.+2.*Gm/i)/(1.+6.*Gm/i)**0.5
          else:
              if (i >= 4.*Gm):
                  emin = (1.-2.*Gm/i)/(1.-3.*Gm/i)**0.5
              else:
                  print "wrong range of r, please reset."
                  return 0
          term2 = itg.quad(intepsi,emin,emax,args=(i))[0]
          rho.append(term1*term2)
      rho = np.asarray(rho)
      return rho
  
  def intepsi(e,r):
      Lc2 = 32.*Gm**2./(36.*e**2.-27.*e**4.-8.+e*(9.*e**2.-8.)**1.5)
      return np.sqrt(e**2.-(1.-2.*Gm/r)*(1.+Lc2/r**2.))*e
  sa = sarho(R)
  sa[np.isnan(sa)] = 0
  rho[np.isnan(rho)] = 0
  print rho
  print J0
  #sa = sa* np.max(rho)/np.max(sa)
  rho = rho * sa[np.argmax(samples)]/rho[np.argmax(samples)]
  #################################################
  plt.loglog(R,rho, label='rho, fname = ' + filename)
  plt.loglog(R,sa, label = 'Analytical')
plt.legend(frameon=False)
plt.ylabel(r"Normalized density")
plt.xlabel(r"Radius ($r_g$)")
plt.show()
