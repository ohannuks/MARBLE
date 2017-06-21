import matplotlib
import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR

import matplotlib
#font = {'weight' : 'normal',
#        'size'   : 22}        # Set font
#matplotlib.rc('font', **font) # Set font
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
#matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rc('text', usetex=True)

print "Usage: python createComparisonPlot.py ProductionData/M10_r1_spin0998 ProductionData/M10_r1_spin0998grid.dat /tmp/test.pdf"

fname = sys.argv[1]
fname2 = sys.argv[2]
fout   = sys.argv[3]


# Calculate density from the data file
data =pp.loadpdata(fname2)
Eps2 = data["Eps_2"]
Lz2  = data["Lz_2"]
K2   = data["K_2"]
M    = data["M_2"]
a    = data["a_2"]
rMin2= data["rMin_2"]
rMax2= data["rMax_2"]
thMin2= data["thMin_2"]
thMax2= data["thMax_2"]
g2   = data["g_2"]

# Get the density from Gyoto
###########################################
r,Density1,Px_sigma,Py_sigma,Pz_sigma,nSamples,nParticles=pp.readAveragedRadialData(fname,0,True, normalizedDensity=False)
r = np.array(r)
nParticles = np.ravel(np.array(nParticles))
# Define proper volume in the FFIO frame
th=np.pi/2.+0.01
detgBL = GR.detg(M/M,a/M,r,th)
pt     = GR.u(Eps=1., Lz=0., K=(a/M)**2, M=M/M, a=a/M, r=r, th=th)[0]
#TODO: Make sure this is correct
propernParticles = -1.*nParticles #/ gtt
properVolume     = -1.*np.sqrt(-1.*detgBL) * pt
Density1         = propernParticles / properVolume
###########################################


if GR.tests.isBound(Eps2,Lz2,K2,M,a) == False:
  raise ValueError("Eps2,Lz2,Q2 not bound; something wrong")

rMin2, rMax2 = GR.rMinMax(Eps2,Lz2,K2,M,a)
thMin2, thMax2 = GR.thMinMax(Eps2,Lz2,K2,M,a)


# Calculate the weighted density using FDensity
f = np.ravel(1./g2)
f[f>1.e20] = 0
if np.sum(f<0) != 0 or np.sum(np.isnan(f)) != 0:
  plt.plot(f)
  plt.show()
  print g[np.isnan(f)]
  raise ValueError("Bad f")
r2 = r
th=np.pi/2.+0.01
Density2 = np.array([GR.FDensity(f, Eps2, Lz2, K2, M, a, rMin2, rMax2, thMin2, thMax2, ri, th) for ri in r2])
Density2 = Density2 * Density1[len(Density2)-40]/Density2[len(Density2)-40]




# Calculate the unweighted density using FDensity
f = np.ones(len(f))*1.
r3 = r
th=np.pi/2.+0.01
Density3 = np.array([GR.FDensity(f, Eps2, Lz2, K2, M, a, rMin2, rMax2, thMin2, thMax2, ri, th) for ri in r3])
Density3 = Density3 * Density1[len(Density3)-40]/Density3[len(Density3)-40]


## Calculate the density using Sadeghian's analytical formula
#############
#import scipy.integrate as itg
#emax=0.999
#f0 = 1.
#Gm = 1.
#def intepsi(e,r):
#    Lc2 = 32.*Gm**2./(36.*e**2.-27.*e**4.-8.+e*(9.*e**2.-8.)**1.5)
#    return np.sqrt(e**2.-(1.-2.*Gm/r)*(1.+Lc2/r**2.))*e
#def sarho(r):
#  rho = []
#  emax=0.999
#  f0 = 1.
#  Gm = 1.
#  for i in r:
#    term1 = 4.*np.pi*f0/(1.-2.*Gm/i)**1.5
#    if (i >= 6.*Gm):
#      emin = (1.+2.*Gm/i)/(1.+6.*Gm/i)**0.5
#    else:
#      if (i >= 4.*Gm):
#        emin = (1.-2.*Gm/i)/(1.-3.*Gm/i)**0.5
#      else:
#        rho.append(0)
#        continue
#    term2 = itg.quad(intepsi,emin,emax,args=(i))[0]
#    rho.append(term1*term2)
#  rho = np.asarray(rho)
#  return rho
#analyticalDensity = sarho(r)
#analyticalDensity[np.isnan(analyticalDensity)] = 0
#analyticalDensity = analyticalDensity * Density2[len(Density3)-40]/analyticalDensity[len(Density3)-40]
#############

# Calculate the old way
#############
f = np.ravel(1./g2)
f[f>1.e20] = 0
def integrand2(f,Eps,L,Lz,rmin,rmax,r):
  cond = ((rmin<r) & (r<rmax))
  Eps = Eps[cond]
  L   = L[cond]
  Lz  = Lz[cond]
  L   = np.abs(L)
  rmin = rmin[cond]
  rmax = rmax[cond]
  f = f[cond]
  sth2 = 1.
  V   = Eps**2-(1.-2./r)*(1.+L**2/r**2)
  K = L**2
  if len(Eps) == 0:
    return Eps, L, Lz, rmin, rmax, r, np.array([0])
  if np.sum(V*(L**2*sth2-Lz**2) < 0)>0:
    raise ValueError("Bad denominator " + str(np.min(V)) + " " + str(np.min(L**2*sth2-Lz**2)))
  else:
    return Eps, L, Lz, rmin, rmax, r, (-2./r**2)*Eps*L*f / np.sqrt(V*(L**2*sth2-Lz**2))
#R = []
#J0 = []
#samples=[]
#for r in r2:
#  L2 = np.sqrt(K2)
#  Eps, L, Lz, rmin, rmax, rnull, result = integrand2(f,Eps2,L2,Lz2,rMin2,rMax2, r)
#  R.append(r)
#  J0.append(np.sum(result[(np.isnan(result)==False)&(np.abs(result)<1.e20)]))
#  samples.append(len(result))
#R = np.array(R)
#J0= np.array(J0)
#samples = np.array(samples)
#r = np.copy(R)
#Delta = r*(r-2.)
#gpp   = r**2*np.sin(th)
#rho   = -J0 * np.sqrt(gpp/Delta)
#rho = rho * Density1[len(rho)-40]/rho[len(rho)-40]

#############

#plt.figure(figsize=(16,9))
plt.loglog(r, Density1,label='Sampled with Gyoto', lw=3)
plt.loglog(r2,Density2,label='with the degeneracy factor', lw=3)
plt.loglog(r3,Density3,label='Without the degeneracy factor',lw=3)
#plt.loglog(r, rho,     label="old density")
#plt.loglog(r,analyticalDensity,label='Analytical Density')
plt.legend()
plt.xlabel(r"r ($r_g$)")
plt.ylabel(r"Density")
plt.tight_layout()
plt.savefig(fout,bbox_tight=True)
plt.show()

plt.plot(r,pt)
plt.ylim([-2,2])
plt.show()



R = np.copy(r)
for r in R:
  success, Eps, Lz, K, rMin, rMax, thMin, thMax, f0 = GR.filterBound(f,Eps2,Lz2,K2,rMin2, rMax2, thMin2, thMax2, r, th)
  if success == False:
    continue
  dth = thMax-thMin
  Pmu = GR.Nmu(f0, Eps, Lz, K, M, a, rMin, rMax, thMin, thMax, r, th)
  g = GR.gmunu(M,a,r,th)
  u_d = GR.u_d(Eps,Lz,K,M,a,r,th)
  print r, len(Eps), np.max(np.sqrt(-1.*np.dot(Pmu,np.dot(g,Pmu)))), dth[np.argmin(u_d[2])]


dth = thMax-thMin
i = np.argmin(dth)























