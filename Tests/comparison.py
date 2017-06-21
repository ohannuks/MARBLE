import numpy as np
import pylab as plt
import sys
import paris as pp
from paris import GR
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

if len(sys.argv) != 5:
  raise ValueError("Bad input to comparison.py")

fname1 = sys.argv[1]
fname2 = sys.argv[2]
fout   = sys.argv[3]

rsmall = 4
rbig   = 100

# Read first file
##########################################
# Get the first 3 lines of the data
# Get the first 3 lines of the data
attr = pp.loadattributes(fname1)
Ntotal            = attr[0]
Nin               = attr[1]
a                 = attr[2]
M                 = attr[3]
a_2               = attr[4]
M_2               = attr[5]
if M != 1:
  raise ValueError("Expecing M = 1 and a = 0 values for comparison! This is because the old testlimits has always assumed M=1")
# Load data
data = pp.loaddata(fname1)
Eps1  = data[:,0]
Lz1   = data[:,1]
Q1    = data[:,2]
# Calculate K constant of motion
K1  = Q1+(Eps1*a-Lz1)**2

rMin1, rMax1 = GR.rMinMax(Eps1,Lz1,K1,M,a)
thMin1,thMax1= GR.thMinMax(Eps1,Lz1,K1,M,a)

#if GR.tests.isBound(Eps1,Lz1,K1,M,a) == False:
#  raise ValueError("Eps1,Lz1,Q1 not bound; something wrong")
##########################################

# Calculate from the 4-current J^\mu J_\nu = - \rho^2
#################################################
r1 = np.logspace(np.log10(rsmall*M), np.log10(rbig*M),150)
f = np.sqrt(K1) # We are not sampling homogenously in E,L,K space
th   = np.pi/2.
Jmu  = np.array([GR.Nmu(f, Eps1, Lz1, K1, M, a, rMin1, rMax1, thMin1, thMax1, ri, th) for ri in r1])
rho2 = np.zeros(len(Jmu))
rhosimple2 = np.zeros(len(Jmu))
for i in range(len(Jmu)):
  gmunu = GR.gmunu(M, a, r1[i], th)
  for mu in xrange(4):
    for nu in xrange(4):
      rho2[i] = rho2[i] -1.*gmunu[mu][nu] * Jmu[i][mu] * Jmu[i][nu]
  rhosimple2[i] = -1.*gmunu[0][0] * Jmu[i][0]**2
Density0 = np.sqrt(rho2)
Density4 = np.sqrt(rhosimple2)
# Error:
nSamples = np.array([np.sum((rMin1<ri)&(ri<rMax1)) for ri in r1])
err0     = Density0 / np.sqrt(nSamples)
err4     = Density0 / np.sqrt(nSamples)
#################################################



# Calculate the density via monte carlo:
# TODO: This is incorrect
print "Warning: FFIODensity is incorrect by comparing with the analytical formulae"
f = np.ones(len(Eps1))*1.
Density1 = np.array([GR.FFIOdensity2(f,Eps1, Lz1, K1, M, a, rMin1, rMax1, ri) for ri in r1])
nSamples = np.array([np.sum((rMin1<ri)&(ri<rMax1)) for ri in r1])
err1     = Density1 / np.sqrt(nSamples)



# Read second file
##########################################
th = np.pi/2.

Eps, Lz, Q, rMin, rMax = np.transpose( np.loadtxt(fname2) )
K = Q+(Eps*a-Lz)**2
L = np.sqrt(K)

if GR.tests.isBound(Eps,Lz,K,M,a) == False:
  raise ValueError("Eps,Lz,K not bound; something wrong")


def integrand(Eps, L, Lz, rmin, rmax, r):
  cond = ((rmin+0.005<r)&(r<rmax-0.005))
  Eps = Eps[cond]
  L   = L[cond]
  Lz  = Lz[cond]
  L   = np.abs(L)
  rmin = rmin[cond]
  rmax = rmax[cond]
  sth2 = 1.
  f=1
  return Eps, L, Lz, rmin, rmax, r, (-4.*np.pi/r**2)*Eps*L*f / np.sqrt(Eps**2-(1-2./r)*(1.+L**2/r**2))

def integrand2(Eps,L,Lz,rmin,rmax,r):
  cond = ((rmin<r) & (r<rmax))
  Eps = Eps[cond]
  L   = L[cond]
  Lz  = Lz[cond]
  L   = np.abs(L)
  rmin = rmin[cond]
  rmax = rmax[cond]
  sth2 = 1.
  f=1
  V   = Eps**2-(1.-2./r)*(1.+L**2/r**2)
  K = L**2
  if np.sum(V*(L**2*sth2-Lz**2) <= 0)>0:
    raise ValueError("Bad denominator " + str(np.min(V)) + " " + str(np.min(L**2*sth2-Lz**2)) + " " + str(np.min(GR.R(Eps, Lz, K, M, a, r))))
  #print np.sum((L**2*sth2-Lz**2)<=0)
  #print np.sum(cond==True)
  #print len(Eps)
  #print (-2./r**2)*Eps*L*f / np.sqrt(V*(L**2*sth2-Lz**2))
  #exit(1)
  if len(Eps) == 0:
    return Eps, L, Lz, rmin, rmax, r, np.array([0])
  else:
    return Eps, L, Lz, rmin, rmax, r, (-2./r**2)*Eps*L*f / np.sqrt(V*(L**2*sth2-Lz**2))

R = []
J0 = []
samples=[]
for r in r1:
  Eps2, L2, Lz2, rmin2, rmax2, r2, result = integrand2(Eps,L,Lz,rMin,rMax, r)
  R.append(r2)
  J0.append(np.sum(result[np.isnan(result)==False]))
  samples.append(len(result))
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
r2 = R
Density2 = rho
err2     = Density2 / np.sqrt(samples)


# Normalize density:
##########################################
Density2 = Density2 / Density2[70] * Density1[70]
sa = sa / sa[70] * Density1[70]
Density0 = Density0 / Density0[70] * Density1[70]
##########################################
fig = plt.figure(figsize=(16,9))
if sys.argv[4] == "production":
  plt.plot(r1,sa,'-', label="Analytical result",markersize=3,color='black',lw=3)
  plt.errorbar(r1,Density1, fmt='o', label=r'Simulated result',yerr=err1,markersize=3,color='blue')
  fig.get_axes()[0].set_yscale("log",nonposx='clip')
  fig.get_axes()[0].set_xscale("log",nonposx='clip')
elif sys.argv[4] == "log":
  #plt.loglog(r1,sa, label="Sadeghian's analytical result",lw=3)
  #plt.loglog(r1,Density1, label=r'rho new',lw=3)
  #plt.loglog(r2,Density2, label=r'rho old',lw=3)
  #plt.loglog(r1,Density0, label=r'Calculated from 4-current',lw=3)
  plt.plot(r1,sa, label="Sadeghian's analytical result",markersize=3)
  plt.errorbar(r1,Density1, label=r'rho new',yerr=err1,markersize=3)
  plt.errorbar(r2,Density2, label=r'rho old',yerr=err2,markersize=3)
  plt.errorbar(r1,Density0, label=r'Calculated from 4-current',yerr=err0,markersize=3)
  fig.get_axes()[0].set_yscale("log",nonposx='clip')
  fig.get_axes()[0].set_xscale("log",nonposx='clip')
else:
  plt.plot(r1,sa, label="Sadeghian's analytical result",lw=3)
  plt.errorbar(r1,Density1, label=r'rho new',yerr=err1,lw=3)
  plt.errorbar(r2,Density2, label=r'rho old',yerr=err2,lw=3)
  plt.errorbar(r1,Density0, label=r'Calculated from 4-current',yerr=err0,lw=3)
plt.legend(frameon=False)
plt.ylabel(r"Relative density")
plt.xlabel(r"Radius ($r_g$)")
plt.tight_layout()
plt.grid()
plt.savefig(fout,bbox_tight=True)
plt.show()
