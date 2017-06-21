import sys
import os
sys.path.insert(0, "..\Tools")
from lmfit import models
import numpy as np
import pylab as plt
import paris as pp
import matplotlib
from scipy.interpolate import spline
import matplotlib
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 24}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

global N
N=2
def loadProcessedVariables(M,n,spin):
  return np.loadtxt("ProductionData/M%d_r%d_spin%04d.radial"%(M,n,spin))


def getNormalizationFactor(thmaxspin0=np.pi, thmax=np.pi/8.):
    ''' Calculates the normalization factor for all simulations, defined as the minimum peak density of all simulations

    '''
    n=N
    masses = [10,11,12,13,14,15,16,17,18,19]
    peakdensities = []
    for M in masses:
        R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt("ProductionData/M%d_r%d_spin0000.radial"%(M,n))
        Densityspin0=np.copy(FFIODensity)
        R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt("ProductionData/M%d_r%d_spin0998.radial"%(M,n))
        Densityspin1=np.copy(FFIODensity)
        #print Densityspin0, Densityspin1
        peakdensities.append(min(max(Densityspin0[(np.isnan(Densityspin0)==False)&(Densityspin0>0)]),max(Densityspin1[(np.isnan(Densityspin1)==False)&(Densityspin1>0)])))
    return min(peakdensities)


normalizationFactor = getNormalizationFactor()

from scipy.special import erf
from scipy import stats
from scipy.optimize import curve_fit

def SkewGaussian(logr, G, A, mu, sigma, gamma):
    A0 = A(G)
    mu0= mu(G)
    sigma0= sigma(G)
    gamma0= gamma(G)
    factor = A0/(sigma0*np.sqrt(2.*np.pi))
    value = np.exp(-(logr-mu0)**2/(2.*sigma0**2))*(1.+erf(gamma0*(logr-mu0)/(sigma0*np.sqrt(2))))
    return factor*value

def SkewGaussian2(logr, A, mu, sigma, gamma):
    A0 = A
    mu0= mu
    sigma0= sigma
    gamma0= gamma
    factor = A0/(sigma0*np.sqrt(2.*np.pi))
    value = np.exp(-(logr-mu0)**2/(2.*sigma0**2))*(1.+erf(gamma0*(logr-mu0)/(sigma0*np.sqrt(2))))
    return factor*value

def calculateError(r, y, nSamples):
  yerr=1.96*np.ravel(y)*np.sqrt(1./np.ravel(nSamples));yerr[np.isnan(yerr)] = np.mean(yerr[np.isnan(yerr) == False])
  return yerr

def linefit(x, y, yerr):
    def line(x, a, b):
      return a*x+b
    popt,pcov=curve_fit(line,x,y,sigma=np.ones(len(x)))#,sigma=yerr)
    p=np.poly1d([popt[0],popt[1]])
    perr = np.sqrt(np.diag(pcov))
    return p, perr


def skewGaussianFitting(thmaxspin0=np.pi,thmax=np.pi/8., growthFactor=True):
  import pylab as plt
  normalizationFactor = getNormalizationFactor()
  #cond = Density !=0
  #R = R[cond]
  #Density = Density[cond]
  growthFactor = False
  from scipy.optimize import curve_fit
  n=N
  masses = [10,11,12,13,14,15,16,17,18,19]
  growthDensities = []
  growthDensitiesError = []
  for M in masses:
    import numpy as np
    R,Density,ZAMODensity,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt("ProductionData/M%d_r%d_spin0998.radial"%(M,n));# plt.plot(R,Density);plt.show()
    #print normalizationFactor
    Density = Density / normalizationFactor
    yerr=calculateError(R,Density,nSamples) 
    p0 = [np.mean(Density),8.5,1,0]
    p0 = [np.mean(Density),8,2,-1]
    logr = np.log(R)
    popt,pcov=curve_fit(SkewGaussian2,logr,Density,p0)#,sigma=yerr)
    A = popt[0]; mu = popt[1]; sigma=popt[2]; gamma=popt[3]
    amplitude=A 
    center=mu 
    sigma =sigma 
    gamma=gamma 
    # Error: 
    perr = np.sqrt(np.diag(pcov)) 
    growthDensities.append([M,amplitude,center,sigma,gamma]) 
    growthDensitiesError.append([M,perr[0],perr[1],perr[2],perr[3]])
    plt.figure()
    plt.plot(R,Density)
    plt.plot(R,SkewGaussian2(logr,A,mu,sigma,gamma))
    plt.show()
  #plt.figure()
  import matplotlib.pyplot as plt
  import numpy as np
  growthDensities = np.array(growthDensities)
  growthDensitiesError = np.array(growthDensitiesError)
  masses       =growthDensities[:,0]
  amplitudes   =growthDensities[:,1]
  amplitudesErr=growthDensitiesError[:,1]
  centers      =growthDensities[:,2]
  centersErr   =growthDensitiesError[:,2]
  sigmas       =growthDensities[:,3]
  sigmasErr    =growthDensitiesError[:,3]
  gammas       =growthDensities[:,4]
  gammasErr    =growthDensitiesError[:,4]
  #residuals=growthDensities[:,5]
  if growthFactor == True:
    G = 1./masses
  else:
    G = masses
  plt.figure(figsize=(16,9))
  ax1 = plt.subplot(511)
  plt.semilogy(G, amplitudes,'.',ms=10, mew=3, color='r')
  err = (np.log10(amplitudes+amplitudesErr)-np.log10(amplitudes-amplitudesErr))/2.
  p, err = linefit(masses,np.log10(amplitudes),err)
  k0 = p[1]; c0=p[0]
  print "log10(A)="+str(k0)+"+-"+str(err[0])+" G + " + str(c0) + "+-"+str(err[1])
  plt.semilogy(G,10.**p(masses), lw=3, color='b')
  plt.setp(ax1.get_xticklabels(), visible=False)
  # share x only
  ax2 = plt.subplot(512, sharex=ax1)
  #plt.plot(G,centers,'.',ms=10, mew=3, color='r')#,yerr=centersErr)
  plt.errorbar(G,centers,fmt='.',ms=10,mew=3,yerr=centersErr,color='r')
  plt.ylim([min(centers)*0.95,max(centers)*1.05])
  #plt.yticks(0.8+np.arange(3)*0.4)
  err = (np.exp(centers+centersErr)-np.exp(centers-centersErr))/2.
  p, err = linefit(masses, np.exp(centers), err)
  k1=p[1]; c1=p[0]
  print "exp(mu)="+str(k1)+"+-"+str(err[0])+" G + " + str(c1)+"+-"+str(err[1])
  plt.plot(G,np.log(p(masses)), lw=3, color='b')
  #plt.ylim([0.9*np.min(np.log(p(G))),1.1*np.max(np.log(p(G)))])
  # make these tick labels invisible
  plt.setp(ax2.get_xticklabels(), visible=False)
  # share x and y
  ax3 = plt.subplot(513, sharex=ax1)
  plt.errorbar(G,sigmas,fmt='.',yerr=sigmasErr,ms=10, mew=3, color='r')#,yerr=sigmasErr)
  #plt.yticks(np.arange(3)*0.05+0.85)
  #plt.ylim([0.9*np.min(sigmas),1.1*np.min(sigmas)])
  err=np.copy(sigmasErr)
  p, err = linefit(masses, sigmas, err)
  k2=p[1]; c2=p[0]
  print "sigma="+str(k2) + "+-" + str(err[0])+" G + " + str(c2) +"+-"+str(err[1])
  plt.plot(G,p(masses), lw=3, color='b')
  plt.ylim([min(sigmas)*0.95,max(sigmas)*1.05])
  plt.setp(ax3.get_xticklabels(), visible=False)
  #plt.yticks([0.9,0.95,1.0])
  #plt.ylim([0.85,1.05])
  # share x and y
  ax4 = plt.subplot(514, sharex=ax1)
  plt.errorbar(G,gammas,fmt='.',yerr=gammasErr,ms=10, mew=3, color='r')#,yerr=gammasErr)
  err=np.copy(gammasErr)
  p,err = linefit(masses, gammas,err)
  k3=p[1]; c3=p[0]
  print "gamma="+str(k3)+"+-"+str(err[0])+" G + " + str(c3)+"+-"+str(err[1])
  plt.plot(G, p(masses), lw=3, color='b')
  plt.ylim([min(gammas)*0.95,max(gammas)*1.05])
  #plt.yticks([6,12,18])
  #plt.ylim([0.9*np.min(gammas),1.1*np.max(gammas)])
  plt.setp(ax4.get_xticklabels(), visible=False)
  
  # Finally, do the residual fitting:
  global A, mu, sigma, gamma
  A = lambda G: 10**(k0*G+c0)
  mu= lambda G: np.log(k1*G+c1)
  sigma= lambda G: k2*G+c2
  gamma= lambda G: k3*G+c3
  
  residuals = []
  for i in xrange(len(masses)):
      M = masses[i]
      n = 1
      R,Density,ZAMODensity,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt("ProductionData/M%d_r%d_spin0998.radial"%(M,n));# plt.plot(R,Density);plt.show()
      # Get radial plot for spin 0, M10
      #R = spin0[M]["R"]; #logr = np.log(R)
      #G = M # By our definition, growth factor = m
      # set initial parameter values
      # set initial parameter values
      Density = Density/normalizationFactor
      #params = model.make_params(amplitude=np.max(Density)*10**8, center=10, sigma=1, gamma=0)
      #print Density
      numberOfSamples=nSamples
      yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples))
      #print yerr
      A = amplitudes[i]
      mu= centers[i]
      sigma=sigmas[i]
      gamma=gammas[i]
      fit = SkewGaussian2(logr, A, mu, sigma, gamma)
      #result = model.fit(Density, params, x=np.log(np.ravel(spin1[M]["R"])),yerr=yerr)
      #fit = result.best_fit
      Density = Density[fit>0.001*max(fit)] # Trimming!
      fit = fit[fit>0.01*max(fit)]
      yerr= yerr[fit>0.01*max(fit)]
      check0 = (fit>0.01*max(fit)) & (yerr!=0) & (np.isnan(yerr)==0)
      w = 1./yerr[check0]
      print w, np.sum(w)
      residual = np.sum((w/np.sum(w)) * (np.abs(Density[check0] - fit[check0])/fit[check0]))
      residuals.append(residual)
  
  # Share x and y
  ax5 = plt.subplot(515, sharex=ax1)
  #plt.xlabel(r"Growth factor $\left(500 \cdot \frac{10 \cdot M_{seed}}{M_{final}}\right)$")
  plt.plot(G,100.*np.asarray(residuals),'.',ms=10, mew=3, color='r')
  #print residuals
  #plt.yticks([2,3,4,5])
  if growthFactor == True:
    plt.xlabel(r"Growth factor $G$")
  else:
    plt.xlabel(r"$5000R_{IF}$")
  
  #ax1.set_ylim([10**(-12),10**(-7)])
  #ax1.set_yticks([-11,-10,-9,-8])
  ax1.set_ylabel(r"$A$")
  ax2.set_ylabel(r"$\mu$")
  ax3.set_ylabel(r"$\widetilde{\sigma}$")
  #ax3.set_ylim([np.min(sigmas)*0.95,np.max(sigmas)*1.05])
  ax4.set_ylabel(r"$\gamma$")
  #ax4.set_ylim([np.min(gammas)*0.6,np.max(gammas)*1.05])
  #ax5.set_ylabel(r"Residual")
  ax5.set_ylabel(r"MAPE(\%)")
  plt.xlim([np.min(G)-0.2,np.max(G)+0.2])
  #plt.tight_layout()
  plt.show()
  return k0,k1,k2,k3,c0,c1,c2,c3


from paris import GR
def gridAnnihilation(fname):
  ''' Calculates the rate of annihilation within 15 M on (r, theta) grid

      :returns: r, th, annihilationRates, sqrt(dispersions)
  '''
  # Calculate density from the data file
  data =pp.loadpdata(fname)
  Ntotal = data["Ntotal"]
  Eps_1  = data["Eps_1"]
  Lz_1   = data["Lz_1"]
  K_1    = data["K_1"]
  M_1    = data["M_1"]
  a_1    = data["a_1"]
  rMin_1 = data["rMin_1"]
  rMax_1 = data["rMax_1"]
  g_1    = np.ravel(data["g_1"])
  
  Eps_2  = data["Eps_2"]
  Lz_2   = data["Lz_2"]
  K_2    = data["K_2"]
  M_2    = data["M_2"]
  a_2    = data["a_2"]
  rMin_2 = data["rMin_2"]
  rMax_2 = data["rMax_2"]
  thMin_2 = data["thMin_2"]
  thMax_2 = data["thMax_2"]
  g_2    = np.ravel(data["g_2"])
  
  N_1 = np.ones(len(Eps_1))*1. / Ntotal
  f_1 = N_1 / g_1
  f_2 = f_1 # Constancy of phase space
  
  R = np.linspace(M_2,15.*M_2,50)
  TH= np.linspace(np.pi/3.,np.pi*(1.-1./3.),10)
  TH= [np.pi/2.+0.01]
  rV                = []
  thV               = []
  annihilationRates = []
  sqrtdispersions       = []
  nSamples          = []
  # Calculate annihilation rate for each "Grid" point
  for r in R:
    for th in TH:
      cond = ((rMin_2<r)&(r<rMax_2)&(thMin_2<th)&(th<thMax_2))
      annihilationRate = GR.ZAMOAnnihilation( f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r, th)
      annihilationRates.append(annihilationRate* 1.e45)
      sqrtdispersions.append(np.sqrt(GR.ZAMODispersion(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r, th)))
      rV.append(r)
      thV.append(th)
      nSamples.append(np.sum(cond))
      print "r, th: %f %f"%(r,th)
  print "Max annihilation rate %lf"%np.max(annihilationRates)
  print "Mean value of degeneracy: %lf"%np.mean(g_2)
  print "Minimum value of degeneracy: %lf"%np.min(g_2)
  print "Maximum value of degeneracy: %lf"%np.max(g_2)
  return rV, thV, annihilationRates, sqrtdispersions, nSamples, M_2, a_2

def gridDensity(fname):
  ''' Calculates the density squared within 15 M on (r, theta) grid

      :returns: r, th, annihilationRates, sqrt(dispersions)
  '''
  # Calculate density from the data file
  data =pp.loadpdata(fname)
  Ntotal = data["Ntotal"]
  Eps_1  = data["Eps_1"]
  Lz_1   = data["Lz_1"]
  K_1    = data["K_1"]
  M_1    = data["M_1"]
  a_1    = data["a_1"]
  rMin_1 = data["rMin_1"]
  rMax_1 = data["rMax_1"]
  g_1    = np.ravel(data["g_1"])
  
  Eps_2  = data["Eps_2"]
  Lz_2   = data["Lz_2"]
  K_2    = data["K_2"]
  M_2    = data["M_2"]
  a_2    = data["a_2"]
  rMin_2 = data["rMin_2"]
  rMax_2 = data["rMax_2"]
  thMin_2 = data["thMin_2"]
  thMax_2 = data["thMax_2"]
  g_2    = np.ravel(data["g_2"])
  
  N_1 = np.ones(len(Eps_1))*1. / Ntotal
  f_1 = N_1 / g_1
  f_2 = f_1 # Constancy of phase space
  
  R = np.linspace(M_2,15.*M_2,80)
  TH= np.linspace(np.pi/3.,np.pi*(1.-1./3.),10)
  TH= [np.pi/2.+0.02]
  rV                = []
  thV               = []
  Densitys = []
  sqrtdispersions       = []
  nSamples          = []
  # Calculate annihilation rate for each "Grid" point
  for r in R:
    for th in TH:
      cond = ((rMin_2<r)&(r<rMax_2)&(thMin_2<th)&(th<thMax_2))
      Density = GR.ZAMODensity( f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r, th)**2
      Densitys.append(Density* 1.e45)
      sqrtdispersions.append(np.sqrt(GR.ZAMODispersion(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r, th)))
      rV.append(r)
      thV.append(th)
      nSamples.append(np.sum(cond))
      print "r, th: %f %f"%(r,th)
  print "Max annihilation rate %lf"%np.max(Densitys)
  print "Mean value of degeneracy: %lf"%np.mean(g_2)
  return rV, thV, Densitys, sqrtdispersions, nSamples, M_2, a_2

