import sys
import os
sys.path.insert(0, "G:\Data\git\DAMAHAO\Tools")
import numpy as np
import pylab as plt
import paris as pp
import matplotlib
from scipy.interpolate import spline
if not os.path.exists("Figs"): os.mkdir("Figs")
############## Setting ##########################
output_path='./Figs/PeakGrowth.pdf'
#################################################
font = {'family' : 'normal',
    'weight' : 'normal',
    'size'   : 22}

matplotlib.rc('font', **font)
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

# Calculate annihilation for every value
def calculateAnnihilation(Px, Py, Pz, Density, Px_sigma, Py_sigma, Pz_sigma, crossSection, threshold=0 ):
  if len(np.shape(Density)) > 1:
   print "BAD LENGTH"
   Density = np.ravel(Density)
   Px = np.ravel(Px)
   Py = np.ravel(Py)
   Pz = np.ravel(Pz)
   Px_sigma = np.ravel(Px_sigma)
   Py_sigma = np.ravel(Py_sigma)
   Pz_sigma = np.ravel(Pz_sigma)

  annihilations = []
  for i in xrange(len(np.ravel(Px))):
    P = np.array([Px[i], Py[i], Pz[i]])
    rho = Density[i]
    dispersion =np.array([Px_sigma[i], Py_sigma[i], Pz_sigma[i]]) 
    if rho == 0 or rho**2 < threshold*np.max(rho)**2:
      annihilation = 0
    else:
      annihilation = pp.rateOfAnnihilationCuba( P, 1.*rho, dispersion, crossSection )
    annihilations.append(annihilation)
  return annihilations



def annihilationMomentumInterplay( thmaxspin1=np.pi/8., peaksOnly=False ):
  from lmfit import models
  # Plot masses
  # Read mass files
  spin0={}
  spin1={}
  masses = [6,7,9,10,11,12,13,14,15,16,17,18]
  for M in masses:
    R,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples,Px,Py,Pz=pp.readAveragedRadialData("M%02d_r2_spin0000.dat"%M,thmax=np.pi,nSamples=True, moments=True)
    threshold=0.05
    annihilations = calculateAnnihilation(Px, Py, Pz, Density, Px_sigma, Py_sigma, Pz_sigma, 1.,threshold)
    spin0[M]={"M":M,"R":R,"Density":Density,"Px_sigma":Px_sigma,"Py_sigma":Py_sigma,"Pz_sigma":Pz_sigma,"nSamples":nSamples,"Px":Px,"Py":Py,"Pz":Pz, "annihilations": annihilations}
    R,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples,Px,Py,Pz=pp.readAveragedRadialData("M%02d_r2_spin0998.dat"%M,thmax=thmaxspin1,nSamples=True,moments=True)
    annihilations = calculateAnnihilation(Px, Py, Pz, Density, Px_sigma, Py_sigma, Pz_sigma, 1.,threshold)
    spin1[M]={"M":M,"R":R,"Density":Density,"Px_sigma":Px_sigma,"Py_sigma":Py_sigma,"Pz_sigma":Pz_sigma,"nSamples":nSamples,"Px":Px,"Py":Py,"Pz":Pz, "annihilations": annihilations}
  fig, ax1 = plt.subplots(figsize=(16, 9))
  growthAnnihilations=[]
  growthPxsigmas =[]
  growthPysigmas =[]
  growthPzsigmas =[]
  peakPxsigmas =[]
  peakPysigmas =[]
  peakPzsigmas =[]
  for M in masses:
    # Get variables
    R = spin0[M]["R"]
    annSpin0 = spin0[M]["annihilations"]
    annSpin1 = spin1[M]["annihilations"]
    densitySpin0 = spin0[M]["Density"]
    densitySpin1 = spin1[M]["Density"]
    sigmaxSpin0 = spin0[M]["Px_sigma"]
    sigmaxSpin1 = spin1[M]["Px_sigma"]
    sigmaySpin0 = spin0[M]["Py_sigma"]
    sigmaySpin1 = spin1[M]["Py_sigma"]
    sigmazSpin0 = spin0[M]["Pz_sigma"]
    sigmazSpin1 = spin1[M]["Pz_sigma"]
    samplesSpin0 = spin0[M]["nSamples"]
    samplesSpin1 = spin1[M]["nSamples"]

    maxsigmaPxspin0 = sigmaxSpin0[np.argmax(annSpin0)] 
    maxsigmaPxspin1 = sigmaxSpin1[np.argmax(annSpin1)]
    maxsigmaPyspin0 = sigmaySpin0[np.argmax(annSpin0)]
    maxsigmaPyspin1 = sigmaySpin1[np.argmax(annSpin1)]
    maxsigmaPzspin0 = sigmazSpin0[np.argmax(annSpin0)]
    maxsigmaPzspin1 = sigmazSpin1[np.argmax(annSpin1)]

    sigmaPxsamplesSpin0 = samplesSpin0[np.argmax(annSpin0)]
    sigmaPxsamplesSpin1 = samplesSpin1[np.argmax(annSpin1)]
    sigmaPysamplesSpin0 = samplesSpin0[np.argmax(annSpin0)]
    sigmaPysamplesSpin1 = samplesSpin1[np.argmax(annSpin1)]
    sigmaPzsamplesSpin0 = samplesSpin0[np.argmax(annSpin0)]
    sigmaPzsamplesSpin1 = samplesSpin1[np.argmax(annSpin1)]

    growthAnnihilations.append([M,np.max(annSpin0), np.max(annSpin1), samplesSpin0[np.argmax(annSpin0)], samplesSpin1[np.argmax(annSpin1)]])
    growthPxsigmas.append( [M,maxsigmaPxspin0, maxsigmaPxspin1, sigmaPxsamplesSpin0, sigmaPxsamplesSpin1])
    growthPysigmas.append( [M,maxsigmaPyspin0, maxsigmaPyspin1, sigmaPysamplesSpin0, sigmaPysamplesSpin1])
    growthPzsigmas.append( [M,maxsigmaPzspin0, maxsigmaPzspin1, sigmaPzsamplesSpin0, sigmaPzsamplesSpin1])

  growthAnnihilations = np.array(growthAnnihilations)
  growthPxsigmas  = np.array(growthPxsigmas)
  growthPysigmas  = np.array(growthPysigmas)
  growthPzsigmas  = np.array(growthPzsigmas)
  peakPxsigmas = growthPxsigmas
  peakPysigmas = growthPysigmas
  peakPzsigmas = growthPzsigmas

  M=growthAnnihilations[:,0]
  #M=500./growthAnnihilations[:,0]
  global normalizationFactor
  normalizationFactor=np.min(growthAnnihilations[:,1])
  annihilationsSpin0=growthAnnihilations[:,1]/np.min(growthAnnihilations[:,1])
  annihilationsSpin1=growthAnnihilations[:,2]/np.min(growthAnnihilations[:,1])
  ax1.plot(M, annihilationsSpin0,"-x",mew=3,label=r'$\rho_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
  ax1.plot(M, annihilationsSpin1,"-*",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
  #ax1.semilogy(M, 10**(p(M)),"-",label=r"$\log(\rho_{peak})=k\cdot G + c$ best fit",color='black',lw=1)
  #ax1.semilogy(M, annihilationsSpinfit,"--o",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3)

  #plt.ylim([1,4e3])
  #ax1.plot(M, annihilationsSpin0/np.min(annihilationsSpin0),"--o",label=r'$\rho_{peak}$; $a=0.0$',color='black',lw=1)
  #ax1.plot(M, annihilationsSpin1/np.min(annihilationsSpin0),"--o",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3)
  ax1.set_ylabel(r"Peak Relative Density")

  plt.grid()
  ax2 = ax1.twinx()
  #ax2.plot([],[],"-*",label=r'$\rho_{peak}$, $a=0.0$',color='red',lw=3,markersize=17)
  ax2.plot([],[],"-x",mew=3,label=r'$\rho_{peak}$, $a=0.0$',color='red',lw=3,markersize=17)
  ax2.plot([],[],"-*",label=r"$\rho_{peak}$, $a=0.998$",color='black',lw=3,markersize=17)
  #ax2.plot(M, growthPxsigmas[:,1],"-*",label=r'$\sigma_{r, peak}^2$; $a=0.0$',color='black',lw=1)
  if peaksOnly==False:
    ax2.plot(M, np.sqrt(growthPysigmas[:,1]),"--o",label=r'$\sigma_{\theta, max}$, $a=0.0$',color='red',markersize=7,lw=3)
  #ax2.plot(M, np.sqrt(peakPysigmas[:,1]),"--*",label=r'$\sigma_{\theta, peak}$, $a=0.0$',color='red',markersize=17,lw=3)
  #ax2.plot(M, np.sqrt(peakPysigmas[:,1]),"--x",mew=3,label=r'$\sigma_{\theta, peak}$, $a=0.0$',color='red',markersize=17,lw=3)
  #ax2.errorbar(M, np.sqrt(peakPysigmas[:,1]),color='red',lw=0,yerr=2.*peakPysigmas[:,1]/peakPysigmas[:,3])
  ax2.errorbar(M, np.sqrt(peakPysigmas[:,1]),fmt='--o',color='red',yerr=4.*np.sqrt(peakPysigmas[:,1])/peakPysigmas[:,3],label=r'$\sigma_{\theta,peak}$, $a=0.998$',markersize=7,lw=3)
  #ax2.errorbar(M, np.sqrt(growthPysigmas[:,1]),fmt="o",color='red',lw=1,yerr=2.*growthPysigmas[:,1]/growthPysigmas[:,3])
  #ax2.plot(M, growthPzsigmas[:,1],"-*",label=r'$\sigma_{\phi, peak}$; $a=0.0$',color='red',lw=1)
  #ax2.errorbar(M, growthPzsigmas[:,1],fmt="o",color='red',lw=1,yerr=4.*growthPzsigmas[:,1]/growthPzsigmas[:,3])
  #ax2.plot(M, growthPxsigmas[:,2],"-o",label=r'$\sigma_{r, peak}^2$; $a=0.998$',color='black',lw=3)
  #ax2.errorbar(M, growthPxsigmas[:,2],fmt="o",color='black',lw=1,yerr=4.*growthPxsigmas[:,2]/growthPxsigmas[:,4])
  if peaksOnly==False:
    ax2.plot(M, np.sqrt(growthPysigmas[:,2]),"--o",label=r'$\sigma_{\theta, max}$, $a=0.998$',color='black',markersize=7,lw=3)
  #ax2.plot(M, np.sqrt(peakPysigmas[:,2]),"--*",label=r'$\sigma_{\theta, peak}$, $a=0.998$',color='black',markersize=17,lw=3)
  #ax2.plot(M, np.sqrt(peakPysigmas[:,2]),"--x",mew=3,label=r'$\sigma_{\theta, peak}$, $a=0.998$',color='black',markersize=17,lw=3)
  ax2.errorbar(M, np.sqrt(peakPysigmas[:,2]),fmt='--o',color='black',yerr=4.*np.sqrt(peakPysigmas[:,2])/peakPysigmas[:,4],label=r'$\sigma_{\theta,peak}$, $a=0.998$',markersize=7,lw=3)
  #ax2.errorbar(M, growthPysigmas[:,2],fmt="o",color='blue',lw=1,yerr=4.*growthPysigmas[:,2]/growthPysigmas[:,4])
  #ax2.plot([], [],"-",label=r"$\log(\rho_{peak})=k\cdot G + c$",color='red',lw=1)
  #ax2.plot(M, growthPzsigmas[:,2],"-o",label=r'$\sigma_{\phi, peak}^2$; $a=0.998$',color='red',lw=3)
  #ax2.errorbar(M, growthPzsigmas[:,2],fmt="o",color='red',lw=1,yerr=4.*growthPzsigmas[:,2]/growthPzsigmas[:,4])
  #plt.ylim([0.1,0.7])
  plt.xlim([5.8,18.2])
  #plt.legend(loc="upper left")
  plt.legend(loc="upper right",frameon=False,fontsize=22)
  ax2.set_ylabel(r"Peak Momentum Dispersion $(mc)$")
  #ax1.set_xlabel(r"Growth factor $\left(500 \cdot \frac{10 \cdot M_{seed}}{M_{final}}\right)$")
  ax1.set_xlabel(r"Growth Factor $G$")
  #ax1.set_xlabel(r"Growth factor $\left(M_{seed}\right)$")
  #plt.xlim([45,63])
  plt.savefig(output_path)
  plt.show()
  return annihilationsSpin1
#annihilationMomentumInterplay( thmaxspin1=np.pi/8. )
#annihilationMomentumInterplay( thmaxspin1=0, peaksOnly=False)
annihilations = annihilationMomentumInterplay( thmaxspin1=0, peaksOnly=True)

