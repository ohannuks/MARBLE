# coding: utf-8
import matplotlib
font = {'weight' : 'normal',
        'size'   : 22}        # Set font
matplotlib.rc('font', **font) # Set font
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)
import paris as pp; import numpy as np; import pylab as plt
from helperFunctions import getNormalizationFactor, gridAnnihilation, gridDensity

def Csc(x):
  return 2.*np.sin(x)/(1.-np.cos(2.*x))

from paris import GR

global N
N=2
print "ProductionFigureGrowthOfAnnihilation.py False log"
masses = [10,11,12,13,14,15,16,17,18,19,20]
def findPeak(spin, calculateDensitySquared=False):
  ''' Calculates the peak density etc 
      :returns: peaks

      Note: Peaks is a 4-dimensional array with each entry being [density, Py_sigma, index, samples]
  '''
  from paris import GR
  # Go through all masses
  peaks = []
  peaks2= []
  for M in masses:
    n=N
    if spin==0:
      fname = "ProductionData/M%d_r%d_spin0000"%(M,n)
    else:
      fname = "ProductionData/M%d_r%d_spin0998"%(M,n)
    normalization = getNormalizationFactor()
    # Get all the annihilation rates:
    if calculateDensitySquared == False:
      r, th, annihilationRates, sqrtdispersions, nSamples, M_2, a_2 = gridAnnihilation(fname)
    else:
      r, th, annihilationRates, sqrtdispersions, nSamples, M_2, a_2 = gridDensity(fname)
    # Get the peak annihilation:
    i = np.argmax(annihilationRates)
    peaks.append( [annihilationRates[i]*GR.dtaudt(r[i], th[i], a_2, M_2), sqrtdispersions[i][1], sqrtdispersions[i][2], sqrtdispersions[i][3], i, nSamples[i]] )
    print peaks
  peaks = np.array(peaks)
  print "Peaks:"
  print peaks
  return masses, peaks


# Plot:
def dualPlot(peaks0, peaks1, masses, growthFactor=True, annihilation=False):
  fig, ax1 = plt.subplots(figsize=(16, 9))
  # Set growth factor
  if growthFactor == True:
    G = 1./np.array(masses)
  else:
    G = np.array(masses)
  if annihilation == False:
    var="\rho"
  else:
    var="R"
  if annihilation==False:
    if sys.argv[2] == "log":
      ax1.semilogy(G, peaks0[:,0],"-x",mew=3,label=r'$\rho_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
      ax1.semilogy(G, peaks1[:,0],"-*",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
    else:
      ax1.plot(G, peaks0[:,0],"-x",mew=3,label=r'$\rho_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
      ax1.plot(G, peaks1[:,0],"-*",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
    ax1.set_ylabel(r"Peak Relative Density")
  else:
    if sys.argv[1] == "log":
      ax1.semilogy(G, peaks0[:,0],"-x",mew=3,label=r'$\Gamma_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
      ax1.semilogy(G, peaks1[:,0],"-*",label=r"$\Gamma_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
    else:
      ax1.plot(G, peaks0[:,0],"-x",mew=3,label=r'$\Gamma_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
      ax1.plot(G, peaks1[:,0],"-*",label=r"$\Gamma_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
    ax1.set_ylabel(r"Peak Relative Annihilation Rate")
  plt.grid()
  ax2 = ax1.twinx()
  ax2.plot([],[],"-x",mew=3,label=r'$'+var+'_{peak}$, $a=0.0$',color='red',lw=3,markersize=17)
  ax2.plot([],[],"-*",label=r"$"+var+"_{peak}$, $a=0.998$",color='black',lw=3,markersize=17)
  ax2.errorbar(G, peaks0[:,2],fmt='--o',color='red',yerr=2.*peaks0[:,2]/peaks0[:,5],label=r'$\Delta P^\theta_{\rm peak}$, $a=0.0$',markersize=7,lw=3)
  ax2.errorbar(G, peaks1[:,2],fmt='--o',color='black',yerr=2.*peaks1[:,2]/peaks1[:,5],label=r'$\Delta P^\theta_{\rm peak}$, $a=0.998$',markersize=7,lw=3)
  #ax2.errorbar(G, peaks0[:,1],fmt='-.o',color='red',yerr=2.*peaks0[:,1]/peaks0[:,5],label=r'$\Delta P^r_{\rm peak}$, $a=0.0$',markersize=7,lw=3)
  #ax2.errorbar(G, peaks1[:,1],fmt='-.o',color='black',yerr=2.*peaks1[:,1]/peaks1[:,5],label=r'$\Delta P^r_{\rm peak}$, $a=0.998$',markersize=7,lw=3)
  #ax2.errorbar(masses, peaks0[:,3],fmt='-.o',color='blue',yerr=2.*peaks0[:,3]/peaks0[:,5],label=r'$\sigma_{r,peak}$, $a=0.0$',markersize=7,lw=3)
  #ax2.errorbar(masses, peaks1[:,3],fmt='-.o',color='green',yerr=2.*peaks1[:,3]/peaks1[:,5],label=r'$\sigma_{r,peak}$, $a=0.998$',markersize=7,lw=3)
  ax2.set_ylabel(r"Momentum Dispersion (mc)")
  ax2.set_xlabel("Growth factor")
  if growthFactor == True:
    ax1.set_xlabel(r"Growth factor $G$")
  else:
    ax1.set_xlabel(r"$5000 R_{\rm IF}$")
  return ax1, ax2


import pathlib
if pathlib.Path('./dttest0').is_file():
  peaks0=np.loadtxt("dttest0");
else:
  masses, peaks0 = findPeak(0,False)
if pathlib.Path('./dttest1').is_file():
  peaks1=np.loadtxt("dttest1");
else:
  masses, peaks1 = findPeak(0.998,False)

normalization = np.min(peaks0[:,0])
peaks0[:,0] = peaks0[:,0]
peaks1[:,0] = peaks1[:,0]

np.savetxt("dttest0",peaks0)
np.savetxt("dttest1",peaks1)


import sys
if len(sys.argv) <= 1 or sys.argv[1] == "false":
  growthFactor = False
else:
  growthFactor = True

annihilation=True
ax1,ax2 = dualPlot(peaks0,peaks1,masses, growthFactor=growthFactor, annihilation=annihilation )
if growthFactor == True:
  plt.legend(frameon=False,loc='upper left')
  #ax1.set_ylim([1,6*10**8])
  #ax1.set_xlim([0.054,0.17])
else:
  plt.legend(frameon=False)

if annihilation == True:
  print "Annihilation set to true"
  #ax1.set_ylim([1e-4,1e13])

#plt.ylim([0,0.6])
#plt.xlim([5,19])
plt.savefig("./Plots/ProductionFigureGrowthOfAnnihilation.pdf")
plt.show()

