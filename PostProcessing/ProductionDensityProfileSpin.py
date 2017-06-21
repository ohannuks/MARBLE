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
output_path='./Plots/ProductionDensityProfileSpin.pdf'
#################################################

from helperFunctions import loadProcessedVariables

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



def createRadialPlot(thmax=np.pi/8., thmaxspin0=np.pi, Mass=10):
  plotAll = False
  import numpy as np
  import pylab as plt
  import paris as pp
  import matplotlib
  from scipy.interpolate import spline
  fig, ax1 = plt.subplots(figsize=(16, 9))
  n = 1
  # Get radial plot for spin 0
  R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=loadProcessedVariables(Mass,n,0)
  numberOfSamples = nSamples
  Density = FFIODensity
  print "Rmax: " + str(np.ravel(R)[np.argmax(np.ravel(Density))])
  relativeDensityNormalization = np.max(Density)
  Density = Density / relativeDensityNormalization
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Density, x_smooth)
  y_smooth2=y_smooth
  #plt.plot(x_smooth,y_smooth,'-o',label=r'Density, $a=0$',color='red',lw=3)
  plt.plot(x_smooth,y_smooth,'-',label=r'Density, $a=0$',color='red',lw=3)
  print np.shape(Density)
  print np.shape(numberOfSamples)
  #ax1.errorbar(R,Density,yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples)),fmt='o',color='red')
  ax1.errorbar(R,Density,yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples)),color='red')
  plt.plot([],[],'--',label=r'Momentum Dispersion in $\theta$ direction, $a=0$',color='red',lw=3)
  R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=loadProcessedVariables(Mass,n,998)
  Density = FFIODensity
  numberOfSamples = nSamples
  #plt.axvline(x_smooth[np.argmax(y_smooth)],color='black',ls='--',lw=2)
  # Spin 1

  Density = Density / relativeDensityNormalization
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Density, x_smooth)
  x_smooth=np.append(1,x_smooth)
  y_smooth = np.append(0,y_smooth)
  #ax1.plot(x_smooth,y_smooth,'-o',label=r'Density, $a=0.998$',color='black',lw=3)
  ax1.plot(x_smooth,y_smooth,'-',label=r'Density, $a=0.998$',color='black',lw=3)
  #ax1.errorbar(R,Density,yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples)),fmt='o',color='black')
  ax1.errorbar(R,Density,yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples)),color='black')
  plt.plot([],[],'--',label=r'Momentum Dispersion in $\theta$ direction, $a=0.998$',color='black',lw=3)
  #ax1.set_yticks([0,1,5,10,15,20])

  #ax1.errorbar(R,Density,yerr=1.96*np.sqrt(np.array(Density)/np.array(numberOfSamples)),ls='o')

  ax1.set_ylabel("Relative Density")
  #ax1.set_xlabel(r'Radius $(\frac{GM}{c^2})$')
  ax1.set_xlabel(r'Radius $(GM/c^2)$')
  plt.grid()
  #ax1.set_ylim([0.,20])
  #ax1.annotate(r'a=0.998',xy=(x_smooth[np.argmax(y_smooth)],np.max(y_smooth)))
  #plt.axvline(x_smooth[np.argmax(y_smooth)],color='black',ls='--',lw=2)

  plt.legend(frameon=False,fontsize=22)
  ax2 = ax1.twinx()
  # Get radial plot for spin 0
  import scipy.ndimage.filters as filter
  R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=loadProcessedVariables(Mass,n,0)
  Density =FFIODensity
  Py_sigma = np.sqrt(Py_sigma); Px_sigma = np.sqrt(Px_sigma); Pz_sigma = np.sqrt(Pz_sigma)
  #x_smooth = np.linspace(min(R), max(R), 300)
  #y_smooth = spline(R, Py_sigma, x_smooth)
  x_smooth = R
  y_smooth = np.append(Py_sigma[Py_sigma<0.01],filter.gaussian_filter1d(Py_sigma[Py_sigma>0.01],3))
  ax2.plot(x_smooth,y_smooth,'--',label=r'$\theta$-direction for $a=0$',color='red',lw=3)
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Px_sigma, x_smooth)
  print "a=0.0, Px_sigma max: " + str(max(y_smooth))
  if plotAll == True:
    ax2.plot(x_smooth,y_smooth,'--',label=r'$r$-direction for $a=0$',color='red',lw=3)
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Pz_sigma, x_smooth)
  print "a=0.0, Pz_sigma max: " + str(max(y_smooth))
  if plotAll == True:
    ax2.plot(x_smooth,y_smooth,'--',label=r'$\phi$-direction for $a=0$',color='red',lw=3)
  # Spin 1
  R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=loadProcessedVariables(Mass,n,998)
  Density= FFIODensity
  Py_sigma = np.sqrt(Py_sigma); Px_sigma = np.sqrt(Px_sigma); Pz_sigma = np.sqrt(Pz_sigma)
  #x_smooth = np.linspace(min(R), max(R), 300)
  #y_smooth = spline(R, Py_sigma, x_smooth)
  x_smooth = R
  y_smooth = np.append(Py_sigma[Py_sigma<0.01],filter.gaussian_filter1d(Py_sigma[Py_sigma>0.01],3))
  x_smooth=np.append(1,x_smooth)
  y_smooth = np.append(0,y_smooth)
  ax2.plot(x_smooth,y_smooth,'--',label=r'$\theta$-direction for $a=0.998$',color='black',lw=3)
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Px_sigma, x_smooth)
  print "a=0.998, Px_sigma max: " + str(max(y_smooth))
  #ax2.plot(x_smooth,y_smooth,label=r'$r$-direction for $a=0.998$',color='blue')
  x_smooth = np.linspace(min(R), max(R), 300)
  y_smooth = spline(R, Pz_sigma, x_smooth)
  print "a=0.998, Pz_sigma max: " + str(max(y_smooth))
  if plotAll==True:
    ax2.plot(x_smooth,y_smooth,label=r'$\phi$-direction for $a=0.998$',color='blue')
  ax2.set_xticks(np.arange(0,20,2))
  ax2.set_ylabel(r"Momentum Dispersion $(mc)$")
  #ax2.set_ylim([0,0.9])
  #ax1.set_xlim([0,20])
  #ax2.set_xlim([0,20])
  #plt.legend()
  plt.savefig(output_path,bbox_tight=True)
  plt.show()
createRadialPlot(thmax=np.pi*0.07, thmaxspin0=np.pi, Mass=10)
#createRadialPlot(thmax=0, thmaxspin0=np.pi, Mass=10)
#createRadialPlot(thmax=np.pi/8., thmaxspin0=np.pi, Mass=14)
plt.show()

