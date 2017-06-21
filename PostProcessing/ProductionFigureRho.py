import sys
import os
import numpy as np
import pylab as plt
import paris as pp
import matplotlib
from scipy.interpolate import spline
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

def calculateError(r, y, nSamples):
  yerr=1.96*np.ravel(y)*np.sqrt(1./np.ravel(nSamples));yerr[np.isnan(yerr)] = np.mean(yerr[np.isnan(yerr) == False])
  return yerr
from helperFunctions import skewGaussianFitting, SkewGaussian, SkewGaussian2, getNormalizationFactor

normalizationFactor = getNormalizationFactor()
from scipy.special import erf
from scipy import stats
k0,k1,k2,k3,c0,c1,c2,c3 = skewGaussianFitting(thmaxspin0=np.pi, growthFactor = False)
x0=np.array([k0,k1,k2,k3,c0,c1,c2,c3])
plt.close()

def plotGaussianWithError(r, density, Px_sigma, Py_sigma, Pz_sigma, nSamples):
  yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples))
  plt.errorbar(r,density,yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(nSamples)))
  return

def bestfit(r, density, Px_sigma, Py_sigma, Pz_sigma, nSamples):
  yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples))

def tryGaussianFitting(x,thmax=np.pi/8.,mass=10,var="density",showResidual=False):
    if len(x)==9:
      k0,k1,k2,k3,c0,c1,c2,c3,k4=tuple(x)
    else:
      k0,k1,k2,k3,c0,c1,c2,c3=tuple(x)
      k4=0
    global A, mu, sigma, gamma
    A = lambda G: 10**(k0*G+c0)
    mu= lambda G: np.log(k1*G+c1)
    sigma= lambda G: k2*G+c2
    gamma= lambda G: k3*G+c3
    G = mass
    n = 2
    r, FFIODensity, density, Px_sigma, Py_sigma, Pz_sigma, nSamples=np.loadtxt("ProductionData/M%d_r%d_spin0998.radial"%(mass,n))
    r = np.ravel(r)#/ 500.
    logr = np.log(np.ravel(r))
    if var=="density":
        variable = np.ravel(FFIODensity)/normalizationFactor
    elif var=="Px_sigma":
        variable = np.ravel(Px_sigma)
    elif var=="Py_sigma":
        variable = np.ravel(Py_sigma)
    elif var=="Pz_sigma":
        variable = np.ravel(Pz_sigma)
    else:
        print "Invalid variable!"
    Pz_sigma = np.ravel(Pz_sigma)
    Py_sigma = np.ravel(Py_sigma)
    Px_sigma = np.ravel(Px_sigma)
    xvals,yvals = logr,variable
    Density = variable
    numberOfSamples=np.ravel(nSamples)
    yerr=1.96*np.ravel(Density)*np.sqrt(1./np.ravel(numberOfSamples))
    yerr[np.isnan(yerr)]=0
    plt.figure(figsize=(16, 9))
    plt.errorbar(np.exp(xvals)/500., yvals,yerr=yerr,fmt='o',color='blue',label='Relative Density (simulation data)',markersize=3)
    plt.ylim([0,np.max(yvals)])
    plt.xlim([np.min(np.exp(xvals))/500.,np.max(np.exp(xvals))/500.])
    print SkewGaussian(logr,G,A, mu, sigma, gamma)
    plt.plot(np.exp(xvals)/500., SkewGaussian(logr,G,A, mu, sigma, gamma),'-',color='black',lw=3,label='Relative Density (fit)')
    plt.grid()
    plt.ylabel("Relative Density")
    #plt.ylim([0,1.2*np.max(yvals)])
    plt.xticks(np.arange(1,np.max(np.exp(xvals)/500.)+0.1,2))
    #plt.xticks([])
    fit = SkewGaussian(logr,G,A, mu, sigma, gamma)
    #plot arrows and letters here
    #plt.annotate("",
    #        xy=(1, np.max(fit)), xycoords='data',
    #        xytext=(1, 0), textcoords='data',
    #        arrowprops=dict(arrowstyle="<->",
    #                        connectionstyle="arc3",lw=3),
    #        )
    #plt.text(1.1,np.max(fit), r'A',fontsize=25)
    #
    #plt.annotate("",
    #        xy=(3.5, np.max(fit)*0.5), xycoords='data',
    #        xytext=(11.1, np.max(fit)*0.5), textcoords='data',
    #        arrowprops=dict(arrowstyle="<->",
    #                        connectionstyle="arc3",lw=3),
    #        )
    #
    #plt.text(7,np.max(fit)*0.52, r'$\sigma$',fontsize=35)
    #plt.annotate("",
    #        xy=(r[np.argmax(fit)], 0), xycoords='data',
    #        xytext=(r[np.argmax(fit)], np.max(fit)), textcoords='data',
    #        arrowprops=dict(arrowstyle="->",
    #                        connectionstyle="arc3",lw=3),
    #        )
    #plt.text(r[np.argmax(fit)],np.max(fit)*1.02, r'$\mu$',fontsize=35)
    #plt.annotate("",
    #        xy=(r[np.argmax(fit)]*1.1, np.max(fit)), xycoords='data',
    #        xytext=(r[np.argmax(fit)]*1.02+5, np.max(fit)), textcoords='data',
    #        arrowprops=dict(arrowstyle="->",
    #                        connectionstyle="arc3",lw=3),
    #        )
    #plt.text(10,np.max(fit)*1.02, r'$\gamma$',fontsize=35)
    plt.legend(frameon=False,fontsize=22)
    plt.ylim([0,np.max(yvals)])
#    if showResidual == True:
#      plt.subplot(313)
#      plt.grid()
#      w = 1./(2.*yerr)
#      check0 = (np.isnan(w)==False) & (w<1e30)
#      w = w[check0]
#      fit=SkewGaussian(logr,G,A, mu, sigma, gamma)[check0]
#      Residual=np.abs((yvals[check0]-SkewGaussian(logr,G,A, mu, sigma, gamma)[check0]))
#      plt.plot(np.exp(xvals)[check0], Residual,'-',color='black',lw=3,label='Residual')
#      print w
#      plt.ylabel("Residual")
#      observed=yvals[20:]
#      expected=SkewGaussian(logr,G,A, mu, sigma, gamma)[20:]
#      chi_squared_stat = (((observed-expected)**2)/expected).sum().sum()
#      print(chi_squared_stat)
#      crit = stats.chi2.ppf(q = 0.6, # Find the critical value for 95% confidence*
#                            df = len(observed)-1)   # *
#      #plt.ylim([-20,20])
#      #plt.yticks([-20,-10,0,10,20])
#      #plt.xticks(np.arange(1,30))
#      #plt.xlim([0,25])
    plt.xlabel(r'Radius $(GM/c^2)$')
    #plt.show()
tryGaussianFitting(mass=10,var="density",thmax=np.pi*0.07,showResidual=False, x=x0)
plt.savefig("Plots/ProductionFigureRho.pdf");
plt.show()


