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


def densityMomentumInterplay( thmaxspin1=np.pi/8., peaksOnly=False ):
    from lmfit import models
    # Plot masses
    # Read mass files
    spin0={}
    spin1={}
    masses = [6,7,9,10,11,12,13,14,15,16,17,18]
    for M in masses:
        R,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=pp.readAveragedRadialData("M%02d_r2_spin0000.dat"%M,thmax=np.pi,nSamples=True)
        spin0[M]={"M":M,"R":R,"Density":Density,"Px_sigma":Px_sigma,"Py_sigma":Py_sigma,"Pz_sigma":Pz_sigma,"nSamples":nSamples}
        R,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=pp.readAveragedRadialData("M%02d_r2_spin0998.dat"%M,thmax=thmaxspin1,nSamples=True)
        spin1[M]={"M":M,"R":R,"Density":Density,"Px_sigma":Px_sigma,"Py_sigma":Py_sigma,"Pz_sigma":Pz_sigma,"nSamples":nSamples}

    fig, ax1 = plt.subplots(figsize=(16, 9))

    growthDensities=[]
    growthPxsigmas =[]
    growthPysigmas =[]
    growthPzsigmas =[]
    peakPxsigmas =[]
    peakPysigmas =[]
    peakPzsigmas =[]
    for M in masses:
        # Get radial plot for spin 0, M10
        R = spin0[M]["R"]
        maxdensityspin0 = np.max(spin0[M]["Density"])
        densitysamplesspin0=spin1[M]["nSamples"][np.argmax(spin0[M]["Density"])]
        model = models.SkewedGaussianModel()
        # set initial parameter values
        params = model.make_params(amplitude=0.1, center=10, sigma=1, gamma=0)
        result = model.fit(np.ravel(spin1[M]["Density"]), params, x=np.ravel(spin1[M]["R"]))
        maxdensityspin1 = np.max(spin1[M]["Density"])
        densitysamplesspinfit = np.max(result.best_fit)
        densitysamplesspin1=spin1[M]["nSamples"][np.argmax(spin1[M]["Density"])]
        maxsigmaPxspin0 = np.max(spin0[M]["Px_sigma"])
        sigmaPxsamplesspin0=spin0[M]["nSamples"][np.argmax(spin0[M]["Px_sigma"])]
        maxsigmaPyspin0 = np.max(spin0[M]["Py_sigma"])
        sigmaPysamplesspin0=spin0[M]["nSamples"][np.argmax(spin0[M]["Py_sigma"])]
        maxsigmaPzspin0 = np.max(spin0[M]["Pz_sigma"])
        sigmaPzsamplesspin0=spin0[M]["nSamples"][np.argmax(spin0[M]["Pz_sigma"])]
        maxsigmaPxspin1 = np.max(spin1[M]["Px_sigma"])
        sigmaPxsamplesspin1=spin1[M]["nSamples"][np.argmax(spin1[M]["Px_sigma"])]
        maxsigmaPyspin1 = np.max(spin1[M]["Py_sigma"])
        sigmaPysamplesspin1=spin1[M]["nSamples"][np.argmax(spin1[M]["Py_sigma"])]
        maxsigmaPzspin1 = np.max(spin1[M]["Pz_sigma"])
        sigmaPzsamplesspin1=spin1[M]["nSamples"][np.argmax(spin1[M]["Pz_sigma"])]
        growthDensities.append([M,maxdensityspin0, maxdensityspin1, densitysamplesspin0, densitysamplesspin1, densitysamplesspinfit])
        growthPxsigmas.append( [M,maxsigmaPxspin0, maxsigmaPxspin1, sigmaPxsamplesspin0, sigmaPxsamplesspin1])
        growthPysigmas.append( [M,maxsigmaPyspin0, maxsigmaPyspin1, sigmaPysamplesspin0, sigmaPysamplesspin1])
        growthPzsigmas.append( [M,maxsigmaPzspin0, maxsigmaPzspin1, sigmaPzsamplesspin0, sigmaPzsamplesspin1])
        def getPeakSigma( variable ):
            maxDensityIndex = np.argmax(spin0[M]["Density"])
            peaksigmaPxspin0 = spin0[M][variable][maxDensityIndex]
            peaksigmaPxsamplesspin0=densitysamplesspin0
            maxDensityIndex = np.argmax(spin1[M]["Density"])
            peaksigmaPxspin1 = spin1[M][variable][maxDensityIndex]
            peaksigmaPxsamplesspin1=densitysamplesspin1
            return (peaksigmaPxspin0, peaksigmaPxsamplesspin0, peaksigmaPxspin1, peaksigmaPxsamplesspin1)
        peakSigmaPxspin0, samplesSpin0, peaksigmaPxspin1, samplesSpin1 = getPeakSigma( "Px_sigma" )
        peakPxsigmas.append([M, peakSigmaPxspin0, peaksigmaPxspin1, samplesSpin0, samplesSpin1])
        peakSigmaPxspin0, samplesSpin0, peaksigmaPxspin1, samplesSpin1 = getPeakSigma( "Py_sigma" )
        peakPysigmas.append([M, peakSigmaPxspin0, peaksigmaPxspin1, samplesSpin0, samplesSpin1])
        peakSigmaPxspin0, samplesSpin0, peaksigmaPxspin1, samplesSpin1 = getPeakSigma( "Pz_sigma" )
        peakPzsigmas.append([M, peakSigmaPxspin0, peaksigmaPxspin1, samplesSpin0, samplesSpin1])

    growthDensities = np.array(growthDensities)
    growthPxsigmas  = np.array(growthPxsigmas)
    growthPysigmas  = np.array(growthPysigmas)
    growthPzsigmas  = np.array(growthPzsigmas)

    peakPxsigmas  = np.array(peakPxsigmas)
    peakPysigmas  = np.array(peakPysigmas)
    peakPzsigmas  = np.array(peakPzsigmas)

    M=growthDensities[:,0]
    #M=500./growthDensities[:,0]
    global normalizationFactor
    normalizationFactor=np.min(growthDensities[:,1])
    densitiesSpin0=growthDensities[:,1]/np.min(growthDensities[:,1])
    densitiesSpin1=growthDensities[:,2]/np.min(growthDensities[:,1])
    densitiesSpinfit=growthDensities[:,5]/np.min(growthDensities[:,1])
    z = np.polyfit(M, np.log10(densitiesSpin0), 1)
    print z
    p = np.poly1d(z)
    #ax1.semilogy(M, 10**(p(M)),"-",label=r"$\log(\rho_{peak})=k\cdot G + c$",color='black',lw=1)
    #ax1.semilogy(M, densitiesSpin0,"-*",label=r'$\rho_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
    ax1.semilogy(M, densitiesSpin0,"-x",mew=3,label=r'$\rho_{peak}$; $a=0.0$',color='red',lw=3,markersize=17)
    z = np.polyfit(M, np.log10(densitiesSpin1), 1)
    print z
    p = np.poly1d(z)
    ax1.semilogy(M, densitiesSpin1,"-*",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3,markersize=17)
    #ax1.semilogy(M, 10**(p(M)),"-",label=r"$\log(\rho_{peak})=k\cdot G + c$ best fit",color='black',lw=1)
    #ax1.semilogy(M, densitiesSpinfit,"--o",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3)

    #plt.ylim([1,4e3])
    #ax1.plot(M, densitiesSpin0/np.min(densitiesSpin0),"--o",label=r'$\rho_{peak}$; $a=0.0$',color='black',lw=1)
    #ax1.plot(M, densitiesSpin1/np.min(densitiesSpin0),"--o",label=r"$\rho_{peak}$; $a=0.998$",color='black',lw=3)
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
    plt.ylim([0.1,0.7])
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
#densityMomentumInterplay( thmaxspin1=np.pi/8. )
#densityMomentumInterplay( thmaxspin1=0, peaksOnly=False)
densityMomentumInterplay( thmaxspin1=0, peaksOnly=True)

