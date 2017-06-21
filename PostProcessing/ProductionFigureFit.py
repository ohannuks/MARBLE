import sys
import os
sys.path.insert(0, "..\Tools")
from lmfit import models
import numpy as np
import pylab as plt
import paris as pp
import matplotlib
from scipy.interpolate import spline
if not os.path.exists("Figs"): os.mkdir("Figs")
############## Setting ##########################
output_path='./Figs/SkewGaussianFit.pdf'
output_path='../Plots/SkewGaussianFit.pdf'
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

from helperFunctions import getNormalizationFactor, SkewGaussian, SkewGaussian2, calculateError, linefit, skewGaussianFitting

normalizationFactor = getNormalizationFactor()

from scipy.special import erf
from scipy import stats
from scipy.optimize import curve_fit
import sys
growthFactor = False

k0,k1,k2,k3,c0,c1,c2,c3 = skewGaussianFitting(thmaxspin0=np.pi,thmax=np.pi*0.07, growthFactor = growthFactor)
plt.savefig("Plots/ProductionFigureFit.pdf")
plt.show()
print "Done!"


