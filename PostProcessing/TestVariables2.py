''' Plots density and other variables for testing purposes
    
    Example usage: 
    python ProductionDataVariables.py Data/limits.radial # Note: the Data/limits must be an npz file processed by the processData.py
    
'''
import sys
import numpy as np
import paris as pp
from paris import GR
import pylab as plt

masses = [10,11,12,13,14,15,16,17,18]

for M in masses:
  R,FFIODensity,Density,Px_sigma,Py_sigma,Pz_sigma,nSamples=np.loadtxt("ProductionData/M%02d_r2_spin0000.radial"%M)
  plt.plot(R/500.,FFIODensity,label="M=%d"%M)
plt.legend()
plt.plot()
