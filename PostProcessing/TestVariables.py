''' Plots density and other variables for testing purposes
    
    Example usage: 
    python ProductionDataVariables.py Data/limits.radial # Note: the Data/limits must be an npz file processed by the processData.py
    
'''
import sys
import numpy as np
import paris as pp
from paris import GR
import pylab as plt

fname = sys.argv[1]

R, FFIODensity, ZAMODensity, ZAMODispersionx, ZAMODispersiony, ZAMODispersionz, nSamples = np.loadtxt(fname)

plt.plot(R/500.,ZAMODensity,label='ZAMO')
print FFIODensity
plt.plot(R/500.,FFIODensity,label='FFIO')
plt.legend()
plt.show()
plt.plot(R/500.,ZAMODispersionx)
plt.show()
plt.plot(R/500.,ZAMODispersiony)
plt.show()
plt.plot(R/500.,ZAMODispersionz)
plt.show()




