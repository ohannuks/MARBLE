import numpy as np
import sys
from scipy.interpolate import RegularGridInterpolator

def readTrajectory(fin):
  ''' Reads a trajectory file

      :param fin: File name of the input file
      :returns: (t,x,y,z) list
  '''
  data=np.loadtxt(fin)
  t=data[:,0]
  x=data[:,1]
  y=data[:,2]
  z=data[:,3]
  return t,x,y,z


