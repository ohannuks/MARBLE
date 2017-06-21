''' Creates all the variables r, rho, Px_sigma, Py_sigma, Pz_sigma, nSamples to make everything in a similar format as the readAverageRadialData
    
    Example usage: 
    python ProductionDataVariables.py Data/limits Data/limits.radial # Note: the Data/limits must be an npz file processed by the processData.py
    
'''
import sys
import matplotlib
matplotlib.use("Agg")
import numpy as np
import paris as pp
from paris import GR

fname = sys.argv[1]
fout  = sys.argv[2]


def calculateVariables(fname, average=False):
  # Load the data files
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
  thMin_2= data["thMin_2"]
  thMax_2= data["thMax_2"]
  g_2    = np.ravel(data["g_2"])
  #if GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1) == False:
  #  raise ValueError("Eps_1,Lz_1,Q_1 not bound; something wrong")
  #if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  #  raise ValueError("Eps_2,Lz_2,Q_2 not bound; something wrong")
  # Data loaded; now calculate the 4-current
  N_1 = np.ones(len(Eps_1)) * 1. / Ntotal
  f_1 = N_1 / g_1
  f_2 = f_1
  th = np.pi/2.+0.01
  FFIODensity = []
  ZAMODensity = []
  ZAMODispersion = []
  nSamples    = []
  R = np.linspace(1.*M_2,50*M_2,1000)
  router = 1.*M_2 + np.sqrt(M_2**2-a_2**2)
  if np.sum(rMin_2<=router) != 0:
    raise ValueError("Bad limits; we are inside the inner horizon with rmin router %lf %lf"%(np.min(rMin_2),router))
  for r in R:
    if average == True:
      TH = np.linspace(np.pi*(1./2.-1./7.), np.pi*(1./2.+1./7.), 10)
    else:
      TH = [np.pi/2.+0.01]
    # Average over TH:
    NTH = len(TH)
    samples    = 0.
    dispersion = np.zeros(4)*1.
    FDensity   = 0.
    ZDensity   = 0.
    for i in np.arange(len(TH)):
      th = TH[i]
      samples    = samples + np.sum(((r<rMax_2)&(rMin_2<r)&(th<thMax_2)&(thMin_2<th)))
      # Get the ZAMO Dispersion
      dispersion = dispersion + GR.ZAMODispersion(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r, th)
      FDensity   = FDensity   + GR.FDensity(f_2,Eps_2,Lz_2,K_2,M_2,a_2,rMin_2, rMax_2, thMin_2, thMax_2, r, th)
      ZDensity   = ZDensity   + GR.ZAMODensity(f_2,Eps_2,Lz_2,K_2,M_2,a_2,rMin_2, rMax_2, thMin_2, thMax_2, r, th)
    # Calculate mean:
    samples    = samples    / float(NTH)
    dispersion = dispersion / float(NTH)
    FDensity   = FDensity   / float(NTH)
    ZDensity   = ZDensity   / float(NTH)
    # Insert:
    FFIODensity.append(FDensity)
    # Note that J^0 = \rho_0 U^0 = \rho = (density at rest) in the ZAMO frame
    ZAMODensity.append(ZDensity)
    ZAMODispersion.append(dispersion)
    nSamples.append(samples)
  FFIODensity    = np.array(FFIODensity)
  ZAMODensity    = np.array(ZAMODensity)
  ZAMODispersion = np.array(ZAMODispersion)
  nSamples       = np.array(nSamples)
  return R, FFIODensity, ZAMODensity, ZAMODispersion[:,1], ZAMODispersion[:,2], ZAMODispersion[:,3], nSamples


R, FFIODensity, ZAMODensity, ZAMODispersionx, ZAMODispersiony, ZAMODispersionz, nSamples = calculateVariables(fname)
condition = np.isnan(ZAMODensity)
FFIODensity[np.isnan(FFIODensity)] = 0
ZAMODensity[condition]    =0
ZAMODispersionx[condition]=0
ZAMODispersiony[condition]=0
ZAMODispersionz[condition]=0
nSamples[condition]       =0
np.savetxt(fout, np.array([R, FFIODensity, ZAMODensity, np.sqrt(ZAMODispersionx), np.sqrt(ZAMODispersiony), np.sqrt(ZAMODispersionz), nSamples]))


R, FFIODensity, ZAMODensity, ZAMODispersionx, ZAMODispersiony, ZAMODispersionz, nSamples = calculateVariables(fname,average=True)
condition = np.isnan(ZAMODensity)
FFIODensity[np.isnan(FFIODensity)] = 0
ZAMODensity[condition]    =0
ZAMODispersionx[condition]=0
ZAMODispersiony[condition]=0
ZAMODispersionz[condition]=0
nSamples[condition]       =0
np.savetxt(fout+".averaged", np.array([R, FFIODensity, ZAMODensity, np.sqrt(ZAMODispersionx), np.sqrt(ZAMODispersiony), np.sqrt(ZAMODispersionz), nSamples]))




