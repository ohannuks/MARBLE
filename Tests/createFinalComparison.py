# Plots/finalcomparisonbetweengyotoandanalytical.pdf: createFinalComparison.py Data/degeneracygrid.dat Data/degeneracydistribution.npz
#        python $< Data/degeneracydistribution Data/degeneracygrid.dat $@
import matplotlib
#font = {'weight' : 'normal',
#        'size'   : 22}        # Set font
#matplotlib.rc('font', **font) # Set font
#matplotlib.rcParams['mathtext.fontset'] = 'custom'
#matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
#matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
#matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rc('text', usetex=True)

import paris as pp
import numpy as np
import pylab as plt
import sys
from paris import GR

fname1    = sys.argv[1]
fname2    = sys.argv[2]
fout      = sys.argv[3]

data =pp.loadpdata(fname2)
Eps_1  = data["Eps_1"]
Lz_1   = data["Lz_1"]
K_1    = data["K_1"]
M_1    = data["M_1"]
a_1    = data["a_1"]
rMin_1 = data["rMin_1"]
rMax_1 = data["rMax_1"]
g_1    = data["g_1"]

Eps_2  = data["Eps_2"]
Lz_2   = data["Lz_2"]
K_2    = data["K_2"]
M_2    = data["M_2"]
a_2    = data["a_2"]
rMin_2 = data["rMin_2"]
rMax_2 = data["rMax_2"]
thMin_2 = data["thMin_2"]
thMax_2 = data["thMax_2"]
g_2    = data["g_2"]


# Get the density via Gyoto
##############################################
r,Density1,Px_sigma,Py_sigma,Pz_sigma,nSamples,nParticles=pp.readAveragedRadialData(fname1,0,True, normalizedDensity=False)
r = np.array(r)
nParticles = np.ravel(np.array(nParticles))
# Define proper volume in the FFIO frame
#  gtt    = -1. + 2./R
#  detgBL = -1.*R**4
#  nParticles = variable
#  propernParticles = nParticles * -1.*gtt
#  properVolume     = -1.*np.sqrt(-1.*detgBL) / gtt
#  Density          = propernParticles / properVolume
#TODO: Make sure this is correct
M = M_2/M_2
a = a_2/M_2
th=np.pi/2.+0.01
detgBL = GR.detg(M,a,r,th)
pt     = GR.u(Eps=1., Lz=0., K=a**2, M=M, a=a, r=r, th=th)[0]
propernParticles = nParticles #* -1.*gtt
properVolume     = -1.*np.sqrt(-1.*detgBL) * pt
Density1         = propernParticles / properVolume

# Give the density out
rGyoto, DensityGyoto = r, Density1
nSamplesGyoto = np.ravel(nSamples)
##############################################

# Calculate density from the data file
##############################################

print "a is " + str(a_2) + " " + str(data["a"])

if GR.tests.isBound(Eps_1,Lz_1,K_1,M_1,a_1) == False:
  raise ValueError("Eps_1,Lz_1,Q_1 not bound; something wrong")
if GR.tests.isBound(Eps_2,Lz_2,K_2,M_2,a_2) == False:
  raise ValueError("Eps_2,Lz_2,Q_2 not bound; something wrong")


# Calculate the density using weighted f
r   = rGyoto * M_2
f_1 = 1. / g_1
f_2 = f_1 # Constancy of phase space
if np.sum(f_2 < 0) != 0:
  raise ValueError("Bad f_2!")
Density = []
nSamples=[]
for ri in r:
  f_2 = np.ravel(f_2)
  Density.append(GR.FDensity(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, ri,th))
  nSamples.append(np.sum( (ri<rMax_2)&(rMin_2<ri)&(thMin_2<th)&(th<thMax_2) ))
Density = np.array(Density)
nSamples= np.array(nSamples)

# Convert density into r_g units:
rCalculated, DensityCalculated, errCalculated = r/M_2, Density*M_2**3, Density*M_2**3/np.sqrt(nSamples)
##############################################

# Normalize:
DensityGyoto = DensityGyoto*DensityCalculated[len(DensityCalculated)-1]/DensityGyoto[len(DensityGyoto)-1] 

#DensityGyoto = DensityGyoto / np.max(DensityGyoto)
#DensityCalculated= DensityCalculated / np.max(DensityCalculated)

plt.figure(figsize=(16,9))
plt.errorbar(rGyoto,DensityGyoto,yerr=2.*DensityGyoto/np.sqrt(nSamplesGyoto),fmt='o',label='Sampled with Gyoto',lw=3,color='blue')
plt.errorbar(rCalculated,DensityCalculated,yerr=2.*DensityCalculated/np.sqrt(nSamples),fmt='o',label='Sampled Directly from Phase Space',color='black')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r"Radius ($GM/c^2$)")
plt.ylabel(r"Relative Density")
plt.tight_layout()
plt.savefig(fout,bbox_tight=True)
plt.show()



