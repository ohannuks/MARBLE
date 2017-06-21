import numpy as np
import sys
import paris
import paris as pp
from pyevtk.hl import gridToVTK
def generateVTK(fin, fout):
  ''' Generates a VTK file out of a given file

      :param fin: File name of the input file
      :param fout: File name of the output file
  '''
  # Read in general info
  count=11
  fname=fin
  from pyevtk.hl import gridToVTK
  r,th,ph,x,y,z,pointData =pp.readDataFile(fname)
  fOut=fout
  gridToVTK(fOut, x, y, z, pointData = pointData)




def generate2DVTK( fin, fout, calculateAnnihilation = False ):
  # Read in general info
  count=11
  fname=fin
  attributes = np.fromfile(fname,dtype=np.float32,count=count)
  rMin = attributes[0]
  thMin= attributes[1]
  phMin= attributes[2]
  rMax = attributes[3]
  thMax= attributes[4]
  phMax= attributes[5]
  rDimensions = int(attributes[6])
  thDimensions= int(attributes[7])
  phDimensions= int(attributes[8])
  nElements = int(attributes[9])
  nVariables= attributes[10]
  # The variables are saved in this order:
  #r, th, ph, density
  rawdata = np.fromfile(fname,dtype=np.float32)
  counter = int(count)
  phDimensions2=1
  from pyevtk.hl import gridToVTK
  r,th,ph,x,y,z,pointData =pp.readDataFile(fname)
  #pointData = {"density" : values, "Vr" : v_r, "Vth" : v_th, "Vph" : v_ph, "Vx" : Vx, "Vy" : Vy, "Vz" : Vz, "sigma_x" : Vx_sigma, "sigma_y" : Vy_sigma, "sigma_z" : Vz_sigma, "dT" : timeVals, "Px" : Px, "Py" : Py, "Pz" : Pz, "sigma_Px" : Px_sigma, "sigma_Py" : Py_sigma, "sigma_Pz" : Pz_sigma, "densityNormalized" : np.array(densityNormalized), "numberOfSamples" : np.array(numberOfSamples), "numberOfParticles" : np.array(numberOfParticles)}

  # Project all variables:
  def project( vals ):
    return np.mean(vals,axis=2)
  pointDataProjected = {}
  for i in pointData.iteritems():
    name   = i[0]
    value  = i[1]
    pointDataProjected[name] = project(value).reshape((rDimensions,thDimensions,phDimensions2))
  # TODO: Fix this
  ## Calculate variables
  #annihilation = []
  #if calculateAnnihilation == True:
  #  print "shape"
  #  print np.shape(Px)
  #  for i in xrange(len(np.ravel(Px))):
  #    rho = np.ravel(values)[i]
  #    P = np.array([np.ravel(Px)[i], np.ravel(Py)[i], np.ravel(Pz)[i]])
  #    dispersion= np.array([np.ravel(Px_sigma)[i], np.ravel(Py_sigma)[i], np.ravel(Pz_sigma)[i]])
  #    dispersion = np.array(dispersion)
  #    P = np.array(P)
  #    print dispersion
  #    print rho
  #    crossSection = 1.
  #    #if rho != 0:
  #    #  #annihilation.append(paris.rateOfAnnihilationCuba( P, 1.*rho, dispersion, crossSection ))
  #    #else:
  #    annihilation.append(0.)
  #    #annihilation.append(rho)
  #  annihilation = np.array(annihilation).reshape((rDimensions,thDimensions,phDimensions2))
  fOut=fout
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions2)]
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  # Save to grid
  #if calculateAnnihilation == True:
  #  pointData["Gamma"] = np.array(annihilation)
  gridToVTK(fOut, x, y, z, pointData = pointData)
  return rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions2, pointData

def generate3DFrom2DVTK( fin, fout ):
  rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions, pointData = generate2DVTK(fin, fout)
  phMin = 0
  phMax = 2.*np.pi
  density = pointData["densityNormalized"]
  phDimensions2 = 300
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions2)]
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  density = pointData["densityNormalized"]
  density3d=np.reshape(np.ones(rDimensions*thDimensions*phDimensions2),(rDimensions,thDimensions,phDimensions2))
  for i in np.arange(rDimensions):
    for j in np.arange(thDimensions):
      for k in np.arange(phDimensions2):
        density3d[i,j,k] = density[i,j,0]
  pointData = {"density" : density3d}
  from pyevtk.hl import gridToVTK
  gridToVTK(fout, x, y, z, pointData = pointData)
  return rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions, pointData


import GR

def generate2DVTKNew( fin, fout):
  # Read in data:
  data =pp.loadpdata(fin)
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
  thMin_2 = data["thMin_2"]
  thMax_2 = data["thMax_2"]
  g_2    = np.ravel(data["g_2"])
  N_1 = np.ones(len(Eps_1))*1. / Ntotal
  f_1 = N_1 / g_1
  f_2 = f_1 # Constancy of phase space
  fOut=fout
  rMin= 1.*M_2
  rMax= 25.*M_2
  rDimensions = 200
  thMin= 0
  thMax= np.pi
  thDimensions= 60
  phMin= 0
  phMax= 2.*np.pi
  phDimensions2 = 1
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions2)]
  th = theta
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  FFIODensity = np.zeros(np.shape(r))
  ZAMODensity = np.zeros(np.shape(r))
  sqrtdisp    = np.zeros((np.shape(r)[0],np.shape(r)[1],np.shape(r)[2],4))
  for i in range(len(r)):
    for j in range(len(r[0])):
      for k in range(len(r[0][0])):
        FFIODensity[i][j][k] =  GR.FDensity(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k],th[i][j][k])
        ZAMODensity[i][j][k] = GR.ZAMODensity(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k], th[i][j][k])
        sqrtdisp[i][j][k]    = np.sqrt(GR.ZAMODispersion(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k],th[i][j][k]))
  FFIODensity = np.array(FFIODensity)
  ZAMODensity = np.array(ZAMODensity)
  sqrtdisp    = np.array(sqrtdisp)
  pointData = {"FFIODensity": FFIODensity, "ZAMODensity": ZAMODensity}#, "sigma_Px": sqrtdisp[:,1], "sigma_Py": sqrtdisp[:,2], "sigma_Pz": sqrtdisp[:,3]}
  gridToVTK(fOut, x, y, z, pointData = pointData)
  return rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions2, pointData

def generate3DFrom2DVTKNew( fin, fout ):
  rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions, pointData = generate2DVTKNew(fin, fout)
  phMin = 0
  phMax = 2.*np.pi
  density = pointData["FFIODensity"]
  phDimensions2 = 300
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions2)]
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  density = pointData["FFIODensity"]
  density3d=np.reshape(np.ones(rDimensions*thDimensions*phDimensions2),(rDimensions,thDimensions,phDimensions2))
  for i in np.arange(rDimensions):
    for j in np.arange(thDimensions):
      for k in np.arange(phDimensions2):
        density3d[i,j,k] = density[i,j,0]
  pointData = {"FFIODensity" : density3d}
  from pyevtk.hl import gridToVTK
  gridToVTK(fout, x, y, z, pointData = pointData)
  return rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions, pointData



def generate3DVTKUltra( fin, fout):
  # Read in data:
  data =pp.loadpdata(fin)
  Ntotal = data["Ntotal"]
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
  N_1 = np.ones(len(Eps_1))*1. / Ntotal
  f_1 = N_1 / g_1
  f_2 = f_1 # Constancy of phase space
  fOut=fout
  rMin= 1.*M_2
  rMax= 25.*M_2
  rDimensions = 200
  thMin= 0
  thMax= np.pi
  thDimensions= 30
  phMin= 0
  phMax= 2.*np.pi
  phDimensions2 = 30
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions2)]
  th = theta
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  FFIODensity = np.zeros(np.shape(r))
  ZAMODensity = np.zeros(np.shape(r))
  sqrtdisp    = np.zeros((np.shape(r)[0],np.shape(r)[1],np.shape(r)[2],4))
  Px          = np.zeros(np.shape(r))
  Py          = np.zeros(np.shape(r))
  Pz          = np.zeros(np.shape(r))
  for i in range(len(r)):
    for j in range(len(r[0])):
      for k in range(len(r[0][0])):
        FFIODensity[i][j][k] =  GR.FDensity(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k],th[i][j][k])
        ZAMODensity[i][j][k] = GR.ZAMODensity(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k], th[i][j][k])
        sqrtdisp[i][j][k]    = np.sqrt(GR.ZAMODispersion(f_2, Eps_2, Lz_2, K_2, M_2, a_2, rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k],th[i][j][k]))
        Nmu                  = GR.Nmu(f_2,Eps_2,Lz_2,K_2,M_2,a_2,rMin_2, rMax_2, thMin_2, thMax_2, r[i][j][k],th[i][j][k])
        sinth = np.sin(th[i][j][k])
        costh = np.cos(th[i][j][k])
        sinphi = np.sin(phi[i][j][k])
        cosphi = np.cos(phi[i][j][k])
        Px[i][j][k] = Nmu[1]*sinth*cosphi + Nmu[2]*costh*cosphi + Nmu[3]*(-1.*sinphi)
        Py[i][j][k] = Nmu[1]*sinth*sinphi + Nmu[2]*costh*sinphi + Nmu[3]*cosphi
        Px[i][j][k] = Nmu[1]*costh        + Nmu[2]*(-1.*sinth)  + Nmu[3]*0.
  FFIODensity = np.array(FFIODensity)
  ZAMODensity = np.array(ZAMODensity)
  sqrtdisp    = np.array(sqrtdisp)
  Px          = np.array(Px)
  Py          = np.array(Py)
  Pz          = np.array(Pz)
  pointData = {"FFIODensity": FFIODensity, "ZAMODensity": ZAMODensity, "Px":Px, "Py":Py, "Pz": Pz}#, "sigma_Px": sqrtdisp[:,1], "sigma_Py": sqrtdisp[:,2], "sigma_Pz": sqrtdisp[:,3]}
  gridToVTK(fOut, x, y, z, pointData = pointData)
  return rMin, rMax, thMin, thMax, phMin, phMax, rDimensions, thDimensions, phDimensions2, pointData











