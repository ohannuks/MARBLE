import numpy as np
import sys
import paris

def savepdata(fname, data):
  np.savez(fname, **data)
  return

def loadpdata(fname):
  return np.load(fname+".npz")

def loaddata(fname):
  return np.loadtxt(fname,skiprows=6)

def loadattributes(fname):
  with open(fname) as myfile:
      firstNlines=myfile.readlines()[0:6] #put here the interval you want
  return np.array(firstNlines).astype(dtype=np.float64) # Convert to floating points


def readDataFile( fin, calculateAnnihilation = False ):
  ''' Reads in data from file fin

      :param fin: Filename
      :param calculateAnnihilation: Calculates the annihilation rate

      :returns (r,theta,phi,x,y,z,pointData)
  '''
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
  rawdata = np.fromfile(fname,dtype=np.float32)
  counter = int(count)
  r, theta, phi = np.mgrid[rMin:rMax:np.complex(0,rDimensions),thMin:thMax:complex(0,thDimensions),phMin:phMax:complex(0,phDimensions)]
  x = r * np.sin(theta) * np.cos(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(theta)
  counter = counter + 3*nElements
  values = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  v_r = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  v_th = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  v_ph = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vx = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vy = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vz = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vx_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vy_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Vz_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  timeVals = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Px = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Py = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Pz = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Px_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Py_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  Pz_sigma = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
  counter = counter + nElements
  recordDensity = False
  if len(rawdata) > counter:
    recordDensity = True
    densityNormalized = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
    counter = counter + nElements
  recordSamples = False
  if len(rawdata) > counter:
    recordSamples = True
    numberOfSamples = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
    counter = counter + nElements
  else:
    numberOfSamples = values
  recordNumberOfParticles = False
  if len(rawdata) > counter:
    numberOfParticles = rawdata[counter:counter+nElements].reshape( (rDimensions, thDimensions, phDimensions))
    counter = counter + nElements
    recordNumberOfParticles = True
  else:
    numberOfParticles = values
  # Calculate variables
  annihilation = []
  if calculateAnnihilation == True:
    print "shape"
    print np.shape(Px)
    for i in xrange(len(np.ravel(Px))):
      rho = np.ravel(values)[i]
      P = np.array([np.ravel(Px)[i], np.ravel(Py)[i], np.ravel(Pz)[i]])
      dispersion= np.array([np.ravel(Px_sigma)[i], np.ravel(Py_sigma)[i], np.ravel(Pz_sigma)[i]])
      dispersion = np.array(dispersion)
      P = np.array(P)
      print dispersion
      print rho
      crossSection = 1.
      if rho != 0:
        annihilation.append(paris.rateOfAnnihilationCuba( P, 1.*rho, dispersion, crossSection ))
      else:
        annihilation.append(0.)
      #annihilation.append(rho)
    annihilation = np.array(annihilation).reshape((rDimensions,thDimensions,phDimensions))
  # Save to grid
  pointData = {"density" : values, "Vr" : v_r, "Vth" : v_th, "Vph" : v_ph, "Vx" : Vx, "Vy" : Vy, "Vz" : Vz, "sigma_x" : Vx_sigma, "sigma_y" : Vy_sigma, "sigma_z" : Vz_sigma, "dT" : timeVals, "Px" : Px, "Py" : Py, "Pz" : Pz, "sigma_Px" : Px_sigma, "sigma_Py" : Py_sigma, "sigma_Pz" : Pz_sigma, "densityNormalized" : np.array(densityNormalized), "numberOfSamples" : np.array(numberOfSamples), "numberOfParticles" : np.array(numberOfParticles)}
  if calculateAnnihilation == True:
    pointData["Gamma"] = np.array(annihilation)
  return (r,theta,phi,x,y,z,pointData)



def readRadialData(fin, nSamples=False, moments = False):
  ''' Reads in data fin

      :param fin: filename to be read in

      :returns: (R,Density,Px_sigma,Py_sigma,Pz_sigma)
  '''
  r,th,ph,x,y,z,pointData = readDataFile(fin)
  R = np.ravel(r[:,(len(th[0])-1)/2])
  Density = pointData["densityNormalized"][:,(len(th[0])-1)/2]
  Px_sigma= pointData["sigma_Px"][:,(len(th[0])-1)/2]
  Py_sigma= pointData["sigma_Py"][:,(len(th[0])-1)/2]
  Pz_sigma= pointData["sigma_Pz"][:,(len(th[0])-1)/2]
  numberOfSamples = pointData["numberOfSamples"][:,(len(th[0])-1)/2]
  Px= pointData["Px"][:,(len(th[0])-1)/2]
  Py= pointData["Py"][:,(len(th[0])-1)/2]
  Pz= pointData["Pz"][:,(len(th[0])-1)/2]
  numberOfParticles = pointData["numberOfParticles"][:,(len(th[0])-1)/2]
  if nSamples == True:
    if moments == True:
      return (R,Density,Px_sigma,Py_sigma,Pz_sigma,numberOfSamples, Px, Py, Pz, numberOfParticles)
    return (R,Density,Px_sigma,Py_sigma,Pz_sigma,numberOfSamples, numberOfParticles)
  return (R,Density,Px_sigma,Py_sigma,Pz_sigma, numberOfParticles)

def readAveragedRadialData(fin, thmax, nSamples=False, moments=False, normalizedDensity=True):
  ''' Reads in data fin

      :param fin: filename to be read in

      :returns: (R,Density,Px_sigma,Py_sigma,Pz_sigma)
  '''
  if thmax == 0:
      return readRadialData(fin,nSamples, moments)
  r,th,ph,x,y,z,pointData = readDataFile(fin)
  R = []
  Density  = []
  Px_sigma = []
  Py_sigma = []
  Pz_sigma = []
  numberOfSamples=[]
  Px = []
  Py = []
  Pz = []
  numberOfParticles=[]
  # Average out according to some theta:
  for i in np.arange(len(r)):
    densityArray = []
    Px_sigmaArray = []
    Py_sigmaArray = []
    Pz_sigmaArray = []
    numberOfSamplesArray = []
    PxArray = []
    PyArray = []
    PzArray = []
    numberOfParticlesArray = []
    for j in np.arange(len(r[0])):
      if abs(th[i][j][0]-np.pi/2.) < thmax:
        if( normalizedDensity ):
          densityArray.append(pointData["densityNormalized"][i][j][0])
        else:
          densityArray.append(pointData["density"][i][j][0])
        Px_sigmaArray.append(pointData["sigma_Px"][i][j][0])
        Py_sigmaArray.append(pointData["sigma_Py"][i][j][0])
        Pz_sigmaArray.append(pointData["sigma_Pz"][i][j][0])
        numberOfSamplesArray.append(pointData["numberOfSamples"][i][j][0])
        PxArray.append(pointData["Px"][i][j][0])
        PyArray.append(pointData["Py"][i][j][0])
        PzArray.append(pointData["Pz"][i][j][0])
        numberOfParticlesArray.append(pointData["numberOfParticles"][i][j][0])
    Density.append(np.mean(densityArray))
    Px_sigma.append(np.mean(Px_sigmaArray))
    Py_sigma.append(np.mean(Py_sigmaArray))
    Pz_sigma.append(np.mean(Pz_sigmaArray))
    R.append(r[i][j][0])
    numberOfSamples.append(np.sum(numberOfSamplesArray))
    Px.append(np.mean(PxArray))
    Py.append(np.mean(PyArray))
    Pz.append(np.mean(PzArray))
    numberOfParticles.append(np.mean(numberOfParticlesArray))
  if nSamples == True:
    if moments == True:
      print "Woow"
      return (R,Density,Px_sigma,Py_sigma,Pz_sigma, numberOfSamples, Px, Py, Pz, numberOfParticles)
    print "Here"
    return (R,Density,Px_sigma,Py_sigma,Pz_sigma, numberOfSamples, numberOfParticles)
  print "No samples :O "
  return (R,Density,Px_sigma,Py_sigma,Pz_sigma, numberOfParticles)

