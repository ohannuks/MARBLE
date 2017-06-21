import numpy as np

def peak( attributes ):
  ''' Calculates the peaks of densities for given attributes
  
      :returns: {"M":M, "Density":Density, "Px_sigma":Px_sigma, "Py_sigma":Py_sigma, "Pz_sigma":Pz_sigma, "DensityErr":DensityError, "Px_sigmaErr":Px_sigmaError, "Py_sigmaErr":Py_sigmaError, "Pz_sigmaErr":Pz_sigmaError}

      :note: Assumes attributes[M]={"M":M,"R":R,"Density":Density,"Px_sigma":Px_sigma,"Py_sigma":Py_sigma,"Pz_sigma":Pz_sigma,"nSamples":nSamples}

      :note: Example: 

      peakValues = peak( spin0 )
      print peakValues["Density"]
  '''
  densities=[];masses=[];Px_sigmas=[];Py_sigmas=[];Pz_sigmas=[]
  densitiesErr=[]; Px_sigmasErr=[]; Py_sigmasErr=[]; Pz_sigmasErr=[]
  # Loop through all possible values
  for M,v in sorted(attributes.items()):
    maxDensity=np.max(attributes[M]["Density"])
    densities.append(maxDensity)
    maxSigmaR =np.max(attributes[M]["Px_sigma"])
    Px_sigmas.append(maxSigmaR)
    maxSigmaTh=np.max(attributes[M]["Py_sigma"])
    Py_sigmas.append(maxSigmaTh)
    maxSigmaPh=np.max(attributes[M]["Pz_sigma"])
    Pz_sigmas.append(maxSigmaPh)
    samples=attributes[M]["nSamples"][np.argmax(attributes[M]["Density"])]
    maxDensityErr=1.96*maxDensity/np.sqrt(samples)
    densitiesErr.append(maxDensityErr)
    samples=attributes[M]["nSamples"][np.argmax(attributes[M]["Px_sigma"])]
    maxSigmaRErr=4.*maxSigmaR/samples
    Px_sigmasErr.append(maxSigmaRErr)
    samples=attributes[M]["nSamples"][np.argmax(attributes[M]["Py_sigma"])]
    maxSigmaThErr=4.*maxSigmaTh/samples
    Py_sigmasErr.append(maxSigmaThErr)
    samples=attributes[M]["nSamples"][np.argmax(attributes[M]["Pz_sigma"])]
    maxSigmaPhErr=4.*maxSigmaPh/samples
    Pz_sigmas.append(maxSigmaPhErr)
    masses.append(M)
  return {"M":masses, "Density":densities, "Px_sigma":Px_sigmas, "Py_sigma":Py_sigmas, "Pz_sigma":Pz_sigmas, "DensityErr":densitiesErr, "Px_sigmaErr":Px_sigmasErr, "Py_sigmaErr":Py_sigmasErr, "Pz_sigmaErr":Pz_sigmasErr}
