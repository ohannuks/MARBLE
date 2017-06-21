# Filepaths for python
############
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/" + "FileGeneration")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/" + "Read")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/" + "Analysis")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/" + "Development")
############

from generateVTK import generateVTK, generate2DVTK, generate3DFrom2DVTK, generate2DVTKNew, generate3DFrom2DVTKNew
from plotTrajectory import plotTrajectory
from cuba import calculateAnnihilationPython
from ReadData import readDataFile, readRadialData, readAveragedRadialData, loadattributes, loaddata, loadpdata, savepdata
from peaks import peak
import GR
