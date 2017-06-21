import sys
import matplotlib
matplotlib.use("Agg")
import paris as pp
from paris import GR
from multiprocessing import Pool
import numpy as np
print "Usage: python $< Data/degeneracydistribution.dat Data/degeneracydistribution"
fname = sys.argv[1]

print "Reading data.."

# Get the first 3 lines of the data
attr = pp.loadattributes(fname)
Ntotal            = attr[0]
Nin               = attr[1]
a                 = attr[2]
M                 = attr[3]
a_2               = attr[4]
M_2               = attr[5]

print "Ntotal: " + str(Ntotal)

