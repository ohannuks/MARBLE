from multiprocessing import Pool
import numpy as np

def ab(a,b):
  return a+b

a = np.arange(10)
b = np.arange(10)
print ab(a,b)

def ab_pool(args):
  a = args[0]
  b = args[1]
  return a+b

pool = Pool(processes=4)
x = np.array(pool.map(ab_pool, zip(a,b),chunksize=4))
pool.close()
pool.join()

print x
