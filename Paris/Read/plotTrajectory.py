# Get data
import numpy as np

def plotTrajectory( fname, spin=1, fig='', plothorizon=True, color='black' ):
  ''' Plots the trajectory of a given xyz file

  '''
  x = []
  y = []
  z = []
  data = np.loadtxt(fname)
  x.append(data[:,1])
  y.append(data[:,2])
  z.append(data[:,3])
  x = np.array(x)
  y = np.array(y)
  z = np.array(z)
  import matplotlib.pyplot as plt
  # Run the matploblib script :)
  import matplotlib as mpl
  from mpl_toolkits.mplot3d import Axes3D
  import matplotlib.pyplot as plt
  
  mpl.rcParams['legend.fontsize'] = 10
  
  if fig == '':
    fig = plt.figure()
  for i in range(0,len(x)):
    ax = fig.gca(projection='3d')
    ax.plot(x[i], y[i], z[i],'-', color=color)
  
  
  # Draw horizon:
  ######################################################
  if plothorizon == True:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    # Correct radius:
    Rs=2.
    r_inner = (Rs+np.sqrt(Rs**2-4.*(spin)**2))/2.
    r_outer = Rs
    
    # Inner:
    #######################
    coefs = (1./r_inner, 1./r_inner, 1./r_inner)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1 
    # Radii corresponding to the coefficients:
    rx, ry, rz = 1/np.sqrt(coefs)
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(v))
    y = ry * np.outer(np.sin(u), np.sin(v))
    z = rz * np.outer(np.ones_like(u), np.cos(v))
    # Plot:
    #ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r', alpha=0.1)
    ax.plot_surface(x, y, z,  rstride=6, cstride=6, color='red', edgecolors='red', alpha=0.5)
    # Adjustment of the axes, so that they all have the same span:
    max_radius = max(rx, ry, rz)
    for axis in 'xyz':
        getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
    if spin == 1:
      #######################
      # Outer:
      #######################
      coefs = (1./r_outer, 1./r_outer, 1./r_inner)  # Coefficients in a0/c x**2 + a1/c y**2 + a2/c z**2 = 1 
      # Radii corresponding to the coefficients:
      rx, ry, rz = 1/np.sqrt(coefs)
      # Set of all spherical angles:
      u = np.linspace(0, 2 * np.pi, 100)
      v = np.linspace(0, np.pi, 100)
      # Cartesian coordinates that correspond to the spherical angles:
      # (this is the equation of an ellipsoid):
      x = rx * np.outer(np.cos(u), np.sin(v))
      y = ry * np.outer(np.sin(u), np.sin(v))
      z = rz * np.outer(np.ones_like(u), np.cos(v))
      # Plot:
      #ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r', alpha=0.1)
      ax.plot_surface(x, y, z,  rstride=6, cstride=6, color='b', edgecolors='b', alpha=0.1)
      # Adjustment of the axes, so that they all have the same span:
      max_radius = max(rx, ry, rz)
      for axis in 'xyz':
          getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
      #######################
  return fig
  
  

