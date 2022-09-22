#!/home/vsriv/miniconda3/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

optim_iter = int(sys.argv[1])
A = np.loadtxt("%04d/gamma_max_params.dat" % (optim_iter))
A = np.reshape(A, [31, 11, 11])

levels = np.linspace(0,1,101)

for i in range(31):
  x1 = np.linspace(0,1,11)
  x2 = np.linspace(0,1,251)
  f = interp2d(x1, x1, A[i,:,:], kind='linear')
  z = f(x2, x2)
  print(z.shape)
  plt.contourf(z, levels=levels, cmap='viridis')
  plt.axis('equal')
  plt.title("Feature 1 = %lf" % (0.1*i))
  plt.colorbar()
  plt.show()
