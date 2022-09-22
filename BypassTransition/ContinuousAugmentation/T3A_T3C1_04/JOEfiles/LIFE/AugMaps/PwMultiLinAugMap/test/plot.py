#!/home/vsriv/miniconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt("plot.dat")
A = np.reshape(A, [101, 101, 3])

B = np.loadtxt("plot_check.dat")
B = np.reshape(B, [101, 101, 3])

plt.contourf(A[:,:,0], A[:,:,1], A[:,:,2], levels=100)
plt.axis('equal')
plt.show()

plt.contourf(B[:,:,0], B[:,:,1], B[:,:,2], levels=100)
plt.axis('equal')
plt.show()
