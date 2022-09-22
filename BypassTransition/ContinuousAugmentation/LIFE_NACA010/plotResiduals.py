import numpy as np
import matplotlib.pyplot as plt

maxIter = 33

def getColor(i):
    return [1-i/maxIter, i/maxIter, 0]

plt.rcParams.update({
    'font.size' : 18,
    'figure.figsize' : [8,6]
})

for i in range(2, maxIter+1):
    if i!=28:
        res = np.loadtxt('%04d/DIRECT/010/residuals.dat' % (i))[:,[0,5]]
        plt.semilogy(res[:,0], res[:,1], 'k', alpha=0.2)

res = np.loadtxt('%04d/DIRECT/010/residuals.dat' % (1))[:,[0,5]]
plt.semilogy(res[:,0], res[:,1], 'k', alpha=0.2, label='Iterate')

res = np.loadtxt('%04d/DIRECT/010/residuals.dat' % (0))[:,[0,5]]
plt.semilogy(res[:,0], res[:,1], 'r', label='Baseline')

res = np.loadtxt('%04d/DIRECT/010/residuals.dat' % (28))[:,[0,5]]
plt.semilogy(res[:,0], res[:,1], 'g', label='Optimal')

plt.xlabel('Forward solver iterations')
plt.ylabel('Energy residuals')
plt.legend(ncol=3, loc='lower right', bbox_to_anchor=[1, 1.02])
plt.tight_layout(pad=0.5)
plt.savefig('figures/NACAresidual.png')
plt.show()

'''

for i in range(maxIter+1):
    res = np.loadtxt('%04d/ADJOINT/010/residualsAD.dat' % (i), skiprows=1)[:,[0,5]]
    plt.semilogy(res[:,0], res[:,1], 'k', alpha=0.3, color=getColor(i))

plt.xlabel('Adjoint solver iterations')
plt.ylabel('Energy adjoint residuals')
plt.show()

#'''
