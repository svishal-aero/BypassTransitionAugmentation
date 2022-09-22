#!/home/svishal/miniconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size':22, 'figure.figsize':[10,8]})
mpl.rc('text', usetex=True)

def plotData():
    plt.figure()
    C = np.loadtxt("markerObjective.out")
    C = C[np.argsort(C[:,0])]
    xmin = C[ 0, 0]
    xmax = C[-1, 0]
    plt.plot([0,1], [0,0], '--', color='grey', label='Zero')
    plt.plot((C[:,0]-xmin)/(xmax-xmin), C[:,-2], '.k', label='Data', alpha=0.3)
    return [xmin, xmax]

def plotProfile(filenames, style, label, bounds):
    A = np.loadtxt(filenames[0]); A = A[np.argsort(A[:,0])]; A = A[15:-15]
    B = np.loadtxt(filenames[1]); B = B[np.argsort(B[:,0])]; B = B[25:-15]
    coordsA = (A[:,0]-bounds[0])/(bounds[1]-bounds[0])
    coordsB = (B[:,0]-bounds[0])/(bounds[1]-bounds[0])
    tauA = np.sign(A[:,7])*(A[:,-7]**2+A[:,-8]**2)**0.5
    tauB = np.sign(B[:,7])*(B[:,-7]**2+B[:,-8]**2)**0.5
    plt.plot(coordsA, tauA, style, lw=2, label=label)
    plt.plot(coordsB, tauB, style, lw=2)

def formatShowAndSavePlot():
    plt.xlabel('Percent Fractional Chord')
    plt.ylabel('Wall shear stress ($\\tau$)')
    plt.title('NACA 65-010')
    plt.ylim(bottom=-3)
    plt.legend()
    plt.tight_layout(pad=0.4)
    #plt.savefig('../../../figures/%s_%s' % (sys.argv[1], fileExt))
    plt.show()

#plotData()
#plotProfile(['bladelower.dat.base', 'bladeupper.dat.base'], '-r', 'Baseline')
#plotProfile(['bladelower.dat.beta1', 'bladeupper.dat.beta1'], '-g', '$$\\beta_1$$')
#formatShowAndSavePlot('1c')

#bounds = plotData()
#plotProfile(['bladelower.dat.beta1', 'bladeupper.dat.beta1'], '-r', '$$\\beta_1$$', bounds)
#plotProfile(['bladelower.dat.beta1beta2', 'bladeupper.dat.beta1beta2'], '-g', '$$\\beta_1\\beta_2$$', bounds)
#formatShowAndSavePlot('2n')

#bounds = plotData()
#plotProfile(['bladelower.dat.beta1', 'bladeupper.dat.beta1'], '-r', '$$\\beta_1$$', bounds)
#plotProfile(['bladelower.dat.beta1beta2', 'bladeupper.dat.beta1beta2'], '-b', '$$\\beta_1\\beta_2$$', bounds)
#plotProfile(['bladelower.dat.beta1beta2sigma', 'bladeupper.dat.beta1beta2sigma'], '-g', '$$\\beta_1\\beta_2^\\sigma$$', bounds)
#formatShowAndSavePlot('2')

bounds = plotData()
plotProfile(['bladelower.dat', 'bladeupper.dat'], '-g', '$$\\beta_1\\beta_2^\\sigma$$', bounds)
formatShowAndSavePlot()

#D = np.loadtxt("pressure.ref", skiprows=1)
#plt.plot((D[:,0]-xmin)/(xmax-xmin), D[:,-1], '.k', label='Data', alpha=0.3)
#plt.plot((A[:,0]-xmin)/(xmax-xmin), A[:,8] - 25, '-g', LineWidth=2, label='Augmented')
#plt.plot((B[:,0]-xmin)/(xmax-xmin), B[:,8] - 25, '-g', LineWidth=2)
#plt.xlabel('Fractional Chord')
#plt.ylabel('Pressure')
#plt.title('NACA 65-010 (Pressure distribution)')
#plt.legend()
#plt.xlim(left=0)
#plt.xlim(right=1)
#plt.tight_layout(pad=1.01)
#plt.show()
