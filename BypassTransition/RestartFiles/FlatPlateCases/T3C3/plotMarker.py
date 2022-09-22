#!/home/vsriv/miniconda3/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({'font.size':22, 'figure.figsize':[10,8]})
mpl.rc('text', usetex=True)

def plotData():
    plt.figure()
    A  = np.loadtxt("Cf.ref", skiprows=1)
    Cf = A[:,-1]
    x  = A[:,0]
    Cf = Cf[np.argsort(x)]
    x  = x[np.argsort(x)]
    plt.scatter(x[::5],Cf[::5],s=10,c='b',label="Data")

def plotProfile(filename, style, label):
    A  = np.loadtxt(filename)
    Cf = A[:,-2]
    x  = A[:,0]
    Cf = Cf[np.argsort(x)]
    x  = x[np.argsort(x)]
    plt.plot(x,Cf,style,label=label)

def formatSaveAndShowPlot(fileExt):
    plt.ylim(top=0.02)
    plt.ylim(bottom=0.0)
    plt.xlabel("x-coordinate")
    plt.ylabel("Skin friction coefficient")
    plt.title('T3C3')
    plt.legend()
    plt.tight_layout(pad=0.4)
    plt.savefig("../../../figures/%s_%s" % (sys.argv[1],fileExt))
    plt.show()

plotData()
plotProfile('FLATPLATE.dat.base', '-r', '$$\\gamma_{max}=1$$')
plotProfile('FLATPLATE.dat.beta1', '-g', '$$\\gamma_{max}=\\beta_1$$')
formatSaveAndShowPlot('1c')

plotData()
plotProfile('FLATPLATE.dat.beta1', '-r', '$$\\gamma_{max}=\\beta_1$$')
plotProfile('FLATPLATE.dat.beta1beta2', '-g', '$$\\gamma_{max}=\\beta_1\\beta_2$$')
formatSaveAndShowPlot('2n')

plotData()
plotProfile('FLATPLATE.dat.beta1beta2', '-r', '$$\\gamma_{max}=\\beta_1\\beta_2$$')
plotProfile('FLATPLATE.dat.beta1', '--k', '$$\\gamma_{max}=\\beta_1$$')
plotProfile('FLATPLATE.dat.beta1beta2sigma', '-g', '$$\\gamma_{max}=\\beta_1\\beta_2^\\sigma$$')
formatSaveAndShowPlot('2')
