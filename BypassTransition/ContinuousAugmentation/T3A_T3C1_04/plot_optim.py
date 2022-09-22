#!/home/vsriv/miniconda3/bin/python

import sys
from pyLIFE import plot_objective_history

if __name__=="__main__":

  plot_objective_history(upto_optim_iter=int(sys.argv[1]))
