param_file = "gamma_max_params.dat"

sens_file  = "gamma_max_sens.dat"

caseDict = {

  "T3A": {
    "fileList"      : ["Joe.in", "restart.out", "adjoint_values.txt", "Cf.ref", "plotMarker.py", "FLATPLATE.dat.base"],
    "nNodes"        : 1,
    "nCoresPerNode" : 20,
    "DirectCmd"     : "srun -n 20 ../../../JOEfiles/joe_direct",
    "AdjointCmd"    : "srun -n 20 ../../../JOEfiles/joe_adjoint",
    "objFile"       : "obj.dat",
    "plotScript"    : "plotMarker.py",
  },
  
  "T3C1": {
    "fileList"      : ["Joe.in", "restart.out", "adjoint_values.txt", "Cf.ref", "plotMarker.py", "FLATPLATE.dat.base"],
    "nNodes"        : 2,
    "nCoresPerNode" : 20,
    "DirectCmd"     : "srun -n 40 ../../../JOEfiles/joe_direct",
    "AdjointCmd"    : "srun -n 40 ../../../JOEfiles/joe_adjoint",
    "objFile"       : "obj.dat",
    "plotScript"    : "plotMarker.py",
  },

}

'''

  "Sep": {
    "fileList"      : ["gamma_max_params.dat", "Joe.in", "restart.out", "adjoint_values.txt", "Cf.ref", "Cf.dat", "Cp.ref", "Cp.dat", "plotMarker.py"],
    "nNodes"        : 2,
    "nCoresPerNode" : 10,
    "DirectCmd"     : "srun -n 20 ../../../JOEfiles/joe_direct",
    "AdjointCmd"    : "srun -n 20 ../../../JOEfiles/joe_adjoint",
    "objFile"       : "obj.dat",
    "plotScript"    : "plotMarker.py",
  },

}

'''
