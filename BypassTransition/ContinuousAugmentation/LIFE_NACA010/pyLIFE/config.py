param_file = "gamma_sep_params.dat"

sens_file  = "gamma_sep_sens.dat"

caseDict = {

  "010": {
    "fileList"      : ["Joe.in", "restart.out", "wallShearStress.ref", "pressure.ref", "tau.ref", "plotMarker.py", "gamma_max_params.dat"],
    "nNodes"        : 6,
    "nCoresPerNode" : 20,
    "DirectCmd"     : "mpirun -np 120 ../../../JOEfiles/joe_direct > log",
    "AdjointCmd"    : "mpirun -np 120 ../../../JOEfiles/joe_adjoint > log",
    "objFile"       : "obj.dat",
    "plotScript"    : "plotMarker.py",
  },

}
