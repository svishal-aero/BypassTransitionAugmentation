import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from subprocess import call
from .config import *

mpl.rc("text", usetex=True)
mpl.rcParams.update({"font.size":22, "figure.figsize":[8,6]})


#------------------------------------------------------------------------------
# Constrain the augmentation function parameters, if needed
#------------------------------------------------------------------------------
def constrain_params(params):
  params[params<0.0] = 0.0
  params[params>1.0] = 1.0


#------------------------------------------------------------------------------
# Change the step length w.r.t. iterations, if desired
#------------------------------------------------------------------------------
def get_step_length(optim_iter):
  return 0.2 #0.1 + 0.1 * np.exp(-optim_iter/5.0)


#------------------------------------------------------------------------------
# Normalize the sensitivities, as desired
#------------------------------------------------------------------------------
def get_normalized_sens(optim_iter):
  sens = 0
  for casename in caseDict.keys():
    case_sens_file = "%04d/ADJOINT/%s/%s" % (optim_iter, casename, sens_file)
    case_sens = np.loadtxt(case_sens_file)
    sens = sens + case_sens
  return sens / np.max(np.abs(sens))


#------------------------------------------------------------------------------
# Given sensitivity vector and step_length, update parameters
#------------------------------------------------------------------------------
def write_updated_parameters(optim_iter):
  params         = np.loadtxt("%04d/%s" % (optim_iter-1, param_file))
  step_length    = get_step_length(optim_iter)
  norm_sens      = get_normalized_sens(optim_iter-1)
  updated_params = params - step_length * norm_sens
  constrain_params(updated_params)
  np.savetxt("%04d/%s" % (optim_iter, param_file), updated_params)


#------------------------------------------------------------------------------
# Check if a direct/adjoint run is complete
#------------------------------------------------------------------------------
def run_is_complete(optim_iter, mode, casename):
  return os.path.exists("%04d/%s/%s/completed" % (optim_iter, mode, casename))


#------------------------------------------------------------------------------
# Check if all direct/adjoint runs are complete
#------------------------------------------------------------------------------
############################# DEPRECATED ######################################
def all_cases_are_complete(optim_iter, mode):
  all_complete = True
  for casename in caseDict.keys():
    all_complete = all_complete and run_is_complete(optim_iter, mode, casename)
  return all_complete


#------------------------------------------------------------------------------
# Do not proceed until all direct/adjoint runs have completed
#------------------------------------------------------------------------------
############################# DEPRECATED ######################################
def wait_until_all_cases_are_complete(optim_iter, mode):
  while True:
    time.sleep(3)
    if all_cases_are_complete(optim_iter, mode): break


#------------------------------------------------------------------------------
# Any callback function (file deletion, etc.) after a direct/adjoint run
#------------------------------------------------------------------------------
def callback_fn(optim_iter, mode, casename):
  if run_is_complete(optim_iter, mode, casename):
    sys.stdout.write("Case \"%s\" completed\n" % (casename))
    sys.stdout.flush()
    if mode=="ADJOINT":
      dirname = "%04d/%s/%s" % (optim_iter, mode, casename)
      sys.stdout.write("Deleting ADOLC files for case %s... "%(casename))
      sys.stdout.flush()
      call("rm -f %s/ADOLC-*" % (dirname), shell=True)
      sys.stdout.write("DONE\n")
      sys.stdout.flush()
      sys.stdout.write("Aborting job for case %s... "%(casename))
      sys.stdout.flush()
      call("bkill %d" % (int(np.loadtxt("%s/job_id" % (dirname)))), shell=True)
      '''
      call("scancel %d" % (int(np.loadtxt("%s/job_id" % (dirname)))), shell=True)
      '''
      sys.stdout.write("DONE\n")
      sys.stdout.flush()
    return True
  else:
    return False


#------------------------------------------------------------------------------
# Call the callback function for all cases after direct/adjoint runs
#------------------------------------------------------------------------------
def callback_for_all_cases(optim_iter, mode):
  pendingCaseList = list(caseDict.keys())
  done = [False for casename in pendingCaseList]
  while True:
    all_done = True
    time.sleep(3)
    for i in range(len(pendingCaseList)):
      if done[i]==False:
        done[i] = callback_fn(optim_iter, mode, pendingCaseList[i])
      all_done = all_done and done[i]
    if all_done:
      break


#------------------------------------------------------------------------------
# Create corresponding directories for direct/adjoint containing copied files
#------------------------------------------------------------------------------
def create_directories(optim_iter, mode):
  if mode=="DIRECT":
    call("mkdir -p %04d" % (optim_iter), shell=True)
    write_updated_parameters(optim_iter)
    for casename in caseDict.keys():
      call("mkdir -p %04d/DIRECT/%s" % (optim_iter, casename), shell=True)
      src_file = "%04d/%s" % (optim_iter, param_file)
      dest_dir = "%04d/DIRECT/%s" % (optim_iter, casename)
      call("cp %s %s" % (src_file, dest_dir), shell=True)
      for filename in caseDict[casename]["fileList"]:
        src_file = "%04d/ADJOINT/%s/%s" % (optim_iter-1, casename, filename)
        dest_dir = "%04d/DIRECT/%s" % (optim_iter, casename)
        call("cp %s %s" % (src_file, dest_dir), shell=True)
  if mode=="ADJOINT":
    for casename in caseDict.keys():
      call("mkdir -p %04d/ADJOINT/%s" % (optim_iter, casename), shell=True)
      src_file = "%04d/%s" % (optim_iter, param_file)
      dest_dir = "%04d/ADJOINT/%s" % (optim_iter, casename)
      call("cp %s %s" % (src_file, dest_dir), shell=True)
      for filename in caseDict[casename]["fileList"]:
        src_file = "%04d/DIRECT/%s/%s" % (optim_iter, casename, filename)
        dest_dir = "%04d/ADJOINT/%s" % (optim_iter, casename)
        call("cp %s %s" % (src_file, dest_dir), shell=True)


#------------------------------------------------------------------------------
# Run direct/adjoint simulations
#------------------------------------------------------------------------------
def run_case(optim_iter, mode, casename):
  os.chdir("%04d/%s/%s" % (optim_iter, mode, casename))
  cfg = caseDict[casename]
  with open("script.sh","w") as f:
    f.write("#BSUB -e job.errors\n")
    f.write("#BSUB -o job.output\n")
    f.write("#BSUB -J %s%s%03d\n" % (mode[0], casename, optim_iter))
    f.write("#BSUB -n %d\n" % (cfg["nNodes"]*cfg["nCoresPerNode"]))
    f.write("#BSUB -R span[ptile=%d]\n" % (cfg["nCoresPerNode"]))
    f.write("#BSUB -W 01:00\n")
    f.write("#BSUB -q normal\n")
    f.write("#BSUB -R rusage[mem=16000]\n")
    f.write("#BSUB -R affinity[thread(1):cpubind=thread:distribute=balance]\n")
    f.write("\n")
    f.write("echo $LSB_JOBID > job_id\n")
    f.write("\n")
    f.write("source ~/.bashrc\n")
    '''
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=%s%s%03d\n" % (mode[0], casename, optim_iter))
    f.write("#SBATCH --nodes=%d\n" % (cfg["nNodes"]))
    f.write("#SBATCH --ntasks-per-node=%d\n" % (cfg["nCoresPerNode"]))
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --mem-per-cpu=1g\n")
    f.write("#SBATCH --time=02:00:00\n")
    f.write("#SBATCH --partition=kdur\n")
    f.write("#SBATCH --mail-type=FAIL\n")
    f.write("#SBATCH --mail-user=vsriv@umich.edu\n")
    f.write("#SBATCH --get-user-env\n\n")
    f.write("source ~/.bashrc\n")
    f.write("echo $SLURM_JOB_ID > job_id\n")
    '''
    if mode=="DIRECT" : f.write("%s\n" % (cfg["DirectCmd"]))
    if mode=="ADJOINT": f.write("%s\n" % (cfg["AdjointCmd"]))

  call("bsub < script.sh", shell=True)
  '''
  call("sbatch script.sh", shell=True)
  '''
  os.chdir("../../..")


#------------------------------------------------------------------------------
# Run all direct/adjoint simulations
#------------------------------------------------------------------------------
def run_all_cases(optim_iter, mode):
  create_directories(optim_iter, mode)
  for casename in caseDict.keys():
    run_case(optim_iter, mode, casename)
  #wait_until_all_cases_are_complete(optim_iter, mode)
  callback_for_all_cases(optim_iter, mode)
  time.sleep(3)


#------------------------------------------------------------------------------
# Optimize augmentation function parameters
#------------------------------------------------------------------------------
def optimize_params(first_optim_iter=1, num_optim_iters=100):
  if first_optim_iter==1:
    call("cp -r Baseline 0000", shell=True)
  for optim_iter in range(first_optim_iter, num_optim_iters):
    sys.stdout.write("Iteration %04d Running Direct...\n" % (optim_iter))
    sys.stdout.flush()
    run_all_cases(optim_iter, "DIRECT")
    time.sleep(10)
    sys.stdout.write("Iteration %04d Running Adjoint...\n" % (optim_iter))
    sys.stdout.flush()
    run_all_cases(optim_iter, "ADJOINT")


#------------------------------------------------------------------------------
# Plot objective function relative to original value across all iterations
#------------------------------------------------------------------------------
def plot_objective_history(upto_optim_iter=0):
  call("mkdir -p figures", shell=True)
  obj_hist_dict = {}
  for casename in caseDict.keys():
    obj_hist_dict[casename] = []
  optim_iter = 0
  while os.path.exists("%04d/ADJOINT" % (optim_iter)):
    for casename in caseDict.keys():
      obj_filename = "%04d/DIRECT/%s/%s" % (optim_iter, casename,\
                                            caseDict[casename]["objFile"])
      obj = float(np.loadtxt(obj_filename))
      obj_hist_dict[casename].append(obj)
    if upto_optim_iter>0 and upto_optim_iter==optim_iter: break
    optim_iter += 1
  plt.figure()
  for casename, obj_hist in obj_hist_dict.items():
    plt.plot(np.array(obj_hist) / obj_hist[0], '-o', label=casename)
    plt.xlabel("Optimization iterations")
    plt.ylabel("Relative objective values")
    plt.title("Objective reduction")
    plt.legend()
  plt.tight_layout(pad=0.4)
  plt.savefig("figures/obj_hist.png")
  plt.show()
