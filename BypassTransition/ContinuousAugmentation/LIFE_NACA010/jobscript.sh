#!/bin/bash
#SBATCH --job-name=FIML_NACA
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1g
#SBATCH --time=96:00:00
#SBATCH --partition=kdur
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vsriv@umich.edu
#SBATCH --get-user-env

source ~/.bashrc
#rm slurm*
./optim.py > optim_output
