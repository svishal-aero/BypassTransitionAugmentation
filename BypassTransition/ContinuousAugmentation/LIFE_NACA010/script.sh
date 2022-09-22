#BSUB -e job.errors
#BSUB -o job.output
#BSUB -J N00
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -W 48:00
#BSUB -q normal
#BSUB -R rusage[mem=16000]
#BSUB -R affinity[thread(1):cpubind=thread:distribute=balance]

echo $LSB_JOBID > job_id

source ~/.bashrc
python ./optim.py > optim_output
