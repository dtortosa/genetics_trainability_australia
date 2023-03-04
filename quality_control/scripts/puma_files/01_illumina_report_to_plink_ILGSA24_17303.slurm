#!/bin/bash
 
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=inputs_prep
### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
### REQUIRED. Specify the PI group for this job
#SBATCH --account=denard
### Optional. Request email when job begins and ends
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=dftortosa@email.arizona.edu
### REQUIRED. Set the partition for your job.
#SBATCH --partition=high_priority
### OPTIONAL. Select the cores bought by David.
#SBATCH --qos=user_qos_denard
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=50
### OPTIONAL. You can set the number of cores per task. This can be useful for scripts in which inside the parallelized process there are other subprocesses, like the iHS calculation with hapbin per population. 
#SBATCH --cpus-per-task=1
### REQUIRED. Set the memory required for this job. I will set 40GC per each of the 100 cores=4000GB; https://public.confluence.arizona.edu/display/UAHPC/Allocation+and+Limits
###SBATCH --mem=400gb
### REQUIRED. Set the gb per core. YOU HAVE TO SELECT --mem or --mem-per-cpu but NOT BOTH.
#SBATCH --mem-per-cpu=10gb
### set the constraint for high memory nodes in case you use a lot of memory per node
###SBATCH --constraint=hi_mem
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=240:00:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### Load required modules/libraries if needed. Recommendation from Puma is use "module --ignore-cache load "singularity"" instead of "module load singularity". It is also possible your cache file is out-of-date.
module load singularity
### change to your script’s directory
cd /xdisk/denard/dftortosa/march_2023/
### Run your work
singularity run 01_illumina_report_to_plink.sif
