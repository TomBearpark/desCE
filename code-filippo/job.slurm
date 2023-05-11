#!/bin/bash
#SBATCH --job-name=fesearch      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=32G                # memory per node (per cpu-core use mem-per-cpu)
#SBATCH --time=4:00:00           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=fpalomba@princeton.edu

Rscript /scratch/network/fpalomba/BP_2023_fesearch/code/main.R
