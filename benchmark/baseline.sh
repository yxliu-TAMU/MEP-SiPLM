#!/bin/bash

##NESSARY JOB SPECIFICATIONS
#SBATCH --job-name=mutation_baseline       #Set the job name to "JobExample4"
#SBATCH --time=20:00:00              #Set the wall clock limit to 1hr and 30min
#SBATCH --ntasks=1                   #Request 1 task
#SBATCH --mem=180G                  #Request 2560MB (2.5GB) per node
#SBATCH --output=mutation_baseline.%1     #Send stdout/err to "Example4Out.[jobID]"
#SBATCH --gres=gpu:a100:1                #Request 1 GPU per node can be 1 or 2
#SBATCH --partition=gpu              #Request the GPU partition/queue

##OPTIONAL JOB SPECIFICATIONS
##SBATCH --account=132715540063             #Set billing account to 123456
##SBATCH --mail-type=ALL              #Send email on all job events
##SBATCH --mail-user=yxliu@tamu.edu    #Send all emails to email_address 

#First Executable Line
conda env list
source /home/yxliu/.bashrc
source activate protein-graph
python baseline.py