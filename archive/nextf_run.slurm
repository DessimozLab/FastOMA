#!/bin/bash


#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 3G
#SBATCH --time 09:00:00


# # SBATCH --tmp=1G now it is a comment

export PYTHONUNBUFFERED=1

export TMPDIR="/scratch/smajidi1/mytmp"

echo "started"
#python main.py
nextflow run nextf_script.nf

echo "finished"
