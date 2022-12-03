#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 3G
#SBATCH --time 05:00:00

export PYTHONUNBUFFERED=2

echo "started"
#python main.py
nextflow run nextf_script.nf

echo "finished"
