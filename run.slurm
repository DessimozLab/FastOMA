#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 100G
#SBATCH --time 00:30:00


export PYTHONUNBUFFERED=1

echo "started"

folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird/"

python ${folder}/gethog3_3Aug4/main_bird.py


echo "finished"