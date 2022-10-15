#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 150G
#SBATCH --time 05:00:00


export PYTHONUNBUFFERED=2

echo "started"
#folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget"
# ${folder}/qfo2/gethog3_3Aug4/main.py

python main.py


echo "finished"