

"""

dask ntoes

- IO issue open same file with slurm
-  threads_per_worker=1 causes error with LocalCluster
-  N_WORKERS = 3 njobs = 2 is not working, but njobs=10 is


"""
# import dask

from dask.distributed import Client
from dask.distributed import LocalCluster
import math
from dask_jobqueue import SLURMCluster


machine = "slurm"  # slurm local

memory_slurm = "2GB"
time_slurm = "00:20:00"

if machine == "local":
    n_workers = 3
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
    n_jobs= 3

elif machine == "slurm":
    n_core = 1
    n_proc = n_core
    n_jobs = 10
    #n_jobs = 3
    cluster = SLURMCluster(cores=n_core, processes=n_proc, memory=str(memory_slurm), walltime=time_slurm)

    # n_jobs = 3
cluster.scale(n_jobs)  # num_sbatch_job = njobs / ncore

 # n_core = n_proc avoid using multiple threads



client_dask = Client(cluster)

print(client_dask)
print(client_dask.dashboard_link)

if machine == "slurm":
    print(cluster.job_script())

