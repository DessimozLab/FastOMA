

"""
dask ntoes
- IO issue open same file with slurm
-  threads_per_worker=1 causes error with LocalCluster
-  N_WORKERS = 3 njobs = 2 is not working, but njobs=10 is
"""

from dask.distributed import Client
from dask.distributed import LocalCluster
from dask_jobqueue import SLURMCluster

import _config

if _config.dask_machine == "local":
    n_workers = _config.dask_n_core
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)

elif _config.dask_machine == "slurm":
    # n_core = n_proc avoid using multiple threads
    cluster = SLURMCluster(cores=_config.dask_n_core, processes=_config.dask_n_core, memory=str(_config.dask_memory_slurm),
                           walltime=_config.dask_time_slurm)

cluster.scale(_config.dask_n_jobs)  # num_sbatch_job = njobs / ncore


client_dask = Client(cluster)

print(client_dask)
print(client_dask.dashboard_link)

if _config.dask_machine == "slurm":
    print(cluster.job_script())
