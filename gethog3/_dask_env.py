

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

memory_slurm = "1GB"
time_slurm = "00:04:00"

if machine == "local":
    n_workers = 3
    n_jobs = 3
    cluster = LocalCluster(n_workers=n_workers)
elif machine == "slurm":
    n_core = 1
    n_proc = n_core
    n_jobs = 60
    cluster = SLURMCluster(cores=n_core, processes=n_proc, memory=str(memory_slurm), walltime=time_slurm)

 # n_core = n_proc avoid using multiple threads


cluster.scale(n_jobs)  # num_sbatch_job = njobs / ncore
client_dask = Client(cluster)

print(client_dask)
print(client_dask.dashboard_link)

if machine == "slurm":
    print(cluster.job_script())

a=24
# cluster.adapt(maximum_jobs=20)  #  This automatically launches and kill workers based on load.



#
# def init_dask(machine, ncore=1, njobs=1, memory_slurm="3GB", time_slurm="00:10:00"):
#     # print("*** client **** ", cluster.dashboard_link, cluster.get_logs())
#     # ncore = 1  # Total number of cores per job
#     # njobs = 1  # Cut the job up into this many processes.
#     # # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
#     nproc = ncore
#     if machine == "local":
#         cluster = LocalCluster()
#     elif machine == "slurm":
#         cluster = SLURMCluster(cores=ncore, processes=nproc, memory=str(memory_slurm), walltime=time_slurm)
#
#     cluster.scale(njobs)  # # ask for one jobs
#     client_dask = Client(cluster)
#
#     return client_dask
#
#
# machine = "slurm"
# client_dask = init_dask(machine, ncore=1, njobs=1, memory_slurm="2GB", time_slurm="00:10:00")
#





#
#
# class working_dask:
#
#     def __init__(self, machine, ncore=1, njobs=1, memory_slurm="3GB", time_slurm="00:10:00"):
#         self.machine = machine
#         # print("*** client **** ", cluster.dashboard_link, cluster.get_logs())
#         # ncore = 1  # Total number of cores per job
#         # njobs = 1  # Cut the job up into this many processes.
#         # # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
#         nproc = ncore
#         if machine == "local":
#             cluster = LocalCluster()
#         elif machine == "slurm":
#             cluster = SLURMCluster(cores=ncore, processes=nproc, memory=str(memory_slurm), walltime=time_slurm)
#
#         cluster.scale(njobs)  # # ask for one jobs
#         client = Client(cluster)
#
#     return client
#
# machine = "slurm"
# client = init_dask(machine, ncore=1, njobs=1, memory_slurm="2GB", time_slurm="00:10:00")





#
# machine = "local"  # slurm local
# ncore = 3
# njobs = 1
# memory_slurm = "1GB"
# time_slurm = "00:05:00"
# # print("*** client **** ", cluster.dashboard_link, cluster.get_logs())
# # ncore = 1  # Total number of cores per job
# # njobs = 1  # Cut the job up into this many processes.
# # # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
# nproc = ncore
# if machine == "local":
#     cluster = LocalCluster(cores=ncore, memory=str(memory_slurm))
# elif machine == "slurm":
#     cluster = SLURMCluster(cores=ncore, processes=nproc, memory=str(memory_slurm), walltime=time_slurm)
#
# cluster.scale(njobs)  # # ask for one jobs
# client_dask = Client(cluster)