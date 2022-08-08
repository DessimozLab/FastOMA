


# import dask

from dask.distributed import Client
from dask.distributed import LocalCluster
from dask_jobqueue import SLURMCluster


machine = "slurm"
ncore = 1
njobs = 20
memory_slurm = "1GB"
time_slurm = "00:15:00"
# print("*** client **** ", cluster.dashboard_link, cluster.get_logs())
# ncore = 1  # Total number of cores per job
# njobs = 1  # Cut the job up into this many processes.
# # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
nproc = ncore
if machine == "local":
    cluster = LocalCluster()
elif machine == "slurm":
    cluster = SLURMCluster(cores=ncore, processes=nproc, memory=str(memory_slurm), walltime=time_slurm)

cluster.scale(njobs)  # # ask for one jobs
client_dask = Client(cluster)

print(client_dask)
print(client_dask.dashboard_link)







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