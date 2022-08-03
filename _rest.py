## the following are needed when we start from a rootHOG fasta file.


# import pickle
import dill as pickle


def infer_HOG_rhog3(rhogid_num_list, gene_id_name):  # , address_rhogs_folder, species_tree_address):
    """
    The prot sequences of a rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
    we want to infer all subhogs of this rootHOG for different taxanomic levels.

    output: a python dict (HOG_thisLevel):  key=taxanomic level, value= a list of subhogs.
    """

    HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
    logger_hog.info("subhogs in thisLevel are " + ' '.join(["[" + str(i) + "]" for i in HOG_thisLevel]) + " .")

    for hog_i in HOG_thisLevel:
        print(hog_i)
        if len(hog_i._members) > 1:
            # could be improved
            HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
            HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
            # groups_xml.append(HOG_thisLevel_xml)
            # print(hog_i._members)
    # HOG_thisLevel_list.append(HOG_thisLevel)




# project="project1",, queue="normal"
# cluster = SLURMCluster(walltime='00:20:00', n_workers = NCORE, cores=NCORE,processes = NCORE,interface='ib0', memory="20GB",scheduler_options={'interface': 'ens2f0' })
#  env_extra=['source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda/etc/profile.d/conda.sh','conda activate ML2'],
print(cluster.job_script())
print(cluster.dashboard_link)

cluster.scale(jobs=njobs)  # # ask for one jobs

import time
time.sleep(5)
print(cluster)
# cluster.adapt(minimum=10, maximum=30)

client = Client(cluster, timeout='1000s', set_as_default=True)



len_HOG_thisLevel_all = []

number_roothog = len(rhogid_num_list_temp)

num_per_parralel = 18
parralel_num = int(number_roothog / num_per_parralel)

for list_idx in range(parralel_num + 1):

if list_idx == parralel_num:
    rhogid_num_list = rhogid_num_list_temp[list_idx * num_per_parralel:]
else:
    rhogid_num_list = rhogid_num_list_temp[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]

(out_len) = dask.delayed(infer_HOG_rhog3)(rhogid_num_list, gene_id_name)
len_HOG_thisLevel_all.append(out_len)

print("before computation", len(len_HOG_thisLevel_all), len_HOG_thisLevel_all[:2])

output_computed = dask.compute(*len_HOG_thisLevel_all)
print(" computation done ")





# import dill as pickle
# from os import listdir
# import xml.etree.ElementTree as ET
# from xml.dom import minidom
# import os

# address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"

# address_rhogs_folder  =address_working_folder + "/rhog_size_g2_s1k/"


# ## create a list of rootHOG IDs  stored in the folder of rHOG .
# rhog_files = listdir(address_rhogs_folder)
# rhogid_num_list= []
# for rhog_file in rhog_files:
#     if rhog_file.split(".")[-1] == "fa":
#         rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
#         rhogid_num_list.append(rhogid_num)
# print(len(rhogid_num_list)," .")

# rhogid_num_list_temp = rhogid_num_list#[:200]



