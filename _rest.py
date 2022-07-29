## the following are needed when we start from a rootHOG fasta file.


# import pickle
import dill as pickle


def infer_HOG_rhog3(rhogid_num_list, gene_id_name):  # , address_rhogs_folder, species_tree_address):
    """
    The prot sequences of a rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
    we want to infer all subHOGs of this rootHOG for different taxanomic levels.

    output: a python dict (HOG_thisLevel):  key=taxanomic level, value= a list of subHOGs.
    """

    HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
    logger_hog.info("subHOGs in thisLevel are " + ' '.join(["[" + str(i) + "]" for i in HOG_thisLevel]) + " .")

    for hog_i in HOG_thisLevel:
        print(hog_i)
        if len(hog_i._members) > 1:
            # could be improved
            HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
            HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
            # groups_xml.append(HOG_thisLevel_xml)
            # print(hog_i._members)
    # HOG_thisLevel_list.append(HOG_thisLevel)


del dic_sub_hogs
del HOG_thisLevel
gc.collect()

with open('/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/pickle_7ka/file_' + str(rhogid_num) + '.pickle',
          'wb') as handle:
    dill_pickle.dump(HOG_thisLevel_xml_all, handle, protocol=pickle.HIGHEST_PROTOCOL)

num_hog = len(HOG_thisLevel_xml_all)

del HOG_thisLevel_xml_all
gc.collect()

return (num_hog)


# def infer_HOG_rhog3(rhogid_num_list, gene_id_name):# , address_rhogs_folder, species_tree_address):
#     """
#     The prot sequences of a rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
#     we want to infer all subHOGs of this rootHOG for different taxanomic levels.

#     output: a python dict (HOG_thisLevel):  key=taxanomic level, value= a list of subHOGs.
#     """

#     HOG_thisLevel_list = []
#     len_HOG_thisLevel_list = []
#     HOG_thisLevel_xml_all = []
#     rhogid_num = 0
#     for rhogid_num in rhogid_num_list:
#         # rhogid_num = rhogid_num_list[rhogid_num_i]
#         logger_hog.info("\n"+"="*50+"\n"+"Working on root hog: "+str(rhogid_num)+". \n")  # +", ",rhogid_num_i,"-th. \n"
#         prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
#         rhog_i = list(SeqIO.parse(prot_address, "fasta"))
#         logger_hog.info("number of proteins in the rHOG is "+str(len(rhog_i))+".")

#         (species_tree) = read_species_tree(species_tree_address)
#         (species_tree, species_names_rhog, prot_names_rhog) = prepare_species_tree(rhog_i, species_tree)
#         #species_tree.write()  print(species_tree.write())

#         dic_sub_hogs = {}
#         # finding hogs at each level of species tree (from leaves to root, bottom up)
#         for node_species_tree in species_tree.traverse(strategy = "postorder"):
#             if node_species_tree.is_leaf() :
#                 # each leaf itself is a subhog
#                 continue
#             logger_hog.info("\n"+"*"*15+"\n"+"Finding hogs for the taxonomic level:"+ str(node_species_tree.name)+ "\n"+str(node_species_tree.write())+"\n")
#             dic_sub_msas = []
#             (dic_sub_hogs) = infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num)
#             HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
#             logger_hog.info("subHOGs in thisLevel are "+' '.join(["["+str(i)+"]" for i in HOG_thisLevel])+" .")

#             for hog_i in HOG_thisLevel:
#                 print(hog_i)
#                 if len(hog_i._members)>1:
#                     # could be improved
#                     HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
#                     HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
#                     #groups_xml.append(HOG_thisLevel_xml)
#                     #print(hog_i._members)
#         #HOG_thisLevel_list.append(HOG_thisLevel)
#         del dic_sub_hogs
#         del HOG_thisLevel
#         gc.collect()

#     with open('/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/pickle_7ka/file_'+str(rhogid_num)+'.pickle', 'wb') as handle:
#         pickle.dump(HOG_thisLevel_xml_all, handle, protocol=pickle.HIGHEST_PROTOCOL)


#     num_hog = len(HOG_thisLevel_xml_all)

#     del HOG_thisLevel_xml_all
#     gc.collect()

#     return (num_hog)


Ncore = 1  # Total number of cores per job
njobs = 20  # Cut the job up into this many processes.
# By default, process ~= sqrt(cores) so that the number of processes and the number of threads per process is roughly the same.
Nproc = Ncore;

cluster = SLURMCluster(cores=Ncore, processes=Nproc, memory="50GB",
                       walltime="00:45:00")
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


# with open(address_working_folder + "/filename.pickle", 'rb') as handle:
#     (groups_xml, gene_id_name,orthoxml_file) = pickle.load(handle)
# len(gene_id_name)

# # print("adding to xml started.")
# # #rhogid_num = 564451
# # for rhogid_num in rhogid_num_list_temp:
# #     pickl_file_path = address_working_folder + '/7i_1/pickle_7ha/file_'+str(rhogid_num)+'.pickle'
# #     if os.path.exists(pickle_file):
# #     else:
# #         print(pickl_file_path.split("/")[-1],"not exist")



