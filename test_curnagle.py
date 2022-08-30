
import xml.etree.ElementTree as ET
import dill as dill_pickle
from os import listdir
from xml.dom import minidom
import os

print("started ")
working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# gene_trees_folder = ""  # working_folder + "/gene_trees_/"
# check gene_trees_folder exist otherwise mkdir this

#address_rhogs_folder = working_folder + "/rhog_g501_done/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
#species_tree_address = working_folder + "/archive/lineage_tree_qfo.phyloxml"
pickle_folder = working_folder + "/pickle_folder_29aug/"
# add warning when pickle folder is not empty
output_xml_name = "out_29aug_11pm.xml"
gene_id_pickle_file = working_folder + "gene_id_28aug.pickle"




orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                               "originVersion": "Nov 2021", "version": "0.3"})  #

if os.path.getsize(gene_id_pickle_file) > 0:
    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = dill_pickle.load(handle)
    # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
print("gene_id_name read ")

for query_species_name, list_prots in gene_id_name.items():

    species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
    database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
    genes_xml = ET.SubElement(database_xml, "genes")

    for (gene_idx_integer, query_prot_name) in list_prots:
        query_prot_name_pure = query_prot_name.split("||")[0].strip().split("|")[1]
        gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})

print("gene_xml created ")
pickle_files_adress = listdir(pickle_folder)
print("number of hog is", len(pickle_files_adress))
hogs_a_rhog_xml_all = []
for pickle_file_adress in pickle_files_adress:
    if os.path.getsize(pickle_folder + pickle_file_adress) > 0:
        with open(pickle_folder + pickle_file_adress, 'rb') as handle:
            hogs_a_rhog_xml_batch = dill_pickle.load(handle)  # hogs_a_rhog_xml_batch is a list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
        # hogs_rhogs_xml_all is a list of hog object.

print("after reading number of hogs is  ", len(hogs_a_rhog_xml_all))

groups_xml = ET.SubElement(orthoxml_file, "groups")

for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
    groups_xml.append(hogs_a_rhog_xml)

xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# print(xml_str[:-1000])

with open(working_folder +output_xml_name, "w") as file_xml:
    file_xml.write(xml_str)
file_xml.close()

print("orthoxml is written in  "+ working_folder +output_xml_name)






# from xml.dom import minidom
# import xml.etree.ElementTree as ET
# import _utils
# import _inferhog
#
# # from _utils import logger_hog
# import gethog3._utils_rhog as _utils_rhog
#
# # from distributed import get_client
# # from dask.distributed import rejoin, secede
#
# if __name__ == '__main__':
#     working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
#     gene_trees_folder = "" # working_folder + "/gene_trees_/"
#     # check gene_trees_folder exist otherwise mkdir this
#
#     address_rhogs_folder = working_folder + "/rhog_g10k/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
#     file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)
#
#
#     oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
#     # working_folder+"omamer_database/oma_path/OmaServer.h5"
#     print("rHOG inferece has started. The oma database address is in ", oma_database_address)
#     (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
#     (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
#     hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)
#
#     (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
#     prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
#
#     for prot_i, prot in enumerate(query_prot_names_species_mapped):
#
#         a=1
#
#
#
#
#     rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
#
#
#
#


# from distributed import Client
# from dask.distributed import worker_client
# from dask.distributed import LocalCluster
# # from dask_jobqueue import SLURMCluster
# import time
#
# def fib(n):
#     if n < 2:
#         return n
#     with worker_client() as client:
#         a_future = client.submit(fib, n - 1)
#         b_future = client.submit(fib, n - 2)
#         a, b = client.gather([a_future, b_future])
#     return a + b
#
# if __name__ == "__main__":
#     # cluster = SLURMCluster(cores=1, processes=1, memory="1G", walltime="00:05:00")
#     cluster = LocalCluster(n_workers=3, threads_per_worker=1)
#     # n_jobs = 2
#     # cluster.scale(n_jobs)
#     client_dask = Client(cluster)
#     time.sleep(5)
#     lst = [20, 30, 22, 55]
#     for i in range(5):
#         future = client_dask.submit(fib, lst[i])
#         result = future.result()
#         print(result)
#
#
# # ###### collect xml files
# #
# # import dill as dill_pickle
# # from os import listdir
# # import xml.etree.ElementTree as ET
# #
# # from xml.dom import minidom
# #
# # def collect_write_xml(working_folder, pickle_folder, output_xml_name):
# #
# #     orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
# #                                                    "originVersion": "Nov 2021", "version": "0.3"})  #
# #
# #     with open(working_folder + '/file_gene_id_name.pickle', 'rb') as handle:
# #         gene_id_name = dill_pickle.load(handle)
# #         # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
# #
# #     for query_species_name, list_prots in gene_id_name.items():
# #
# #         species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
# #         database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
# #         genes_xml = ET.SubElement(database_xml, "genes")
# #
# #         for (gene_idx_integer, query_prot_name) in list_prots:
# #             query_prot_name_pure = query_prot_name.split("||")[0].strip()
# #             gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
# #
# #         #groups_xml = ET.SubElement(orthoxml_file, "groups")
# #
# #     pickle_files_adress = listdir(pickle_folder)
# #
# #     hogs_a_rhog_xml_all = []
# #     for pickle_file_adress in pickle_files_adress:
# #         with open(pickle_folder + pickle_file_adress, 'rb') as handle:
# #             hogs_a_rhog_xml = dill_pickle.load(handle)
# #             hogs_a_rhog_xml_all += hogs_a_rhog_xml
# #
# #     print(len(hogs_a_rhog_xml_all))
# #
# #     groups_xml = ET.SubElement(orthoxml_file, "groups")
# #
# #     for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
# #         groups_xml.append(hogs_a_rhog_xml)
# #
# #     xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# #     print(xml_str)
# #
# #     with open(working_folder+output_xml_name, "w") as file_xml:
# #         file_xml.write(xml_str)
# #     file_xml.close()
# #
# #     print("orthoxml is written in "+ working_folder+output_xml_name)
# #     return 1
# #
# #
# #
# # output_xml_name = "out12a.xml"
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # pickle_folder = working_folder + "/pickle_folder/"
# #
# # collect_write_xml(working_folder, pickle_folder, output_xml_name)
#
#
# # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# #
# #
# # species_prot_dic = {}
# # for rhogid_num in rhogid_num_list:
# #     prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
# #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# #
# #     for prot_i in rhog_i:
# #         if format_prot_name == 1:  # qfo dataset
# #             prot_name = prot_i.name  # 'tr|E3JPS4|E3JPS4_PUCGT
# #             species_i = prot_name.split("|")[-1].split("_")[-1].strip()
# #             if species_i == 'RAT': species_i = "RATNO"
# #         elif format_prot_name == 0:  # bird dataset
# #             # rec.name  CLIRXF_R07389
# #             # prot_name = prot_i.name
# #             prot_descrip = prot_i.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
# #             species_i = prot_descrip.split(" ")[1].split("|")[-1]
# #             # species_name = prot_name.split("_")[0].strip()
# #
# #         # species_i = prot_i.id.split("|")[-1].split("_")[-1]
# #         if species_i in species_prot_dic:
# #             species_prot_dic[species_i].append(prot_i.id)
# #         else:
# #             species_prot_dic[species_i] = [prot_i.id]
# #         # all_prot_temp_list.append(prot_i.id)
# #
# # print("there are species ", len(species_prot_dic))
# # orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
# #                                                "originVersion": "Nov 2021", "version": "0.3"})  #
# #
# # number_roothog = len(rhogid_num_list)
# # num_per_parralel = 10
# # parralel_num = int(number_roothog / num_per_parralel)
# # if number_roothog != parralel_num * num_per_parralel: parralel_num += 1
# # rhogid_batch_list = []
# # for list_idx in range(parralel_num):
# #     if list_idx == parralel_num:
# #         rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
# #     else:
# #         rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
# #     rhogid_batch_list.append(rhogid_num_list_portion)
# #
# #
# #
# # for rhogid_batch_idx in range(len(rhogid_batch_list)):
# #     rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
# #
# #     gene_counter = 1000000 + rhogid_batch * 10000
# #     gene_id_name = {}
# #     query_species_names_rhogs = list(species_prot_dic.keys())
# #     for species_name in query_species_names_rhogs:
# #
# #         species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
# #         database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
# #         genes_xml = ET.SubElement(database_xml, "genes")
# #
# #     prot_list = species_prot_dic[species_name]
# #     for prot_itr in range(len(prot_list)):  # [12:15]
# #         prot_i_name = prot_list[prot_itr]
# #         gene_id_name[prot_i_name] = gene_counter
# #         if "|" in prot_i_name:
# #             prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
# #         else:
# #             prot_i_name_short = prot_i_name
# #
# #         gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
# #         gene_counter += 1
# #
# # groups_xml = ET.SubElement(orthoxml_file, "groups")
# #
#
#
#
#
# #
# # # Proteins in each file belong to the same species.
# #
# # # change the name of each file based on the species name inside each prot id
# #
# #
# # from os import listdir
# # from Bio import SeqIO
# # import os
# #
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # prot_folder = working_folder + "/omamer_search_old/proteome/"
# # project_files = listdir(prot_folder)
# # query_species_names_old = []
# # query_species_names_new = []
# # for file in project_files:
# #     if file.split(".")[-1] == "fa":
# #         file_name_split = file.split(".")[:-1]
# #         query_species_name_old = '.'.join(file_name_split)
# #         prot_address = prot_folder + query_species_name_old + ".fa"
# #         prots_record = list(SeqIO.parse(prot_address, "fasta"))
# #         prot_record = prots_record[0]
# #         prot_name = prot_record.name  # 'tr|E3JPS4|E3JPS4_PUCGT
# #         query_species_name_new = prot_name.split("|")[-1].split("_")[-1].strip()
# #         # if query_species_name_new == 'RAT': query_species_name_new = "RATNO"
# #         query_species_names_old.append(query_species_name_old)
# #         query_species_names_new.append(query_species_name_new)
# #
# # os.mkdir(working_folder+"/omamer_search")
# # os.mkdir(working_folder+"/omamer_search/proteome/")
# # os.mkdir(working_folder+"/omamer_search/hogmap")
# #
# #
# # for idx, query_species_name_old in enumerate(query_species_names_old):
# #     query_species_name_new = query_species_names_new[idx]
# #
# #     prot_address_old = working_folder + "omamer_search_old/proteome/" + query_species_name_old + ".fa"
# #     prot_address_new = working_folder + "omamer_search/proteome/" + query_species_name_new + "_.fa"
# #     os.system('cp ' + prot_address_old + ' ' + prot_address_new)
# #
# #     hogmap_address_old = working_folder + "omamer_search_old/hogmap/" + query_species_name_old + ".hogmap"
# #     hogmap_address_new = working_folder + "omamer_search/hogmap/" + query_species_name_new + "_.hogmap"
# #     os.system('cp ' + hogmap_address_old + ' ' + hogmap_address_new)
# #
# #
# # # 13:54:16 - the species DANRE  already exists in the oma database, remove them first
# #
# #
# #
# # print("done")
# #
#
#
#
#
#
#
# #
# # # Proteins in each file belong to the same species.
# #
# # from os import listdir
# # from Bio import SeqIO
# # import os
# #
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # prot_folder = working_folder + "/omamer_search_old/proteome/"
# # project_files = listdir(prot_folder)
# # query_species_names_old = []
# # query_species_names_new = []
# # for file in project_files:
# #     if file.split(".")[-1] == "fa":
# #         file_name_split = file.split(".")[:-1]
# #         query_species_name_old = '.'.join(file_name_split)
# #         prot_address = prot_folder + query_species_name_old + ".fa"
# #         prots_record = list(SeqIO.parse(prot_address, "fasta"))
# #         prot_record = prots_record[0]
# #         prot_name = prot_record.name  # 'tr|E3JPS4|E3JPS4_PUCGT
# #         query_species_name_new = prot_name.split("|")[-1].split("_")[-1].strip()
# #         if query_species_name_new == 'RAT': query_species_name_new = "RATNO"
# #         query_species_names_old.append(query_species_name_old)
# #         query_species_names_new.append(query_species_name_new)
# #
# # os.mkdir(working_folder+"/omamer_search")
# # os.mkdir(working_folder+"/omamer_search/proteome/")
# # os.mkdir(working_folder+"/omamer_search/hogmap")
# #
# #
# # for idx, query_species_name_old in enumerate(query_species_names_old):
# #     query_species_name_new = query_species_names_new[idx]
# #
# #     prot_address_old = working_folder + "omamer_search_old/proteome/" + query_species_name_old + ".fa"
# #     prot_address_new = working_folder + "omamer_search/proteome/" + query_species_name_new + ".fa"
# #     os.system('cp ' + prot_address_old + ' ' + prot_address_new)
# #
# #     hogmap_address_old = working_folder + "omamer_search_old/hogmap/" + query_species_name_old + ".hogmap"
# #     hogmap_address_new = working_folder + "omamer_search/hogmap/" + query_species_name_new + ".hogmap"
# #     os.system('cp ' + hogmap_address_old + ' ' + hogmap_address_new)
# #
# #
# # print("done")
#
#
#
#
# # # #
# # # # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# # # # rhogid_num_list_input = rhogid_num_list[9:13]
# # # #
# # # # (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_num_list_input,
# # # #                                                                                 address_rhogs_folder)
# # # # for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
# # # #     groups_xml.append(hogs_a_rhog_xml)
# # # # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# # # # print(xml_str)
# # # #
# # #
# # #
# # # #
# # # # def incc(x):
# # # #     return x + 1
# # # #
# # # # from dask.distributed import Client
# # # #
# # # # client = Client()  # start local workers as processes
# # # #
# # # # futures = client.map(incc, range(4))
# # #
# # # # a=2
# # # #
# # # # results = client.gather(futures)
# # # #
# # # # print(results)
# # # #
# # # # print("here")
# # # import time
# # #
# # # # def slow_pow(x,y):
# # # #     time.sleep(1)
# # # #     return x ** y
# # # #
# # # # print("s1\n",slow_pow(3,5))
# # #
# # # from dask.distributed import Client
# # # client = Client(processes=False)  # n_workers=2, threads_per_worker=2
# # # print("s0\n",client)
# # #
# # # # res = client.submit(slow_pow, 2,3)
# # # # print("s3\n", res)
# # # # print("s4\n", res.result())
# # #
# # #
# # # # powers_of_10 = []
# # # # for i in range(1,11):
# # # #     future = client.submit(slow_pow, i, 10)
# # # #     powers_of_10.append(future)
# # # # print("s3\n")
# # # # print("s4\n", [future.result() for future in powers_of_10])
# # #
# # #
# # #
# # # # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5)
# # # # print("s3\n")
# # # #
# # # # out1 = [future.result() for future in futures]
# # # #
# # # # print("s4\n", out1)
# # #
# # #
# # # # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5)
# # # #
# # # # [future.result() for future in futures]
# # # #
# # #
# # #
# # #
# # # # futures = client.map(incc, range(4))
# # #
# # #
# # # def slow_pow(x,y,z):
# # #     time.sleep(1)
# # #     print("s1", z)
# # #     return x * y
# # #
# # # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5, [23]*5)
# # # print("s2\n")
# # #
# # # out1 = [future.result() for future in futures]
# # #
# # # print("s3 ", out1)
# #
# #
# #
# #
# # # def square(n):
# # #     return n*n
# # # my_list = [2,3,4,5,6,7,8,9]
# # # updated_list = map(square, my_list)
# # # print(updated_list)
# # # print(list(updated_list))
# #
# #
# # # def myMapFunc(list1, list2):
# # #     return list1+list2
# # #
# # # my_list1 =list(range(5))
# # # my_list2 = [10]
# # #
# # # updated_list = map(myMapFunc, my_list1,my_list2)
# # # print(updated_list)
# # # print(list(updated_list))
# #
# #
# #
# #
# # # def cal(a,b,c,d):
# # #
# # #     print(a+b+c+d)
# # #     return 1
# # # list1=(2,3,4)
# # # print(cal(1,list1))
# # #
# # # from _utils import logger_hog
# # # import _utils
#
# #
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
# #
# #
# # # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# #
# # # print("num", len(rhogid_num_list))
# # #
# # # rhogid_len_list = []
# # # for rhogid_num in rhogid_num_list:
# # #     prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
# # #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# # #     rhogid_len_list.append(len(rhog_i))
# # #
# # # print(len(rhogid_len_list), rhogid_len_list[:2])
# # #
# #
# #
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # address_rhogs_folder = working_folder + "/old3/rhog_all/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
# # file = open(address_rhogs_folder+"size.txt")
# # rhogid_len_list =[]
# # for f in file:
# #     rhogid_len_list.append(int(f.strip()))
# # # print(rhogid_len_list)
# # import matplotlib.pyplot as plt
# #
# # plt.hist(rhogid_len_list, bins=100)  # , density=True
# # plt.yscale('log', nonposy='clip')
# # plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/size_qfo_all2.png")
# #
# # print("here2")
# #
# #
# #
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird/"
# # address_rhogs_folder = working_folder + "/rhogs_all/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
# # file = open(address_rhogs_folder+"size.txt")
# # rhogid_len_list =[]
# # for f in file:
# #     rhogid_len_list.append(int(f.strip()))
# # # print(rhogid_len_list)
# #
# #
# # plt.figure()
# # plt.hist(rhogid_len_list, bins=100)  # , density=True
# # plt.yscale('log', nonposy='clip')
# # plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/size_bird_all2.png")
# #
# # print("here")
#
# #
# # from datetime import datetime
# # import time
# #
# # current_time = datetime.now().strftime("%H:%M:%S")
# # #print(current_time)
# # # current_time += "sina"
# # #print(current_time)
# #
# # def aa(a):
# #     print("here2")
# #     time.sleep(5)
# #     return a*100
# #
# # from _dask_env import client_dask
# #
# #
# # futures = client_dask.map(aa, [200,100,23,42])
# # print([i.result() for i in futures])
# #
# #
# # futures = client_dask.map(aa, [222,22,11,432])
# # print([i.result() for i in futures])
# #
# #
# # futures = client_dask.map(aa, [2400,1400,323,423])
# # print([i.result() for i in futures])
# #
# #
# #
#
# # futures= []
# # print("here1")
# # future = client_dask.submit(aa, 200)
# # futures.append(future)
# # print("here3")
# # print(future)
# # print("here4")
# # #print(future.result())
# # print("here5")
# #
# #
# #
# #
# #
# # print("here10")
# # future = client_dask.submit(aa, 123)
# # futures.append(future)
# # print("here3")
# # print(future)
# # print("here4")
# # # print(future.result())
# # print("here5")
# #
# #
# # print("here10")
# # future = client_dask.submit(, 123)
# # print("here3")
# # print(future)
# # print("here4")
# # futures.append(future)
# # # futures.append(future)print(future.result())
# # print("here5")
# #
# # print([i.result() for i in futures])
# #
# # print(future)
# # #
# # # print("here2")
# # # future = client_dask.submit(aa, 200)
# # # print(future.result())
# # #
# # #
# # # print("here3")
# # # future = client_dask.submit(aa, 200)
# # # print(future.result())
# # #
# # def Fibonacci(n):
# #     # Check if input is 0 then it will
# #     # print incorrect input
# #     if n < 0:
# #         print("Incorrect input")
# #
# #     # Check if n is 0
# #     # then it will return 0
# #     elif n == 0:
# #         return 0
# #
# #     # Check if n is 1,2
# #     # it will return 1
# #     elif n == 1 or n == 2:
# #         return 1
# #
# #     else:
# #         return Fibonacci(n - 1) + Fibonacci(n - 2)
# #
# #
# # # Driver Program
# # print(Fibonacci(19))
#
# #
# #
# # from distributed import Client, get_client
# #
# # def fib(n):
# #     if n < 2:
# #         return n
# #     client = get_client()
# #     a_future = client.submit(fib, n - 1)
# #     b_future = client.submit(fib, n - 2)
# #     a, b = client.gather([a_future, b_future])
# #     return a + b
# #
# # if __name__ == "__main__":
# #     client = Client()
# #     future = client.submit(fib, 10)
# #     result = future.result()
# #     print(result)  # prints "55"
# #
