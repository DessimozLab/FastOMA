# # #
import dill as dill_pickle
import os
import _utils
from os import listdir
import xml.etree.ElementTree as ET

from Bio import SeqIO
from xml.dom import minidom


working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# gene_trees_folder = working_folder + "/gene_trees_test/"
address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # "/rhog_size_g2_s500/" sample_rootHOG
# species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"
pickle_address = working_folder + "/pickle_folder_6aug/"
format_prot_name = 1



pickle_files_adress = listdir(pickle_address)[:10]

hogs_a_rhog_xml_all = []
for pickle_file_adress in pickle_files_adress:
    with open(pickle_address+ pickle_file_adress, 'rb') as handle:

        hogs_a_rhog_xml_all += dill_pickle.load(handle)

print(len(hogs_a_rhog_xml_all))


print(2)
rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)


species_prot_dic = {}
for rhogid_num in rhogid_num_list:
    prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
    rhog_i = list(SeqIO.parse(prot_address, "fasta"))

    for prot_i in rhog_i:
        if format_prot_name == 1:  # qfo dataset
            prot_name = prot_i.name  # 'tr|E3JPS4|E3JPS4_PUCGT
            species_i = prot_name.split("|")[-1].split("_")[-1].strip()
            if species_i == 'RAT': species_i = "RATNO"
        elif format_prot_name == 0:  # bird dataset
            # rec.name  CLIRXF_R07389
            # prot_name = prot_i.name
            prot_descrip = prot_i.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
            species_i = prot_descrip.split(" ")[1].split("|")[-1]
            # species_name = prot_name.split("_")[0].strip()

        # species_i = prot_i.id.split("|")[-1].split("_")[-1]
        if species_i in species_prot_dic:
            species_prot_dic[species_i].append(prot_i.id)
        else:
            species_prot_dic[species_i] = [prot_i.id]
        # all_prot_temp_list.append(prot_i.id)

print("there are species ", len(species_prot_dic))
orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                               "originVersion": "Nov 2021", "version": "0.3"})  #

number_roothog = len(rhogid_num_list)
num_per_parralel = 10
parralel_num = int(number_roothog / num_per_parralel)
if number_roothog != parralel_num * num_per_parralel: parralel_num += 1
rhogid_batch_list = []
for list_idx in range(parralel_num):
    if list_idx == parralel_num:
        rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
    else:
        rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
    rhogid_batch_list.append(rhogid_num_list_portion)



for rhogid_batch_idx in range(len(rhogid_batch_list)):
    rhogid_batch = rhogid_batch_list[rhogid_batch_idx]

    gene_counter = 1000000 + rhogid_batch * 10000
    gene_id_name = {}
    query_species_names_rhogs = list(species_prot_dic.keys())
    for species_name in query_species_names_rhogs:

        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")

    prot_list = species_prot_dic[species_name]
    for prot_itr in range(len(prot_list)):  # [12:15]
        prot_i_name = prot_list[prot_itr]
        gene_id_name[prot_i_name] = gene_counter
        if "|" in prot_i_name:
            prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
        else:
            prot_i_name_short = prot_i_name

        gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
        gene_counter += 1

groups_xml = ET.SubElement(orthoxml_file, "groups")




# # #
# # # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# # # rhogid_num_list_input = rhogid_num_list[9:13]
# # #
# # # (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_num_list_input,
# # #                                                                                 address_rhogs_folder)
# # # for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
# # #     groups_xml.append(hogs_a_rhog_xml)
# # # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# # # print(xml_str)
# # #
# #
# #
# # #
# # # def incc(x):
# # #     return x + 1
# # #
# # # from dask.distributed import Client
# # #
# # # client = Client()  # start local workers as processes
# # #
# # # futures = client.map(incc, range(4))
# #
# # # a=2
# # #
# # # results = client.gather(futures)
# # #
# # # print(results)
# # #
# # # print("here")
# # import time
# #
# # # def slow_pow(x,y):
# # #     time.sleep(1)
# # #     return x ** y
# # #
# # # print("s1\n",slow_pow(3,5))
# #
# # from dask.distributed import Client
# # client = Client(processes=False)  # n_workers=2, threads_per_worker=2
# # print("s0\n",client)
# #
# # # res = client.submit(slow_pow, 2,3)
# # # print("s3\n", res)
# # # print("s4\n", res.result())
# #
# #
# # # powers_of_10 = []
# # # for i in range(1,11):
# # #     future = client.submit(slow_pow, i, 10)
# # #     powers_of_10.append(future)
# # # print("s3\n")
# # # print("s4\n", [future.result() for future in powers_of_10])
# #
# #
# #
# # # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5)
# # # print("s3\n")
# # #
# # # out1 = [future.result() for future in futures]
# # #
# # # print("s4\n", out1)
# #
# #
# # # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5)
# # #
# # # [future.result() for future in futures]
# # #
# #
# #
# #
# # # futures = client.map(incc, range(4))
# #
# #
# # def slow_pow(x,y,z):
# #     time.sleep(1)
# #     print("s1", z)
# #     return x * y
# #
# # futures = client.map(slow_pow, [1,2,3,4,5], [10]*5, [23]*5)
# # print("s2\n")
# #
# # out1 = [future.result() for future in futures]
# #
# # print("s3 ", out1)
#
#
#
#
# # def square(n):
# #     return n*n
# # my_list = [2,3,4,5,6,7,8,9]
# # updated_list = map(square, my_list)
# # print(updated_list)
# # print(list(updated_list))
#
#
# # def myMapFunc(list1, list2):
# #     return list1+list2
# #
# # my_list1 =list(range(5))
# # my_list2 = [10]
# #
# # updated_list = map(myMapFunc, my_list1,my_list2)
# # print(updated_list)
# # print(list(updated_list))
#
#
#
#
# # def cal(a,b,c,d):
# #
# #     print(a+b+c+d)
# #     return 1
# # list1=(2,3,4)
# # print(cal(1,list1))
# #
# # from _utils import logger_hog
# # import _utils

#
# working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
#
#
# # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
#
# # print("num", len(rhogid_num_list))
# #
# # rhogid_len_list = []
# # for rhogid_num in rhogid_num_list:
# #     prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
# #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# #     rhogid_len_list.append(len(rhog_i))
# #
# # print(len(rhogid_len_list), rhogid_len_list[:2])
# #
#
#
# working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# address_rhogs_folder = working_folder + "/old3/rhog_all/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
# file = open(address_rhogs_folder+"size.txt")
# rhogid_len_list =[]
# for f in file:
#     rhogid_len_list.append(int(f.strip()))
# # print(rhogid_len_list)
# import matplotlib.pyplot as plt
#
# plt.hist(rhogid_len_list, bins=100)  # , density=True
# plt.yscale('log', nonposy='clip')
# plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/size_qfo_all2.png")
#
# print("here2")
#
#
#
# working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird/"
# address_rhogs_folder = working_folder + "/rhogs_all/"  # rhogs "/rhog_size_g2_s500/" sample_rootHOG
# file = open(address_rhogs_folder+"size.txt")
# rhogid_len_list =[]
# for f in file:
#     rhogid_len_list.append(int(f.strip()))
# # print(rhogid_len_list)
#
#
# plt.figure()
# plt.hist(rhogid_len_list, bins=100)  # , density=True
# plt.yscale('log', nonposy='clip')
# plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/size_bird_all2.png")
#
# print("here")

#
# from datetime import datetime
# import time
#
# current_time = datetime.now().strftime("%H:%M:%S")
# #print(current_time)
# # current_time += "sina"
# #print(current_time)
#
# def aa(a):
#     print("here2")
#     time.sleep(5)
#     return a*100
#
# from _dask_env import client_dask
#
#
# futures = client_dask.map(aa, [200,100,23,42])
# print([i.result() for i in futures])
#
#
# futures = client_dask.map(aa, [222,22,11,432])
# print([i.result() for i in futures])
#
#
# futures = client_dask.map(aa, [2400,1400,323,423])
# print([i.result() for i in futures])
#
#
#

# futures= []
# print("here1")
# future = client_dask.submit(aa, 200)
# futures.append(future)
# print("here3")
# print(future)
# print("here4")
# #print(future.result())
# print("here5")
#
#
#
#
#
# print("here10")
# future = client_dask.submit(aa, 123)
# futures.append(future)
# print("here3")
# print(future)
# print("here4")
# # print(future.result())
# print("here5")
#
#
# print("here10")
# future = client_dask.submit(, 123)
# print("here3")
# print(future)
# print("here4")
# futures.append(future)
# # futures.append(future)print(future.result())
# print("here5")
#
# print([i.result() for i in futures])
#
# print(future)
# #
# # print("here2")
# # future = client_dask.submit(aa, 200)
# # print(future.result())
# #
# #
# # print("here3")
# # future = client_dask.submit(aa, 200)
# # print(future.result())
# #
# def Fibonacci(n):
#     # Check if input is 0 then it will
#     # print incorrect input
#     if n < 0:
#         print("Incorrect input")
#
#     # Check if n is 0
#     # then it will return 0
#     elif n == 0:
#         return 0
#
#     # Check if n is 1,2
#     # it will return 1
#     elif n == 1 or n == 2:
#         return 1
#
#     else:
#         return Fibonacci(n - 1) + Fibonacci(n - 2)
#
#
# # Driver Program
# print(Fibonacci(19))

#
#
# from distributed import Client, get_client
#
# def fib(n):
#     if n < 2:
#         return n
#     client = get_client()
#     a_future = client.submit(fib, n - 1)
#     b_future = client.submit(fib, n - 2)
#     a, b = client.gather([a_future, b_future])
#     return a + b
#
# if __name__ == "__main__":
#     client = Client()
#     future = client.submit(fib, 10)
#     result = future.result()
#     print(result)  # prints "55"
#
