# # #
import dill as dill_pickle
import os
import _utils
from os import listdir
import xml.etree.ElementTree as ET

from xml.dom import minidom


working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
gene_trees_folder = working_folder + "/gene_trees_test/"
address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # "/rhog_size_g2_s500/" sample_rootHOG
species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"
pickle_address = working_folder + "/pickle_folder/"

pickle_files_adress = listdir(pickle_address)

hogs_a_rhog_xml_all = []
for pickle_file_adress in pickle_files_adress:
    with open(pickle_address+ pickle_file_adress, 'rb') as handle:

        hogs_a_rhog_xml_all += dill_pickle.load(handle)

print(hogs_a_rhog_xml_all)
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
# # from os import listdir
# # from Bio import SeqIO
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
