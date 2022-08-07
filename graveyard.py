
# if format_prot_name == 1:  # qfo dataset
#     prot_name = prot_i.name  # 'tr|E3JPS4|E3JPS4_PUCGT
#     species_i = prot_name.split("|")[-1].split("_")[-1].strip()
#     if species_i == 'RAT': species_i = "RATNO"
# elif format_prot_name == 0:  # bird dataset      # rec.name  CLIRXF_R07389
#     prot_descrip = prot_i.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
#     species_i = prot_descrip.split(" ")[1].split("|")[-1]


def gene_num_convertor_old(rhogid_num_list_input, address_rhogs_folder, format_prot_name, rhogid_batch=1):
    species_prot_dic = {}
    rhogid_len_list = []
    for rhogid_num in rhogid_num_list_input:
        prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        rhogid_len_list.append(len(rhog_i))
        for prot_i in rhog_i:
            # qfo : >tr|A0A0N7KF21|A0A0N7KF21_ORYSJ=#ORYSJ_=#1000000344 tr|A0A0N7KF21|A0A0N7KF21_ORYSJ Os02g0264501 protein OS=Oryza sativa subsp. japonica (Rice) OX=39947 GN=Os02g0264501 PE=4 SV=1
            prot_id = prot_i.id.split("=#")
            prot_name = prot_id[2]   # for debugging  prot_id[0] readable prot name,  for xml prot_id[2]
            species_i = prot_id[1][:-1]

            if species_i in species_prot_dic:
                species_prot_dic[species_i].append(prot_i.id)
            else:
                species_prot_dic[species_i] = [prot_i.id]
            # all_prot_temp_list.append(prot_i.id)

    print("Number of species in the batch is ", len(species_prot_dic))
    gene_counter = 1000000 + rhogid_batch * 10000
    gene_id_name = {}
    query_species_names_rhogs = list(species_prot_dic.keys())
    for species_name in query_species_names_rhogs:
        prot_list = species_prot_dic[species_name]
        # for prot_itr in range(len(prot_list)):
        #     prot_i_name = prot_list[prot_itr]
        for prot_itr, prot_i_name in enumerate(prot_list):
            gene_id_name[prot_i_name] = gene_counter
            gene_counter += 1

    return gene_id_name

    # for hogs in hogs_a_rhog:
    #     if isinstance(hogs, list): # bwecause of dask future, sometimes it is  a list of list, sometimes list of hogs. need to be imporved
    #         for hog_i in hogs:
    #             print(hog_i)
    #             if len(hog_i._members) > 1:
    #                 # could be improved
    #                 hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
    #                 hogs_rhogs_xml.append(hogs_a_rhog_xml)
    #     # else:  # bwecause of dask future, sometimes it is  a list of list, sometimes list of hogs. need to be imporved
    #     #     if len(hogs._members) > 1:
    #     #         hogs_a_rhog_xml = hogs.to_orthoxml(**gene_id_name)
    #     #         hogs_rhogs_xml.append(hogs_a_rhog_xml)


# def prepare_xml_old(rhogid_num_list_input, address_rhogs_folder, format_prot_name, rhogid_batch = 1):
#     species_prot_dic = {}
#     # all_prot_temp_list= []
#     rhogid_len_list = [ ]
#     for rhogid_num in rhogid_num_list_input:
#         prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
#         rhog_i = list(SeqIO.parse(prot_address, "fasta"))
#         rhogid_len_list.append(len(rhog_i))
#
#         for prot_i in rhog_i:
#             if format_prot_name == 1:  # qfo dataset
#                 prot_name = prot_i.name  # 'tr|E3JPS4|E3JPS4_PUCGT
#                 species_i = prot_name.split("|")[-1].split("_")[-1].strip()
#                 if species_i == 'RAT': species_i = "RATNO"
#             elif format_prot_name == 0:  # bird dataset
#                 # rec.name  CLIRXF_R07389
#                 # prot_name = prot_i.name
#                 prot_descrip = prot_i.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
#                 species_i = prot_descrip.split(" ")[1].split("|")[-1]
#                 # species_name = prot_name.split("_")[0].strip()
#
#             # species_i = prot_i.id.split("|")[-1].split("_")[-1]
#             if species_i in species_prot_dic:
#                 species_prot_dic[species_i].append(prot_i.id)
#             else:
#                 species_prot_dic[species_i] = [prot_i.id]
#             # all_prot_temp_list.append(prot_i.id)
#
#     print("there are species ", len(species_prot_dic))
#     orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
#                                                    "originVersion": "Nov 2021", "version": "0.3"})  #
#     gene_counter = 1000000 + rhogid_batch * 10000
#     gene_id_name = {}
#     query_species_names_rhogs = list(species_prot_dic.keys())
#     for species_name in query_species_names_rhogs:
#         no_gene_species = True  # for code develop ment
#         species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
#         database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
#         genes_xml = ET.SubElement(database_xml, "genes")
#
#         prot_list = species_prot_dic[species_name]
#         # for prot_itr in range(len(prot_list)):  # [12:15]
#         #     prot_i_name = prot_list[prot_itr]
#         for prot_itr, prot_i_name in enumerate(prot_list):
#             gene_id_name[prot_i_name] = gene_counter
#             if "|" in prot_i_name:
#                 prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
#             else:
#                 prot_i_name_short = prot_i_name
#
#             gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
#             gene_counter += 1
#
#     groups_xml = ET.SubElement(orthoxml_file, "groups")
#
#
#     return (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list)
#


# ## the following are needed when we start from a rootHOG fasta file.
#
#
# # import pickle
# import dill as pickle
#
#
# def infer_HOG_rhog3(rhogid_num_list, gene_id_name):  # , address_rhogs_folder, species_tree_address):
#     """
#     The prot sequences of a rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
#     we want to infer all subhogs of this rootHOG for different taxanomic levels.
#
#     output: a python dict (HOG_thisLevel):  key=taxanomic level, value= a list of subhogs.
#     """
#
#     HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
#     logger_hog.info("subhogs in thisLevel are " + ' '.join(["[" + str(i) + "]" for i in HOG_thisLevel]) + " .")
#
#     for hog_i in HOG_thisLevel:
#         print(hog_i)
#         if len(hog_i._members) > 1:
#             # could be improved
#             HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
#             HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
#             # groups_xml.append(HOG_thisLevel_xml)
#             # print(hog_i._members)
#     # HOG_thisLevel_list.append(HOG_thisLevel)
#
#
#
#
# # project="project1",, queue="normal"
# # cluster = SLURMCluster(walltime='00:20:00', n_workers = NCORE, cores=NCORE,processes = NCORE,interface='ib0', memory="20GB",scheduler_options={'interface': 'ens2f0' })
# #  env_extra=['source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda/etc/profile.d/conda.sh','conda activate ML2'],
# print(cluster.job_script())
# print(cluster.dashboard_link)
#
# cluster.scale(jobs=njobs)  # # ask for one jobs
#
# import time
# time.sleep(5)
# print(cluster)
# # cluster.adapt(minimum=10, maximum=30)
#
# client = Client(cluster, timeout='1000s', set_as_default=True)
#
#
#
# len_HOG_thisLevel_all = []
#
# number_roothog = len(rhogid_num_list_temp)
#
# num_per_parralel = 18
# parralel_num = int(number_roothog / num_per_parralel)
#
# for list_idx in range(parralel_num + 1):
#
# if list_idx == parralel_num:
#     rhogid_num_list = rhogid_num_list_temp[list_idx * num_per_parralel:]
# else:
#     rhogid_num_list = rhogid_num_list_temp[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
#
# (out_len) = dask.delayed(infer_HOG_rhog3)(rhogid_num_list, gene_id_name)
# len_HOG_thisLevel_all.append(out_len)
#
# print("before computation", len(len_HOG_thisLevel_all), len_HOG_thisLevel_all[:2])
#
# output_computed = dask.compute(*len_HOG_thisLevel_all)
# print(" computation done ")
#
#
#
#
#
# # import dill as pickle
# # from os import listdir
# # import xml.etree.ElementTree as ET
# # from xml.dom import minidom
# # import os
#
# # address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
#
# # address_rhogs_folder  =address_working_folder + "/rhog_size_g2_s1k/"
#
#
# # ## create a list of rootHOG IDs  stored in the folder of rHOG .
# # rhog_files = listdir(address_rhogs_folder)
# # rhogid_num_list= []
# # for rhog_file in rhog_files:
# #     if rhog_file.split(".")[-1] == "fa":
# #         rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
# #         rhogid_num_list.append(rhogid_num)
# # print(len(rhogid_num_list)," .")
#
# # rhogid_num_list_temp = rhogid_num_list#[:200]
#
#
#

    # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
    # logger_hog.info("Number of root hogs is "+str(len(rhogid_num_list))+".")
    # print(rhogid_num_list[:2])
    #
    # rhogid_num_list_input = rhogid_num_list[:100]
    # logger_hog.info("Number of working root hog is " + str(len(rhogid_num_list_input)) + ".")
    # (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_num_list_input, address_rhogs_folder, format_prot_name )
    # # # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
    # # #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
    # len(gene_id_name)
    #
    # dask_future = False
    #
    # if dask_future:
    #     # print("*** client **** ")
    #
    #     # print(cluster.dashboard_link)
    #     # print(cluster.get_logs())
    #     ncore = 1 # Total number of cores per job
    #     njobs = 2  # Cut the job up into this many processes.
    #     # # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
    #     nproc = ncore
    #
    #     cluster = LocalCluster()
    #     # cluster = SLURMCluster(cores=ncore, processes=nproc, memory="20GB", walltime="01:00:00")
    #     cluster.scale(njobs)  # # ask for one jobs
    #     client_dask = Client(cluster)
    #
    #     # futures = client.map(score, x_values)
    #     # results = client.gather(futures)
    #     # hogs_a_rhog_xml_all = results
    #
    #     len_tresh = 1000
    #     dask_out_list =[]
    #     for rhogid_num_i in range(len(rhogid_num_list_input)):
    #         rhogid_num = rhogid_num_list_input[rhogid_num_i]
    #         rhogid_len = rhogid_len_list[rhogid_num_i]
    #         if rhogid_len < len_tresh:
    #
    #             dask_future_taxon = False
    #             vars_input = (gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future, dask_future_taxon, format_prot_name)
    #             vars_input_future = client.scatter(vars_input)
    #             dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, vars_input_future)
    #             dask_out_list.append(dask_out)
    #
    #             # dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, vars_input)
    #             # dask_out_list.append(dask_out)
    #             print("*a*" * 100)
    #         else:
    #             print("*b*" * 100)
    #             dask_future_taxon = True  # second level of parralelizion
    #             print(rhogid_num_i, rhogid_num, rhogid_len)
    #             # dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name,
    #             #                          address_rhogs_folder, species_tree_address, gene_trees_folder,
    #             #                          pickle_address, dask_future, dask_future_taxon)
    #             # dask_out_list.append(dask_out)
    #             print("here")
    #
    #     # for dask_out in dask_out_list :
    #     #     hogs_a_rhog_xml_all = dask_out.result()
    #     #     print(hogs_a_rhog_xml_all)
    #
    # else:
    #     print("*d*" * 100)
    #     dask_future_taxon = False
    #     for rhogid_num_i in range(len(rhogid_num_list_input)):
    #         rhogid_num = rhogid_num_list_input[rhogid_num_i]
    #         rhogid_len = rhogid_len_list[rhogid_num_i]
    #         if rhogid_len > 50:
    #             vars_input = (
    #             gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future,
    #             dask_future_taxon, format_prot_name)
    #             hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, vars_input)
    #             exit
    #
    #
    # # for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
    # #     groups_xml.append(hogs_a_rhog_xml)
    # # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # # print(xml_str)
    #
    #
    # print("test")

    # dask_working.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')
    # print("visualized.")
    # dask_result = dask_working.compute()
    # print(dask_result)  # prints "55"

    # dask_a = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address,
    #  gene_trees_folder)
    # print("this is dask_a: \n ",dask_a)
    # dask_a.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')

    # dask.compute()

    # hogs_rhogs_xml_all = []
    # for rhogid_num in  rhogid_num_list_input:
    #     hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder,
    #     species_tree_address, gene_trees_folder)
    #     hogs_rhogs_xml_all.append(hogs_a_rhog_xml_all)
    # print("here")

    #client = Client(processes=False)  # start local workers as processes
    # future_1 = client.submit(infer_hogs_for_a_rhog, species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    # rhogid_num, gene_trees_folder)
    # future = client.scatter(parameters)

    #dask.visualize(gene_id_name)



    #out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder)

    # (dic_sub_hogs)= future_1.result()

    # futures = client.map(inc, range(1000))
    # as completed
    # future.cancel()

    # rhogid_num_list_temp = [836500]  # rhogid_num_list[23]  # [833732]

    #     (dic_sub_hogs) = infer_hogs_for_a_rhog(species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    #                                                        rhogid_num, gene_trees_folder)


# # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
# #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
    # orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA", "originVersion": "Nov 2021", "version": "0.3"})  #
    # groups_xml = ET.SubElement(orthoxml_file, "groups")
    # for hog_xml in hogs_a_rhog_xml_all:
    #     groups_xml.append(hog_xml)
    # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str)
