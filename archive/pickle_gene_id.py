

import xml.etree.ElementTree as ET
import dill as dill_pickle
from os import listdir
from xml.dom import minidom
import os
from Bio import SeqIO
# import dill as dill_pickle
import pickle


def prepare_xml(rhogid_num_list_temp):

    species_prot_dic = {}
    # all_prot_temp_list= []
    for rhogid_num in rhogid_num_list_temp:
        prot_address = address_rhogs_folder +"HOG_B" +str(rhogid_num).zfill(7 ) +".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        for prot_i in rhog_i:
            species_i = prot_i.id.split("||")[1][:-1]  # prot_i.id.split("|")[-1].split("_")[-1]
            if species_i in species_prot_dic:
                species_prot_dic[species_i].append(prot_i.id)
            else:
                species_prot_dic[species_i] = [prot_i.id]
            # all_prot_temp_list.append(prot_i.id)

    print("there are species " ,len(species_prot_dic))
    orthoxml_file  = ET.Element("orthoXML", attrib={"xmlns":"http://orthoXML.org/2011/", "origin" :"OMA", "originVersion" :"Nov 2021", "version" :"0.3"} ) #
    gene_counter =100000
    gene_id_name = {}
    query_species_names_rHOGs = list(species_prot_dic.keys())
    for species_name in query_species_names_rHOGs:
        no_gene_species = True  # for code develop ment
        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")

        prot_list = species_prot_dic[species_name]
        for prot_itr in range(len(prot_list)):  # [12:15]
            prot_i_name = prot_list[prot_itr]
            gene_id_name[prot_i_name] = gene_counter
            prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})

    groups_xml = ET.SubElement(orthoxml_file, "groups")

    return (groups_xml, gene_id_name, orthoxml_file)

address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/ali_code_31aug/"

address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/rhog_all_v3_g2_s500/"
#
# address_pickles_folder = address_working_folder + "pickle/"
#
# species_tree_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/lineage_tree_qfo.phyloxml"
#
# gene_trees_folder = address_working_folder + "/gene_trees_test_mid/"
#
# address_logs_folder = address_working_folder + "logs/"
# address_group_xml_ortho = address_working_folder + "group_xml_ortho_adjusted_family_40.pickle"

rhog_files = listdir(address_rhogs_folder)[:]

rhog_files = listdir(address_rhogs_folder)
rhogid_num_list = []
for rhog_file in rhog_files:
    if rhog_file.split(".")[-1] == "fa":
        rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
        rhogid_num_list.append(rhogid_num)

rhogid_num_list_temp = rhogid_num_list

(groups_xml, gene_id_name, orthoxml_file) = prepare_xml(rhogid_num_list_temp)
address_group_xml_ortho = address_working_folder+"group_xml_ortho_adjusted_family_40_newrhg.pickle"


with open(address_group_xml_ortho, 'wb') as handle:
    # dill_pickle.dump(gene_id_name, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
    pickle.dump((groups_xml, gene_id_name, orthoxml_file) , handle, protocol=pickle.HIGHEST_PROTOCOL)
