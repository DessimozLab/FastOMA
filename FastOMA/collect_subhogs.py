
import xml.etree.ElementTree as ET
import pickle
from os import listdir
from xml.dom import minidom
# from . import _config
from ._config import logger_hog

def collect_subhogs():

    logger_hog.info("started collecting pickle files ")

    # todo as input argument/option in nextflow
    protein_format_qfo_dataset_before2022 = False
    # in benchamrk dataset the output prot names should be short
    # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ
    # for qfo benchamrk, the middle should be wirtten in the file

    pickle_folder = "./" # pickle_rhogs
    output_xml_name = "./output_hog_.orthoxml"
    gene_id_pickle_file = "./gene_id_dic_xml.pickle"


    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
    logger_hog.debug("We read the gene_id_name dictionary with"+str(len(gene_id_name))+"items")
    logger_hog.debug("Now creating the header of orthoxml")

   #  #### create the header of orthoxml ####
    for query_species_name, list_prots in gene_id_name.items():
        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "database ", "version": "2023"})
        genes_xml = ET.SubElement(database_xml, "genes")
        for (gene_idx_integer, query_prot_name) in list_prots:
            if protein_format_qfo_dataset_before2022:
                # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ   for qfo benchamrk, the middle should be wirtten in the file
                query_prot_name_pure = query_prot_name.split("|")[1]
            else:
                query_prot_name_pure = query_prot_name
            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
    logger_hog.debug("gene_xml is created.")

    #  #### create the groups of orthoxml   ####
    pickle_files_adress_raw = listdir(pickle_folder)
    pickle_files_adress = [i for i in pickle_files_adress_raw if i.endswith(".pickle") and i.startswith("file_")]

    logger_hog.info("number of pickle files is "+str(len(pickle_files_adress))+".")
    logger_hog.info("pickle files are " + str(pickle_files_adress) + ".")
    hogs_a_rhog_xml_all = []
    for pickle_file_adress in pickle_files_adress:
        with open(pickle_folder + pickle_file_adress, 'rb') as handle:
            hogs_a_rhog_xml_batch = pickle.load(handle)  #  list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
            # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.

    logger_hog.debug("number of hogs in all batches is "+str(len(hogs_a_rhog_xml_all))+" .")
    groups_xml = ET.SubElement(orthoxml_file, "groups")

    for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
        groups_xml.append(hogs_a_rhog_xml)
    logger_hog.debug("converting the xml object to string.")

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str[:-1000])
    logger_hog.debug("writing orthoxml started")
    with open(output_xml_name, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    logger_hog.info("orthoxml is written in " + output_xml_name)

