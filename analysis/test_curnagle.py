
import xml.etree.ElementTree as ET
# import dill as dill_pickle
import pickle
from os import listdir
from xml.dom import minidom
import os

print("started ")
working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
pickle_folder = working_folder + "/pickle_folder_3sep10pm_allv3_new/"


output_xml_name = "out_3sep_all_new_6.xml"
gene_id_pickle_file = working_folder + "gene_id_28aug.pickle"


if os.path.getsize(gene_id_pickle_file) > 0:
    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)   # (groups_xml, gene_id_name, orthoxml_file)

print("gene_id_name read ")


orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                "originVersion": "Nov 2021", "version": "0.3"})  #


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
print("Following files are with 0 byte (if any, so hog inferecnce should be repeated after deleting the files):")
for pickle_file_adress in pickle_files_adress:
    if os.path.getsize(pickle_folder + pickle_file_adress) > 0:
        with open(pickle_folder + pickle_file_adress, 'rb') as handle:
            hogs_a_rhog_xml_batch = pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
    else:
        print(pickle_file_adress)

print("after reading number of hogs is  ", len(hogs_a_rhog_xml_all))



groups_xml = ET.SubElement(orthoxml_file, "groups")

for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
    groups_xml.append(hogs_a_rhog_xml)

xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")


with open(working_folder +output_xml_name, "w") as file_xml:
    file_xml.write(xml_str)
file_xml.close()

print("orthoxml is written in  "+ working_folder +output_xml_name)



#


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
#         orthoxml_to_newick.py=1
#
#
#
#
#     rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
#
#
#
#

