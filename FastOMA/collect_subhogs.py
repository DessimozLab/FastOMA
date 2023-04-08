
import xml.etree.ElementTree as ET
import pickle
from os import listdir
from xml.dom import minidom
# from . import _config



def collect_subhogs():


    print("started ")

    qfo_analysis = False
    if qfo_analysis:
        qfo_bird = "qfo3"
        # in benchamrk dataset the output prot names should be short
        # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ   for qfo benchamrk, the middle should be wirtten in the file
    else:
        qfo_bird = "bird_hog" # for bird dataset

    # in_folder = _config.in_folder

    # python ${FastOMA} / collect_subhogs.py ${pickle_rhogs} ${gene_id_dic_xml}


    if qfo_bird == "bird_hog":
        protein_format_qfo_dataset = False
    else:
        protein_format_qfo_dataset = True

    pickle_folder = "./pickle_rhogs/"   #  "pi_rest_rhog/" #"/pickle_b_0.5_3000/" pickles_rhog
    # add warning when pickle folder is not empty
    output_xml_name = "./output_hog_.orthoxml"
    gene_id_pickle_file = "./gene_id_dic_xml.pickle"




    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
    print("gene_id_name read ")


    for query_species_name, list_prots in gene_id_name.items():

        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")


        if protein_format_qfo_dataset:
            for (gene_idx_integer, query_prot_name) in list_prots:
                # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ   for qfo benchamrk, the middle should be wirtten in the file
                query_prot_name_pure = query_prot_name.split("|")[1]
                gene_xml = ET.SubElement(genes_xml, "gene",
                                         attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
        else:
            for (gene_idx_integer, query_prot_name) in list_prots:
                query_prot_name_pure = query_prot_name
                gene_xml = ET.SubElement(genes_xml, "gene",
                                         attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})

    print("gene_xml created ")
    pickle_files_adress_raw = listdir(pickle_folder)
    pickle_files_adress = [i for i in  pickle_files_adress_raw if i.endswith(".pickle")]

    print(len(pickle_files_adress))
    hogs_a_rhog_xml_all = []
    for pickle_file_adress in pickle_files_adress:
        with open(pickle_folder + pickle_file_adress, 'rb') as handle:
            hogs_a_rhog_xml_batch = pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
            # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.

    print("number of hogs in all batches is ", len(hogs_a_rhog_xml_all))

    groups_xml = ET.SubElement(orthoxml_file, "groups")

    for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
        groups_xml.append(hogs_a_rhog_xml)
    print("convert to string")

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str[:-1000])
    print("writing")
    with open(output_xml_name, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    print("orthoxml is written in  "+ output_xml_name)


