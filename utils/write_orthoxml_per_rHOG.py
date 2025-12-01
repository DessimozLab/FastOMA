import xml.etree.ElementTree as ET

from os import listdir
from xml.dom import minidom

import pickle



folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run_1june/"

# create this folder /out_folder/orthoxml_out/

gene_id_pickle_file = folder + "/out_folder/gene_id_dic_xml.pickle"

with open(gene_id_pickle_file, 'rb') as handle:
    gene_id_name = pickle.load(handle)
    # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
print("gene_id_name read ", len(gene_id_name))

pickle_folder = folder + "/out_folder/pickle_rhogs_/"
pickle_files_adress = listdir(pickle_folder)

orthoxml_out_folder = folder + "/out_folder/orthoxml_out/"
check = listdir(orthoxml_out_folder)


print("gene_xml created ")
# hogs_a_rhog_xml_all = []
for idx, pickle_file_adress in enumerate(pickle_files_adress):

    if idx % 100 == 0: print(idx)
    with open(pickle_folder + pickle_file_adress, 'rb') as handle:
        hogs_a_rhog_xml_batch = pickle.load(
            handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
    handle.close()
    # hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
    # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.
    # print("number of hogs is batch is ", len(hogs_a_rhog_xml_batch))

    xml_str = ""
    for i in hogs_a_rhog_xml_batch:
        xml_str += minidom.parseString(ET.tostring(i)).toprettyxml(indent="   ")
    xs = xml_str.split("\n")
    list_geneid = []
    for x in xs:
        if "geneRef id" in x:
            list_geneid.append(int(x.split("\"")[1]))
    print(len(list_geneid))

    query_species_name_list = []
    for query_species_name, list_prots in gene_id_name.items():
        for gene in list_prots:
            if gene.numeric_id in list_geneid:
                query_species_name_list.append(query_species_name)
                break  # early abort as we found one protein for this species

    query_species_name_set = list(set(query_species_name_list))

    output_xml_name = orthoxml_out_folder + pickle_files_adress[0] + "_.orthoxml"
    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    for query_species_name, list_prots in gene_id_name.items():
        if query_species_name in query_species_name_set:
            species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
            database_xml = ET.SubElement(species_xml, "database", attrib={"name": " database ", "version": "2020"})
            genes_xml = ET.SubElement(database_xml, "genes")

            for gene in list_prots:
                if gene.numeric_id in list_geneid:  # +[1007003758]
                    attribs = {"id": str(gene.numeric_id), "protId": gene.prot_id}
                    if gene.main_isoform is not None:
                        attribs["main_isoform"] = str(gene.main_isoform)
                    gene_xml = ET.SubElement(genes_xml, "gene", attrib=attribs)

    groups_xml = ET.SubElement(orthoxml_file, "groups")

    for hogs_a_rhog_xml in hogs_a_rhog_xml_batch:
        groups_xml.append(hogs_a_rhog_xml)
    # print("convert to string")

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")

    with open(output_xml_name, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    print("orthoxml is written in  " + output_xml_name)
