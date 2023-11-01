
import xml.etree.ElementTree as ET
import pickle
from os import listdir
from xml.dom import minidom
# from . import _config
from ._config import logger_hog

from FastOMA.zoo.hog import extract_flat_groups_at_level

from ete3 import Tree
# import sys
import os
from FastOMA.zoo.hog.convert import orthoxml_to_newick
from Bio import SeqIO

# This code collect subhogs and writes outputs.


def collect_subhogs():

    logger_hog.info("started collecting pickle files ")

    # todo as input argument/option in nextflow
    protein_format_qfo_dataset_before2022 = True #False
    # in benchamrk dataset the output prot names should be short
    # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ
    # for qfo benchmark, the middle should be written in the file

    pickle_folder = "./pickles_temp/" #pickle_rhogs
    output_xml_name = "./output_hog.orthoxml"
    gene_id_pickle_file = "./gene_id_dic_xml.pickle" # this file includes integer mapping of each protein


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
    logger_hog.debug("pickle files are " + str(len(pickle_files_adress)) + ".")
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







    logger_hog.info("Now writing OG fasta files ")


    def max_og_tree(tree):
        for node in tree.traverse("preorder"):
            # for node in xml_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "attriremoved") and n.attriremoved==True):
            if not node.is_leaf() and hasattr(node, "Ev") and node.Ev == 'duplication':  # node.name[:3] == "dup"
                dup_node = node
                children = dup_node.get_children()
                list_num_species = []
                for child in children:
                    child_name_leaves = child.get_leaves()
                    species_list = []
                    for leaf in child_name_leaves:
                        name = leaf.name
                        if name[:3] == "no_":
                            name = leaf.name.split("_")[-1]
                        if name in species_dic:
                            species_name = species_dic[name]
                            species_list.append(species_name)
                        else:
                            print("species not in the dic ", name)
                    species_set = set(species_list)
                    list_num_species.append(len(species_set))
                index_max_species = list_num_species.index(max(list_num_species))
                # if there are few children with identical number of species, the case would be not a polytomi but two children with one species
                # num_occurence = [1 for i in list_num_species if i == max(list_num_species)]
                # if len(num_occurence) > 1:
                #    print("please check this case with the developer the tool. The tree has polytomy.")
                child_max_species = children[index_max_species]
                children_to_remove = [i for i in children if i != child_max_species]
                for child_to_remove in children_to_remove:
                    for i in child_to_remove.get_leaves():
                        i.in_og = "no"

        og_prot_list = []
        for node in tree.traverse("preorder"):
            if node.is_leaf():
                if hasattr(node, "in_og") and node.in_og == "no":
                    pass  # print(node.name)
                else:
                    og_prot_list.append(node.name)

        return og_prot_list

    input_orthoxml = output_xml_name # sys.argv[1]  # "out_folder/output_hog_.orthoxml"
    rhog_all_folder = "./omamer_rhogs/" #sys.argv[2] + "/"  # "out_folder/rhogs_all/"
    fasta_format = "fa"  # of the rhogs

    output_file_og_tsv = "OrthologousGroups.tsv"

    trees, species_dic = orthoxml_to_newick(input_orthoxml, return_gene_to_species=True)  # encode_levels_as_nhx=False,  xref_tag="protId",
    print("We extracted " + str(len(trees)) + " trees  in NHX format from the input HOG orthoxml" + input_orthoxml)

    OGs = {}
    for hog_id, tree_string in trees.items():
        tree = Tree(tree_string, format=1)
        og_prot_list = max_og_tree(tree)
        if len(og_prot_list) >= 2: # a group should have at least 1 member/protein
            OGs[hog_id] = og_prot_list

    print("done")

    with open(output_file_og_tsv, 'w') as handle:
        for hog_id, og_prot_list in OGs.items():
            line_text = str(hog_id) + "\t" + str(og_prot_list)[1:-1] + "\n"
            handle.write(line_text)
    handle.close()

    print("We wrote the protein families information in the file " + output_file_og_tsv)

    out_folder_ogs = "OrthologousGroupsFasta/"
    os.makedirs(out_folder_ogs)

    print("start writing " + str(len(OGs)) + " OGs as fasta files in folder " + out_folder_ogs)
    for hog_id, og_prot_list in OGs.items():  # hog_id="HOG_0667494_sub10524"
        rhog_id = "_".join(hog_id.split("_")[:2])

        omamer_rhogs_all_address = rhog_all_folder + rhog_id + "." + fasta_format
        omamer_rhogs_all_prots = list(SeqIO.parse(omamer_rhogs_all_address, "fasta"))

        og_prots = []
        og_prot_list = OGs[hog_id]
        for rhogs_prot in omamer_rhogs_all_prots:
            if rhogs_prot.id.split("||")[0] in og_prot_list:
                sp = rhogs_prot.id.split("||")[1]
                rhogs_prot.description += " [" + sp + "]"
                og_prots.append(rhogs_prot)

        og_id = "OG_" + hog_id  # one OG per rootHOG      # "/HOG_"+ rhogid
        SeqIO.write(og_prots, out_folder_ogs + og_id + ".fa", "fasta")
    print("writing done")



    # import sys
    # input_orthoxml = output_xml_name
    output_file = "rootHOGs.tsv"

    toplevel_groups = []
    for grp in extract_flat_groups_at_level(input_orthoxml):
        toplevel_groups.append(set(g.xref for g in grp))

    # toplevel_groups is a list of sets

    print("We extracted "+str(len(toplevel_groups))+" protein families from the input HOG orthoxml"+input_orthoxml)
    print("The first one contain "+str(len(toplevel_groups[0]))+" proteins.")

    with open(output_file, 'w') as handle:
        for toplevel_group_idx, toplevel_group in enumerate(toplevel_groups):
            line_text = str(toplevel_group_idx)+"\t"+str(toplevel_group)[1:-1]+"\n"
            handle.write(line_text)
    handle.close()

    print("We wrote the protein families information in the file "+output_file)


