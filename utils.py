

from os import listdir

from ete3 import Tree


import logging

logging.basicConfig()
logger_hog = logging.getLogger("hog")
logger_hog.setLevel(logging.INFO)  # WARN


from  hog_class import HOG
import wrappers


def list_rhog_fastas(address_rhogs_folder):
    """
     create a list of rootHOG IDs  stored in the folder of rHOG .
     input: folder address
     output: list of rhog Id (integer)
    """
    rhog_files = listdir(address_rhogs_folder)
    rhogid_num_list= []
    for rhog_file in rhog_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)

    return rhogid_num_list



import logging
from ete3 import Phyloxml


def read_species_tree(species_tree_address):
    """
    reading a species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
    project = Phyloxml()
    project.build_from_file(species_tree_address)
    # Each tree contains the same methods as a PhyloTree object
    for species_tree in project.get_phylogeny():
        species_tree = species_tree

    for node_species_tree in species_tree.traverse(strategy="postorder"):
        if node_species_tree.is_leaf():
            temp1 = node_species_tree.phyloxml_clade.get_taxonomy()[0]
            # print(temp1.get_code())
            node_species_tree.name = temp1.get_code()
    # print(len(species_tree)); print(species_tree)


    return (species_tree)


def prepare_species_tree(rhog_i, species_tree):
    """
    a function for extracting a subtree from the input species tree  a.k.a pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    species_names_rhog = []
    prot_names_rhog = []
    for rec in rhog_i:
        prot_name = rec.name  # 'tr|E3JPS4|E3JPS4_PUCGT
        # prot_name = prot_name_full.split("|")[1].strip() # # 'tr|E3JPS4|E3JPS4_PUCGT
        species_name = prot_name.split("|")[-1].split("_")[-1]
        if species_name == 'RAT': species_name = "RATNO"
        species_names_rhog.append(species_name)
        prot_names_rhog.append(prot_name)

    species_names_uniqe = set(species_names_rhog)
    logger_hog.info("The number of unique species in the rHOG is " + str(len(species_names_uniqe)) + ".")
    species_tree.prune(species_names_uniqe, preserve_branch_length=True)
    # species_tree.write()
    for node in species_tree.traverse(strategy="postorder"):
        node_name = node.name
        num_leaves_no_name=0
        if len(node_name) < 1:
            if node.is_leaf():
                node.name = "leaf_" + str(num_leaves_no_name)
            else:
                node_children = node.children
                list_children_names = [node_child.name for node_child in node_children]
                node.name = '_'.join(list_children_names)
    print("Working on the following species tree.")
    print(species_tree)

    return (species_tree, species_names_rhog, prot_names_rhog)


def lable_SD_internal_nodes(tree_out):
    """
    for the input gene tree, run the species overlap method
    and label internal nodes of the gene tree

    output: labeled gene tree
    """
    species_name_dic = {}
    counter_S = 0
    counter_D = 0

    for node in tree_out.traverse(strategy="postorder"):
        # print("** now working on node ",node.name) # node_children
        if node.is_leaf():
            prot_i = node.name
            species_name_dic[node] = {str(prot_i).split("|")[-1].split("_")[-1]}
        else:
            node.name = "S/D"
            leaves_list = node.get_leaves()  # print("leaves_list", leaves_list)
            species_name_set = set([str(prot_i).split("|")[-1].split("_")[-1] for prot_i in leaves_list])
            # print("species_name_set", species_name_set)
            species_name_dic[node] = species_name_set

            node_children = node.children  # print(node_children)
            node_children_species_list = [species_name_dic[node_child] for node_child in node_children]  # list of sets
            # print("node_children_species_list", node_children_species_list)
            node_children_species_intersection = set.intersection(*node_children_species_list)

            if node_children_species_intersection:  # print("node_children_species_list",node_children_species_list)
                counter_D += 1
                node.name = "D" + str(counter_D)
            else:
                counter_S += 1
                node.name = "S" + str(counter_S)
    return tree_out


def infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num, gene_trees_folder):

    # during parralleization, there will be a problem, few times it wants to creat the folder
    # if not os.path.exists(gene_trees_folder) :
    #    os.mkdir(gene_trees_folder)
    #  File "code7d_4.py", line 470, in infer_HOG_thisLevel
    # os.mkdir(gene_trees_folder)
    # FileExistsError: [Errno 17] File exists: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2//gene_trees_test_7d_3/'

    logger_hog.info(
        "\n" + "*" * 15 + "\n" + "Finding hogs for the taxonomic level:" + str(node_species_tree.name) + "\n" + str(
            node_species_tree.write()) + "\n")

    if len(rhog_i) == 0:
        logger_hog.warning('There is no protein in the rHOG: ' + str(rhogid_num))
        dic_sub_hogs[node_species_tree.name] = []
        return (dic_sub_hogs)

    elif len(rhog_i) == 1:
        logger_hog.warning('There is only one protein in the rHOG: ' + str(rhogid_num))
        node_species_name = node_species_tree.children[0].name  # there is only one species (for the one protein)
        prot = rhog_i[0]
        sub_hog_leaf = HOG(prot, node_species_name, rhogid_num)
        subHOGs_children = [sub_hog_leaf]
        HOG_this_level = subHOGs_children
        dic_sub_hogs[node_species_tree.name] = HOG_this_level
        return (dic_sub_hogs)

    sub_msa_list_lowerLevel = []  # including subHOGS of lower level
    subHOGs_children = []

    # print("working on node", node_species_tree.name,"with",len(node_species_tree.children),"children.")
    for node_child in node_species_tree.children:
        if node_child.is_leaf():
            node_species_name = node_child.name
            # extracting those proteins of the rHOG that belongs to this species (child node of species tree)
            interest_list = [idx for idx in range(len(species_names_rhog)) if
                             species_names_rhog[idx] == node_species_name]
            rhog_part = [rhog_i[i] for i in interest_list]
            # sub_msa = [MultipleSeqAlignment([i]) for i in rhog_part]             #print("len",len(rhog_part))

            for prot in rhog_part:
                sub_hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
                # list_all_hogs_ever.append(sub_hog_leaf)
                subHOGs_children.append(sub_hog_leaf)
        else:  # the child node is an internal node, subHOGs are inferred till now during traversing.
            # print("sub msa for internal node", node_child.name,"is read from dic.")
            if node_child.name in dic_sub_hogs:
                sub_hogs_child = dic_sub_hogs[node_child.name]
                subHOGs_children += sub_hogs_child
            else:
                logger_hog.error("Error 131, no sub msa for the internal node ", node_child.name, node_child, "\n",
                                 dic_sub_hogs)
                assert 2 == 1
    temp11 = []
    for temp in [i._members for i in subHOGs_children]:
        temp11.append([prot.split('|')[2] for prot in temp])
    # print("there are ",len(subHOGs_children),"subHOGs lower of this level:",[i._hogid for i in subHOGs_children],temp11)
    # print("We want to infer subHOGs at this level,i.e. merge few of them.")
    subHOG_to_be_merged_set_other_Snodes = []

    if len(subHOGs_children) == 0:
        logger_hog.error('Error 139, There is no protein in this subhog, for rhog' + str(rhogid_num))

    elif len(subHOGs_children) == 1:
        HOG_this_level = subHOGs_children
        # print("**** error 134 *** ", len(subHOGs_children),subHOGs_children) #return (-1,-1,-1)

    else:
        sub_msa_list_lowerLevel_ready = [hog._msa for hog in subHOGs_children]
        merged_msa = wrappers.merge_msa(sub_msa_list_lowerLevel_ready)
        logger_hog.info("All subHOGs are merged, merged msa is with length of" + str(len(merged_msa)) + " " + str(
            len(merged_msa[0])) + ".")

        gene_tree_file_addr = gene_trees_folder + "/tree_" + str(rhogid_num) + "_" + str(
            node_species_tree.name) + ".nwk"
        gene_tree_raw = wrappers.infer_gene_tree(merged_msa, gene_tree_file_addr)
        gene_tree = Tree(gene_tree_raw + ";", format=0)
        logger_hog.info("Gene tree is infered with length of " + str(len(gene_tree)) + ".")
        # gene_tree_i +=1
        R = gene_tree.get_midpoint_outgroup()
        gene_tree.set_outgroup(R)  # print("Midpoint rooting is done for gene tree.")
        gene_tree = lable_SD_internal_nodes(gene_tree)
        # print("Overlap speciation is done for internal nodes of gene tree, as following:")
        print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")

        tree_leaves = [i.name for i in gene_tree.get_leaves()]
        # assigned_leaves_to_hog = []        #sub_msas_list_this_level = []
        subHOGs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
        HOG_this_level = []
        subHOG_to_be_merged_set_other_Snodes = []
        subHOG_to_be_merged_set_other_Snodes_flattned_temp = []
        for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n,
                                                                                         "processed") and n.processed == True):  # start from root
            # print("Leaves assigned to hog are ", assigned_leaves_to_hog)   #print("Traversing gene tree. Now at node", node.name)
            if not node.is_leaf():
                node_leaves_name = [i.name for i in node.get_leaves()]
                # print(node_leaves_name)

                if node.name[0] == "S":  # this is a sub-hog.
                    subHOG_to_be_merged = []
                    for node_leave_name in node_leaves_name:  # print(node_leave_name)
                        for subHOG in subHOGs_children:
                            subHOG_members = subHOG._members
                            if node_leave_name in subHOG_members:  # could be improved
                                if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp:
                                    subHOG_to_be_merged.append(subHOG)
                                    subHOGs_id_children_assigned.append(subHOG._hogid)
                                else:
                                    print("issue 184", node.name, subHOG._hogid, node_leave_name)
                                    if "processed" in node:
                                        print(node.name)
                                    else:
                                        print("processed not in ",
                                              node.name)  # print(node_leave_name,"is in ",subHOG._hogid)
                    if subHOG_to_be_merged:
                        subHOG_to_be_merged_set = set(subHOG_to_be_merged)
                        taxnomic_range = node_species_tree.name
                        HOG_this_node = HOG(subHOG_to_be_merged_set, taxnomic_range, rhogid_num, msa=merged_msa)
                        HOG_this_level.append(HOG_this_node)
                        subHOG_to_be_merged_set_other_Snodes.append([i._hogid for i in subHOG_to_be_merged_set])
                        subHOG_to_be_merged_set_other_Snodes_flattned_temp = [item for items in
                                                                              subHOG_to_be_merged_set_other_Snodes for
                                                                              item in items]
                        #  I don't need to traverse deeper in this clade
                    node.processed = True  # print("?*?*  ", node.name)

            subHOG_to_be_merged_set_other_Snodes_flattned = [item for items in subHOG_to_be_merged_set_other_Snodes for
                                                             item in items]
            if [i._hogid for i in subHOGs_children] == subHOG_to_be_merged_set_other_Snodes_flattned:
                break
        for subHOG in subHOGs_children:  # for the single branch  ( D include a  subhog and a S node. )
            if subHOG._hogid not in subHOGs_id_children_assigned:  # print("here", subHOG)
                HOG_this_level.append(subHOG)
        prot_list_sbuhog = [i._members for i in HOG_this_level]
        prot_list_sbuhog_short = []
        for prot_sub_list_sbuhog in prot_list_sbuhog:
            prot_list_sbuhog_short.append([prot.split('|')[2] for prot in prot_sub_list_sbuhog])
        logger_hog.info("- " + str(
            len(prot_list_sbuhog_short)) + "HOGs are inferred at the level " + node_species_tree.name + ": " + " ".join(
            [str(i) for i in prot_list_sbuhog_short]))
    # print("By merging ",subHOG_to_be_merged_set_other_Snodes)

    # check for conflicts in merging
    #     for i in range(subHOG_to_be_merged_set_other_Snodes):
    #         if
    #         for i in range(subHOG_to_be_merged_set_other_Snodes):
    # print("*&*& ",node_species_tree.name)
    dic_sub_hogs[node_species_tree.name] = HOG_this_level
    return (dic_sub_hogs)
