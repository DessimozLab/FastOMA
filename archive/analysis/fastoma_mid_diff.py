


def prepare_species_tree(rhog_i, species_tree):
    for rec in rhog_i:
        prot_name = rec.name  # 'tr|E3JPS4|E3JPS4_PUCGT
        species_name = prot_name.split("|")[-1].split("_")[-1]




def infer_HOG_rhog3(rhogid_num_list, gene_id_name):  # , address_rhogs_folder, species_tree_address):
    """
    The prot sequences of orthoxml_to_newick.py rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
    we want to infer all subHOGs of this rootHOG for different taxanomic levels.

    output: orthoxml_to_newick.py python dict (HOG_thisLevel):  key=taxanomic level, value= orthoxml_to_newick.py list of subHOGs.
    """

    HOG_thisLevel_list = []
    len_HOG_thisLevel_list = []
    HOG_thisLevel_xml_all = []
    rhogid_num = 0
    for rhogid_num in rhogid_num_list:
        # rhogid_num = rhogid_num_list[rhogid_num_i]
        logger_hog.info(
            "\n" + "=" * 50 + "\n" + "Working on root hog: " + str(rhogid_num) + ". \n")  # +", ",rhogid_num_i,"-th. \n"
        prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        logger_hog.info("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")

        (species_tree) = read_species_tree(species_tree_address)
        (species_tree, species_names_rhog, prot_names_rhog) = prepare_species_tree(rhog_i, species_tree)
        # species_tree.write()  print(species_tree.write())

        dic_sub_hogs = {}
        # finding hogs at each level of species tree (from leaves to root, bottom up)
        for node_species_tree in species_tree.traverse(strategy="postorder"):
            if node_species_tree.is_leaf():
                # each leaf itself is orthoxml_to_newick.py subhog
                continue
            logger_hog.info("\n" + "*" * 15 + "\n" + "Finding hogs for the taxonomic level:" + str(
                node_species_tree.name) + "\n" + str(node_species_tree.write()) + "\n")
            dic_sub_msas = []
            (dic_sub_hogs) = infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
                                                 rhogid_num)
            HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
            logger_hog.info("subHOGs in thisLevel are " + ' '.join(["[" + str(i) + "]" for i in HOG_thisLevel]) + " .")

            for hog_i in HOG_thisLevel:
                print(hog_i)
                if len(hog_i._members) > 1:
                    # could be improved
                    HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
                    HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
                    # groups_xml.append(HOG_thisLevel_xml)
                    # print(hog_i._members)
        # HOG_thisLevel_list.append(HOG_thisLevel)
        del dic_sub_hogs
        del HOG_thisLevel
        gc.collect()

    with open(f'{address_pickles_folder}file_' + str(rhogid_num) + '.pickle', 'wb') as handle:
        pickle.dump(HOG_thisLevel_xml_all, handle, protocol=pickle.HIGHEST_PROTOCOL)

    num_hog = len(HOG_thisLevel_xml_all)

    del HOG_thisLevel_xml_all
    gc.collect()

    return (num_hog)


def infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num):
    # during parralleization, there will be orthoxml_to_newick.py problem, few times it wants to creat the folder
    # if not os.path.exists(gene_trees_folder) :
    #    os.mkdir(gene_trees_folder)
    #  File "code7d_4.py", line 470, in infer_HOG_thisLevel
    # os.mkdir(gene_trees_folder)
    # FileExistsError: [Errno 17] File exists: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2//gene_trees_test_7d_3/'

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
        merged_msa = merge_msa(sub_msa_list_lowerLevel_ready)
        logger_hog.info("All subHOGs are merged, merged msa is with length of" + str(len(merged_msa)) + " " + str(
            len(merged_msa[0])) + ".")

        gene_tree_file_addr = gene_trees_folder + "/tree_" + str(rhogid_num) + "_" + str(
            node_species_tree.name) + ".nwk"
        gene_tree_raw = infer_gene_tree(merged_msa, gene_tree_file_addr)
        gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
        logger_hog.info("Gene tree is infered with length of " + str(len(gene_tree)) + ".")
        # gene_tree_i +=1

        R = gene_tree.get_midpoint_outgroup()
        gene_tree.set_outgroup(R)

        gene_tree = lable_SD_internal_nodes(gene_tree)
        # print("Overlap speciation is done for internal nodes of gene tree, as following:")
        print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")

        # merge_subhogs

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


class HOG:
    _hogid_iter = 1000

            # print("we are here   ********???--??? ",self._hogid)
            list_member_first = list(self._members)[0]
            geneRef_elemnt = ET.Element('geneRef', attrib={
                'id': str(gene_id_name[list_member_first])})  # # gene_id_name[query_prot_record.id]
            # hog_elemnt.append(geneRef_elemnt)
            # could be improved when the rhog contains only one protein
            return geneRef_elemnt  # hog_elemnt



def prepare_xml(rhogid_num_list_temp):
    species_prot_dic = {}
    # all_prot_temp_list= []
    for rhogid_num in rhogid_num_list_temp:
        prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        for prot_i in rhog_i:
            species_i = prot_i.id.split("|")[-1].split("_")[-1]
            if species_i in species_prot_dic:
                species_prot_dic[species_i].append(prot_i.id)
            else:
                species_prot_dic[species_i] = [prot_i.id]
            # all_prot_temp_list.append(prot_i.id)

    print("there are species ", len(species_prot_dic))
    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #
    gene_counter = 100000
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
            gene_counter += 1

    groups_xml = ET.SubElement(orthoxml_file, "groups")

    return (groups_xml, gene_id_name, orthoxml_file)


def distribute_rhogs(rhogs: List[Tuple[str, int]], start_index: int, n_workers: int) -> List[Tuple[str, int]]:
    """Distribute rhogs to workers in orthoxml_to_newick.py way that each worker gets both big and small rhogs.
    inputs:
        rhogs: orthoxml_to_newick.py list of tuples of (rhogs name, rhogs size in bytes)
        start_index: the index of the worker
        n_workers: the number of workers
    outputs:
        orthoxml_to_newick.py list of tuples of (rhogs name, rhogs size in bytes)"""
    if start_index >= n_workers:
        raise ValueError("index out of range")
    else:
        temp = []
        for i in range(start_index, len(rhogs), n_workers):
            temp.append(rhogs[i])
        return temp


if __name__ == "__main__":

    logging.basicConfig()
    logger_hog = logging.getLogger("hog")
    logger_hog.setLevel(logging.INFO)  # WARN
    # make sure addresses end with "/"
    address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/ayazdiza/fastoma-dask/"
    address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/ayazdiza/family-score-adjusted/AdjustedFamilyScore_All_rHOGs/"
    address_pickles_folder = "/work/FAC/FBM/DBC/cdessim2/default/ayazdiza/fastoma_repo/temp_results/pickles/mid_adjustedfamily_40_1/"
    species_tree_address = address_working_folder + "lineage_tree_qfo.phyloxml"
    gene_trees_folder = address_working_folder + "/gene_trees_test_mid/"
    address_logs_folder = "/work/FAC/FBM/DBC/cdessim2/default/ayazdiza/fastoma_repo/logs/"
    address_group_xml_ortho = "/work/FAC/FBM/DBC/cdessim2/default/ayazdiza/family-score-adjusted/group_xml_ortho_adjusted_family_40.pickle"
    this_worker_index = int(sys.argv[1])
    n_workers = int(sys.argv[2])

    with open(address_group_xml_ortho, 'rb') as handle:
        (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)

    ## create orthoxml_to_newick.py list of rootHOG IDs  stored in the folder of rHOG .
    rhog_files = listdir(address_rhogs_folder)[:]
    print("#", rhog_files, len(rhog_files))

    file_sizes = [os.path.getsize(f'{address_rhogs_folder}{i}') for i in rhog_files]
    file_name_size_dict = dict(zip(rhog_files, file_sizes))
    file_name_size_list_sorted = sorted(file_name_size_dict.items(), key=lambda x: x[1])
    this_worker_files = distribute_rhogs(file_name_size_list_sorted, this_worker_index, n_workers)
    print("# This worker total rhog size (bytes):", sum([i[1] for i in this_worker_files]))
    this_worker_files = [i[0] for i in this_worker_files]

    rhogid_num_list = []
    for rhog_file in this_worker_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)
    print("##", rhogid_num_list, len(rhogid_num_list))

    print(infer_HOG_rhog3(rhogid_num_list, gene_id_name))