

import numpy as np
from Bio.Seq import Seq  # , UnknownSeq
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from ete3 import Tree


from . import _wrappers
from . import _utils_subhog
from ._wrappers import logger

fragment_detection_msa_merge = True
fragment_detection = True
iteration_sp_overlap= 10
#gene_trees_write_all = True


"""
This code is for detecting fragments in MSA due to poor gene annotation and used a merged version in orthology relationship "_|_" in the gen trees.
reported in orthoxml either as DubiousMergedfragment or Dubiousfragment

We also remove proteins with low SO (species_overlap) score. 
30 Apr 2022
"""

def read_msa(input_msa):
    ids = []
    seqs = []
    species = []
    coords = []
    for rec in input_msa:
        ids.append(rec.id)
        seqs.append(np.frombuffer(str(rec.seq).upper().encode('ascii'), dtype='S1'))
        species.append(rec.id.split('||')[1][:-1])

        # compute coordinates in matrix
        # todo  if all the element is --, this will be empty and causing error
        ii = np.argwhere(seqs[-1] != b'-')
        if len(ii):
            coords.append((np.min(ii), np.max(ii)))
        else:
            logger.warning("issue 1321230 all of the seq is gap.  "+str(rec))
            coords.append((0, 0))

    ids = np.array(ids)
    seqs = np.vstack(seqs)
    species = np.array(species)
    coords = np.array(coords)
    input_msa_np = (ids, seqs, species, coords)
    return input_msa_np



def compute_identity(s1, s2):
    def compute(s1, s2):
        # ignore gap columns in first sequence
        f = (s1 != b'-')
        n = f.sum()
        return ((s1[f] == s2[f]).sum() / n) if n > 0 else 1

    return compute(s1, s2), compute(s2, s1)


def split_candidates(input_msa_np, margin=0):


    (ids_all, seqs_all, species_all, coords_all) = input_msa_np
    '''
    Count candidates where it is unambiguous which fragments to merge.
    by Alex Warwick Vesztrocy
    '''
    assert 0 <= margin <= 1, "Margin must be in [0,1]"
    candidates = []

    # get species counts
    c = np.unique(species_all, return_counts=True)
    total_seqs = c[1].sum()

    # filter to those with more than 1 seq
    f = (c[1] > 1)
    for sp in c[0][f]:
        ii = np.argwhere(species_all == sp).flatten()
        # print('species', sp, margin, len(ii), total_seqs)
        x = coords_all[ii]

        order = np.argsort(x[:, 0])
        coords = x[order]
        ids = ids_all[ii[order]]
        seqs = seqs_all[ii[order]]

        #  find whether to merge for unambiguous cases
        merge = True
        ident = []
        for i in range(len(coords) - 1):
            (s1, e1) = coords[i]
            (s2, e2) = coords[i + 1]
            i_margin = min(margin * (e1 - s1),
                           margin * (e2 - s2))
            # s1-------e1
            #       s2------e2
            #  (e1-s2) < 0 no overlap; (e1-s2) > 0 overlap.
            overlap = (e1 - s2)
            # print(overlap)
            # todo  we could improve finding overlap between fragments
            if overlap < 0:
                ident.append((1, 1))
            elif overlap <= i_margin:
                # minimal overlap, compute identities
                ident.append(
                    compute_identity(
                        seqs[i][s2:e1],
                        seqs[i + 1][s2:e1],
                    )
                )
            else:
                # overlapping too much
                merge = False
                break
        if merge:
            candidates.append((tuple(ids), ident))
    return candidates


def find_prot_dubious_msa(input_msa, conf_infer_subhhogs):

    input_msa_np = read_msa(input_msa)
    candidates = split_candidates(input_msa_np, conf_infer_subhhogs.overlap_fragments)

    # candidates.append((tuple(ids), ident))
    #prot_dubious_msa_list = [(str(i[0][0]), str(i[0][1])) for i in candidates]  # [i[0] for i in candidates]
    # todo are there only two always? few fragments
    prot_dubious_msa_list = []
    seq_dubious_msa_list = []

    if candidates :
        for candidate in candidates:
            prot_dubious_msa_list.append(list(candidate[0]))

        for prots_dubious_msa in prot_dubious_msa_list:
            seqs_dubious_msa = [i for i in input_msa if i.name in prots_dubious_msa]
            seq_dubious_msa_list.append(seqs_dubious_msa)

    return prot_dubious_msa_list,  seq_dubious_msa_list  # list of pairs


def insert_dubious_prots_hog_hierarchy_toleaves(hog_host, fragment_host, fragments_list_nothost):
    for subhog in hog_host._subhogs:

        insert_dubious_prots_hog_hierarchy_toleaves(subhog, fragment_host, fragments_list_nothost)

    result_insersion = hog_host.insert_dubious_prots(fragment_host, fragments_list_nothost)

    return 1



def remove_prot_hog_hierarchy_toleaves(hog_ii, prot_to_remove):
    # import copy; a= copy.deepcopy(hog_ii)
    for subhog in hog_ii._subhogs:
        if prot_to_remove in subhog._members:
            result_removing0 = remove_prot_hog_hierarchy_toleaves(subhog, prot_to_remove)

    result_removing = hog_ii.remove_prot_from_hog(prot_to_remove)
    subhogs_old = hog_ii._subhogs
    hog_ii._subhogs = [i for i in subhogs_old if len(i._members)]  # for getting rid of  empty subhogs #  hogID=HOG_0026884_sub12868,len=0, tax_least=Ardeidae

    return result_removing



def find_prot_dubious_sd_remove(gene_tree, all_species_dubious_sd_dic):
    # todo this function need to double check with cases of with and without dubious

    #prot_dubious_sd_allspecies = []
    prot_dubious_sd_remove_list = []
    # todo not sure postorder or preorder
    for node in gene_tree.traverse(strategy="postorder"):
        # print("** now working on node ",node.name) # node_children
        if not node.is_leaf() and 'D' in node.name:
            node_name = node.name #d, intersection, union = node_name.split("_")  # if int(intersection) / int(union) < threshold_dubious_sd:
            if node_name in all_species_dubious_sd_dic: # a duplication node with low score,
                node_children = node.children
                all_species_dubious_sd = all_species_dubious_sd_dic[node_name]
                # prot_dubious_sd_ = []
                for species_dubious_sd in all_species_dubious_sd:
                    child_size = []  # gene tree is binary for fasttree
                    prot_dubious_list = []
                    for node_child in node_children:
                        list_leaves = [i.name for i in node_child.get_leaves()]
                        child_size.append(len(list_leaves))
                        dubious_children = []
                        for prot_name in list_leaves:
                            if prot_name.split("||")[1] == species_dubious_sd:
                                #prot_dubious_list.append(prot_name)
                                dubious_children.append(prot_name)
                        prot_dubious_list.append(dubious_children)

                    #try:
                    # this need change, prot_dubious_list is now a list of list  subhogs_list = [i.split("|_|")[1][:-1] for i in prot_dubious_list]  # subhog id at child level
                    #  ["'tr|Q4V8S5|Q4V8S5_DANRE||DANRE||1000020519|_|sub12962'"]
                    # the last char is ', watch out!

                    # todo check ! is it safe or not!
                    # if len(set(subhogs_list)) > 1:
                    #     # we are removing all sequences of this species on the the side of internal node (gene tree), with least leaves
                    #     child_size_min_indx = child_size.index(min(child_size))
                    #     prot_dubious_sd_remove_list.append(prot_dubious_list[child_size_min_indx])
                    #
                    # else:
                    #     logger.debug( "This species (protein from the same subhog) is safe to keep "+ str(node_name)+" "+str(species_dubious_sd))
                    #     #all of them are from the same subhog, so it doesn't matter, a duplication event doesn't affect when all are from the same subhog at children level

                    child_size_min_indx = child_size.index(min(child_size))
                    #prot_dubious_sd_remove_list.append(prot_dubious_list[child_size_min_indx]) # ["'sp|Q9PRL8|ACBP_CHICK||CHICK||1020017457|_|sub10101'"]
                    prot_dubious_sd_remove_list += prot_dubious_list[child_size_min_indx]

                    # except:
                    #    logger.warning("issue 2495869: prot_dubious_list doesnt include the hog id . so we'll  keep it" + str(prot_dubious_list)+ " " +str(gene_tree.write(format=1, format_root_node=True) ) )


    return prot_dubious_sd_remove_list




def handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic, hogs_children_level_list, conf_infer_subhhogs):
    #  prot_dubious_sd_list, node_species_tree, genetree_msa_file_addr, hogs_children_level_list

    all_species_dubious_sd_dic_updated = all_species_dubious_sd_dic
    itr_so = 1 # iteration of species overlap (so)
    while all_species_dubious_sd_dic_updated and itr_so< iteration_sp_overlap:
        logger.debug("These are found with low SO score all_species_dubious_sd_dic " + str(all_species_dubious_sd_dic_updated)+" which are now being handled itr"+str(itr_so)+" .")
        prot_dubious_sd_remove_list = find_prot_dubious_sd_remove(gene_tree, all_species_dubious_sd_dic_updated)

        if prot_dubious_sd_remove_list:
            rest_leaves = set([i.name for i in gene_tree.get_leaves()]) - set(prot_dubious_sd_remove_list)
            gene_tree.prune(rest_leaves, preserve_branch_length=True)
            if conf_infer_subhhogs.gene_trees_write or conf_infer_subhhogs.gene_rooting_method=="mad":
                genetree_msa_file_addr = genetree_msa_file_addr[:-1]+str(int(genetree_msa_file_addr[-1])+itr_so)
                gene_tree.write(format=1, outfile=genetree_msa_file_addr + "_dubiousSD.nwk")

            (gene_tree, all_species_dubious_sd_dic_updated, genetree_msa_file_addr) = _utils_subhog.genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr, conf_infer_subhhogs)

            hogs_children_level_list_raw = hogs_children_level_list
            for prot_dubious_sd_remove in prot_dubious_sd_remove_list:
                logger.debug("** we are removing the sequence "+str(prot_dubious_sd_remove)+"due to low species overlap score")

            for prot_dubious_sd_remove in prot_dubious_sd_remove_list:
                for hog in hogs_children_level_list_raw:
                    if prot_dubious_sd_remove in hog._members:
                        result_removing = remove_prot_hog_hierarchy_toleaves(hog, prot_dubious_sd_remove)
                        if result_removing == 0:  # the hog is empty
                            hogs_children_level_list.remove(hog)

        else:
            all_species_dubious_sd_dic_updated = []

        itr_so += 1

    return (gene_tree, hogs_children_level_list, genetree_msa_file_addr)


def merge_prots_name_hierarchy_toleaves(hog_host, fragment_name_host, merged_fragment_name):

    # todo do I really need this ? as I don't check the members of those subhogs close to the leaves anymore.. to keep it consistent
    for subhog in hog_host._subhogs:
        if fragment_name_host in subhog._members:
            merge_prots_name_hierarchy_toleaves(subhog, fragment_name_host, merged_fragment_name)

    result_merging = hog_host.merge_prots_name_hog(fragment_name_host, merged_fragment_name)

    return hog_host


def merge_fragments_hogclass(fragments_set, seq_dubious_msa, hogs_children_level_list_raw, merged_msa):
    # a bit of redundancy fragments_set is the names which also are stored in  seq_dubious_msa
    hog_host = ""
    fragments_list = list(fragments_set)
    fragment_name_host = fragments_list[0]
    for hog in hogs_children_level_list_raw:
        if fragment_name_host in hog._members:
            hog_host = hog
    merged_fragment_name = "_|_".join(fragments_list)  # fragment_name_host + "_|_" + fragments_list
    fragments_list_remove = fragments_list[1:]
    hogs_to_remove_list =[]
    for fragment_idx in range(1, len(fragments_list)):  # the 0 element is the host

        fragment_name_remove = fragments_list[fragment_idx]
        # fragment_seq_remove = seq_dubious_msa[fragment_idx]
        # remove
        for hog in hogs_children_level_list_raw:
            if fragment_name_remove in hog._members:
                result_removing = remove_prot_hog_hierarchy_toleaves(hog, fragment_name_remove)
                if result_removing == 0:  # the hog is empty
                    hogs_to_remove_list.append(hog)

    hogs_children_level_list = [i for i in hogs_children_level_list_raw if i not in hogs_to_remove_list]

    merged_sequence = ""
    aa_consensus = ''
    for column_idx in range(len(seq_dubious_msa[0])):
        # amino-acid
        aa_col_list = [str(rec.seq)[column_idx] for rec in seq_dubious_msa]
        aa_col_set = set(aa_col_list)
        if aa_col_set == {'-'}:
            aa_consensus='-'
        elif len(aa_col_set-{'-'}) == 1:
            aa_consensus = next(iter(aa_col_set-{'-'}))
        elif len(aa_col_set-{'-'}) > 1:
            aa_consensus = 'X'
        else:
            logger.WARNING("issue 123124124"+str(aa_col_set))

        merged_sequence += aa_consensus
    # seq0 = str(seq_dubious_msa[0].seq)
    # seq1 = str(seq_dubious_msa[1].seq)
    # assert len(seq0) == len(seq1)
    # merged_sequence = ""
    # for aa_idx in range(len(seq1)):
    #     if seq0[aa_idx] == '-':
    #         merged_sequence += seq1[aa_idx]
    #     elif seq1[aa_idx] == '-':
    #         merged_sequence += seq0[aa_idx]
    #     elif seq0[aa_idx] == seq1[aa_idx]:   #  and seq0 != '-' and seq1 != '-'
    #         merged_sequence += seq0[aa_idx]
    #     else:
    #         merged_sequence += "X"
    # assert len(merged_sequence) == len(seq0)
    # merged_fragment_name = "_|_" .join(fragments_list) # fragment_name_host + "_|_" + fragments_list
    if len(merged_fragment_name) > 220:
        logger.warning("The length of sequence id which now being merged is getting very long > 220, we should make sure that it won't cause any issues with fasttree nd biopython and Mafft"+str(merged_fragment_name))
    merged_msa_new_list = []
    if merged_msa and merged_msa[0]:
        assert len(merged_msa[0]) == len(merged_sequence), str(fragment_name_host)
    else:
        logger.warning("issue 15723 merged_msa is empty ?")
    for seq_rec in merged_msa:
        if seq_rec.id == fragment_name_host:
            seq_rec_edited = SeqRecord(Seq(merged_sequence), id=merged_fragment_name, name=merged_fragment_name)
            merged_msa_new_list.append(seq_rec_edited)
        elif seq_rec.id in fragments_list_remove:
            pass
        else:
            merged_msa_new_list.append(seq_rec)
    merged_msa_new = MultipleSeqAlignment(merged_msa_new_list)


    # fragments_list_remove
    result_merging1 = merge_prots_name_hierarchy_toleaves(hog_host, fragment_name_host, merged_fragment_name)

    #hog_host.merge_prots_msa(fragment_name_host, fragment_name_remove, merged_sequence, merged_msa_new)
    hog_host.merge_prots_msa(merged_fragment_name, merged_msa_new)
    logger.debug("these proteins fragments are merged  into one "+str(merged_fragment_name)+" but reported seperately in orthoxml")


    return fragments_list_remove, hogs_children_level_list, merged_msa_new



def handle_fragment_msa(prot_dubious_msa_list, seq_dubious_msa_list, gene_tree, node_species_tree, genetree_msa_file_addr, hogs_children_level_list, merged_msa, conf_infer_subhhogs):
    merged_msa_new = merged_msa
    if not prot_dubious_msa_list:  # empty list
        return gene_tree, hogs_children_level_list, merged_msa_new

    logger.debug("** these are found prot_dubious_msa_list " + str(prot_dubious_msa_list))
    fragments_set_list, seq_dubious_msa_list_checked = check_prot_dubious_msa(prot_dubious_msa_list, gene_tree, seq_dubious_msa_list)
    fragments_remove_list = [] # this list include only one of the fragments for each set, the other one is merged version afterwards
    msa_filt_row_col_new = ""
    if fragments_set_list and len(fragments_set_list[0]) > 1:
        # logger.debug("** these are found fragments_set_list " + str(fragments_set_list))
        # remove fragments from gene tree
        if fragment_detection and fragment_detection_msa_merge:
            merged_msa_new = merged_msa
            for fragments_set_idx, fragments_set in enumerate(fragments_set_list):
                fragments_set = fragments_set_list[fragments_set_idx]
                seq_dubious_msa = seq_dubious_msa_list_checked[fragments_set_idx]
                fragments_list_remove, hogs_children_level_list, merged_msa_new = merge_fragments_hogclass(fragments_set, seq_dubious_msa, hogs_children_level_list, merged_msa_new)
                # merged_msa_new  contains the merged_seq
                # note that merged_msa is the merging of different subhog childrend and is not related to merge fragments  of two seqeuncs with bad annotation

            if merged_msa_new:
                genetree_msa_file_addr = genetree_msa_file_addr[:-1]+str(int(genetree_msa_file_addr[-1])+1)
                if conf_infer_subhhogs.msa_write:
                    SeqIO.write(merged_msa_new, genetree_msa_file_addr + "_fragmentsMerged.fa", "fasta")
                (msa_filt_row_col_new, msa_filt_col, hogs_children_level_list, genetree_msa_file_addr) = _utils_subhog.filter_msa(merged_msa_new, genetree_msa_file_addr, hogs_children_level_list, conf_infer_subhhogs)

            if len(msa_filt_row_col_new) > 1 and len(msa_filt_row_col_new[0]) > 3:
                gene_tree_raw = _wrappers.infer_gene_tree(msa_filt_row_col_new, genetree_msa_file_addr, conf_infer_subhhogs)
                gene_tree = Tree(gene_tree_raw , format=0, quoted_node_names=True) #+ ";"
            else:
                logger.warning("** issue 861956")
                gene_tree = ""
                # fragments_remove_list += fragments_list_remove # for now fragments_list_remove include 1 prots
        elif fragment_detection and (not fragment_detection_msa_merge):
            fragments_remove_set = set.union(*fragments_set_list)
            rest_leaves = set([i.name for i in gene_tree.get_leaves()]) - fragments_remove_set
            if len(rest_leaves) < 2:
                # todo
                logger.warning("** issue 86194")
                hogs_children_level_list = []
                gene_tree = ""
                return gene_tree, hogs_children_level_list, merged_msa_new
            else:
                gene_tree.prune(rest_leaves, preserve_branch_length=True)

        try:
            if  len(gene_tree) and (conf_infer_subhhogs.gene_trees_write or conf_infer_subhhogs.gene_rooting_method == "mad"): # len(gene_tree) > 1 and
                genetree_msa_file_addr = genetree_msa_file_addr[:-1]+str(int(genetree_msa_file_addr[-1])+1)
                gene_tree.write(outfile=genetree_msa_file_addr+"_dubiousMSA.nwk",format=1)
        except:
            logger.warning("issue 24966013  couldn't write the file "+genetree_msa_file_addr+"_dubiousMSA.nwk")


        if len(gene_tree) > 1:
            (gene_tree, all_species_dubious_sd_dic2, genetree_msa_file_addr) = _utils_subhog.genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr, conf_infer_subhhogs, [])


        # todo check the following is needed
        # if all_species_dubious_sd_dic2:
        #     # logger.debug("these are  found after removing with msa , all_species_dubious_sd_dic2 "+str(all_species_dubious_sd_dic2))
        #     (gene_tree, hogs_children_level_list) = handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic2, hogs_children_level_list)
        # the handle_fragment_sd might have the issue as the gene tree might not contain the subHOG ids

        # todo the following was part of this if condition. are they needed
        # for fragments_set in fragments_set_list:
        #     fragments_list = list(fragments_set)
        #     fragment_host = fragments_list[0]  # host of the new small hog consisting few fragments
        #     for hog in hogs_children_level_list:
        #         if fragment_host in hog._members:
        #             hog_host = hog
        #             break  # once we found, we don't need to continue searching in hog
        #     fragments_list_nothost = fragments_list[1:]
        #     for fragment in fragments_list_nothost:
        #         hogs_children_level_list_raw = hogs_children_level_list
        #         for hog in hogs_children_level_list_raw:
        #             if fragment in hog._members:
        #                 result_removing = remove_prot_hog_hierarchy_toleaves(hog, fragment)
        #                 if result_removing == 0:  # the hog is empty
        #                     hogs_children_level_list.remove(hog)
        #                     # print(hogs_children_level_list)
        #     insert_dubious_prots_hog_hierarchy_toleaves(hog_host, fragment_host, fragments_list_nothost)

    return gene_tree, hogs_children_level_list, merged_msa_new

def check_prot_dubious_msa(prot_dubious_msa_list, gene_tree, seq_dubious_msa_list):

    farthest, max_dist_numNodes = gene_tree.get_farthest_node(topology_only=True)  # furthest from the node
    farthest, max_dist_length = gene_tree.get_farthest_node()  # furthest from the node
    print("max_dist_numNodes, max_dist_length ", max_dist_numNodes, max_dist_length)
    fragments_set_list = []
    # some genes may not be in the merged_msa because of trimming rows of msa, we are not using merged_msa
    gene_tree_leaves_name = set([i.name for i in gene_tree.get_leaves()])


    for prot_dubious_msa_idx , prot_dubious_msa_set in enumerate(prot_dubious_msa_list):
        # print(prot_dubious_msa_set)
        fragments = []
        prot_dubious_msa_set_edited = [i for i in prot_dubious_msa_set if i in gene_tree_leaves_name]
        if len(prot_dubious_msa_set_edited)>1:
            prot_host = prot_dubious_msa_set_edited[0]
            for prot in prot_dubious_msa_set_edited[1:]:
                # todo following could be imporved, during filtering row/col msa, a fragments could be removed and not in gene tree anymore,
                # there might be few fragments, checking the distance of the first one with the rest
                # todo this could be improved check all vs all
                assert len(gene_tree.get_leaves_by_name(prot_host)) == 1, "the prot name is not in gene tree or there are more than one" + str(prot_host)
                assert len(gene_tree.get_leaves_by_name(prot)) == 1, "the prot name is not in gene tree or there are more than one" + str(prot)

                # dist_numNodes = gene_tree.get_distance(prot_host, prot, topology_only=True)

                node_prot = gene_tree.get_leaves_by_name(prot)[0]
                node_prot_host = gene_tree.get_leaves_by_name(prot_host)[0]
                dist_length = gene_tree.get_distance(node_prot_host, node_prot)

                dist_length_corrected = dist_length - abs(node_prot.dist- node_prot_host.dist)

                print("check_prot_dubious_msa dist_numNodes, dist_length ", dist_length)
                if dist_length_corrected < max(0.005, max_dist_length / 5):   # or (dist_length_corrected - 2*max_dist_length)< 0.001
                    # todo important hazard to test
                    # dist_numNodes < max(max_dist_numNodes * 1 / 5, 3) or
                    fragments.append(prot)

            if fragments:
                fragments.append(prot_host)
                fragments_set_list.append(set(fragments))

    seq_dubious_msa_list_checked = []
    for fragments_set in fragments_set_list:
        seq_dubious_msa = []
        for seq_original in seq_dubious_msa_list:
            seq_dubious_msa = [i for i in seq_original if i.id in fragments_set]
            if seq_dubious_msa:
                seq_dubious_msa_list_checked.append(seq_dubious_msa)

    logger.debug("These are fragments found " + str(fragments_set_list))

    return fragments_set_list, seq_dubious_msa_list_checked

