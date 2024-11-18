import collections
import csv
import itertools
import sys
from pathlib import Path
import gzip

import dendropy
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import concurrent.futures
import time
import os
import shutil
import pickle
import gc
import random
from ete3 import Tree, TreeNode
import xml.etree.ElementTree as ET
from typing import List

# import networkx as nx
# import matplotlib.pyplot as plt
from . import _wrappers
from . import _utils_subhog
from . import _utils_frag_SO_detection
from ._hog_class import HOG, Representative, split_hog
from ._utils_subhog import MSAFilter, MSAFilterElbow, MSAFilterTrimAL

from ._wrappers import logger
from .zoo.utils import unique

low_so_detection = True # detection of proteins with low species overlap score in gene tree
fragment_detection = True  # this also need to be consistent in _hog_class.py
keep_subhog_each_pickle = False
inferhog_resume_subhog = True
inferhog_max_workers_num = 6 # how parallel to work in a species tree, (won't help near the root)
inferhog_min_hog_size_xml = 2     # by setting this as 1, pyham won't work on xml output.
orthoxml_v03 = True


def read_infer_xml_rhogs_batch(rhogid_batch_list, inferhog_concurrent_on, pickles_rhog_folder, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a list of rootHOGs
    """
    logger.debug("Inferring subHOGs for batch of %d rootHOGs started.", len(rhogid_batch_list))
    hogs_rhog_xml_len_batch = []
    for rhogid in rhogid_batch_list:
        hogs_rhogs_xml_len = read_infer_xml_rhog(rhogid, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
        hogs_rhog_xml_len_batch.append(hogs_rhogs_xml_len)

    return hogs_rhog_xml_len_batch


def read_infer_xml_rhog(rhogid, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a  rootHOGs
    """

    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    if not os.path.exists(pickles_subhog_folder):
        os.makedirs(pickles_subhog_folder)

    # if (_config.gene_trees_write or _config.msa_write) and not os.path.exists("./genetrees"):
    #     os.makedirs("./genetrees")

    logger.info("\n" + "==" * 10 + "\n Start working on root hog: " + rhogid + ". \n")
    rhog_i_prot_address = rhogs_fa_folder + "/HOG_" + rhogid+ ".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    logger.info("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")
    # the file "species_tree_checked.nwk" is created by the check_input.py
    (species_tree) = _utils_subhog.read_species_tree(conf_infer_subhhogs.species_tree)

    (species_tree, species_names_rhog, prot_names_rhog) = _utils_subhog.prepare_species_tree(rhog_i, species_tree, rhogid)
    species_names_rhog = list(set(species_names_rhog))
    logger.info("Number of unique species in rHOG " + rhogid + " is " + str(len(species_names_rhog)) + ".")

    if inferhog_concurrent_on:  # for big HOG we use parallelization at the level taxonomic level using concurrent
        infer_hogs_concurrent(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
    else:
        infer_hogs_for_rhog_levels_recursively(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)

    #####  Now read the final pickle file for this rootHOG
    root_node_name = species_tree.name
    pickle_subhog_file = pickles_subhog_folder + str(root_node_name) + ".pickle"
    with open(pickle_subhog_file, 'rb') as handle:
        hogs_a_rhog = pickle.load(handle)

    if not keep_subhog_each_pickle:
        shutil.rmtree(pickles_subhog_folder)

    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        if len(hog_i._members) >= inferhog_min_hog_size_xml:
            # could be improved   # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_a_rhog_xml_raw = hog_i.to_orthoxml()    # <generef  >      <paralg object >
            if orthoxml_v03 and 'paralogGroup' in str(hogs_a_rhog_xml_raw) :
                # in version v0.3 of orthoxml, there shouldn't be any paralogGroup at root level. Let's put them inside an orthogroup should be in
                hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(hog_i._hogid)})
                property_element = ET.SubElement(hog_elemnt, "property", attrib={"name": "TaxRange", "value": str(hog_i._tax_now.name)})
                hog_elemnt.append(hogs_a_rhog_xml_raw)
                hogs_a_rhog_xml = hog_elemnt
            else:
                hogs_a_rhog_xml = hogs_a_rhog_xml_raw
            hogs_rhogs_xml.append(hogs_a_rhog_xml)
        else:
            logger.debug("we are not reporting due to fastoma singleton hog |*|  " + str(list(hog_i._members)[0]))
            if len(hog_i._members) > 1:
                logger.warning("issue 166312309 this is not a singleton "+str(hog_i._members))


    pickles_rhog_file = pickles_rhog_folder + '/file_' + rhogid + '.pickle'
    with open(pickles_rhog_file, 'wb') as handle:
        # dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("All subHOGs for the rootHOG %s as OrthoXML format is written in %s", rhogid, pickles_rhog_file) # todo report how many genes are in the group for this rootHOG (out of how many were in the initial grouping from omamer)
    # to see orthoxml as string, you might need to do it for different idx
    # idx=0; from xml.dom import minidom; import xml.etree.ElementTree as ET; minidom.parseString(ET.tostring(hogs_rhogs_xml[idx])).toprettyxml(indent="   ")
    del hogs_a_rhog  # to be memory efficient
    gc.collect()
    hogs_rhogs_xml_len = len(hogs_rhogs_xml)
    return hogs_rhogs_xml_len


def infer_hogs_concurrent(species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG using multi-threading (in parallel) on different taxanomic levels of species tree
    """

    pending_futures = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers= inferhog_max_workers_num) as executor:
        for node in species_tree.traverse(strategy="preorder"):
            node.dependencies_fulfilled = set()  # a set
            # node.infer_submitted = False
            if node.is_leaf():
                future_id = executor.submit(singletone_hog_, node, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
                # singletone_hog_(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
                # {<Future at 0x7f1b48d9afa0 state=finished raised TypeError>: 'KORCO_'}
                pending_futures[future_id] = node.name
                # node.infer_submitted = True
        # to keep, have a family
        while len(pending_futures) > 0:
            time.sleep(0.01)
            future_id_list = list(pending_futures.keys())
            for future_id in future_id_list:
                if future_id.done(): # done means finished, but may be unsecussful.

                    species_node_name = pending_futures[future_id]
                    logger.debug("checking for rootHOG id "+str(rhogid)+" future object is done for node " +str(species_node_name))
                    future_id.result()
                    #logger.debug("the result of future is  " + str(future_id.result()))

                    del pending_futures[future_id]
                    species_node = species_tree.search_nodes(name=species_node_name)[0]
                    if species_node == species_tree:  # we reach the root
                        # assert len(pending_futures) == 0, str(species_node_name)+" "+rhogid_
                        assert species_node.name == species_tree.name
                        assert len(pending_futures) == 0
                        break
                    parent_node = species_node.up
                    parent_node.dependencies_fulfilled.add(species_node_name)  # a set
                    childrend_parent_nodes = set(node.name for node in parent_node.get_children())
                    if parent_node.dependencies_fulfilled == childrend_parent_nodes:
                        #  if not parent_node.infer_submitted:
                        future_id_parent = executor.submit(infer_hogs_this_level, parent_node, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs)
                        # parent_node.infer_submitted = True
                        # future_id_parent= parent_node.name+"aaa"
                        pending_futures[future_id_parent] = parent_node.name
                        # for future_id:  del pending_futures[future_id] i need another dictionary the other way arround to removes this futures


def infer_hogs_for_rhog_levels_recursively(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG using recursive function to traverse species tree (different taxanomic levels)
    """

    if sub_species_tree.is_leaf():
        singletone_hog_(sub_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder)
        return

    children_nodes = sub_species_tree.children
    for node_species_tree_child in children_nodes:
        infer_hogs_for_rhog_levels_recursively(node_species_tree_child, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)
    infer_hogs_this_level(sub_species_tree, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs)


def singletone_hog_(node_species_tree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    create subHOGs for leaves of tree (species level). Each protein is a subHOG.
    """

    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    this_level_node_name = node_species_name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    # logger.debug(" ** inferhog_resume_subhog is " + str(_config.inferhog_resume_subhog))
    if inferhog_resume_subhog:
        # logger.debug("inferhog_resume_subhog is " + str(_config.inferhog_resume_subhog) + " so, we are reading from pickles.")
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
        # open already calculated subhogs , but not completed till root in previous run
        if os.path.exists(pickle_subhog_file):
            if os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
                with open(pickle_subhog_file, 'rb') as handle:
                    # i don't even need to open this even
                    # is output of pickle.load(handle) is chlired or this level ?
                    # todo I think I don't need to read the pickle file
                    hogs_this_level_list = pickle.load(handle) #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                    if hogs_this_level_list:
                        logger.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                        return len(hogs_this_level_list)
                    else:
                        logger.debug(" Issue  1238510: the pickle file for single tone is empty "+ str(hogs_this_level_list)+" " + rhogid)

    # logger.debug("reading protien / singletone HOG of  " + str(this_level_node_name))
    hogs_this_level_list = []
    rhog_i_prot_address = rhogs_fa_folder +"/HOG_"+rhogid+".fa"
    with open(rhog_i_prot_address, 'rt') as fasta:
        species_seq_generator = (rec for rec in SeqIO.parse(fasta, format="fasta")
                                 if rec.id.split('||')[1] == node_species_name)
        for prot in species_seq_generator:
            hog_leaf = HOG(prot, node_species_tree, rhogid)  # node_species_tree.name
            hogs_this_level_list.append(hog_leaf)
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)+".pickle"
    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs is written in pickle file.")

    return len(hogs_this_level_list)


def read_children_hogs(node_species_tree, rhogid, pickles_subhog_folder_all):
    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    # pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
    # if _config.inferhog_resume_subhog:  # open already calculated subhogs , but may not be completed in previous run
    #     if os.path.exists(pickle_subhog_file):
    #         if os.path.getsize(pickle_subhog_file) > 3:  # bytes
    #             with open(pickle_subhog_file, 'rb') as handle:
    #                 # i don't even need to open this even
    #                 hogs_children_level_list = pickle.load(handle)
    #                 if hogs_children_level_list:
    #                     return hogs_children_level_list

    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list = []
    for child_name in children_name:
        pickle_subhog_file = pickles_subhog_folder + str(child_name) + ".pickle"
        with open(pickle_subhog_file, 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))  # when there is an error somewhere else, probably with paralelization, for big rootHOG, the root reason of the error won't shown up. but you stope here at worst which is late.
            #  todo, check file exist, how to handle if not
    logger.debug("Finding hogs for rhogid: " + rhogid+ ", for taxonomic level:" + str(
        node_species_tree.name) + " for species sub-tree:\n  " + str(node_species_tree.write(format=1, format_root_node=True)) + "\n")
    return hogs_children_level_list


RepLookup = collections.namedtuple("RepLookup", ["hog", "representative"])
class LevelHOGProcessor:
    def __init__(self, node_species_tree:TreeNode, subhogs:List[HOG], rhogid:str, conf):
        self.node_species_tree = node_species_tree
        self._frozen_hogs_before = set([])
        self.subhogs = {}
        for hog in subhogs:
            if hog.active:
                self.subhogs[hog.hogid] = hog
            else:
                self._frozen_hogs_before.add(hog)
        self.rhogid = rhogid
        self.conf = conf
        self._base_msa_tree_filename = "HOG_"+rhogid+"_"+str(node_species_tree.name)
        self._outname_step = 0
        self._removed_rep = set([])
        self._rep_lookup = self._prepare_lookups()
        self._msa_filter = self._instantiate_msa_filter()

    def _instantiate_msa_filter(self):
        if self.conf.msa_filter_method == "col-row-threshold":
            return MSAFilter(self, self.conf)
        elif self.conf.msa_filter_method == "col-elbow-row-threshold":
            return MSAFilterElbow(self, self.conf)
        elif self.conf.msa_filter_method == "trimal":
            return MSAFilterTrimAL(self, self.conf)
        else:
            raise ValueError("Invalid msa_filter_method value")

    def get_name_of_output(self, is_msa=False, is_tree=False):
        if (is_msa and self.conf.msa_write) or (is_tree and self.conf.gene_trees_write):
            self._outname_step += 1
            return f"{self._base_msa_tree_filename}_it{self._outname_step}"
        return None

    def _prepare_lookups(self):
        lookup = {}
        for hog in self.subhogs.values():
            for rep in hog.get_representatives():
                if rep.get_id() not in self._removed_rep:
                    lookup[rep.get_id()] = RepLookup(hog, rep)
        return lookup

    def find_most_divergent_representatives_from_genetree(self, genetree_subtree:TreeNode):
        """this method returns the N most divergent representatives from a rooted tree. within each cluster
        it selects the one representative sequence that corresponds to the the median distance from the root
        among the sequnces in the cluster.

        :param genetree: the rooted tree (with distances) that should be processed
        :return: the N most divergent representatives for that tree. N is taken from self.conf"""
        genetree = copy_tree(genetree_subtree)

        # we remove the non-enabled representatives from the tree
        keep = [n for n in genetree.iter_leaves() if self._rep_lookup[n.name].representative.enabled]
        if len(keep) < len(genetree):
            genetree.prune(keep, preserve_branch_length=True)
        if len(genetree) <= self.conf.number_of_samples_per_hog:
            return [Representative(self._rep_lookup[n.name].representative) for n in genetree.iter_leaves()]

        # first, annotate the distance from the root
        for n in genetree.traverse():
            n.add_feature('dist_from_root', n.dist + (n.up.dist_from_root if n.up is not None and hasattr(n.up, 'dist_from_root') else 0))
        # now, let's get a list of all the internal nodes ordered with the dist_from_root and find which distance
        # results in number_of_samples_per_hog subtrees that are most divergent
        ordered = sorted([n for n in genetree.traverse() if not n.is_leaf()], key=lambda x: x.dist_from_root)
        idx = min(len(ordered), self.conf.number_of_samples_per_hog) - 1  # in case of multi-fucations, len(ordered) can be shorter than number_of_samples_per_hog
        cutoff = ordered[idx].dist_from_root
        sub_clades = [sub for sub in genetree.iter_leaves(
            is_leaf_fn=lambda n: len(n.children) == 0 or n.dist_from_root >= cutoff
        )]
        # now, let's select for each sub_clade the one leave that is close to the median dist_from_root (no outliers)
        representatives = []
        for sub in sub_clades:
            order_by_dist = sorted((n for n in sub.iter_leaves()), key=lambda n: n.dist_from_root)
            median_elem = order_by_dist[len(order_by_dist)//2]
            new_rep = Representative(self._rep_lookup[median_elem.name].representative,
                                     [self._rep_lookup[z.name].representative for z in order_by_dist])
            representatives.append(new_rep)
        return representatives

    def _get_merge_candidates_with_hogids(self, reconciled_genetree:TreeNode):
        for top_speciation_node in reconciled_genetree.iter_leaves(
                is_leaf_fn=lambda n: len(n.children) == 0 or (hasattr(n, 'evoltype') and n.evoltype == "S")
        ):
            hogids_in_subtree = set([])
            for n in top_speciation_node.iter_leaves():
                n.add_feature('hogid', self._rep_lookup[n.name].hog.hogid)
                hogids_in_subtree.add(n.hogid)
            top_speciation_node.add_feature('hogids', hogids_in_subtree)
            yield top_speciation_node

    def write_msa_or_tree_if_necessary(self, elem, fn_suffix="", features=[]):
        is_tree = isinstance(elem, TreeNode)
        is_msa = isinstance(elem, MultipleSeqAlignment)
        fn = self.get_name_of_output(is_tree=is_tree, is_msa=is_msa)
        if fn is not None:
            fn = Path(fn + fn_suffix)
            with open(fn, 'wt') as fout:
                if is_tree:
                    fout.write(elem.write(format=1, features=features, format_root_node=True))
                elif is_msa:
                    SeqIO.write(elem, fout, format="fasta")
            if is_tree:
                fn2 = fn.stem + ".tsv"
                with open(fn2, 'wt') as fout:
                    writer = csv.writer(fout, dialect="excel-tab")
                    writer.writerow(["id", "hogid", "hog_size", "representative_size", "enabled"])
                    for n in elem.iter_leaves():
                        rep_val = self._rep_lookup[n.name]
                        writer.writerow([n.name, rep_val.hog.hogid, len(rep_val.hog.get_members()), len(rep_val.representative.get_subelements()), rep_val.representative.enabled])
                fn3 = fn.stem + "-rep.tsv.gz"
                with gzip.open(fn3, 'wt') as fout:
                    writer = csv.writer(fout, dialect="excel-tab")
                    writer.writerow(["id", "representative"])
                    for n in elem.iter_leaves():
                        rep_val = self._rep_lookup[n.name]
                        for member in rep_val.representative.get_subelements():
                            writer.writerow([n.name, member])

    def align_subhogs(self):
        sub_msas = [hog.get_msa() for hog in self.subhogs.values() if len(hog.get_msa()) > 0]
        logger.debug(f"Merging {len(sub_msas)} MSAs for rhog: {self.rhogid}, level: {self.node_species_tree.name}")
        if len(sub_msas) == 0:
            logger.info(
                f"Issue 1455, merged_msa is empty {self.rhogid}, for taxonomic level: {self.node_species_tree.name}")
            raise ValueError("Issue 1455, merged_msa is empty")
        if len(sub_msas) > 1:
            merged_msa = _wrappers.merge_msa(sub_msas)
            self.write_msa_or_tree_if_necessary(merged_msa, fn_suffix=".fa")

            if fragment_detection:
                prot_dubious_msa_list, seq_dubious_msa_list = _utils_frag_SO_detection.find_prot_dubious_msa(
                    merged_msa, self.conf)
                if len(prot_dubious_msa_list) > 0:
                    logger.info(f"found fragments: {prot_dubious_msa_list} with {seq_dubious_msa_list}")
                    logger.error("Handling of fragments not yet implemented. Ignoring for now")
        else:
            merged_msa = sub_msas[0]   # when only on  child, the rest msa is empty.
        logger.debug(f"All sub-hogs are merged, merged_msa ({len(merged_msa)},{len(merged_msa[0])}) for "
                     f"rhog: {self.rhogid}, taxonomic level: {self.node_species_tree.name}")
        return merged_msa

    def filter_msa(self, msa):
        if self._msa_filter is None:
            self._msa_filter = MSAFilter(self, self.conf)
        filtered_msa, removed_ids = self._msa_filter.filter_msa(msa)
        if len(removed_ids) > 0:
            self._remove_representatives(ids=removed_ids)
        return filtered_msa

    def _remove_representatives(self, ids):
        self._removed_rep.update(ids)
        hogids_to_remove = []
        for hog_id, hog in self.subhogs.items():
            for rep in hog.get_representatives():
                if rep.get_id() not in self._removed_rep:
                    break
            else:
                # going here if not a break -> hog disappeared. remove it from current level and freeze it at level below
                hogids_to_remove.append(hog_id)
                hog.active = False
                self._frozen_hogs_before.add(hog)
        for h in hogids_to_remove:
            self.subhogs.pop(h)
        self._rep_lookup = self._prepare_lookups()


    def infer_genetree_from_msa(self, msa):
        genetree_nwk = _wrappers.infer_gene_tree(msa)
        try:
            genetree = Tree(genetree_nwk + ";", format=0, quoted_node_names=True)
        except ValueError:
            logger.error(f"cannot load genetree from {genetree_nwk}")
            raise RuntimeError("cannot load genetree from {}".format(genetree_nwk))
        self.write_msa_or_tree_if_necessary(genetree, fn_suffix=".nwk")
        return genetree

    def infer_rooted_genetree(self, gene_tree: TreeNode):
        if self.conf.gene_rooting_method == "midpoint":
            r_outgroup = gene_tree.get_midpoint_outgroup()
            gene_tree.set_outgroup(r_outgroup)  # print("Midpoint rooting is done for gene tree.")

        elif self.conf.gene_rooting_method == "midpoint-dendropy":
            dt = dendropy.Tree.get_from_string(gene_tree.write(format=1), schema="newick")
            dt.reroot_at_midpoint()
            nw = dt.as_string(schema="newick", suppress_rooting=True)
            gene_tree = Tree(nw, format=1, quoted_node_names=True)

        elif self.conf.gene_rooting_method == "Nevers_rooting":
            logger.info("Nevers_rooting started for " + str(gene_tree.write(format=1, format_root_node=True)))
            species = Tree("species_tree.nwk", format=1)
            gene_tree = _utils_subhog.get_score_all_root(gene_tree, species)
            logger.info("Nevers_rooting finished for " + str(gene_tree.write(format=1, format_root_node=True)))

        elif self.conf.gene_rooting_method == "mad":
            gene_tree = _wrappers.mad_rooting(gene_tree)  # todo check with qouted gene tree

        elif self.conf.gene_rooting_method == "outlier":  # not yet implmented completely, todo need check with new gene tree
            outliers = _utils_subhog.find_outlier_leaves(gene_tree)
            r_outgroup = _utils_subhog.midpoint_rooting_outgroup(gene_tree, leaves_to_exclude=outliers)
            gene_tree.set_outgroup(r_outgroup)
        else:
            logger.warning("rooting method not found !!   * * * * *  *")
            raise ValueError("invalid rooting method: {}".format(self.conf.gene_rooting_method))
        return gene_tree

    def infer_reconciliation(self, genetree:TreeNode, sos_threshold=0.0):
        """Annotate each internal node with a 'evoltype' and 'sos' attribute.

        The evoltype attribute will be either 'S' or 'D' (speciation or duplication),
        and the 'sos' attribute will contain the fraction of species that overlap among the
        subtrees. The parameter sos_threshold specifies the minimum fraction of species that
        must overlap among the subtrees in order to be considered a duplication.

        :param genetree: the rooted tree that should be processed
        :param sos_threshold: the minimum fraction of species that must overlap among the
                              subtrees in order to be considered a duplication.
        """
        logger.info("\n" + genetree.get_ascii(show_internal=True, attributes=['name', 'hogid']))
        cnt_D, cnt_S = 0, 0
        for n in genetree.traverse('postorder'):
            if n.is_leaf():
                n.add_feature('species', self._rep_lookup[n.name].representative.get_species())
                n.add_feature('hogid', self._rep_lookup[n.name].hog.hogid)
            else:
                assert len(n.children) >= 2, f"{n.name} has {len(n.children)} children"
                specs = tuple(c.species for c in n.children)
                for c in n.children:
                    c.del_feature('species')  # cleanup to avoid excessive memory
                sp_inter = set.intersection(*specs)
                sp_union = set.union(*specs)
                sos = len(sp_inter) / len(sp_union)
                n.add_feature('species', sp_union)
                n.add_feature('sos', sos)
                n.add_feature('so_tuple', (len(sp_inter), len(sp_union)))
                n.add_feature('evoltype', 'D' if sos > sos_threshold else 'S')
                if sos > sos_threshold:
                    cnt_D += 1
                    n.name = f"D{cnt_D}_{len(sp_inter)}_{len(sp_union)}"
                else:
                    cnt_S += 1
                    n.name = f"S{cnt_S}"
        genetree.del_feature('species')
        self.write_msa_or_tree_if_necessary(genetree, fn_suffix="_rec.nwk", features=['evoltype', 'sos', 'so_tuple', 'hogid'])

    def _compute_between_and_within_distances(self, genetree, subtrees_node, partitions):
        mrca = genetree.get_common_ancestor(*subtrees_node)
        id2part = {id_: i for i, p in enumerate(partitions) for id_ in p}
        N = len(partitions)
        dists = collections.defaultdict(list)
        for l1, l2 in itertools.combinations(mrca.iter_leaves(), 2):
            i1 = id2part.get(l1.name, N)
            i2 = id2part.get(l2.name, N)
            dists[(i1, i2)].append(mrca.get_distance(l1, l2))

        return dists

    def _resolve_conflicts(self, gene_tree: TreeNode, subtrees_node, partitions) -> bool:
        """resolves cases of conflicting subhog splittings where
        the representatives of a subhog are split among different subtrees
        and have been merged even before the previous level.

        the method returns a boolean value indicating whether a new reconciliation
        is necessary (e.g. in case the tree has been modified and hence the labeling
        of speciation/duplication might have changed).
        """
        updated_subtree_nodes = []
        for n in subtrees_node:
            # previous rounds of _resolve_conflicts might have left disconnected subtree_nodes
            # so we identify the first node which is either a leaf or a bifurcating node.
            while len(n.children) == 1:
                n = n.children[0]
            updated_subtree_nodes.append(n)
        subtrees_node = updated_subtree_nodes
        mrca = gene_tree.get_common_ancestor(*subtrees_node)
        min_support = min(n.support for n in subtrees_node)
        min_branch_to_subtree = min(gene_tree.get_distance(mrca, n) for n in subtrees_node)
        partitions = sorted(partitions, key=lambda x: [-sum(len(self._rep_lookup[r].representative.get_subelements()) for r in x), -len(x),])
        if False and (min_support < 0.7 or min_branch_to_subtree < 0.01):
            # we collapse things, i.e. separate nodes in tree based on hogid
            # we mv all leaves from the non-biggest partitions to the mrca of the biggest partion
            if len(partitions[0]) > 1:
                target_node = gene_tree.get_common_ancestor(partitions[0])
            else:
                target_node = gene_tree.search_nodes(name=partitions[0][0])[0].up
            for small_partition in partitions[1:]:
                for name in small_partition:
                    n = gene_tree.search_nodes(name=name)[0]
                    parent = n.up
                    n = n.detach()
                    target_node.add_child(n)
                    if len(parent.children) == 1 and parent.up is None:
                        parent.children[0].delete(prevent_nondicotomic=True, preserve_branch_length=True)
                    elif len(parent.children) == 1:
                        parent.delete(prevent_nondicotomic=True, preserve_branch_length=True)
                    # deactivate representative for next level
                    self._rep_lookup[n.name].representative.enabled = False
            return True
        else:
            # we remove the bogus representatives, but keep the labeling of the
            # speciation/duplication nodes. that way, the two subhogs will not get merged
            # at this level (since they are in different subtrees).
            stats = [{'len': len(p),
                      'rep_size': sum(len(self._rep_lookup[r].representative.get_subelements()) for r in p),
                     } for p in partitions]
            logger.info(f"Resolving conflict for {self._rep_lookup[partitions[0][0]].hog.hogid} by ignoring representatives {partitions[1:]}")
            logger.debug(f"{len(partitions)} partitions: {stats}")
            for small_partition in partitions[1:]:
                for name in small_partition:
                    n = gene_tree.search_nodes(name=name)[0]
                    n.delete(prevent_nondicotomic=True, preserve_branch_length=True)
            return False

    def merge_subhogs(self, reconciled_genetree:TreeNode, msa:MultipleSeqAlignment):
        while True:
            hogids_in_subtrees = collections.defaultdict(list)
            for n in self._get_merge_candidates_with_hogids(reconciled_genetree):
                for hogid in n.hogids:
                    hogids_in_subtrees[hogid].append(n)
            if any(len(z) > 1 for z in hogids_in_subtrees.values()):
                logger.info("At least one subhog is split. here is the full labeled genetree:\n"
                            + reconciled_genetree.get_ascii(show_internal=True, attributes=['evoltype', 'sos', 'hogids', 'hogid']))
            else:
                break
            redo_reconcilation = False
            for hogid, subtrees in hogids_in_subtrees.items():
                if len(subtrees) > 1:
                    logger.info(f"Representaives of {hogid} are split among {len(subtrees)} candidate subtrees.")
                    split_parts = [list(n.name for n in sub.iter_leaves() if n.hogid == hogid) for sub in subtrees]
                    split_hogs = split_hog(self.subhogs[hogid], self.node_species_tree.name, *split_parts)
                    if split_hogs and len(split_hogs) > 1:
                        # we could split the current hog.
                        self.subhogs.pop(hogid)
                        self.subhogs.update({h.hogid: h for h in split_hogs})
                    else:
                        redo_reconcilation |= self._resolve_conflicts(reconciled_genetree, subtrees, split_parts)
            self._rep_lookup = self._prepare_lookups()
            if redo_reconcilation:
                if len(reconciled_genetree.children) == 1:
                    reconciled_genetree.children[0].delete(preserve_branch_length=True)
                self.infer_reconciliation(reconciled_genetree, sos_threshold=self.conf.threshold_dubious_sd)
        new_hogs = []
        processed_nodes = set([])
        for subtrees in hogids_in_subtrees.values():
            for subtree in subtrees:
                if subtree in processed_nodes:
                    continue
                processed_nodes.add(subtree)
                new_repr = self.find_most_divergent_representatives_from_genetree(subtree)
                hog = HOG([self.subhogs[h] for h in subtree.hogids],
                          self.node_species_tree,
                          rhogid=self.rhogid,
                          msa=msa,
                          representatives=new_repr,
                          conf_infer_subhhogs=self.conf)
                new_hogs.append(hog)
        return new_hogs

    def process(self) -> List[HOG]:
        """Process the given roothog at the given taxonomic level.

        This method does all the processessing of the subhogs from the
        previous taxonomic levels. It does

        1. Align the subhogs from the previous levels as a new MSA.
        2. Filter the initial MSA (remove gappy columns and rows)
        3. Infer a genetree for the new filtered MSA
        4. Root (mid-point) and reconcile (with species-overlap) the gene tree
        5. Merge the subhogs from the previous levels into new HOGs
        """
        msa = self.align_subhogs()
        filtered_msa = self.filter_msa(msa)
        if len(filtered_msa) == 0:
            return list(self.subhogs.values())
        gene_tree = self.infer_genetree_from_msa(filtered_msa)
        if len(gene_tree) > 1:
            rooted_gene_tree = self.infer_rooted_genetree(gene_tree)
            # TODO: dealing with fragments is not done yet.
            self.infer_reconciliation(rooted_gene_tree, sos_threshold=self.conf.threshold_dubious_sd)
        else:
            rooted_gene_tree = gene_tree
        new_hogs = self.merge_subhogs(rooted_gene_tree, msa=filtered_msa)
        new_hogs.extend(self._frozen_hogs_before)
        return new_hogs


def infer_hogs_this_level(node_species_tree, rhogid, pickles_subhog_folder_all, conf_infer_subhhogs):
    """
    infer subHOGs for a rootHOG at a taxanomic level
    """

    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + rhogid + "/"
    if node_species_tree.is_leaf():
        logger.warning("issue 1235,it seems that there is only on species in this tree, singleton hogs are treated elsewhere "+ rhogid)

    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"

    # TODO arrage resume with nextflow and also for when read single_tone pickles
    if inferhog_resume_subhog:
        if os.path.exists(pickle_subhog_file) and os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
            with open(pickle_subhog_file, 'rb') as handle:
                # todo : do I really need to read the pickle file
                hogs_this_level_list = pickle.load(handle)  #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                if hogs_this_level_list:
                    logger.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                    return len(hogs_this_level_list)
                else:
                    logger.debug(" Issue  1238510: the pickle file for single tone is empty " + str(hogs_this_level_list) + " " + rhogid)

    hogs_children_level_list = read_children_hogs(node_species_tree, rhogid, pickles_subhog_folder_all)

    if len(hogs_children_level_list) == 1:
        hogs_this_level_list = hogs_children_level_list
        with open(pickle_subhog_file, 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return len(hogs_children_level_list)

    level_processor = LevelHOGProcessor(node_species_tree, hogs_children_level_list, rhogid, conf_infer_subhhogs)
    hogs_this_level_list = level_processor.process()

    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return len(hogs_this_level_list)


def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid, merged_msa, conf_infer_subhhogs):
    """
    merge subhogs based on the gene tree specieciaton node of gene tree by creating inter-HOG graph (implicitley )
    """

    subhogs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
    hogs_this_level_list = []
    subHOG_to_be_merged_set_other_Snodes = []
    subHOG_to_be_merged_set_other_Snodes_flattned_temp = []
    ##  the following if for debugging and visualisation of connected component of inter-HOG graph
    # hoggraph_node_name = [i._hogid.split("_")[1][3:] for i in hogs_children_level_list]
    # hog_size_dic = {}
    # dic_hog = {}
    # for hog in hogs_children_level_list:
    #     hog_id_short = hog._hogid.split("_")[1][3:]
    #     for prot in hog._members:
    #         dic_hog[prot] = hog_id_short + "_len" + str(len(hog._members))
    #     hog_size_dic[hog_id_short] = len(hog._members)
    # for hog in hogs_children_level_list:
    #     hog_id_short = hog._hogid.split("_")[1][3:] + "_len" + str(len(hog._members))
    #     # print(hog_id_short, ":", hog._members)
    # hoggraph_node_name_len = [i+"_len"+str(hog_size_dic[i]) for i in hoggraph_node_name]
    # hoggraph = nx.Graph()
    # hoggraph.add_nodes_from(hoggraph_node_name_len)

    for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "processed") and n.processed==True):

        if not node.is_leaf() and node.name[0] == "S":
            node_leaves_name_raw0 = [i.name for i in node.get_leaves()]
            # leaves names  with subhog id  'HALSEN_R15425||HALSEN|_|1352015793_sub10149',
            node_leaves_name_raw1 = [i.split("|_|")[0] for i in node_leaves_name_raw0]

            if node_leaves_name_raw1[0].startswith("'"):
                node_leaves_name = [i[1:] for i in node_leaves_name_raw1] # since splited , there is no ' at the end. so -1] is not needed.
                # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005", gene tree is quoted, there are both ' and " !
            else:
                node_leaves_name = node_leaves_name_raw1
            # node_leaves_name =[]
            # for name_i in node_leaves_name_:
            #     node_leaves_name += name_i.split("_|_")

            # num_prot = len(s_gene_tree_leaves)
            # for i in range(num_prot):
            #     hog_i = dic_hog[s_gene_tree_leaves[i]]
            #     for j in range(i):
            #         hog_j = dic_hog[s_gene_tree_leaves[j]]
            #         if hoggraph.has_edge(hog_i, hog_j):
            #             hoggraph[hog_i][hog_j]['weight'] += 1
            #         else:
            #             hoggraph.add_edge(hog_i, hog_j, weight=1)
            # print(hoggraph.edges(data=True))
            # todo this piece of code could be neater, for extracting connected component of inter-hog graph
            subHOG_to_be_merged = []
            for node_leave_name in node_leaves_name:  # print(node_leave_name)
                for subHOG in hogs_children_level_list:
                    subHOG_members = subHOG._members
                    if node_leave_name in subHOG_members:  # could be improved
                        if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp: #todo 18march this means that we are not merging a few together
                            subHOG_to_be_merged.append(subHOG)
                            subhogs_id_children_assigned.append(subHOG._hogid)
                        else:  # this hog is already decided to be merged  print(node.name, subHOG._hogid, node_leave_name)
                            if "processed" in node:
                                logger.warning("issue 1863 "+ str(node.name)+str(subHOG._hogid)+ str(node_leave_name)) # print("processed", node.name) #else: #    print("processed not in ", node.name)  # print(node_leave_name,"is in ",subHOG._hogid)
            if len(subHOG_to_be_merged) == 1:
                logger.warning("issue 125568313 "+str(subHOG_to_be_merged)+" "+node.name+" probably the subhog was merged previously" )
                logger.debug("issue 125568313 " + str(subHOG_to_be_merged_set_other_Snodes_flattned_temp))
                logger.debug("issue 125568313 " +str(node.name)+" "+ str(gene_tree.write(format=1,format_root_node=True)))

            elif len(subHOG_to_be_merged)>1:
                subHOG_to_be_merged_set = set(subHOG_to_be_merged)
                HOG_this_node = HOG(subHOG_to_be_merged_set, node_species_tree, rhogid, merged_msa, conf_infer_subhhogs)
                if len(HOG_this_node._msa) == 1:
                    logger.warning("issue 1258313"+str(HOG_this_node)+str(HOG_this_node._msa)+" "+node.name  )
                hogs_this_level_list.append(HOG_this_node)

                subHOG_to_be_merged_set_other_Snodes.append([i._hogid for i in subHOG_to_be_merged_set])
                subHOG_to_be_merged_set_other_Snodes_flattned_temp = [item for items in subHOG_to_be_merged_set_other_Snodes for  item in items]
                #  I don't need to traverse deeper in this clade
            node.processed = True  # print("?*?*  ", node.name)
        subHOG_to_be_merged_set_other_Snodes_flattned = [item for items in subHOG_to_be_merged_set_other_Snodes for item in items]
        if [i._hogid for i in hogs_children_level_list] == subHOG_to_be_merged_set_other_Snodes_flattned:
            break
    ##  the following if fore debugging and visualisation of connected component of inter-HOG graph
    # fig = plt.figure(figsize=(300, 200), dpi=60)
    # pos = nx.spring_layout(hoggraph, k=0.25, iterations=30)  # For better example looking  # smaller k, biger space between
    # nx.draw(hoggraph, pos, with_labels=True, node_color='y', node_size=500, font_size=16) # , alpha=0.4
    # # nx.draw(G, pos,, edge_color="r", font_size=16, with_labels=True)
    # labels = {e: hoggraph.edges[e]['weight'] for e in hoggraph.edges}
    # nx.draw_networkx_edge_labels(hoggraph, pos, edge_labels=labels, font_size=16)
    # num = random.randint(3, 1000000)
    # plt.savefig("./hoggraph/" + hogs_children_level_list[0]._hogid[4:] + "file_rndm"+str(num)+".jpg")
    # plt.show()
    for subHOG in hogs_children_level_list:  # for the single branch  ( D include subhog and S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    for hog_j in hogs_this_level_list:
        hog_j._tax_now = node_species_tree
    ##  the following if fore debugging of connected component of inter-HOG graph
    ## check for conflicts in merging
    #     for i in range(subHOG_to_be_merged_set_other_Snodes):  if
    #         for i in range(subHOG_to_be_merged_set_other_Snodes):  print("*&*& ",node_species_tree.name)
    # prot_list_sbuhog = [i._members for i in hogs_this_level_list]
    # prot_list_sbuhog_short = []
    # for prot_sub_list_sbuhog in prot_list_sbuhog:
    #     if format_prot_name == 0:  # bird dataset TYTALB_R04643
    #         prot_list_sbuhog_short = prot_sub_list_sbuhog
    #     elif format_prot_name == 1:  # qfo dataset  'tr|E3JPS4|E3JPS4_PUCGT
    #         prot_list_sbuhog_short.append([prot.split('|')[2] for prot in prot_sub_list_sbuhog])
    # logger_hog.debug(str(len(hogs_this_level_list))+" hogs are inferred at the level "+node_species_tree.name+": "+' '.join(
    #     [str(i) for i in prot_list_sbuhog_short]))

    return hogs_this_level_list


def copy_tree(t:TreeNode):
    """This function creates a deep copy of an ete3 tree while
    preserving the names

    standard code in ete3:
       new_node = self.__class__(self.write(features=[]))
    """
    features = set([])
    for n in t.traverse():
        features.union(n.features)
    features -= {"name"}
    new_node = t.__class__(t.write(features=features, quoted_node_names=True, format_root_node=True),
                           quoted_node_names=True)
    return new_node


