
import csv
import itertools

import networkx as nx
from Bio import SeqIO
import pickle
from os import listdir
import os
import sys

from .zoo.unionfind import UnionFind
from ._wrappers import logger
import collections

filter_nonchild_rootHOG= False
mmseqs_executable_path ="mmseqs"

HOGMapData = collections.namedtuple("HOGMapData", ("hogid", "score", "seqlen", "subfamily_medianseqlen"))


"""UnionFind.py

Union-find data structure. Based on Josiah Carlson's code,
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
with significant additional changes by D. Eppstein and
Adrian Altenhoff.
"""


class UnionFind(object):
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self, elements=None):
        """Create a new union-find structure.

        If elements is not None, the structure gets initialized
        with each element as a singleton component.

        :param elements: an iterable to initialize the structure.
        """

        self.weights = {}
        self.parents = {}
        if elements is not None:
            for elem in iter(elements):
                self.parents[elem] = elem
                self.weights[elem] = 1

    def __getitem__(self, obj):
        """return the name of set which contains obj.

        :param obj: the query object

        :SeeAlso: :meth:`find`"""
        return self.find(obj)

    def find(self, obj):
        """Find and return the name of the set containing the obj.

        If the object is not found in any set, a new singleton set
        is created that holds only this object until it is further merged."""

        # check for previously unknown obj. If unknown, add it
        # as a new cluster
        if obj not in self.parents:
            self.parents[obj] = obj
            self.weights[obj] = 1
            return obj

        # find path of objects leading to the root
        path = [obj]
        root = self.parents[obj]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def remove(self, obj):
        """Remove an object from the sets.

        Removes an object entirly from the datastructure. The
        containing set will shrink by this one element.

        :Note: If one tries to accessed it afterwards using
            :meth:`find`, it will be created newly and put as a
            singleton.
        """
        if obj not in self.parents:
            return
        comp = self.find(obj)
        self.weights[comp] -= 1
        self.parents.pop(obj)

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them.

        any number of objects can be passed to this method and
        all of them will be merged into one set containing at
        least these objects.

        :param objects: the objects to be merged. they have to be all
            hashable. If they haven't been initialy added to the UnionFind
            datastructre at instantiation time, they are added at this point
            in time.
        """
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r], r) for r in roots], key=lambda x: x[0])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

    def get_components(self):
        """return a list of sets corresponding to the connected
        components of the structure."""
        comp_dict = collections.defaultdict(set)
        for elem in iter(self):
            comp_dict[self[elem]].add(elem)
        comp = list(comp_dict.values())
        return comp


def parse_proteomes(folder=None, min_sequence_length=0):  # list_oma_species
    """
    parsing fasta files of proteins located in /proteome/
    using Bio.SeqIO.parse
    Each fasta file is for one species.  The file name is the species name.
    All proteomes should be with the same extension.
    output: prot_recs_lists: a dic with key as species name and  its value list of Biopython record of species.
    """
    if folder is None:
        folder = "./proteome"

    fasta_format_keep = ""
    project_files = listdir(folder)
    logger.debug("using '%s' as proteome folder, found %d files", folder, len(project_files))
    species_names = []  # query/input species name based on the file name
    for file in project_files:
        species_name, ext = file.rsplit('.', 1)
        logger.debug("%s: %s | %s", file, species_name, ext)
        if ext in ("fa", "fasta"):
            species_names.append(species_name)
            fasta_format_keep = ext  # last one is stored either fa or fasta

    # todo accept all fasta formats in the input prtoeome folder, fasta, fa, fna, ..
    prot_recs_lists = {} # key: species name, value is a dic of query protein Biopython records. # 'MYCGE': [SeqRecord(seq=Seq('MDFDK
    #smallprot_recs_lists ={}
    for species_name in species_names:
        prot_address = os.path.join(folder, species_name + "." + fasta_format_keep)
        prots_record = list(rec for rec in SeqIO.parse(prot_address, "fasta") if len(rec) >= min_sequence_length)
        #prots_record_small = list(rec for rec in SeqIO.parse(prot_address, "fasta") if len(rec) < min_sequence_length)
        # logger.debug(prots_record)
        logger.info(f"{species_name} contains {len(prots_record)} that are at least {min_sequence_length} long.")
        prot_recs_lists[species_name] = prots_record
        #smallprot_recs_lists[species_name]=prots_record_small

    logger.info("There are %d species in the proteome folder.", len(species_names))
    return species_names, prot_recs_lists, fasta_format_keep #, smallprot_recs_lists


def add_species_name_prot_id( prot_recs_lists):
    """
    adding the name of species to each protein record
        - based on file name
    adding protein idx number, integer needed by xml format
    output:  prot_recs_all =  {'MYCGE': {'sp|P47500|RF1_MYCGE||MYCGE||1000000001': SeqRecord(seq=
    """
    prot_idx_name_pickle_file = "./gene_id_dic_xml.pickle"
    start_num_prot = int(1e9)
    start_num_prot_per_sp = int(1e6) #
    prot_recs_all = {} # {'MYCGE': {'sp|P47500|RF1_MYCGE||MYCGE||1000000001': SeqRecord(seq=
    prot_idx_name = {} # {'MYCGE': [(1000000001, 'sp|P47500|RF1_MYCGE'),(1000000002, 'sp|P13927|EFTU_MYCGE'),
    species_idx = -1
    for species_name, prot_recs_list in prot_recs_lists.items():
        species_idx += 1
        prot_idx = start_num_prot + species_idx * start_num_prot_per_sp
        prot_recs_all[species_name]={}
        prot_idx_name[species_name] = []
        for prot_rec in prot_recs_list:
            prot_idx+=1
            prot_name= prot_rec.id

            if len(prot_name) > 230:
                logger.info("We are truncating the prot name as it may be problematic for mafft, " + str(prot_name))
                prot_name = prot_name[:230]

            # todo, this could be a dic
            prot_idx_name[species_name].append((prot_idx, prot_name))

            prot_name_new = prot_name+ "||"+species_name+"||"+str(prot_idx) # orthoxml file needs an integer as
            prot_rec.id = prot_name_new
            prot_recs_all[species_name][prot_name] = prot_rec

    with open(prot_idx_name_pickle_file, 'wb') as handle:
        pickle.dump(prot_idx_name, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return prot_recs_all


def parse_hogmap_omamer(proteomes, fasta_format_keep, folder=None):
    """
     function for parsing output of omamer (hogmap files) located in /hogmap/
    Each hogmap file correspond to one fasta file of species, with the same name.
    Note that some records of fasta may removed in hogmap, probably becuase of being so short.
    hogmap file example:
    # qseqid	hogid	hoglevel	family_p	family_count	family_normcount	subfamily_score	subfamily_count	qseqlen	subfamily_medianseqlen	qseq_overlap
    # sp|O66907|ATPA_AQUAE	HOG:C0886513.1b	Eukaryota	696.5519485850672	187	0.37302594463771466	0.2643067275644273	120	504	551	0.8230616302186878
    output is dic of dic for all species:
    """
    if folder is None:
        folder = "./hogmap"

    hogmaps = {}
    unmapped = collections.defaultdict(list)
    for species_name, prot_reqs in proteomes.items():
        proteome = set(rec.id.split('||', 1)[0] for rec in prot_reqs)
        hogmap_address = os.path.join(folder, species_name + "." + fasta_format_keep + ".hogmap")
        cur_species_hogmap = collections.defaultdict(list)
        with open(hogmap_address, 'rt') as hogmap_file:
            reader = csv.DictReader((line for line in hogmap_file if not line.startswith('!')),
                                    dialect="excel-tab")
            for row in reader:
                if row['hogid'] == "N/A":
                    unmapped[species_name].append(row['qseqid'])
                elif row['qseqid'] not in proteome:
                    logger.debug(f"Ignoring {row} [{species_name}] from hogmap as not in proteome")
                else:
                    cur_species_hogmap[row['qseqid']].append(
                        HOGMapData(row['hogid'], row['family_p'], row['qseqlen'], row['subfamily_medianseqlen']))
        hogmaps[species_name] = cur_species_hogmap
        logger.info("hogmap %s: %d proteins mapped to %d hogs, %d not mapped",
                        species_name, len(cur_species_hogmap), sum(len(x) for x in cur_species_hogmap.values()),
                        len(unmapped[species_name]))

    logger.info("There are %d species in the hogmap folder.", len(hogmaps))
    return hogmaps, unmapped


def group_prots_roothogs(hogmaps):
    """
    function for finding those proteins that are mapped to the same rootHOG.
    output: rhogs_dic
    """

    rhogs_prots = collections.defaultdict(list)
    for species_name, prots_map in hogmaps.items():
        # prot_recs = prot_recs_all[species_name]

        for prot_id, prot_map in prots_map.items():
            # omamer output is sorted based on normcount (not family-p). but that's intended by omamer developers
            # first hit (highest score) -> get roothog
            hogid = prot_map[0].hogid
            rhogid = hogid.split(".")[0].split(":")[1]
            # add to species, prot_id to list for best roothog
            rhogs_prots[rhogid].append((species_name, prot_id))

    logger.info("There are " + str(len(rhogs_prots)) + " rootHOGs.")
    # from old code, to check later: Keep prot seq in hog class. We can write species idx and prot idx to  improve speed of code for omamer thresholding
    return rhogs_prots


def handle_singleton(rhogs_prots, hogmaps, conf_infer_roothogs):
    query_singleton_rhog = []  # they are the only one to hit to this HOG
    for rhog, prot_sp_list in rhogs_prots.items():  # 'C0602129': [('GORGO', 'GORGO00025'), ('GORGO', 'GORGO00026')],
        if len(prot_sp_list) == 1:
            query_singleton_rhog += prot_sp_list
    logger.debug("There are " + str(len(query_singleton_rhog)) + " singleton. The hog on which only one query proteins mapped.")

    # adding singleton to already rhogs with >1 members
    query_singleton_rhog_solved = []
    for (species, prot) in query_singleton_rhog:
        prot_maps = hogmaps[species][prot]
        if len(prot_maps) > 1:
            #prot_maps2 = prot_maps
            # omamer output is sorted based on normcount. but that's intended by omamer developers
            # family_p	family_count	family_normcount
            best, minor_mappings = prot_maps[0], prot_maps[1:]
            # consistency check: assert that best roothog id is indeed singleton
            rhogid0 = best.hogid.split(".")[0].split(":")[-1]
            if len(rhogs_prots[rhogid0]) != 1:
                logger.error(f"singleton roothog is not singleton: {species}, {prot}: rhog: {rhogid0}: {rhogs_prots[rhogid0]}")
                raise RuntimeError("singleton roothog is not singleton")

            for cand_map in minor_mappings:
                rhogid = cand_map.hogid.split(".")[0].split(":")[-1]
                if rhogid in rhogs_prots:
                    if len(rhogs_prots[rhogid]) > 1 and float(cand_map.score) > conf_infer_roothogs.mergHOG_fscore_thresh: # todo  check
                        rhogs_prots[rhogid].append((species, prot))
                        del rhogs_prots[rhogid0]
                        query_singleton_rhog_solved.append((species, prot))
                        break
    logger.debug("We add " + str(
        len(query_singleton_rhog_solved)) + " proteins/singleton to another rHOG based on omamer multi-hit.")

    query_singleton_rhog_solved_set = set(query_singleton_rhog_solved)
    query_singleton_rhog_remained = [i for i in query_singleton_rhog if i not in query_singleton_rhog_solved_set]
    logger.debug("However, " + str(len(query_singleton_rhog_remained)) + " proteins/singletons are remained.")

    # try to create group from multi-hits of remained singleton
    dic_singlton_remained = collections.defaultdict(set)
    for (species, prot) in query_singleton_rhog_remained:
        prot_maps = hogmaps[species][prot]
        for cand_map in prot_maps[1:]:
            if float(cand_map.score) > conf_infer_roothogs.mergHOG_fscore_thresh:
                rhogid = cand_map.hogid.split(".")[0].split(":")[-1]
                dic_singlton_remained[rhogid].add((species, prot))
    logger.debug(f"These are associated to {len(dic_singlton_remained)} HOGs considering all multi-hits.")

    # get lookup from prot -> roothog (to be used for deleting a new cluster assignment)
    prots_rhogs_dic = {sp_prots[0]: rhog for rhog, sp_prots in rhogs_prots.items() if len(sp_prots) == 1}
    clusters = UnionFind()
    for rhog, prot_set in dic_singlton_remained.items():
        clusters.union(*prot_set)

    nr_new_rhogs, nr_prot_in_new_rhogs = 0, 0
    for cc in clusters.get_components():
        if len(cc) > 1:
            nr_new_rhogs += 1
            # check that indeed the proteins are singletons in the rhog_prots dict
            if any(len(rhogs_prots.get(prots_rhogs_dic[prot], [])) > 1 for prot in cc):
                raise RuntimeError("issue 15219325")
            for prot in cc:
                rhog = prots_rhogs_dic[prot]
                if rhog in rhogs_prots:
                    del rhogs_prots[rhog]
            cc = list(cc)
            rhog = prots_rhogs_dic[cc[0]]
            rhogs_prots[rhog] = cc
            nr_prot_in_new_rhogs += len(cc)

    logger.debug(f"We merged {nr_prot_in_new_rhogs} in {nr_new_rhogs} additional groups with >1 members.")

    counter_rhog_singlton_final = 0
    counter_notsingl_final = 0
    for rhogid, spec_prot_list in rhogs_prots.items():
        if len(set(spec_prot_list)) > 1:
            counter_notsingl_final += 1
        else:
            counter_rhog_singlton_final += 1
    logger.debug("Now, we have " + str(counter_notsingl_final) + " rootHOGs with >1 proteins and  and " + str(
        counter_rhog_singlton_final) + " singleton rootHOGs ")

    return rhogs_prots


def filter_big_roothogs(hogmaps, rhogs_prots, conf_infer_roothogs):


    prots_list_big = []  # as backup
    rhogids_big = []
    for rhogid, sp_prot_list in rhogs_prots.items():
        if len(sp_prot_list) > conf_infer_roothogs.big_rhog_size:
            prots_list_big.append(sp_prot_list)  # (species_name,prot_id)
            rhogids_big.append(rhogid)
            logger.info(
                "a big rootHOG was found " + rhogid + " with " + str(len(sp_prot_list)) + " query proteins.")
    logger.info("There are " + str(len(rhogids_big)) + " big rootHOGs.")

    for rhogid in rhogids_big:

        sp_prot_list = rhogs_prots[rhogid]
        sp_prot_list_filt = []
        hogids = []

        # removing proteins with low omamer score from big rootHOG
        for species_name, prot_id in sp_prot_list:
            prot_maps = hogmaps[species_name][prot_id]
            # if len(prot_maps) > 1:  # probably for big rootHOG, there won't be multi-hits
            #     scores = [float(i[1]) for i in prot_maps]  # (hogid,score,seqlen,subfamily_medianseqlen)
            #     hogids = [i[0] for i in prot_maps]
            #     max_score = max(scores)
            #     max_index = scores.index(max_score)
            #     hogid = hogids[max_index]
            # else:
            #     hogid = prot_maps[0][0]
            #     max_score = float(prot_maps[0][1])
            # todo: highest normcount or pavalue , default of omamer ?
            hogid = prot_maps[0][0]
            max_score = float(prot_maps[0][1])

            if max_score > conf_infer_roothogs.big_fscore_thresh:
                sp_prot_list_filt.append((species_name, prot_id))
                hogids.append(hogid)

        logger.info("For big rootHOG " + rhogid + ", " + str(
            len(sp_prot_list_filt)) + " proteins left after filtering with threshold " + str(conf_infer_roothogs.big_fscore_thresh))

        if len(sp_prot_list_filt) < conf_infer_roothogs.big_rhog_size:
            sp_prot_list_filt2 = sp_prot_list_filt

        elif  filter_nonchild_rootHOG :  # removing proteins that are mapped to rootHOG (= HOGC123 not the levelsHOGC123.1a) from big rootHOG
            #hogids2 = []
            sp_prot_list_filt2 = []
            for prot_idx, sp_prot in enumerate(sp_prot_list_filt): #[('UP000192223_224129', 'tr|A0A1W4WAU6|A0A1W4WAU6_AGRPL'), ('UP000192223_224129', 'tr|A0A1W4WU99|A0A1W4WU99_AGRPL'),
                species , protein = sp_prot
                prot_maps =hogmaps[species][protein]
                hogid = prot_maps[0][0]
                if len(hogid.split(".")) > 1: # if mapped to subHOG (not to the rootHOG)
                    sp_prot_list_filt2.append(sp_prot)
                    #hogids2.append(hogid)
        else:
            sp_prot_list_filt2 = sp_prot_list_filt
        logger.info("For big rootHOG " + rhogid + ", " + str(len(sp_prot_list_filt2)) + " proteins left after removing non-child subhogs if activate: "+str(filter_nonchild_rootHOG))

        if len(sp_prot_list_filt2):
            rhogs_prots[rhogid] = sp_prot_list_filt2
        else:
            del rhogs_prots[rhogid]

    return rhogs_prots


def write_rhog(rhogs_prot_records, prot_recs_all, address_rhogs_folder, min_rhog_size=1):
    # max_rhog_size =1e12
    #address_rhogs_folder = folder + "rhog"
    logger.info("Writing Sequences of roothogs are fasta file in " + address_rhogs_folder)
    if not os.path.exists(address_rhogs_folder):
        os.mkdir(address_rhogs_folder)


    rhogid_written_list = []
    for rhogid, rhog_prots in rhogs_prot_records.items():
        rhog_recs = []
        for (species_name, prot_name) in rhog_prots: #
            if prot_name in prot_recs_all[species_name]: # some small prots are removed in the begining min_sequence_length
                prot_rec = prot_recs_all[species_name][prot_name]
                rhog_recs.append(prot_rec)


        if min_rhog_size <= len(rhog_recs):  # <= max_rhog_size:
            # todo add the release id   to file names  rhogids_list[:2] > ['HOG:C0884658', 'HOG:C0709155']
            SeqIO.write(rhog_recs, address_rhogs_folder + "/HOG_" + rhogid + ".fa", "fasta")
            rhogid_written_list.append(rhogid)
        else:
            logger.debug("The roothog " +str(rhogid)+" was too small with size of "+str(len(rhog_recs))+" which is smaller than threshold "+str(min_rhog_size))
            #logger.debug(" issue 12314405 to check" )
            #for prot1 in rhog_recs:
                #logger.debug("we are removing due to omamer signleton hog |*|" + str(prot1.id))

    logger.info("Writing Sequences of " + str(len(rhogid_written_list)) + " roothogs finished.")

    return rhogid_written_list


def find_rhog_candidate_pairs(hogmaps, rhogs_prots, conf_infer_roothogs): # rhogs_prots
    """find pairs of roothogs that could be merged based on the number
    of shared assignments and sizes.
    """
    pair_rhogs_count = collections.defaultdict(int)
    rhogs_size = collections.defaultdict(int)
    for species_name, prt_prot_maps in hogmaps.items():
        for prot, prot_maps in prt_prot_maps.items():
            # prot_maps= [HOGMapData(hogid='HOG:E0315075.2a.2b.13b', score='3115.977589303038', seqlen='574', subfamily_medianseqlen='532'),
            rhogs = []
            for prot_map in prot_maps: #
                # prot_map: HOGMapData record
                if float(prot_map.score) > conf_infer_roothogs.mergHOG_fscore_thresh:
                    rhogid = prot_map.hogid.split(".")[0].split(":")[1]
                    rhogs_size[rhogid] += 1
                    rhogs.append((rhogid, float(prot_map.score)))

            rhogs.sort()  # once sorted, iteration will yield sorted tuples in pairs
            for (hogi, score_i), (hogj, score_j) in itertools.combinations(rhogs, 2):
                # skip pairs of hogs if one of them has not been identified as the best
                # scoring hog in omamer. this could lead to big hairball graphs
                if not (hogi in rhogs_prots and hogj in rhogs_prots):
                    continue
                pair = (hogi, hogj)
                pair_rhogs_count[pair] += 1

    print(len(pair_rhogs_count))
    logger.debug("There are " + str(len(pair_rhogs_count)) + " pairs of rhogs.")

    candidates_pair = []
    for (hogi, hogj), count_shared in pair_rhogs_count.items():
        if hogi in rhogs_size and hogj in rhogs_size:  # during previous functions, we might
            ratioMax = count_shared / max(rhogs_size[hogi], rhogs_size[hogj])
            ratioMin = count_shared / min(rhogs_size[hogi], rhogs_size[hogj])
            # todo: big_rhog_size make it as a ratio, could be problometic when we have many more species -> merging too much
            condition_merge = (max(rhogs_size[hogi], rhogs_size[hogj]) < conf_infer_roothogs.big_rhog_size / 6 and \
                                ((ratioMax > conf_infer_roothogs.mergHOG_ratioMax_thresh or ratioMin > conf_infer_roothogs.mergHOG_ratioMin_thresh) and count_shared > conf_infer_roothogs.mergHOG_shared_thresh))  \
                              and rhogs_size[hogi] < conf_infer_roothogs.big_rhog_size / 6 and rhogs_size[hogj] < conf_infer_roothogs.big_rhog_size / 6
            if condition_merge:
                if rhogs_size[hogi] >= rhogs_size[hogj]:
                    candidates_pair.append((hogi, hogj))  # bigger first
                else:
                    candidates_pair.append((hogj, hogi))
                # print(hogi,"(",rhogs_size[hogi],")",hogj,"(",rhogs_size[hogj],")", count_shared,round(ratioMax,2),round(ratioMin,2))
            # D0651051 ( 59 ) D0658569 ( 180 ) 62 0.34 1.05 # todo why bigger than 1?
            # D0646495 ( 2 ) D0631227 ( 14 ) 25 1.79 12.5

    logger.info("There are %d candidate pairs of rhogs for merging.", len(candidates_pair))
    return candidates_pair


def cluster_rhogs(candidates_pair):
    cluster_rhogs_obj = UnionFind()
    for pair in candidates_pair:
        cluster_rhogs_obj.union(pair[0], pair[1])
    cluster_rhogs_sets = cluster_rhogs_obj.get_components()

    cluster_rhogs_list = [ ]
    for cl in cluster_rhogs_sets:
        cluster_rhogs_list.append(list(cl))

    return cluster_rhogs_list
#
# def cluster_rhogs_old(candidates_pair): # doesnt work correctly
#     # init
#     all_hog_raw = []
#     for pair in candidates_pair:
#         all_hog_raw.append(pair[0])
#         all_hog_raw.append(pair[1])
#     all_hog = list(set(all_hog_raw))
#
#     allcc = []  # connected compoenets
#     for hog in all_hog:
#         allcc.append([hog])
#
#     dic_where = {}
#     for idx, cc in enumerate(allcc):
#         dic_where[cc[0]] = idx  # in the beginning, each inner list has only on element, a unique hog
#
#     # print(len(all_hog_raw),len(all_hog),len(allcc),allcc[:2])
#     logger.debug("There are " + str(len(all_hog)) + " all_hog.")
#
#     # print(dic_where)
#     for pair in candidates_pair:
#         # print(pair)
#         idx_0 = dic_where[pair[0]]
#         idx_1 = dic_where[pair[1]] # this will be problametic for sequence of pairs
#
#         dic_where[pair[1]] = idx_0
#
#         # print(idx_0)
#         allcc[idx_0] += allcc[idx_1]
#         allcc[idx_1] = []
#
#         # print(dic_where,allcc)
#
#     cluster_rhogs_list = [i for i in allcc if len(i) > 1]
#     # print(cluster_rhogs_list)
#     logger.debug("There are " + str(len(cluster_rhogs_list)) + " cluster_rhogs.")
#
#     # allcc_cleaned_len = [len(i) for i in cluster_rhogs]
#     # len(allcc),len(cluster_rhogs),np.sum(allcc_cleaned_len)
#
#     return cluster_rhogs_list



def merge_rhogs2(hogmaps, rhogs_prots, conf_infer_roothogs):
    logger.debug("started merging  ")

    candidates_pair = find_rhog_candidate_pairs(hogmaps, rhogs_prots, conf_infer_roothogs) # rhogs_prots

    print("There are " + str(len(candidates_pair)) + " candidate pairs of rhogs for merging.")
    cluster_rhogs_list = cluster_rhogs(candidates_pair)

    print("There are " + str(len(cluster_rhogs_list)) + " clusters.")
    logger.debug("There are " + str(len(cluster_rhogs_list)) + " clusters.")

    print("** the recursion limit is "+str( sys.getrecursionlimit()))
    logger.debug("** the recursion limit is "+str( sys.getrecursionlimit()))
    cluster_rhogs_list = cluster_rhogs_nx(cluster_rhogs_list, candidates_pair)

    print("There are " + str(len(cluster_rhogs_list)) + " selected clusters.")

    # cluster_rhogs_list = cluster_rhogs(candidates_pair)

    num_rhog_g1 = 0
    for rhg, list_pr in rhogs_prots.items():
        if len(list_pr) > 1:
            num_rhog_g1 += 1

    logger.debug("There are " + str(num_rhog_g1) + " rhogs (size>1) before merging.")
    print("len(rhogs_prots) is ", len(rhogs_prots))

    file_out_merge =  open("merging_rhogs.txt","w")
    file_out_merge.write("#first column is the host hog, the rest will be merged here.\n")
    for cluster in cluster_rhogs_list:
        cluster_host = min(cluster) # write the host rhog as the first element
        file_out_merge.write(cluster_host + "\t")
        for cluster_i in cluster:
            if cluster_i != cluster_host:
                file_out_merge.write(cluster_i+"\t")
        file_out_merge.write("\n")
    file_out_merge.close()

    counter_merged_prots=0
    newhost= ""
    rhogs_host_size = 0
    for cluster in cluster_rhogs_list:
        if len(cluster)>30:
            print(len(cluster))
        #prots = [rhogs_prots[hog] for hog in cluster]
         #  cluster[0] #
        found = False

        while not found and len(cluster)>0:
            host_hog = min(cluster)
            if host_hog in rhogs_prots:
                found = True
                rhogs_host_size = len(rhogs_prots[host_hog])
            else:
                # logger print rhog
                cluster.remove(host_hog)

        all_prots = []
        all_prots2 =[]
        for hog in cluster:
            # todo check, create smallers clusters
            if len(all_prots) < conf_infer_roothogs.big_rhog_size:
                if hog in rhogs_prots:
                    all_prots += rhogs_prots[hog]
                    del rhogs_prots[hog]
            else:
                if newhost:
                    if len(all_prots2) < conf_infer_roothogs.big_rhog_size:
                        if hog in rhogs_prots:
                            all_prots2 += rhogs_prots[hog]
                            del rhogs_prots[hog]
                    else:
                        print("cluster is very big "+str(len(cluster)))
                else:
                    newhost = hog

        if newhost:
            all_prots_uniq2 = list(set(all_prots2))  # there might be repeated prots
            rhogs_prots[newhost] = all_prots_uniq2

        all_prots_uniq = list(set(all_prots))  # there might be repeated prots
        rhogs_prots[host_hog] = all_prots_uniq
        # merging D0562038 and D0559070
        # tr | C3ZG56 | C3ZG56_BRAFL HOG: D0562038, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0562038, tr | C3ZG56 | C3ZG56_BRAFL HOG: D0559070, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0559070
        if len(all_prots_uniq) > rhogs_host_size:
            counter_merged_prots += len(all_prots_uniq)
        # otherwise, merging didn't help

    print(len(rhogs_prots),counter_merged_prots)

    num_rhog_g1 = 0
    for rhg, list_pr in rhogs_prots.items():
        if len(list_pr) > 1:
            num_rhog_g1 += 1

    logger.debug("There are " + str(num_rhog_g1) + " rhogs (size>1)   by merging "+str(counter_merged_prots)+" proteins in total.")

    return rhogs_prots


#
# def merge_rhogs(hogmaps, rhogs_prots, conf_infer_roothogs):
#     logger.debug("started merging  ")
#
#     candidates_pair = find_rhog_candidate_pairs(hogmaps, rhogs_prots)
#
#     cluster_rhogs_list = cluster_rhogs(candidates_pair)
#     num_rhog_g1 = 0
#     for rhg, list_pr in rhogs_prots.items():
#         if len(list_pr) > 1:
#             num_rhog_g1 += 1
#
#     logger.debug("There are " + str(num_rhog_g1) + " rhogs (size>1) before merging.")
#     print(len(rhogs_prots))
#
#     file_out_merge =  open("merging_rhogs.txt","w")
#     file_out_merge.write("#first column is the host hog, the rest will be merged here.\n")
#     for cluster in cluster_rhogs_list:
#         cluster_host = min(cluster) # write the host rhog as the first element
#         file_out_merge.write(cluster_host + "\t")
#         for cluster_i in cluster:
#             if cluster_i != cluster_host:
#                 file_out_merge.write(cluster_i+"\t")
#         file_out_merge.write("\n")
#     file_out_merge.close()
#
#     counter_merged_prots=0
#     newhost= ""
#     for cluster in cluster_rhogs_list:
#         if len(cluster)>30:
#             print(len(cluster))
#         #prots = [rhogs_prots[hog] for hog in cluster]
#          #  cluster[0] #
#         found = False
#         while not found and len(cluster)>0:
#             host_hog = min(cluster)
#             if host_hog in rhogs_prots:
#                 found = True
#                 rhogs_host_size = len(rhogs_prots[host_hog])
#             else:
#                 # logger print rhog
#                 cluster.remove(host_hog)
#
#         all_prots = []
#         all_prots2 =[]
#         for hog in cluster:
#             # todo check, create smallers clusters
#             if len(all_prots) < conf_infer_roothogs.big_rhog_size:
#                 if hog in rhogs_prots:
#                     all_prots += rhogs_prots[hog]
#                     del rhogs_prots[hog]
#             else:
#                 if newhost:
#                     if len(all_prots2) < conf_infer_roothogs.big_rhog_size:
#                         if hog in rhogs_prots:
#                             all_prots2 += rhogs_prots[hog]
#                             del rhogs_prots[hog]
#                     else:
#                         print("cluster is very big "+str(len(cluster)))
#                 else:
#                     newhost = hog
#
#         if newhost:
#             all_prots_uniq2 = list(set(all_prots2))  # there might be repeated prots
#             rhogs_prots[newhost] = all_prots_uniq2
#
#         all_prots_uniq = list(set(all_prots))  # there might be repeated prots
#         rhogs_prots[host_hog] = all_prots_uniq
#         # merging D0562038 and D0559070
#         # tr | C3ZG56 | C3ZG56_BRAFL HOG: D0562038, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0562038, tr | C3ZG56 | C3ZG56_BRAFL HOG: D0559070, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0559070
#         if len(all_prots_uniq) > rhogs_host_size:
#             counter_merged_prots += len(all_prots_uniq)
#         # otherwise, merging didn't help
#
#     print(len(rhogs_prots),counter_merged_prots)
#
#     num_rhog_g1 = 0
#     for rhg, list_pr in rhogs_prots.items():
#         if len(list_pr) > 1:
#             num_rhog_g1 += 1
#
#     logger.debug("There are " + str(num_rhog_g1) + " rhogs (size>1)   by merging "+str(counter_merged_prots)+" proteins in total.")
#
#     return rhogs_prots


def collect_unmapped_singleton(rhogs_prots, unmapped,prot_recs_all,unmapped_singleton_fasta= "singleton_unmapped.fa"):
    unmapped_recs = []
    for species_name, prot_names in unmapped.items():
        for prot_name in prot_names:
            if prot_name in prot_recs_all[species_name]:  # some small prots are removed in the begining min_sequence_length
                prot_rec = prot_recs_all[species_name][prot_name]
                unmapped_recs.append(prot_rec)

    print(len(unmapped_recs))
    singleton_recs = []
    for rhogid, sp_prot_list in rhogs_prots.items():
        if len(sp_prot_list) == 1:
            species_name = sp_prot_list[0][0]
            prot_name = sp_prot_list[0][1]
            prot_rec = prot_recs_all[species_name][prot_name]
            singleton_recs.append(prot_rec)
    print(len(singleton_recs))

    singleton_unmapped_recs = unmapped_recs + singleton_recs

    SeqIO.write(singleton_unmapped_recs, unmapped_singleton_fasta , "fasta")

    return len(singleton_unmapped_recs)


import subprocess


def run_linclust(fasta_to_cluster="singleton_unmapped.fa"):   # todo move run_linclust to _wrapper.py
    num_threads = 10  # todo how to assign more cpu for this step in nextflow
    command_clust = mmseqs_executable_path +" easy-linclust --threads " + str(
        num_threads) + " " + fasta_to_cluster + " singleton_unmapped tmp_linclust"

    logger.debug("linclust rooting started" + command_clust)
    process = subprocess.Popen(command_clust.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # if verbose:
    #    print("output:\n", output)
    #    print("error:\n", error)

    # if "Error analyzing file" in str(output) or error:
    #    try:

    return "done"


def write_clusters(address_rhogs_folder, min_rhog_size):
    cluster_output_address = "singleton_unmapped_all_seqs.fasta"
    cluster_file = open(cluster_output_address, 'r')
    # cluster_dic = {}
    previous_line_start = " "
    # parsing memseqs fasta-like format https://github.com/soedinglab/MMseqs2/wiki#cluster-fasta-like-format
    # >A0A0R4IFW6
    # >tr|A0A0R4IFW6|A0A0R4IFW6_DANRE||DANRE||1000001584 tr|A0A0R
    # MSRKTTSKRHYKPSSEIDDAALARKREYW
    clusters = []
    cluster = []
    # line_strip = " "
    # line_number =1
    for line in cluster_file:
        line_strip = line.strip()
        # print(line_number,previous_line_start,line_strip[0],cluster )# "_",previous_line_start,line_strip[0])
        if previous_line_start == " ":
            clusters = []
            cluster = []
        elif line_strip[0] == ">" and previous_line_start == ">":  # new cluster started at previous line
            if len(cluster) >= 5:  # more than one records [ID1, seq1,ID2, seq2, newcluster_id]
                clusters.append(cluster[:-1])  # last item is the id of new cluster
            cluster = [line_strip]
        elif line_strip[0] == ">" and previous_line_start != ">":
            cluster.append(line_strip)
        elif line_strip[0] != ">":
            cluster.append(line_strip)

        previous_line_start = line_strip[0]
        # line_number+=1

    # for last cluster
    if len(cluster) >= 4:  # more than one records  [ID1, seq1,ID2, seq2]
        clusters.append(cluster)

    logger.debug("Number of linclust clusters raw is " + str(len(clusters)))
    # the record id is parsed by mmseqs very ad hoc. so better use the fasta-like file
    # cluster_output_address = "singleton_unmapped_cluster.tsv"
    # cluster_file = open(cluster_output_address, 'r')
    # cluster_dic = {}
    # for line in cluster_file:
    #     line_strip = line.strip()
    #     rep, prot= line_strip.split()
    #     if rep in cluster_dic:
    #         cluster_dic[rep].append(prot) # the frist line includ (rep,rep)
    #     else:
    #         cluster_dic[rep]=[prot]
    # cluster_list = []
    # for rep, prot_list in cluster_dic.items():
    #     if len(prot_list)>1:
    #         cluster_list.append(prot_list)


    cluster_iter= 1000*10
    for cluster in clusters:
        if len(cluster) >= 2 * min_rhog_size:
            species_names= []     # it seems that genes from same species tend to cluster together, discard such clusters
            for pr_idi in  range(int(len(cluster)/2)):
                species_name = cluster[pr_idi*2].split(" ")[0].split("||")[1]
                species_names.append(species_name)
            if len(set(species_names))>1:
                file_idx = open(address_rhogs_folder + "/HOG_clust" + str(cluster_iter) + ".fa", "w")
                for line in cluster:
                    file_idx.write(line + "\n")
                file_idx.close()
                cluster_iter+=1

    return cluster_iter-1


def parse_isoform_file(species_names, folder=None):
    if folder is None:
        folder = "./splice"
    isoform_by_gene_all = {}

    for species_name in species_names:  # from fasta file
        file_splice_name = os.path.join(folder,  species_name + ".splice")
        if os.path.isfile(file_splice_name):
            # from OMArk
            isoform_by_gene = []
            with open(file_splice_name) as handle:
                for line in handle.readlines():
                    line = line.strip('\n')
                    splice = line.split(";")
                    isoform_by_gene.append(splice)
        else:
            isoform_by_gene = []
        isoform_by_gene_all[species_name] = isoform_by_gene
    return isoform_by_gene_all


def find_nonbest_isoform(species_names, isoform_by_gene_all, hogmaps):
    isoform_selected = {}
    isoform_not_selected = {}

    for species_name in species_names:
        isoform_not_selected[species_name] = []
        isoform_selected[species_name] = []
        isoform_list_sp = isoform_by_gene_all[species_name]
        hogmaps_sp = hogmaps[species_name]
        for isoform_list in isoform_list_sp:
            score_isof_best = 0
            protname_best = isoform_list[0]  # to keep the order if all of the isoforms' omamer score are 'na'
            for isoform_name in isoform_list:
                if isoform_name in hogmaps_sp:

                    prot_maps = hogmaps_sp[isoform_name]
                    if len(prot_maps) > 1:  # for multi-hit omamer output with -n
                        scores = [float(i[1]) for i in prot_maps]  # (hogid,score,seqlen,subfamily_medianseqlen)
                        # hogids = [i[0] for i in prot_maps]
                        # seq_lens = [i[2] for i in prot_maps]
                        # subf_med_lens = [i[3] for i in prot_maps]
                        family_score = max(scores)
                        # idx_max = scores.index(family_score)
                        # seq_len = seq_lens[idx_max]
                        # subf_med_len = subf_med_lens[idx_max]
                        # hogid = hogids[max_index]
                    else:
                        # hogid = prot_maps[0][0]
                        family_score = float(prot_maps[0][1])
                        # seq_len = prot_maps[0][2]
                        # subf_med_len = prot_maps[0][3]
                    score_isof = float(family_score) # for omamer v2 fscore is enough * min(int(seq_len), int(subf_med_len))
                    # print(isoform_name, score_isof,score_isof_best)

                    if score_isof >= score_isof_best:  # when there is a tie, the last one is selected!
                        protname_best = isoform_name
                        score_isof_best = score_isof

            isoform_selected[species_name].append(protname_best)
            isoform_not_selected[species_name] += [i for i in isoform_list if i != protname_best]  # flatten for all

    return isoform_selected, isoform_not_selected


def write_isoform_selected(isoform_by_gene_all, isoform_selected, prot_recs_lists):
    """
    write isoforms

    """
    selected_isoforms_folder = "selected_isoforms/"
    if not os.path.exists(selected_isoforms_folder):
        os.mkdir(selected_isoforms_folder)

    for species, isoform_selected_sp in isoform_selected.items():
        isoform_by_gene = isoform_by_gene_all[species]
        assert len(isoform_by_gene) == len(isoform_selected_sp)

        file_handle = open(selected_isoforms_folder+species + "_selected_isoforms.tsv", "w")
        for gene_idx, isof_list in enumerate(isoform_by_gene):
            isoform_selected_sp_i = isoform_selected_sp[gene_idx]
            file_handle.write(";".join(isof_list) + "\t")
            file_handle.write(isoform_selected_sp_i + "\n")

        all_genes = []
        for list_genes in isoform_by_gene_all[species]:
            all_genes += list_genes
        len(all_genes)

        prot_all = [i.id.split("||")[0] for i in prot_recs_lists[species]]
        len(prot_all), len(all_genes)

        not_in_isoform = set(prot_all) - set(all_genes)

        for isof1 in not_in_isoform:  # for these there is no isoform information, we write them at end of file
            file_handle.write(isof1 + "\t" + isof1 + "\n")

        file_handle.close()

    return 1


def handle_splice(hogmaps,isoform_not_selected):
    hogmaps_selected_isof = {}
    for species, hogmap in hogmaps.items():
        hogmaps_selected_isof[species]={}

        isoform_not_selected_sp = set(isoform_not_selected[species])
        for prot in hogmap:
            if prot not in isoform_not_selected_sp:
                hogmaps_selected_isof[species][prot]=hogmap[prot]

    return hogmaps_selected_isof





def find_outgroup_species(species_tree):
    dic_outgroup_species = {}
    outgroup_species = []

    for node in species_tree.traverse():
        if node.is_leaf() or node.is_root():
            continue

        parent_node = node.up
        outgroup_leaves = []
        if parent_node:
            grandparent_node = parent_node.up
            if grandparent_node:
                # parent_node, grandparent_node,grandparent_node.children
                great_grandparent_node = grandparent_node.up
                if great_grandparent_node:
                    if great_grandparent_node.children:
                        very_distant_cousins = [i for i in great_grandparent_node.children if i != grandparent_node]
                        for i in very_distant_cousins:
                            outgroup_leaves += i.get_leaves()  # a list
                if grandparent_node.children and len(outgroup_leaves) <= 5:
                    distant_cousins = [i for i in grandparent_node.children if i != parent_node]
                    for i in distant_cousins:
                        outgroup_leaves += i.get_leaves()  # a list

            if parent_node.children and len(outgroup_leaves) <= 5:
                cousins = [i for i in parent_node.children if i != node]
                for i in cousins:
                    outgroup_leaves += i.get_leaves()  # a list

        if outgroup_leaves:
            outgroup_species = [i.name for i in outgroup_leaves][:5]

        # print(node.name, outgroup_species)
        dic_outgroup_species[node.name] = outgroup_species

    return dic_outgroup_species


def find_outgroup_prot(dic_outgroup_species, rhog_prot):
    dic_prot_sp = {}
    for sp, prot in rhog_prot:
        if sp in dic_prot_sp:
            dic_prot_sp[sp].append(prot)
        else:
            dic_prot_sp[sp] = [prot]

    dic_outgroup_prot = {}
    for tax, outgroup_species in dic_outgroup_species.items():
        outgroup_prots = []
        for sp in outgroup_species:
            prot_list = dic_prot_sp[sp]
            if len(prot_list) > 1:
                outgroup_prots.append((sp, prot_list[0]))
                outgroup_prots.append((sp, prot_list[1]))
            else:
                outgroup_prots.append((sp, prot_list[0]))

        dic_outgroup_prot[tax] = outgroup_prots

    return dic_outgroup_prot


def write_outgroups(dic_outgroup_prot, rhogid, address_outgorup,prot_recs_all):
    dic_rhog_recs = {}
    for tax, outgroup_prot in dic_outgroup_prot.items():
        rhog_recs = []
        for (species_name, prot_name) in outgroup_prot:
            prot_rec = prot_recs_all[species_name][prot_name]
            rhog_recs.append(prot_rec)
        dic_rhog_recs[tax] = rhog_recs

    for tax, rhog_recs in dic_rhog_recs.items():
        SeqIO.write(rhog_recs, address_outgorup + "/outgroup_" + str(rhogid) + "_" + tax + ".fa", "fasta")

    return 1


def write_outgroups_all(rhogs_prots,prot_recs_all):
    from ete3 import Tree
    sp_folder = "in_folder/"
    address_outgorup = "./outgroup/"
    if not os.path.exists(address_outgorup):
        os.mkdir(address_outgorup)

    print("start working on  outgroup generation for rhogs  "+str(len(rhogs_prots)))
    for rhog, rhog_prot in rhogs_prots.items():
        if len(rhog_prot) > 2:
            # print(rhog)
            species_tree = Tree(sp_folder + "species_tree.nwk", format=1)
            list_species = list(set([i[0] for i in rhog_prot]))
            species_tree.prune(list_species)
            # species_tree.write(format=1)
            # if len(species_tree)>2:
            #    counter+=1
            dic_outgroup_species = find_outgroup_species(species_tree)

            dic_outgroup_prot = find_outgroup_prot(dic_outgroup_species, rhog_prot)

            write_outgroups(dic_outgroup_prot, rhog, address_outgorup,prot_recs_all)
    print("done writing")
    return 1



def cluster_rhogs_nx(cluster_rhogs_list, candidates_pair):
    count=0
    candidates_pair = set(candidates_pair)
    new_cluster = []
    for c in cluster_rhogs_list:
        count+=1

        G = nx.Graph()
        for i in range(len(c)):
            pairA = c[i]

            if not G.has_node(pairA):
                G.add_node(pairA)
            for j in range(i+1, len(c)):

                pairB = c[j]
                if not G.has_node(pairB):
                    G.add_node(pairB)
                if (pairA,pairB) in candidates_pair or (pairB,pairA) in candidates_pair:
                    G.add_edge(pairA,pairB)
        density = nx.density(G) # todo: do we really need this line ?
        if len(G) <1000: # todo make it as a parameter
            # limit number of merging rootHOGs due to O(N^3) of HCS/cc
            newG = HCS(G)
            clusters = [ list(x) for x in nx.connected_components(newG) if len(x)>1]
        else:
            logger.debug("This cluster of rootHOg hit the limit of 1000 rHOGs, we are not merging these rootHOGs.")
            clusters =  [ [i] for i in list(G.nodes)] # do not merge at all

        new_cluster += clusters
    return new_cluster

def highly_connected(G, E):
    """Checks if the graph G is highly connected

    Highly connected means, that splitting the graph G into subgraphs needs more than 0.5*|V| edge deletions
    This definition can be found in Section 2 of the publication.

    :param G: Graph G
    :param E: Edges needed for splitting G
    :return: True if G is highly connected, otherwise False
    """

    return len(E) >= len(G.nodes) / 2


def remove_edges(G, E):
    """Removes all edges E from G

    Iterates over all edges in E and removes them from G
    :param G: Graph to remove edges from
    :param E: One or multiple Edges
    :return: Graph with edges removed
    """

    for edge in E:
        G.remove_edge(*edge)
    return G


def HCS(G):
    """Basic HCS Algorithm
https://github.com/53RT/Highly-Connected-Subgraphs-Clustering-HCS/blob/master/hcs.py
    cluster labels, removed edges are stored in global variables

    :param G: Input graph
    :return: Either the input Graph if it is highly connected, otherwise a Graph composed of
    Subgraphs that build clusters
    """

    E = nx.algorithms.connectivity.cuts.minimum_edge_cut(G)

    if not highly_connected(G, E):
        G = remove_edges(G, E)
        sub_graphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]

        if len(sub_graphs) == 2:
            H = HCS(sub_graphs[0])
            _H = HCS(sub_graphs[1])

            G = nx.compose(H, _H)

    return G
