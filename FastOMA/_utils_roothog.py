
import csv
import itertools
from pathlib import Path
import subprocess
from itertools import chain

import networkx as nx
from Bio import SeqIO
import pickle
from os import listdir
import os
import sys

from .zoo.unionfind import UnionFind
from . import logger
import collections

filter_nonchild_rootHOG= False
mmseqs_executable_path ="mmseqs"

HOGMapData = collections.namedtuple("HOGMapData", ("hogid", "score", "seqlen", "subfamily_medianseqlen"))
Gene = collections.namedtuple("Gene", ("numeric_id", "prot_id", "main_isoform"), defaults=(None,) )
    

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

    species_names.sort()   # sort the species names alphabetically -> reproducible results
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



def create_root_hogs(hogmaps, conf):
    """Create root HOGs from HOG mappings"""
    
    # Initial grouping
    rhogs_prots = group_prots_by_roothogs(hogmaps)
    logger.info(f"Initial grouping: {len(rhogs_prots)} root HOGs")
    
    # Handle singletons
    rhogs_prots = resolve_singletons(rhogs_prots, hogmaps, conf)

    rhogs_prots = merge_rhogs2(hogmaps, rhogs_prots, conf)
    rhogs_prots = filter_big_roothogs(hogmaps, rhogs_prots, conf)
    
    # # Merge similar HOGs
    # if conf.enable_merging:
    #     rhogs_prots = merge_similar_hogs(rhogs_prots, hogmaps, conf)
    #
    # # Filter large HOGs
    # rhogs_prots = filter_large_hogs(rhogs_prots, hogmaps, conf)
    
    return rhogs_prots


def handle_splice_variants(species_names, hogmaps, splice_folder):
    """Handle splice variants and return filtered data"""
    
    # Parse isoform files
    isoform_by_gene = parse_isoform_files(species_names, splice_folder)
    
    # Select best isoforms
    selected_isoforms, non_selected = select_best_isoforms(
        species_names, isoform_by_gene, hogmaps
    )
    
    # Filter hogmaps to keep only selected isoforms
    filtered_hogmaps = filter_hogmaps_by_isoforms(hogmaps, non_selected)
    
    return {
        'isoform_by_gene': isoform_by_gene,
        'selected_isoforms': selected_isoforms,
        'non_selected': non_selected,
        'filtered_hogmaps': filtered_hogmaps
    }


def add_species_name_prot_id(prot_recs_lists):
    """
    adding the name of species to each protein record
        - based on file name
    adding protein idx number, integer needed by xml format
    output:  prot_recs_all =  {'MYCGE': {'sp|P47500|RF1_MYCGE||MYCGE||1000000001': SeqRecord(seq=
    """
    start_num_prot = int(1e9)
    start_num_prot_per_sp = int(1e6) #
    prot_recs_all = {} # {'MYCGE': {'sp|P47500|RF1_MYCGE||MYCGE||1000000001': SeqRecord(seq=
    species_idx = -1
    for species_name, prot_recs_list in prot_recs_lists.items():
        species_idx += 1
        prot_idx = start_num_prot + species_idx * start_num_prot_per_sp
        prot_recs_all[species_name]={}
        for prot_rec in prot_recs_list:
            prot_idx+=1
            prot_name= prot_rec.id

            if len(prot_name) > 230:
                logger.info("We are truncating the prot name as it may be problematic for mafft, " + str(prot_name))
                prot_name = prot_name[:230]

            prot_name_new = prot_name+ "||"+species_name+"||"+str(prot_idx) # orthoxml file needs an integer as
            prot_rec.id = prot_name_new
            prot_recs_all[species_name][prot_name] = prot_rec
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

def group_prots_by_roothogs(hogmaps):
    """Group proteins by their root HOG assignments"""
    
    rhogs_prots = collections.defaultdict(list)
    
    for species_name, prots_map in hogmaps.items():
        for prot_id, prot_mappings in prots_map.items():
            # Get best mapping (first one, as they're sorted by score)
            best_mapping = prot_mappings[0]
            rhogid = extract_root_hog_id(best_mapping.hogid)
            rhogs_prots[rhogid].append((species_name, prot_id))
    
    logger.info(f"Grouped proteins into {len(rhogs_prots)} root HOGs")
    return dict(rhogs_prots)

def extract_root_hog_id(hogid):
    """Extract root HOG ID from full HOG ID"""
    return hogid.split(".")[0].split(":")[1]


def resolve_singletons(rhogs_prots, hogmaps, conf):
    """Resolve singleton HOGs by checking alternative mappings"""
    
    singletons = [(rhog, prots[0]) for rhog, prots in rhogs_prots.items() if len(prots) == 1]
    resolved_count, remained_count = 0, 0
    singletons_remained = collections.defaultdict(set)

    for singleton_rhog, (species, prot) in singletons:
        prot_mappings = hogmaps[species][prot]
        
        # Check alternative mappings
        could_resolve = False
        for alt_mapping in prot_mappings[1:]:
            if float(alt_mapping.score) < conf.mergHOG_fscore_thresh:
                continue
                
            alt_rhog = extract_root_hog_id(alt_mapping.hogid)
            
            # If alternative HOG has multiple members, move protein there
            if alt_rhog in rhogs_prots and len(rhogs_prots[alt_rhog]) > 1:
                rhogs_prots[alt_rhog].append((species, prot))
                del rhogs_prots[singleton_rhog]
                resolved_count += 1
                could_resolve = True
                break
        if not could_resolve:
            remained_count += 1
            for alt_mapping in prot_mappings[1:]:
                if float(alt_mapping.score) < conf.mergHOG_fscore_thresh:
                    continue
                singletons_remained[extract_root_hog_id(alt_mapping.hogid)].add((species, prot))    
    
    logger.info(f"Resolved {resolved_count} singleton HOGs based on omammer multi-hits.")
    logger.info(f"{remained_count} singleton HOGs remained before merging based on minor mappings.")
    
    # get lookup from prot -> roothog (to be used for deleting a new cluster assignment)
    prots_rhogs_dic = {sp_prots[0]: rhog for rhog, sp_prots in rhogs_prots.items() if len(sp_prots) == 1}

    # Merge remaining singletons into groups based on shared alternative mappings of multi-hits
    clusters = UnionFind()
    for rhog, prot_set in singletons_remained.items():
        clusters.union(*prot_set)

    nr_new_rhogs, nr_prot_in_new_rhogs = 0, 0
    for cc in clusters.get_components():
        if len(cc) <= 1:
            continue

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
    logger.info(f"We merged {nr_prot_in_new_rhogs} in {nr_new_rhogs} additional groups with >1 members.")

    counter_rhog_singleton = 0
    counter_not_singleton  = 0
    for rhogid, spec_prot_list in rhogs_prots.items():
        if len(set(spec_prot_list)) > 1:
            counter_not_singleton += 1
        else:
            counter_rhog_singleton += 1
    logger.info(f"Now, we have {counter_not_singleton} rootHOGs with >1 proteins and {counter_rhog_singleton} singleton rootHOGs")
    return rhogs_prots

def save_gene_id_mapping(prot_recs_all, isoform_data=None, out_path="gene_id_dic_xml.pickle"):
    """Save comprehensive gene ID mapping including isoform information"""
    
    gene_mapping = {}
    for species_name, prot_dict in prot_recs_all.items():
        species_data = []

        isoform_to_main = {}
        if isoform_data and species_name in isoform_data['selected_isoforms']:
            isoform_to_main = {gene: main for main, genes in zip(
                    isoform_data['selected_isoforms'][species_name], isoform_data['isoform_by_gene'][species_name]
                ) for gene in genes}
        
        id2numeric = {}
        for prot_name, prot_rec in prot_dict.items():
            parts = prot_rec.id.split('||')
            original_id = parts[0]
            numeric_id = int(parts[2])
            id2numeric[original_id] = numeric_id
        
        species_data = [Gene(num, orig_id, id2numeric.get(isoform_to_main.get(orig_id))) for orig_id, num in id2numeric.items()]
        gene_mapping[species_name] = species_data
    
    # Save to pickle file
    with open(out_path, 'wb') as handle:
        pickle.dump(gene_mapping, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    logger.info(f"Saved comprehensive gene ID mapping to {out_path}")
    return gene_mapping


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
                    rhogid = extract_root_hog_id(prot_map.hogid)
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

    logger.debug(f"There are {len(candidates_pair)} candidate pairs of rhogs for merging.")
    cluster_rhogs_list = cluster_rhogs(candidates_pair)
    logger.debug(f"There are {len(cluster_rhogs_list)} clusters.")
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


def collect_unmapped_singleton(rhogs_prots, unmapped, prot_recs_all, unmapped_singleton_fasta="singleton_unmapped.fa"):
    unmapped_recs = []
    for species_name, prot_names in unmapped.items():
        all_recs_of_species = prot_recs_all[species_name]
        for prot_name in prot_names:
            try:
                unmapped_recs.append(all_recs_of_species[prot_name])
            except KeyError:
                # some small prots are removed in the begining min_sequence_length
                pass
    logger.debug(f"collected sequence records of {len(unmapped_recs)} unmapped proteins.")

    singleton_recs = []
    for rhogid, sp_prot_list in rhogs_prots.items():
        if len(sp_prot_list) == 1:
            species_name = sp_prot_list[0][0]
            prot_name = sp_prot_list[0][1]
            prot_rec = prot_recs_all[species_name][prot_name]
            singleton_recs.append(prot_rec)
    logger.debug(f"collected sequence records of {len(singleton_recs)} singleton proteins.")

    singleton_unmapped_recs = unmapped_recs + singleton_recs
    SeqIO.write(singleton_unmapped_recs, unmapped_singleton_fasta , "fasta")
    logger.debug(f"Wrote {len(singleton_unmapped_recs)} sequences of unmapped and singleton proteins to {unmapped_singleton_fasta}.")
    return len(singleton_unmapped_recs)


def run_linclust(fasta_to_cluster="singleton_unmapped.fa"):   # todo move run_linclust to _wrapper.py .  see easy-cluster below.
    # todo: change FastOMA.nf to assgin more cpus to infer_roothog step.  mmseqs uses all cpus by default  # num_threads = 10  , --threads " + str(num_threads) + "
    command_clust = mmseqs_executable_path +" easy-cluster  " + fasta_to_cluster + " singleton_unmapped tmp_linclust"
    # easy-cluster is much better than easy-linclust but a bit slower. todo: make it as an arugment for user

    logger.debug("clustering rooting started " + command_clust)
    process = subprocess.Popen(command_clust.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    logger.debug(f"clustering rooting finished with returncode {process.returncode}.")
    if process.returncode != 0:
        logger.error(f"mmseqs easy-cluster failed with returncode {process.returncode}.")
        logger.error(f"stdout:\n{output.decode('utf-8')}\n******************")
        logger.error(f"stderr:\n{error.decode('utf-8')}\n******************")
        raise RuntimeError("mmseqs easy-cluster failed.")
    return "singleton_unmapped_all_seqs.fasta"  # the fasta-like file of clusters


def write_clusters(fasta_clustered, address_rhogs_folder, min_rhog_size):
    """parse the fasta-like output of mmseqs linclust and write each cluster

    The format is described here: https://github.com/soedinglab/MMseqs2/wiki#cluster-fasta-like-format
    Example:
        >O67032                                 <- cluster 1
        >sp|O67032|RF1_AQUAE                    <- member1 of cluster 1
        MLKEAYISRLDKLQEKYRKLQEELSKPEVIQD...
        >sp|O84026|RF1_CHLTR                    <- member2 of cluster 1
        MEIKVLECLKRLEEVEKQISDPNIFSNPKEYS...
        >sp|P47500|RF1_MYCGE                    <- member3 of cluster 1
        MDFDKQLFFNVEKIVELTEQLEKDLNKPNLSF...
        >O66429                                 <- cluster 2
        >sp|O66429|EFTU_AQUAE                   <- member1 of cluster 2
        MAKEKFERTKEHVNVGTIGHVDHGKS
        >sp|P0CD71|EFTU_CHLTR                   <- member1 of cluster 1
        MSKETFQRNKPHINIGTIGHVDHGKTTLTAAITRAL
    """

    with open(fasta_clustered, 'rt') as cluster_file:
        previous_line_start = " "
        clusters = []
        cluster = []
        for line in cluster_file:
            line_strip = line.strip()
            if previous_line_start == " ":
                clusters = []
                cluster = []
            elif line_strip[0] == ">" and previous_line_start == ">":  # new cluster started at previous line
                if len(cluster) >= 5:
                    # more than one records [ID1, seq1,ID2, seq2, newcluster_id], i.e. non-trivial cluster
                    clusters.append(cluster[:-1])  # last item is the id of new cluster
                cluster = [line_strip]
            elif line_strip[0] == ">" and previous_line_start != ">":
                cluster.append(line_strip)
            elif line_strip[0] != ">":
                cluster.append(line_strip)
            previous_line_start = line_strip[0]

    # for last cluster
    if len(cluster) >= 4:  # more than one records  [ID1, seq1,ID2, seq2]
        clusters.append(cluster)
    logger.debug(f"Number of non-trivial linclust clusters raw is {len(clusters)}.")

    cluster_iter = 10000
    for cluster in clusters:
        # cluster is like ['>sp|O67032|RF1_AQUAE', 'MLKEAYISRLDKLQEKYRKLQEELSKPEVIQD...', ...], i.e. an id and its sequence.
        # thus, check for double the min_rhog_size to have at least min_rhog_size proteins
        if len(cluster) >= 2 * min_rhog_size:
            species_names = set([])     # it seems that genes from same species tend to cluster together, discard such clusters
            for pr_idi in range(int(len(cluster)/2)):
                species_name = cluster[pr_idi*2].split(" ")[0].split("||")[1]
                species_names.add(species_name)
            if len(species_names) > 1:
                with open(address_rhogs_folder + "/HOG_clust" + str(cluster_iter) + ".fa", "w") as file_idx:
                    for line in cluster:
                        file_idx.write(line + "\n")
                cluster_iter += 1
    return cluster_iter - 10000


def parse_isoform_files(species_names, folder=None):
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


def select_best_isoforms(species_names, isoform_by_gene_all, hogmaps):
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
            isoform_not_selected[species_name].extend([i for i in isoform_list if i != protname_best])  # flatten for all

    return isoform_selected, isoform_not_selected


def write_selected_isoforms(isoform_data, prot_recs_lists, output_folder="selected_isoforms"):
    """
    write isoforms to a tsv file for each species
    """
    selected_isoforms_folder = Path(output_folder)
    selected_isoforms_folder.mkdir(exist_ok=True)

    isoform_selected = isoform_data['selected_isoforms']
    isoform_by_gene_all = isoform_data['isoform_by_gene']
    
    for species, isoform_selected_sp in isoform_selected.items():
        isoform_by_gene = isoform_by_gene_all[species]
        assert len(isoform_by_gene) == len(isoform_selected_sp)

        main_to_minor = {
            main: set(gene_isoforms) - {main}
            for main, gene_isoforms in zip(isoform_selected_sp, isoform_by_gene)
        }
        all_isoforms = set(chain.from_iterable(isoform_by_gene))
        fpath = selected_isoforms_folder / f"{species}_selected_isoforms.fa"
        with open(fpath, "wt") as fout:
            for rec in prot_recs_lists[species]:
                id_ = rec.id.split("||")[0]
                if id_ in main_to_minor or id_ not in all_isoforms:
                    if id_ in main_to_minor:
                        desc_ = "minor_isoforms=" + ",".join(main_to_minor.get(id_, []))
                    else:
                        desc_ = ""
                    fout.write(f">{id_} {desc_}\n{rec.seq}\n")


def filter_hogmaps_by_isoforms(hogmaps, non_selected_isoforms):
    hogmaps_selected_isof = {}
    for species, hogmap in hogmaps.items():
        hogmaps_selected_isof[species]={}

        non_selected_isoforms_sp = set(non_selected_isoforms[species])
        for prot in hogmap:
            if prot not in non_selected_isoforms_sp:
                hogmaps_selected_isof[species][prot]=hogmap[prot]

    return hogmaps_selected_isof


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
