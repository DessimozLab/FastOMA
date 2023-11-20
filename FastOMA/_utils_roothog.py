import collections
import csv

from Bio import SeqIO
import pickle
from os import listdir
import os

from ._utils_subhog import logger_hog
from . import _config
HOGMapData = collections.namedtuple("HOGMapData", ("hogid", "score", "seqlen", "subfamily_medianseqlen"))


def parse_proteomes(folder=None):  # list_oma_species
    """
    parsing fasta files of proteins located in /proteome/
    using Bio.SeqIO.parse
    Each fasta file is for one species.  The file name is the species name.
    All proteoems should be with the same extension.
    output: prot_recs_lists: a dic with key as species name and  its value list of Biopython record of species.
    """
    if folder is None:
        folder = "./proteome"

    fasta_format_keep = ""
    project_files = listdir(folder)
    logger_hog.debug("using '%s' as proteome folder, found %d files", folder, len(project_files))
    species_names = []  # query/input species name based on the file name
    for file in project_files:
        species_name, ext = file.rsplit('.', 1)
        logger_hog.debug("%s: %s | %s", file, species_name, ext)
        if ext in ("fa", "fasta"):
            species_names.append(species_name)
            fasta_format_keep = ext  # last one is stored either fa or fasta

    # todo accept all fasta formats in the input prtoeome folder, fasta, fa, fna, ..
    prot_recs_lists = {} # key: species name, value is a dic of query protein Biopython records. # 'MYCGE': [SeqRecord(seq=Seq('MDFDK
    for species_name in species_names:
        prot_address = os.path.join(folder, species_name + "." + fasta_format_keep)
        prots_record = list(SeqIO.parse(prot_address, "fasta"))
        prot_recs_lists[species_name] = prots_record

    logger_hog.info("There are %d species in the proteome folder.", len(species_names))
    return species_names, prot_recs_lists, fasta_format_keep


def add_species_name_prot_id(species_names, prot_recs_lists):
    """
    adding the name of species to each protein record
        - based on file name
    adding protein idx number, integer needed by xml format
    output:  prot_recs_all =  {'MYCGE': {'sp|P47500|RF1_MYCGE||MYCGE||1000000001': SeqRecord(seq=
    """
    prot_idx_name_pickle_file = "./gene_id_dic_xml.pickle"
    start_num_prot = int(1e9)
    start_num_prot_per_sp = int(1e6)
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
                logger_hog.info("We are truncating the prot name as it may be problamatic for mafft, " + str(prot_name))
                prot_name = prot_name[:230]

            # todo, this could be a dic
            prot_idx_name[species_name].append((prot_idx, prot_name))

            prot_name_new = prot_name+ "||"+species_name+"||"+str(prot_idx) # orthoxml file needs an integer as
            prot_rec.id = prot_name_new
            prot_recs_all[species_name][prot_name] = prot_rec

    with open(prot_idx_name_pickle_file, 'wb') as handle:
        pickle.dump(prot_idx_name, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return prot_recs_all


def parse_hogmap_omamer(species_names, fasta_format_keep,  folder=None):
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
    for species_name in species_names:
        hogmap_address = os.path.join(folder, species_name + "." + fasta_format_keep + ".hogmap")
        cur_species_hogmap = collections.defaultdict(list)
        with open(hogmap_address, 'rt') as hogmap_file:
            reader = csv.DictReader((line for line in hogmap_file if not line.startswith('!')),
                                    dialect="excel-tab")
            for row in reader:
                if row['hogid'] == "N/A":
                    unmapped[species_name].append(row['qseqid'])
                else:
                    cur_species_hogmap[row['qseqid']].append(
                        HOGMapData(row['hogid'], row['family_p'], row['qseqlen'], row['subfamily_medianseqlen']))
        hogmaps[species_name] = cur_species_hogmap
        logger_hog.info("hogmap %s: %d proteins mapped to %d hogs, %d not mapped",
                        species_name, len(cur_species_hogmap), sum(len(x) for x in cur_species_hogmap.values()),
                        len(unmapped[species_name]))

    logger_hog.info("There are %d species in the hogmap folder.", len(hogmaps))
    return hogmaps, unmapped


def group_prots_roothogs(hogmaps):
    """
    function for finding those proteins that are mapped to the same rootHOG.

    output: rhogs_dic
    """

    rhogs_prots = {}
    for species_name, prots_map in hogmaps.items():
        # prot_recs = prot_recs_all[species_name]

        for prot_id, prot_map in prots_map.items():
            # omamer output is sorted based on normcount (not family-p). but that's intended by omamer developers
            # this helps me in other functions like handle_singleton in this
            #  this should be commented
            # if len(prot_map)>1:
            #     scores = [float(i[1]) for i in prot_map] # (hogid,score,seqlen,subfamily_medianseqlen)
            #     rhogids =[i[0] for i in prot_map]
            #     # select hog with the best score
            #     max_index = scores.index(max(scores))
            #     hogid =rhogids[max_index]
            # else:
            #     hogid = prot_map[0][0]
            hogid = prot_map[0][0]
            rhogid = hogid.split(".")[0].split(":")[1]

            if rhogid in rhogs_prots:
                rhogs_prots[rhogid].append((species_name, prot_id))
            else:
                rhogs_prots[rhogid] = [(species_name, prot_id)]
    logger_hog.info("There are " + str(len(rhogs_prots)) + " rootHOGs.")
    # from old code, to check later: Keep prot seq in hog class. We can write species idx and prot idx to  improve speed of code for omamer tresholidng
    return rhogs_prots


def handle_singleton(rhogs_prots,hogmaps):
    query_singleton_rhog = []  # they are the only one to hit to this HOG
    for rhog, prot_sp_list in rhogs_prots.items():  # 'C0602129': [('GORGO', 'GORGO00025'), ('GORGO', 'GORGO00026')],
        if len(prot_sp_list) == 1:
            query_singleton_rhog += prot_sp_list
    logger_hog.debug(
        "There are " + str(len(query_singleton_rhog)) + " singleton. The hog on which only one query proteins mapped.")

    # adding singleton to already rhogs with >1 members
    query_singleton_rhog_solved = []
    for (species, prot) in query_singleton_rhog:
        prot_maps = hogmaps[species][prot]
        if len(prot_maps) > 1:
            prot_maps2 = prot_maps
            # omamer output is sorted based on normcount. but that's intended by omamer developers
            # family_p	family_count	family_normcount
            scores = [float(i[1]) for i in prot_maps]  # (hogid,score,seqlen,subfamily_medianseqlen)
            hogids = [i[0] for i in prot_maps]
            for hogid in hogids[1:]:  # the 0-index is already checked and is a singleton
                rhogid = hogid.split(".")[0].split(":")[1]
                if rhogid in rhogs_prots:
                    if len(rhogs_prots[rhogid]) > 1:
                        rhogs_prots[rhogid].append((species, prot))
                        rhogid0 = hogids[0].split(".")[0].split(":")[1]  # 'HOG:D0903929.1a'
                        del rhogs_prots[rhogid0]
                        query_singleton_rhog_solved.append((species, prot))
                        break
    logger_hog.debug("We add " + str(
        len(query_singleton_rhog_solved)) + " proteins/singleton to another rHOG based on omamer multi-hit.")

    query_singleton_rhog_solved_set = set(query_singleton_rhog_solved)
    query_singleton_rhog_remained = [i for i in query_singleton_rhog if i not in query_singleton_rhog_solved_set]
    logger_hog.debug("However, " + str(len(query_singleton_rhog_remained)) + " proteins/singletons are remained.")

    # try to create group from multi-hits of remained singleton
    dic_singlton_remained = {}
    for (species, prot) in query_singleton_rhog_remained:
        prot_maps = hogmaps[species][prot]
        scores = [float(i[1]) for i in prot_maps]  # (hogid,score,seqlen,subfamily_medianseqlen)
        hogids = [i[0] for i in prot_maps]
        for hogid in hogids[1:]:
            rhogid = hogid.split(".")[0].split(":")[1]
            if rhogid in dic_singlton_remained:
                dic_singlton_remained[rhogid].append((species, prot))
            else:
                dic_singlton_remained[rhogid] = [(species, prot)]
    logger_hog.debug("These are associated to " + str(len(dic_singlton_remained)) + " HOGs considering all multi-hits.")

    count_new_rhogids = 0
    count_updated_rhogids = 0
    already_grouped = []
    for rhogid, spec_prot_list in dic_singlton_remained.items():
        if len(set(spec_prot_list)) > 1:
            already_grouped_set = set(already_grouped)
            val = sum([1 for i in spec_prot_list if i in already_grouped_set])

            if val == 0:  # the prot has not yet grouped
                if rhogid in rhogs_prots:
                    count_updated_rhogids += 1
                else:
                    count_new_rhogids += 1
                rhogs_prots[rhogid] = spec_prot_list
                already_grouped += spec_prot_list
    logger_hog.debug("We updated " + str(count_updated_rhogids) + " and created  " + str(
        count_new_rhogids) + " new HOGs  which are not singleton anymore.")

    counter_rhog_singlton_final = 0
    counter_notsingl_final = 0
    for rhogid, spec_prot_list in rhogs_prots.items():
        if len(set(spec_prot_list)) > 1:
            counter_notsingl_final += 1
        else:
            counter_rhog_singlton_final += 1

    logger_hog.debug("Now, we have " + str(counter_notsingl_final) + " rootHOGs with >1 proteins and  and " + str(
        counter_rhog_singlton_final) + " singleton rootHOGs ")

    return rhogs_prots


def filter_big_roothogs(hogmaps, rhogs_prots):


    prots_list_big = []  # as backup
    rhogids_big = []
    for rhogid, sp_prot_list in rhogs_prots.items():
        if len(sp_prot_list) > _config.big_rhog_size:
            prots_list_big.append(sp_prot_list)  # (species_name,prot_id)
            rhogids_big.append(rhogid)
            logger_hog.info(
                "a big rootHOG was found " + rhogid + " with " + str(len(sp_prot_list)) + " query proteins.")
    logger_hog.info("There are " + str(len(rhogids_big)) + " big rootHOGs.")

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

            if max_score > _config.omamer_family_threshold:
                sp_prot_list_filt.append((species_name, prot_id))
                hogids.append(hogid)

        logger_hog.info("For big rootHOG " + rhogid + ", " + str(
            len(sp_prot_list_filt)) + " proteins left after filtering with threshold " + str(_config.omamer_family_threshold))

        if len(sp_prot_list_filt) < _config.big_rhog_size:
            sp_prot_list_filt2 = sp_prot_list_filt

        else:  # removing proteins that are mapped to rootHOG (= HOGC123 not the levelsHOGC123.1a) from big rootHOG
            #hogids2 = []
            sp_prot_list_filt2 = []
            for prot_idx, sp_prot in enumerate(sp_prot_list_filt): #[('UP000192223_224129', 'tr|A0A1W4WAU6|A0A1W4WAU6_AGRPL'), ('UP000192223_224129', 'tr|A0A1W4WU99|A0A1W4WU99_AGRPL'),
                species , protein = sp_prot
                prot_maps =hogmaps[species][protein]
                hogid = prot_maps[0][0]
                if len(hogid.split(".")) > 1: # if mapped to subHOG (not to the rootHOG)
                    sp_prot_list_filt2.append(sp_prot)
                    #hogids2.append(hogid)

        logger_hog.info("For big rootHOG " + rhogid + ", " + str(
            len(sp_prot_list_filt2)) + " proteins left after removing non-child subhogs")
        if len(sp_prot_list_filt2):
            rhogs_prots[rhogid] = sp_prot_list_filt2
        else:
            del rhogs_prots[rhogid]

    return rhogs_prots


def write_rhog(rhogs_prot_records, prot_recs_all, address_rhogs_folder, min_rhog_size=1):
    # max_rhog_size =1e12
    #address_rhogs_folder = folder + "rhog"
    logger_hog.info("Writing Sequences of roothogs are fasta file in " + address_rhogs_folder)
    if not os.path.exists(address_rhogs_folder):
        os.mkdir(address_rhogs_folder)


    rhogid_written_list = []
    for rhogid, rhog_prots in rhogs_prot_records.items():
        rhog_recs = []
        for (species_name, prot_name) in rhog_prots:
            prot_rec = prot_recs_all[species_name][prot_name]
            rhog_recs.append(prot_rec)

        if min_rhog_size <= len(rhog_recs):  # <= max_rhog_size:
            # todo add the release id   to file names  rhogids_list[:2] > ['HOG:C0884658', 'HOG:C0709155']
            SeqIO.write(rhog_recs, address_rhogs_folder + "/HOG_" + rhogid + ".fa", "fasta")
            rhogid_written_list.append(rhogid)
        else:
            for prot1 in rhog_recs:
                logger_hog.debug("we are removing due to omamer signleton hog |*|" + str(prot1.id))

    logger_hog.info("Writing Sequences of " + str(len(rhogid_written_list)) + " roothogs finished.")

    return rhogid_written_list


def find_rhog_candidate_pairs(hogmaps, rhogs_prots):
    threshod_f_score_merging = 70
    pair_rhogs_count = {}
    for rhog, prt_prot_maps in hogmaps.items():
        for prot, prot_maps in prt_prot_maps.items():
            # [('HOG:D0017631.5a', '1771.7328874713658', '253', '234'), ('HOG:D0863448', '163.60700392437903', '253', '244'),
            rhogids = []
            for prot_map in prot_maps:
                hog, score = prot_map[:2]
                if float(score) > threshod_f_score_merging:
                    rhogid = hog.split(".")[0].split(":")[1]
                    rhogids.append(rhogid)
            for ii in range(len(rhogids)):
                for jj in range(ii + 1, len(rhogids)):
                    hogi, hogj = rhogids[ii], rhogids[jj]
                    if (hogi, hogj) in pair_rhogs_count:
                        pair_rhogs_count[(hogi, hogj)] += 1
                    else:
                        pair_rhogs_count[(hogi, hogj)] = 1

    print(len(pair_rhogs_count))
    logger_hog.debug("There are " + str(len(pair_rhogs_count)) + " pairs of rhogs.")

    rhogs_size = {}
    for rhog, list_prot in rhogs_prots.items():
        rhogs_size[rhog] = len(list_prot)

    candidates_pair = []
    # dic_pair_ratio = {}
    for (hogi, hogj), count_shared in pair_rhogs_count.items():
        if hogi in rhogs_size and hogj in rhogs_size:  # during previous functions, we might
            ratioMax = count_shared / max(rhogs_size[hogi], rhogs_size[hogj])
            ratioMin = count_shared / min(rhogs_size[hogi], rhogs_size[hogj])

            if (ratioMax > _config.mergHOG_ratioMax_thresh or ratioMin >  _config.mergHOG_ratioMin_thresh ) and _config.mergHOG_shared_thresh > 20 and\
                    rhogs_size[hogi] < _config.big_rhog_size / 2 and  rhogs_size[hogj] < _config.big_rhog_size / 2:
                if rhogs_size[hogi] >= rhogs_size[hogj]:
                    candidates_pair.append((hogi, hogj))  # bigger first
                else:
                    candidates_pair.append((hogj, hogi))
                # print(hogi,"(",rhogs_size[hogi],")",hogj,"(",rhogs_size[hogj],")", count_shared,round(ratioMax,2),round(ratioMin,2))
            # D0651051 ( 59 ) D0658569 ( 180 ) 62 0.34 1.05 # todo why bigger than 1?
            # D0646495 ( 2 ) D0631227 ( 14 ) 25 1.79 12.5

    logger_hog.debug("There are " + str(len(candidates_pair)) + " candidate pairs of rhogs for merging.")
    print(len(candidates_pair))

    return candidates_pair


def cluster_rhogs(candidates_pair):
    # init
    all_hog_raw = []
    for pair in candidates_pair:
        all_hog_raw.append(pair[0])
        all_hog_raw.append(pair[1])
    all_hog = list(set(all_hog_raw))

    allcc = []  # connected compoenets
    for hog in all_hog:
        allcc.append([hog])

    dic_where = {}
    for idx, cc in enumerate(allcc):
        dic_where[cc[0]] = idx  # in the beginning, each inner list has only on element, a unique hog

    # print(len(all_hog_raw),len(all_hog),len(allcc),allcc[:2])
    logger_hog.debug("There are " + str(len(all_hog)) + " all_hog.")

    # print(dic_where)
    for pair in candidates_pair:
        # print(pair)
        idx_0 = dic_where[pair[0]]
        idx_1 = dic_where[pair[1]]

        dic_where[pair[1]] = idx_0

        # print(idx_0)
        allcc[idx_0] += allcc[idx_1]
        allcc[idx_1] = []

        # print(dic_where,allcc)

    cluster_rhogs_list = [i for i in allcc if len(i) > 1]
    # print(cluster_rhogs_list)
    logger_hog.debug("There are " + str(len(cluster_rhogs_list)) + " cluster_rhogs.")

    # allcc_cleaned_len = [len(i) for i in cluster_rhogs]
    # len(allcc),len(cluster_rhogs),np.sum(allcc_cleaned_len)

    return cluster_rhogs_list


def merge_rhogs(hogmaps, rhogs_prots):
    logger_hog.debug("started merging  ")

    candidates_pair = find_rhog_candidate_pairs(hogmaps, rhogs_prots)

    cluster_rhogs_list = cluster_rhogs(candidates_pair)

    logger_hog.debug("There are " + str(len(rhogs_prots)) + " rhogs before merging.")
    print(len(rhogs_prots))

    file_out_merge =  open("merging_rhogs.txt","w")
    file_out_merge.write("#first column is the host hog, the rest will be merged here.\n")
    for cluster in cluster_rhogs_list:
        for cluster_i in cluster:
            file_out_merge.write(cluster_i+"\t")
        file_out_merge.write("\n")
    file_out_merge.close()
    counter_merged_prots=0
    for cluster in cluster_rhogs_list:
        #prots = [rhogs_prots[hog] for hog in cluster]
        host_hog = cluster[0]
        rhogs_host_size= len(rhogs_prots[host_hog])
        all_prots = []
        for hog in cluster:
            all_prots += rhogs_prots[hog]
            del rhogs_prots[hog]
        all_prots_uniq = list(set(all_prots))  # there might be repeated prots
        rhogs_prots[host_hog] = all_prots_uniq
        # merging D0562038 and D0559070
        # tr | C3ZG56 | C3ZG56_BRAFL HOG: D0562038, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0562038, tr | C3ZG56 | C3ZG56_BRAFL HOG: D0559070, tr | H2Y1V7 | H2Y1V7_CIOIN HOG: D0559070
        if len(all_prots_uniq) > rhogs_host_size:
            counter_merged_prots += len(all_prots_uniq)
        # otherwise, merging didn't help

    print(len(rhogs_prots),counter_merged_prots)
    logger_hog.debug("There are " + str(len(rhogs_prots)) + " rhogs by merging "+str(counter_merged_prots)+" proteins in total.")

    return rhogs_prots


def collect_unmapped_singleton(rhogs_prots, unmapped,prot_recs_all,unmapped_singleton_fasta= "singleton_unmapped.fa"):
    unmapped_recs = []
    for species_name, prot_names in unmapped.items():
        for prot_name in prot_names:
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


def run_linclust(fasta_to_cluster="singleton_unmapped.fa"):
    num_threads = 10 # todo how to assign more cpu for this step
    command_clust = _config.mmseqs_executable_path +" easy-linclust --threads " + str(
        num_threads) + " " + fasta_to_cluster + " singleton_unmapped tmp_linclust"

    logger_hog.debug("linclust rooting started" + command_clust)
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
    cluster_dic = {}
    previous_line_start = " "
    # parsing memseqs fasta-like format https://github.com/soedinglab/MMseqs2/wiki#cluster-fasta-like-format
    # >A0A0R4IFW6
    # >tr|A0A0R4IFW6|A0A0R4IFW6_DANRE||DANRE||1000001584 tr|A0A0R
    # MSRKTTSKRHYKPSSEIDDAALARKREYW
    clusters = []
    cluster = []
    line_strip = " "
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

    logger_hog.debug("Number of linclust clusters raw is " + str(len(clusters)))
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
                        hogids = [i[0] for i in prot_maps]
                        seq_lens = [i[2] for i in prot_maps]
                        subf_med_lens = [i[3] for i in prot_maps]
                        family_score = max(scores)
                        idx_max = scores.index(family_score)
                        seq_len = seq_lens[idx_max]
                        subf_med_len = subf_med_lens[idx_max]
                        # hogid = hogids[max_index]
                    else:
                        # hogid = prot_maps[0][0]
                        family_score = float(prot_maps[0][1])
                        seq_len = prot_maps[0][2]
                        subf_med_len = prot_maps[0][3]
                    score_isof = float(family_score) * min(int(seq_len), int(subf_med_len))
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
