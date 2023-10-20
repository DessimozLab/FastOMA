
from Bio import SeqIO
import pickle
from os import listdir
import os

from ._utils_subhog import logger_hog
from . import _config



def parse_proteomes(folder="./"):  # list_oma_species
    """
    parsing fasta files of proteins located in /proteome/
    using Bio.SeqIO.parse
    Each fasta file is for one species.  The file name is the species name.
    All proteoems should be with the same extension.
    output: prot_recs_lists: a dic with key as species name and  its value list of Biopython record of species.
    """

    project_files = listdir(folder+"/proteome/")
    species_names = [] # query/input species name based on the file name
    for file in project_files:
        fasta_format = file.split(".")[-1]
        if fasta_format == "fa" or fasta_format == "fasta":
            file_name_split = file.split(".")[:-1]
            species_names.append('.'.join(file_name_split))
            fasta_format_keep = fasta_format # last one is stored either fa or fasta
    prot_recs_lists = {} # key: species name, value is a dic of query protein Biopython records.
    # 'MYCGE': [SeqRecord(seq=Seq('MDFDK

    for species_name in species_names:
        prot_address = folder+"/proteome/" + species_name + "."+fasta_format_keep
        prots_record = list(SeqIO.parse(prot_address, "fasta"))
        prot_recs_lists[species_name]=prots_record

    logger_hog.info("The are "+str(len(species_names))+" species in the proteome folder.")
    return species_names, prot_recs_lists



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


def parse_hogmap_omamer(species_names, folder="./"):
    """
     function for parsing output of omamer (hogmap files) located in /hogmap/
    Each hogmap file correspond to one fasta file of species, with the same name.
    Note that some records of fasta may removed in hogmap, probably becuase of being so short.
    hogmap file example:
    # qseqid	hogid	hoglevel	family_p	family_count	family_normcount	subfamily_score	subfamily_count	qseqlen	subfamily_medianseqlen	qseq_overlap
    # sp|O66907|ATPA_AQUAE	HOG:C0886513.1b	Eukaryota	696.5519485850672	187	0.37302594463771466	0.2643067275644273	120	504	551	0.8230616302186878
    output is dic of dic for all species:
    """
    hogmaps = {}
    for species_name in species_names:
        hogmap_address = folder + "/hogmap/" + species_name + ".fa.hogmap"
        hogmap_file = open(hogmap_address, 'r')

        for line in hogmap_file:
            line_strip = line.strip()
            if not line_strip.startswith('!') and not line_strip.startswith('qs'):
                # qseqid	hogid	hoglevel	family_p	family_count	family_normcount	subfamily_score	subfamily_count	qseqlen	subfamily_medianseqlen	qseq_overlap
                line_split = line_strip.split("\t")
                prot_id = line_split[0]
                hogid = line_split[1]
                score = line_split[3]
                seqlen = line_split[8]
                subfamily_medianseqlen = line_split[9]

                if hogid == "N/A":
                    continue

                if species_name in hogmaps:
                    if prot_id in hogmaps[species_name]:
                        hogmaps[species_name][prot_id].append((hogid, score, seqlen, subfamily_medianseqlen))

                    else:
                        hogmaps[species_name][prot_id] = [(hogid, score, seqlen, subfamily_medianseqlen)]
                else:
                    hogmaps[species_name] = {prot_id: [(hogid, score, seqlen, subfamily_medianseqlen)]}

    logger_hog.info("There are " + str(len(hogmaps)) + " species in the hogmap folder.")
    logger_hog.info("The first species " + species_names[0] + " contains " + str(len(hogmaps[species_names[0]])) + " proteins.")

    return hogmaps


def group_prots_roothogs(hogmaps):
    """
    function for finding those proteins that are mapped to the same rootHOG.

    output: rhogs_dic
    """

    rhogs_prots = {}
    for species_name, prots_map in hogmaps.items():
        # prot_recs = prot_recs_all[species_name]

        for prot_id, prot_map in prots_map.items():
            if len(prot_map)>1:
                scores = [float(i[1]) for i in prot_map] # (hogid,score,seqlen,subfamily_medianseqlen)
                rhogids =[i[0] for i in prot_map]
                # select hog with the best score
                max_index = scores.index(max(scores))
                hogid =rhogids[max_index]
            else:
                hogid = prot_map[0][0]

            rhogid = hogid.split(".")[0].split(":")[1]

            if rhogid in rhogs_prots:
                rhogs_prots[rhogid].append((species_name, prot_id))
            else:
                rhogs_prots[rhogid] = [(species_name, prot_id)]
    logger_hog.info("There are " + str(len(rhogs_prots)) + " rootHOGs.")
    # from old code, to check later: Keep prot seq in hog class. We can write species idx and prot idx to  improve speed of code for omamer tresholidng
    return rhogs_prots


def roothogs_postprocess(hogmaps, rhogs_prots):


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

        # removing protines with low omamer score from big rootHOG
        for species_name, prot_id in sp_prot_list:
            prot_maps = hogmaps[species_name][prot_id]
            if len(prot_maps) > 1:  # probably for big rootHOG, there won't be multi-hits
                scores = [float(i[1]) for i in prot_maps]  # (hogid,score,seqlen,subfamily_medianseqlen)
                hogids = [i[0] for i in prot_maps]
                max_score = max(scores)
                max_index = scores.index(max_score)
                hogid = hogids[max_index]
            else:
                hogid = prot_maps[0][0]
                score = float(prot_maps[0][1])

            if score > _config.omamer_family_threshold:
                sp_prot_list_filt.append((species_name, prot_id))
                hogids.append(hogid)

        logger_hog.info("For big rootHOG " + rhogid + ", " + str(
            len(sp_prot_list_filt)) + " proteins left after filtering with threshold " + str(_config.omamer_family_threshold))

        if len(sp_prot_list_filt) < _config.big_rhog_size:
            sp_prot_list_filt2 = sp_prot_list_filt


        else:  # removing protines that are mapped to rootHOG (=HOGC123 not the levelsHOGC123.1a) from big rootHOG
            hogids2 = []
            sp_prot_list_filt2 = []
            for prot_idx, sp_prot in enumerate(sp_prot_list_filt):
                hogid = hogids[prot_idx]

                if len(hogid.split(".")) > 1:
                    sp_prot_list_filt2.append(sp_prot)
                    hogids2.append(hogid)

        logger_hog.info("For big rootHOG " + rhogid + ", " + str(
            len(sp_prot_list_filt2)) + " proteins left after removing non-child subhogs")

        rhogs_prots[rhogid] = sp_prot_list_filt2

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




# import pyoma.browser.db as db
# def parse_oma_db(oma_database_address):
#     """
#    This is in  new branch
#     orthoxml_to_newick.py function for loading an oma database in hdf5 format using pyoma.browser.db.
#     output: oma_db, list_oma_species
#     """
#     oma_db = db.Database(oma_database_address)
#     logger_hog.info("OMA data is parsed and its release name is:" + oma_db.get_release_name())
#
#     list_oma_species = [z.uniprot_species_code for z in oma_db.tax.genomes.values()]
#     logger_hog.info("There are "+str(len(list_oma_species))+" species in the OMA database.")
#     return oma_db, list_oma_species
#
#
# def parse_proteome():  # list_oma_species
#     """
#     orthoxml_to_newick.py function for parsing fasta files of proteins located in /proteome/
#     using Bio.SeqIO.parse
#     Each fasta file is for one species.  The file name is the species name.
#     output: query_species_names: list of species name, query_prot_recs: list of Biopython record of species
#     """
#     # project_files = listdir(_config.in_folder + "/proteome/")
#     project_files = listdir("./proteome/")
#     query_species_names = []
#     for file in project_files:
#         fasta_format = file.split(".")[-1]
#         if fasta_format == "fa" or fasta_format == "fasta":
#             file_name_split = file.split(".")[:-1]
#             query_species_names.append('.'.join(file_name_split))
#             fasta_format_keep = fasta_format
#
#     query_prot_recs = []
#     for query_species_names_idx, query_species_name in enumerate(query_species_names):
#         prot_address = "./proteome/" + query_species_name + "."+fasta_format_keep
#         prots_record = list(SeqIO.parse(prot_address, "fasta"))
#         query_prot_recs.append(prots_record)
#
#     query_species_num = len(query_species_names)
#     logger_hog.info("The are "+str(query_species_num)+" species in the proteome folder.")
#     # for development
#     # for species_i in range(query_species_num):
#     #    species_name_i = query_species_names[species_i]
#     #    if species_name_i in list_oma_species:
#     #        logger_hog.error("The species"+species_name_i+" already exists in the oma database, remove/rename it first.")
#     #        exit()
#     # The proteins are parsed using  Bio.SeqIO.parse
#     # the first part of the header line before space
#     # >tr|A0A2I3FYY2|A0A2I3FYY2_NOMLE Uncharacterized protein OS=Nomascus leucogenys OX=61853 GN=CLPTM1L PE=3 SV=1
#     # will be ">tr|A0A2I3FYY2|A0A2I3FYY2_NOMLE"
#     # [i.id for i in query_prot_recs[0] if len(i.id)!=30 and len(i.id)!=22 ] #'sp|O47892|CYB_NOMLE',
#     return query_species_names, query_prot_recs
#
#
# def add_species_name_gene_id(query_prot_recs, query_species_names):
#     """
#     adding the name of species to each protein record
#         - based on file name
#     adding gene id number, integer imposed by xml format
#     output: updated version of input
#     """
#     #  _config.in_folder +
#     gene_id_pickle_file = "./gene_id_dic_xml.pickle"
#     max_num_prot = int(1e9)
#     max_num_prot_per_sp = int(1e6)
#     gene_id_name = {}
#     for query_species_idx, query_species_name in enumerate(query_species_names):
#         query_prot_records = query_prot_recs[query_species_idx]
#         gene_counter = max_num_prot + query_species_idx * max_num_prot_per_sp
#         gene_id_name[query_species_name] = []
#         for query_prot_idx, query_prot_record in enumerate(query_prot_records):
#             gene_idx_integer = gene_counter + query_prot_idx
#             query_prot_name = query_prot_record.id
#             if len(query_prot_name) > 230:
#                 logger_hog.info("We are truncating the prot name as it may be problamatic for mafft, " + str(query_prot_name))
#                 query_prot_name = query_prot_name[:230]
#             query_prot_record.id = query_prot_name + "||"+query_species_name+"||"+str(gene_idx_integer)
#             gene_id_name[query_species_name].append((gene_idx_integer, query_prot_name))
#     # this is used to create the first part of xml file.
#     with open(gene_id_pickle_file, 'wb') as handle:
#         pickle.dump(gene_id_name, handle, protocol=pickle.HIGHEST_PROTOCOL)
#
#     return query_prot_recs
#
#
# def parse_hogmap_omamer(query_species_names):
#     """
#     orthoxml_to_newick.py function for parsing output of omamer (hogmap files) located in /hogmap/
#     Each hogmap file correspond to one fasta file of species, with the same name.
#     Note that some records of fasta may removed in hogmap, becuase of being so short.
#     hogmap file example:
#     qseqid hogid overlap family-score subfamily-score qseqlen subfamily-medianseqlen
#     A0A140TAT7_CIOIN HOG:B0833785.1c.8b 1 0.99 0.9 490 503
#     output as list of list for all species:
#     prots_hogmap_name_allspecies, prots_hogmap_hogid_allspecies,
#     prots_hogmap_subfscore_allspecies, prots_hogmap_seqlen_allspecies,
#     prots_hogmap_subfmedseqlen_allspecies
#     The order of species is the same as query_species_names.
#     """
#     prots_hogmap_name_allspecies = []
#     prots_hogmap_hogid_allspecies = []
#     prots_hogmap_overlp_allspecies = []
#     prots_hogmap_fscore_allspecies = []
#     prots_hogmap_seqlen_allspecies = []
#     prots_hogmap_subfmedseqlen_allspecies = []
#     for query_species_name in query_species_names:
#         omamer_output_address =  "./hogmap/" + query_species_name + ".fa.hogmap"
#         omamer_output_file = open(omamer_output_address, 'r')
#         prots_hogmap_name = []
#         prots_hogmap_hogid = []
#         prots_hogmap_seqlen = []
#         prots_hogmap_fscore = []
#         prots_hogmap_overlp = []
#         prots_hogmap_subfmedseqlen = []
#         for line in omamer_output_file:
#             line_strip = line.strip()
#             if not line_strip.startswith('!') and not line_strip.startswith('qs'):
#                 line_split = line_strip.split("\t")
#                 # header of omamer v2
#                 # qseqid	hogid	hoglevel	family_p	family_count	5)family_normcount	subfamily_score	subfamily_count	8)qseqlen	9)subfamily_medianseqlen	10)qseq_overlap
#                 prots_hogmap_name.append(line_split[0])
#                 if len(line_split[1].split(":"))>1:
#                     prots_hogmap_hogid.append(line_split[1].split(":")[1]) # HOG:B0012312
#                 else:
#                     prots_hogmap_hogid.append(line_split[1])  # N/A
#                 prots_hogmap_overlp.append(line_split[10]) # a value > 100 is good
#                 prots_hogmap_fscore.append(line_split[3])
#                 prots_hogmap_seqlen.append(line_split[8])
#                 prots_hogmap_subfmedseqlen.append(line_split[9])
#         prots_hogmap_name_allspecies.append(prots_hogmap_name)
#         prots_hogmap_hogid_allspecies.append(prots_hogmap_hogid)
#         prots_hogmap_overlp_allspecies.append(prots_hogmap_overlp)
#         prots_hogmap_fscore_allspecies.append(prots_hogmap_fscore)
#         prots_hogmap_seqlen_allspecies.append(prots_hogmap_seqlen)
#         prots_hogmap_subfmedseqlen_allspecies.append(prots_hogmap_subfmedseqlen)
#
#     logger_hog.info("There are "+str(len(prots_hogmap_name_allspecies))+" species in the hogmap folder.")
#     logger_hog.info("The first species "+query_species_names[0]+" contains "+str(len(prots_hogmap_hogid_allspecies[0]))+" proteins.")
#     logger_hog.info("The first protein of first species is "+prots_hogmap_name_allspecies[0][0])
#     hogmap_allspecies = (prots_hogmap_name_allspecies, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
#                          prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies)
#     return hogmap_allspecies
#
#
# def filter_prot_mapped(query_species_names, query_prot_recs, query_prot_names_species_mapped):
#     """
#     orthoxml_to_newick.py function for filtering biopython records in query_prot_recs based on hogmaps
#     The reason is that some very short records of fasta are removed in hogmap.
#     So, we may loose track of order comparing hogmap and fasta file.
#     The goal here is to remove those from seq record (of the fasta file).
#     output: query_prot_recs_filt
#     """
#     logger_hog.info("Filtering proteins started.")
#     query_prot_recs_filt = []
#     logger_hog.error("warning: we are reporting protes whose names are truncated. Becuase it is not in hogmap.")
#     for species_idx, query_species_name in enumerate(query_species_names):  # from fasta file
#         query_prot_recs_i = query_prot_recs[species_idx]
#         # we added the species name and the as the fasta record using || but this is not done in the hog map
#         # record.id = 'tr|A0A024FLK4|A0A024FLK4_ORYSJ' || ORYSJ || 1000000
#         query_prot_names_records = [record.id.split("||")[0] for record in query_prot_recs_i]
#         # from hogmap file without proteins that are not mapped on any hogs
#         query_prot_names_species_i = query_prot_names_species_mapped[species_idx]
#         if len(query_prot_names_species_i) != len(query_prot_recs_i):
#             query_prot_records_filterd_sp = []
#             for query_prot_name in query_prot_names_species_i:
#                 if query_prot_name in query_prot_names_records:
#                     prot_record_idx = query_prot_names_records.index(query_prot_name)
#                     prot_record = query_prot_recs_i[prot_record_idx]
#                     query_prot_records_filterd_sp.append(prot_record)
#                 else:
#                     logger_hog.error("Error 1349 " + query_species_name + " " + query_prot_name+". This shouldn't happen many times.")
#             logger_hog.info("For the species "+query_species_name+", some proteins were ignored by omamer, probably cause prot length .")
#             logger_hog.info("Before filtering: in hogmap "+str(len(query_prot_names_species_i))+", in proteome "+str(len(query_prot_recs_i)))
#             logger_hog.info("After filtering: in hogmap "+str(len(query_prot_names_species_i))+" in proteome "+str(len(query_prot_records_filterd_sp)))
#         else:
#             query_prot_records_filterd_sp = query_prot_recs_i
#         query_prot_recs_filt.append(query_prot_records_filterd_sp)
#     logger_hog.info("For the rest of species, all proteins were mapped using OMAmer.")
#     return query_prot_recs_filt
#
#
# def group_prots_roothogs(prots_hogmap_hogid_allspecies, query_species_names, query_prot_recs_filt):
#     """
#     orthoxml_to_newick.py function for finding those proteins that are mapped to the same rootHOG.
#     Then, we write each rootHOG as orthoxml_to_newick.py seprate fasta file in the address_rhogs_folder folder
#     output: rhogid_list, rhogids_prot_records_query
#     """
#     # extract rootHOG ID  "B0833755.5c.10g.24e.16c.18b" ->"B0833755"
#     prots_hogmap_rhogid_allspecies = []
#     for prots_hogmap_hogid in prots_hogmap_hogid_allspecies:
#         prots_hogmap_rhogid = []
#         for prot_hogmap_hogid in prots_hogmap_hogid:
#             prot_hogmap_rhogid = prot_hogmap_hogid.split(".")[0]
#             prots_hogmap_rhogid.append(prot_hogmap_rhogid)
#         prots_hogmap_rhogid_allspecies.append(prots_hogmap_rhogid)
#     # gathering name of prots from all species,  group them based on rHOG that they mapped on
#     rhogid_prot_idx_dic = {}
#     for species_idx in range(len(query_species_names)):
#         prots_hogmap_rhogid = prots_hogmap_rhogid_allspecies[species_idx]
#         for prots_hogmap_idx, prot_hogmap_rhogid in enumerate(prots_hogmap_rhogid):
#             if prot_hogmap_rhogid in rhogid_prot_idx_dic:
#                 rhogid_prot_idx_dic[prot_hogmap_rhogid].append((species_idx, prots_hogmap_idx))
#             else:
#                 rhogid_prot_idx_dic[prot_hogmap_rhogid] = [(species_idx, prots_hogmap_idx)]
#     # extracting prot records for each rootHOG
#     rhogids_prot_records_query = []
#     rhogids_list = []
#     for rhogid in rhogid_prot_idx_dic.keys():
#         rhogid_prot_records = []
#         species_idx_rhogid = []
#         if rhogid != "N/A" and len(rhogid) >= 1:  # ignore un-mapped prots
#             rhogids_list.append(rhogid)
#             rhogid_prot_idx = rhogid_prot_idx_dic[rhogid]
#             for (species_idx, prots_hogmap_idx) in rhogid_prot_idx:
#                 prot_record = query_prot_recs_filt[species_idx][prots_hogmap_idx]
#                 """
#                 Keep prot seq in hog class. We can write species idx and prot idx to  improve speed of code for omamer tresholidng
#                 """
#                 rhogid_prot_records.append(prot_record)
#                 species_idx_rhogid.append(species_idx)
#             rhogids_prot_records_query.append(rhogid_prot_records)
#     logger_hog.info("There are " + str(len(rhogids_list)) + " rhogs, no matter their size.")
#     return rhogids_list, rhogids_prot_records_query
#
#
# def filter_rhog(rhogids_list, rhogids_prot_records_query, prots_hogmap_fscore_allspecies, query_species_names,  prots_hogmap_name_allspecies):
#     """
#     Some of the rhogs are very big. We filter those rhogs where many proteins several tousands are mapped on.
#     The treshold is set in the _config.py file.
#     """
#     logger_hog.info("Filtering rhogs with fscore treshold "+str(_config.omamer_fscore_treshold_big_rhog)+"for rhogs size > "+str(_config.omamer_treshold_big_rhog_szie) )
#
#     rhogids_prot_records_query_filt2 = []
#     rhogids_list_filt = []
#     for rhogid_idx, rhogid in enumerate(rhogids_list):
#         rhogid_prot_record_query = rhogids_prot_records_query[rhogid_idx]
#         if len(rhogid_prot_record_query) < _config.omamer_treshold_big_rhog_szie:
#             rhogid_prot_record_query_filt2 = rhogid_prot_record_query  # without change for small rhogs
#         else:
#             rhogid_prot_record_query_filt = []
#             for i in range(len(rhogid_prot_record_query)):
#                 prot_bio_seq = rhogid_prot_record_query[i]
#                 prot_name_trunc, species_name, prot_idx_xml = prot_bio_seq.id.split("||")
#                 prot_name = prot_bio_seq.name    # .id is the truncated one but .name is full
#                 specis_idx = query_species_names.index(species_name)
#                 prot_list = prots_hogmap_name_allspecies[specis_idx]
#                 prot_idx = prot_list.index(prot_name)
#                 fsore = float(prots_hogmap_fscore_allspecies[specis_idx][prot_idx])
#                 if fsore > _config.omamer_fscore_treshold_big_rhog:
#                     rhogid_prot_record_query_filt.append(prot_bio_seq)
#                 else:
#                     logger_hog.info("we are removing due to filtering rhogs with fscore treshold1 " + str(prot_name))
#
#             if len(rhogid_prot_record_query_filt) < _config.omamer_treshold_big_rhog_szie2: # 40 * 1000
#                     rhogid_prot_record_query_filt2 = rhogid_prot_record_query_filt  # without change for small rhogs
#             else:
#                 logger_hog.info("Second round of filtering rhogs with fscore treshold " + str(_config.omamer_fscore_treshold_big_rhog2) + "for rhogs size > " + str(_config.omamer_treshold_big_rhog_szie2))
#                 rhogid_prot_record_query_filt2 = []
#                 for i in range(len(rhogid_prot_record_query_filt)):
#                     prot_bio_seq = rhogid_prot_record_query_filt[i]
#                     prot_name_trunc, species_name, prot_idx_xml = prot_bio_seq.id.split("||")
#                     prot_name = prot_bio_seq.name  # .id is the truncated one but .name is full
#                     specis_idx = query_species_names.index(species_name)
#                     prot_list = prots_hogmap_name_allspecies[specis_idx]
#                     prot_idx = prot_list.index(prot_name)
#                     fsore = float(prots_hogmap_fscore_allspecies[specis_idx][prot_idx])
#                     if fsore > _config.omamer_fscore_treshold_big_rhog2:
#                         rhogid_prot_record_query_filt2.append(prot_bio_seq)
#                     else:
#                         logger_hog.info("we are removing due to second round of filtering rhogs with fscore treshold2 " + str(prot_name))
#
#         if rhogid_prot_record_query_filt2:  # at least one prot in the rhog
#             rhogids_prot_records_query_filt2.append(rhogid_prot_record_query_filt2)
#             rhogids_list_filt.append(rhogid)
#     return rhogids_list_filt, rhogids_prot_records_query_filt2
#
#

def write_rhog_old(rhogids_list, rhogids_prot_records_query, address_rhogs_folder, min_rhog_size=1, max_rhog_size=1e100,merge_split=0):

    logger_hog.info("Writing Sequences of roothogs are fasta file in " + address_rhogs_folder)
    if not os.path.exists(address_rhogs_folder):
        os.mkdir(address_rhogs_folder)


    ## under development
    #merge_split=0
    rhogid_list = []
    rhogs_merged = []
    if merge_split ==1:
        import pickle
        try:
            with open('/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/allcc_cleaned.pickle', 'rb') as handle:
                allcc_cleaned2 = pickle.load(handle)
        except:
            print("no merge split roothog info.")

        for cc in allcc_cleaned2:
            rhogs_merged += cc

        for cc in allcc_cleaned2:
            rhogid_prot_rec_query_merged =[ ]
            for rhogid in cc:
                if rhogid in rhogids_list:
                    rhogid_idx = rhogids_list.index(rhogid)
                    rhogids_prot_record_query = rhogids_prot_records_query[rhogid_idx]
                    rhogid_prot_rec_query_merged += rhogids_prot_record_query
                else:
                    print("not in the list",rhogid)
            if min_rhog_size <= len(rhogid_prot_rec_query_merged) <= max_rhog_size:
                SeqIO.write(rhogid_prot_rec_query_merged, address_rhogs_folder + "/HOG_" + rhogid + ".fa", "fasta") # name of lastrhogid
                rhogid_list.append(rhogid)

    for rhogid_idx, rhogid in enumerate(rhogids_list):
        if rhogid not in rhogs_merged:
            rhogid_prot_rec_query = rhogids_prot_records_query[rhogid_idx]

            if min_rhog_size <= len(rhogid_prot_rec_query) <= max_rhog_size:
                # todo add the release id   to file names  rhogids_list[:2] > ['HOG:C0884658', 'HOG:C0709155']
                SeqIO.write(rhogid_prot_rec_query, address_rhogs_folder + "/HOG_" + rhogid + ".fa", "fasta")
                rhogid_list.append(rhogid)
            else:
                for prot1 in rhogid_prot_rec_query:
                    logger_hog.debug("we are removing due to omamer signleton hog |*|" + str(prot1.id))


    logger_hog.info("Writing Sequences of roothogs finished." )

    return rhogid_list


import os.path


def parse_isoform_file(query_species_names):
    isoform_by_gene_all = []
    for species_idx, query_species_name in enumerate(query_species_names):  # from fasta file

        file_splice_name = "./splice/" + query_species_name + ".splice"

        if os.path.isfile(file_splice_name):
            # from OMArk
            isoform_by_gene = list()
            with open(file_splice_name) as handle:
                for line in handle.readlines():
                    line = line.strip('\n')
                    splice = line.split(";")
                    isoform_by_gene.append(splice)
        else:
            isoform_by_gene = []

        isoform_by_gene_all.append(isoform_by_gene)
    return isoform_by_gene_all


def find_nonbest_isoform(hogmap_allspecies_elements, isoform_by_gene_all):
    (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
     prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies,
     prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

    not_selected_isofroms_all = []
    for species_idx in range(len(query_prot_names_species_mapped)):

        not_selected_isofroms = []  # when there is no splice file for a species
        query_prot_names = query_prot_names_species_mapped[species_idx]
        prots_hogmap_hogid = prots_hogmap_hogid_allspecies[species_idx]
        prots_hogmap_fscore = prots_hogmap_fscore_allspecies[species_idx]
        prots_hogmap_seqlen = prots_hogmap_seqlen_allspecies[species_idx]
        prots_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen_allspecies[species_idx]
        not_selected_isofroms_list = []
        isoform_list_sp = isoform_by_gene_all[species_idx]
        for isoform_list in isoform_list_sp:
            score_best = 0
            protname_best = isoform_list[0]  # to keep the order if all of the isoforms' omamer score are 'na'
            for isoform_name in isoform_list:
                if isoform_name in query_prot_names:
                    isof_idx = query_prot_names.index(isoform_name)
                    family_score = prots_hogmap_fscore[isof_idx]
                    seq_len = prots_hogmap_seqlen[isof_idx]
                    subf_med_len = prots_hogmap_subfmedseqlen[isof_idx]
                    if family_score != "N/A":
                        score = float(family_score) * min(int(seq_len), int(subf_med_len))
                        # print(isoform_name,score, family_score,seq_len,subf_med_len)
                        if score > score_best:  # when there is a tie, the last one is selected!
                            protname_best = isoform_name
                            score_best = score
            # selected_isofroms.append([protname_best)
            # print("*", protname_best, score_best)
            not_selected_isofroms_list += [i for i in isoform_list if i != protname_best]
        not_selected_isofroms_all.append(not_selected_isofroms_list)

    return not_selected_isofroms_all # those isform that are not selected as best (which will be removed )


# def select_best_isoform(hogmap_allspecies_elements, isoform_by_gene_all):
#     (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
#      prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies,
#      prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
#
#     selected_isofroms_all = []
#     for species_idx in range(len(query_prot_names_species_mapped)):
#
#         selected_isofroms = []
#         query_prot_names = query_prot_names_species_mapped[species_idx]
#         # prots_hogmap_hogid = prots_hogmap_hogid_allspecies[species_idx]
#         prots_hogmap_fscore = prots_hogmap_fscore_allspecies[species_idx]
#         prots_hogmap_seqlen = prots_hogmap_seqlen_allspecies[species_idx]
#         prots_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen_allspecies[species_idx]
#
#         isoform_list_sp = isoform_by_gene_all[species_idx]
#         for isoform_list in isoform_list_sp:
#             score_best = 0
#             protname_best = isoform_list[0]  # to keep the order if all of the isoforms' omamer score are 'na'
#             for isoform_name in isoform_list:
#                 if isoform_name in query_prot_names:
#                     isof_idx = query_prot_names.index(isoform_name)
#                     family_score = prots_hogmap_fscore[isof_idx]
#                     seq_len = prots_hogmap_seqlen[isof_idx]
#                     subf_med_len = prots_hogmap_subfmedseqlen[isof_idx]
#                     if family_score != "na":
#                         score = float(family_score) * min(int(seq_len), int(subf_med_len))
#                         # print(isoform_name,score, family_score,seq_len,subf_med_len)
#                         if score > score_best:  # when there is a tie, the last one is selected!
#                             protname_best = isoform_name
#                             score_best = score
#             selected_isofroms.append(protname_best)
#             # print("*", protname_best, score_best)
#
#         selected_isofroms_all.append(selected_isofroms)
#
#     return selected_isofroms_all
#


def handle_splice(prots_hogmap_hogid_allspecies, query_prot_recs_filt, not_selected_isofroms_all,
                  query_prot_names_species_mapped):

    query_prot_recs_filt_ = []
    prots_hogmap_hogid_allspecies_ = []
    for species_idx, query_prot_rec_filt in enumerate(query_prot_recs_filt):
        not_selected_isofroms = not_selected_isofroms_all[species_idx]
        query_prot_recs_sp = []
        prots_hogmap_hogids_ = []
        for query_prot in query_prot_rec_filt:
            if query_prot.id not in not_selected_isofroms:
                query_prot_recs_sp.append(query_prot)
        query_prot_recs_filt_.append(query_prot_recs_sp)

        prots_hogmap_hogids = prots_hogmap_hogid_allspecies[species_idx]
        query_prot_names = query_prot_names_species_mapped[species_idx]

        for prot_idx, prots_hogmap_hogid in enumerate(prots_hogmap_hogids):
            query_prot_name = query_prot_names[prot_idx]
            if query_prot_name not in not_selected_isofroms:
                prots_hogmap_hogids_.append(prots_hogmap_hogid)

        prots_hogmap_hogid_allspecies_.append(prots_hogmap_hogids_)

    return prots_hogmap_hogid_allspecies_, query_prot_recs_filt_


