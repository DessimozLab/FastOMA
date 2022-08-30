import pyoma.browser.db as db
from datetime import datetime
from Bio import SeqIO
# import dill as dill_pickle
import pickle
from os import listdir
import os


from _utils import logger_hog



def parse_oma_db(oma_database_address):
    """
    a function for loading an oma database in hdf5 format using pyoma.browser.db.
    output: oma_db, list_oma_species
    """
    oma_db = db.Database(oma_database_address)
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- OMA data is parsed and its release name is:", oma_db.get_release_name())
    list_oma_species = [z.uniprot_species_code for z in oma_db.tax.genomes.values()]
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- There are", len(list_oma_species), "species in the OMA database.")
    return oma_db, list_oma_species


def parse_proteome(list_oma_species, working_folder):
    """
    a function for parsing fasta files of proteins located in /omamer_search/proteome/
    using Bio.SeqIO.parse
    Each fasta file is for one species.  The file name is the species name.
    output: query_species_names: list of species name, query_prot_recs: list of Biopython record of species
    """
    project_files = listdir(working_folder + "/omamer_search/proteome/")
    query_species_names = []
    for file in project_files:
        if file.split(".")[-1] == "fa" or file.split(".")[-1] == "fasta":
            file_name_split = file.split(".")[:-1]
            query_species_names.append('.'.join(file_name_split))

    # we may assert existence of query_species_name+".fa/hogmap"
    query_prot_recs = []
    # for query_species_names_idx in range(len(query_species_names)):
    #     query_species_name = query_species_names[query_species_names_idx]
    for query_species_names_idx, query_species_name in enumerate(query_species_names):
        prot_address = working_folder + "omamer_search/proteome/" + query_species_name + ".fa"
        prots_record = list(SeqIO.parse(prot_address, "fasta"))
        query_prot_recs.append(prots_record)

    query_species_num = len(query_species_names)
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- The are", str(query_species_num), "species in the proteome folder.")
    # for development
    for species_i in range(query_species_num):
        # len_prot_record_i = len(query_prot_recs[species_i])
        species_name_i = query_species_names[species_i]
        # print(species_name_i,len_prot_record_i)
        if species_name_i in list_oma_species:
            current_time = datetime.now().strftime("%H:%M:%S")
            print(current_time, "- the species", species_name_i,
                  " already exists in the oma database, remove them first")
            exit()
    # The proteins are parsed using  Bio.SeqIO.parse
    # the first part of the header line before space
    # >tr|A0A2I3FYY2|A0A2I3FYY2_NOMLE Uncharacterized protein OS=Nomascus leucogenys OX=61853 GN=CLPTM1L PE=3 SV=1
    # will be ">tr|A0A2I3FYY2|A0A2I3FYY2_NOMLE"
    # [i.id for i in query_prot_recs[0] if len(i.id)!=30 and len(i.id)!=22 ] #'sp|O47892|CYB_NOMLE',
    return query_species_names, query_prot_recs


def add_species_name_gene_id(query_prot_recs, query_species_names, gene_id_pickle_file):
    """
    adding the name of species to each protein record
        - based on file name
    adding gene id number, integer imposed by xml format
    output: updated version of input
    """
    max_num_prot = int(1e9)
    max_num_prot_per_sp = int(1e6)

    gene_id_name = {}
    for query_species_idx, query_species_name in enumerate(query_species_names):
        query_prot_records = query_prot_recs[query_species_idx]
        gene_counter = max_num_prot + query_species_idx * max_num_prot_per_sp
        gene_id_name[query_species_name] = []
        for query_prot_idx, query_prot_record in enumerate(query_prot_records):
            gene_idx_integer = gene_counter + query_prot_idx
            query_prot_name = query_prot_record.id
            query_prot_record.id += "||"+query_species_name+"||"+str(gene_idx_integer)
            gene_id_name[query_species_name].append((gene_idx_integer, query_prot_name))
    # this is used to creat the first part of xml file, gene name and
    with open(gene_id_pickle_file, 'wb') as handle:
        # dill_pickle.dump(gene_id_name, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(gene_id_name, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return query_prot_recs



def parse_hogmap_omamer(query_species_names, working_folder):
    """
    a function for parsing output of omamer (hogmap files) located in /omamer_search/hogmap/
    Each hogmap file correspond to one fasta file of species, with the same name.
    Note that some records of fasta may removed in hogmap, becuase of being so short.
    hogmap file example:
    qseqid hogid overlap family-score subfamily-score qseqlen subfamily-medianseqlen
    A0A140TAT7_CIOIN HOG:B0833785.1c.8b 1 0.99 0.9 490 503
    output as list of list for all species:
    prots_hogmap_name_allspecies, prots_hogmap_hogid_allspecies,
    prots_hogmap_subfscore_allspecies, prots_hogmap_seqlen_allspecies,
    prots_hogmap_subfmedseqlen_allspecies
    The order of species is the same as query_species_names.
    """
    prots_hogmap_name_allspecies = []
    prots_hogmap_hogid_allspecies = []
    prots_hogmap_subfscore_allspecies = []
    prots_hogmap_seqlen_allspecies = []
    prots_hogmap_subfmedseqlen_allspecies = []
    for query_species_name in query_species_names:
        omamer_output_address = working_folder + "omamer_search/hogmap/" + query_species_name + ".hogmap"
        omamer_output_file = open(omamer_output_address, 'r')
        prots_hogmap_name = []
        prots_hogmap_hogid = []
        prots_hogmap_subfscore = []
        prots_hogmap_seqlen = []
        prots_hogmap_subfmedseqlen = []
        for line in omamer_output_file:
            line_strip = line.strip()
            if not line_strip.startswith('qs'):
                line_split = line_strip.split("\t")
                # if line_split[1]!='na':
                prots_hogmap_name.append(line_split[0])
                prots_hogmap_hogid.append(line_split[1])
                prots_hogmap_subfscore.append(line_split[4])  # subfamily
                prots_hogmap_seqlen.append(line_split[5])
                prots_hogmap_subfmedseqlen.append(line_split[6])
        prots_hogmap_name_allspecies.append(prots_hogmap_name)
        prots_hogmap_hogid_allspecies.append(prots_hogmap_hogid)
        prots_hogmap_subfscore_allspecies.append(prots_hogmap_subfscore)
        prots_hogmap_seqlen_allspecies.append(prots_hogmap_seqlen)
        prots_hogmap_subfmedseqlen_allspecies.append(prots_hogmap_subfmedseqlen)

    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- There are ", len(prots_hogmap_name_allspecies), " species in the hogmap folder.")
    print(current_time, "- The first species", query_species_names[0], " contains ",
          len(prots_hogmap_hogid_allspecies[0]), " proteins.")
    print(current_time, "- The first protein of first species is ", prots_hogmap_name_allspecies[0][0])
    hogmap_allspecies = (prots_hogmap_name_allspecies, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
                         prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies)
    return hogmap_allspecies


def filter_prot_mapped(query_species_names, query_prot_recs, query_prot_names_species_mapped):
    """
    a function for filtering biopython records in query_prot_recs based on hogmaps
    The reason is that some very short records of fasta are removed in hogmap.
    So, we may loose track of order comparing hogmap and fasta file.
    The goal here is to remove those from seq record (of the fasta file).
    output: query_prot_recs_filt
    """
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- Filtering proteins started.")
    query_prot_recs_filt = []
    # for species_idx in range(len(query_species_names)):
    #     query_species_name = query_species_names[species_idx]
    for species_idx, query_species_name in enumerate(query_species_names): # from fasta file
        # print(query_species_name)
        query_prot_recs_i = query_prot_recs[species_idx]
        # we added the species name and the as the fasta record usign ||
        # but this is not done in the hog map
        # I beleive .split("||")[0] solve this.
        query_prot_ids_records = [record.id.split("||")[0] for record in query_prot_recs_i]
        # from hogmap file without proteins that are not mapped on any hogs
        query_prot_names_species_i = query_prot_names_species_mapped[species_idx]
        if len(query_prot_names_species_i) != len(query_prot_recs_i):
            query_prot_records_filterd_sp = []
            for query_prot_name in query_prot_names_species_i:
                if query_prot_name in query_prot_ids_records:
                    prot_record_idx = query_prot_ids_records.index(query_prot_name)
                    prot_record = query_prot_recs_i[prot_record_idx]
                    query_prot_records_filterd_sp.append(prot_record)
                else:
                    current_time = datetime.now().strftime("%H:%M:%S")
                    logger_hog.error(str(current_time) + "- Error 149 " + query_species_name + " " + query_prot_name+". This shouldn't happen many times.")

            current_time = datetime.now().strftime("%H:%M:%S")
            print(current_time, "- For the species", query_species_name, ", few proteins were ignored by omamer.")
            print(current_time, "- before filtering: in hogmap", len(query_prot_names_species_i), "in proteome",
                  len(query_prot_recs_i))
            print(current_time, "- After filtering:  in hogmap", len(query_prot_names_species_i), "in proteome",
                  len(query_prot_records_filterd_sp))
        else:
            query_prot_records_filterd_sp = query_prot_recs_i
        query_prot_recs_filt.append(query_prot_records_filterd_sp)
    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- For the rest of species, all proteins were mapped using OMAmer.")
    return query_prot_recs_filt

def group_prots_roothogs(prots_hogmap_hogid_allspecies,  address_rhogs_folder, query_species_names, query_prot_recs_filt):
    """
    a function for finding those proteins that are mapped to the same rootHOG.
    Then, we write each rootHOG as a seprate fasta file in the address_rhogs_folder folder
    output: rhogid_num_list, rhogids_prot_records_query
    """

    if not os.path.exists(address_rhogs_folder):
        os.mkdir(address_rhogs_folder)

    print_hog_more_than_one_species = False  # if false keep those rootHOG containing only one species

    # extract rootHOG ID  "B0833755.5c.10g.24e.16c.18b" ->"B0833755"
    prots_hogmap_rhogid_allspecies = []
    for prots_hogmap_hogid in prots_hogmap_hogid_allspecies:
        prots_hogmap_rhogid = []
        for prot_hogmap_hogid in prots_hogmap_hogid:
            prot_hogmap_rhogid = prot_hogmap_hogid.split(".")[0]
            prots_hogmap_rhogid.append(prot_hogmap_rhogid)
        prots_hogmap_rhogid_allspecies.append(prots_hogmap_rhogid)

    # gathering name of prots from all species,  group them based on rHOG that they mapped on
    rhogid_prot_idx_dic = {}
    for species_idx in range(len(query_species_names)):
        # species_name = query_species_names[species_idx]
        prots_hogmap_rhogid = prots_hogmap_rhogid_allspecies[species_idx]
        # for prots_hogmap_idx in range(len(prots_hogmap_rhogid)):
        #     prot_hogmap_rhogid = prots_hogmap_rhogid[prots_hogmap_idx]
        for prots_hogmap_idx, prot_hogmap_rhogid in enumerate(prots_hogmap_rhogid):
            if prot_hogmap_rhogid in rhogid_prot_idx_dic:
                rhogid_prot_idx_dic[prot_hogmap_rhogid].append((species_idx, prots_hogmap_idx))
            else:
                rhogid_prot_idx_dic[prot_hogmap_rhogid] = [(species_idx, prots_hogmap_idx)]
    # print(len(rhogid_prot_idx_dic)) #  rhogid_prot_idx_dic['HOG:0018405']

    # extracting prot records for each rootHOG
    rhogids_prot_records_query = []
    rhogids_list = []
    for rhogid in rhogid_prot_idx_dic.keys():
        rhogid_prot_records = []
        species_idx_rhogid = []
        if rhogid != "na" and len(rhogid) > 1:  # ignore un-mapped prots
            rhogids_list.append(rhogid)
            rhogid_prot_idx = rhogid_prot_idx_dic[rhogid]
            for (species_idx, prots_hogmap_idx) in rhogid_prot_idx:
                prot_record = query_prot_recs_filt[species_idx][prots_hogmap_idx]
                rhogid_prot_records.append(prot_record)
                species_idx_rhogid.append(species_idx)

            if print_hog_more_than_one_species:
                if len(set(species_idx_rhogid)) > 1:
                    rhogids_prot_records_query.append(rhogid_prot_records)
            else:  # print roothog even with one species
                rhogids_prot_records_query.append(rhogid_prot_records)
                # else:
        #    print("root hog na / lenght of one ",rhogid)

    # print(len(rhogids_prot_records_query),len(rhogids_prot_records_query[0]))
    rhogid_num_list = []
    # for rhogid_idx in range(len(rhogids_list)):
    #     rhogid = rhogids_list[rhogid_idx]
    for rhogid_idx, rhogid in enumerate(rhogids_list):
        rhogid_prot_rec_query = rhogids_prot_records_query[rhogid_idx]
        rhogid_B = rhogid.split(":")[1]
        rhogid_num = int(rhogid_B[1:])  # # B0613860
        rhogid_num_list.append(rhogid_num)

        if 1 < len(rhogid_prot_rec_query) < 500:
            SeqIO.write(rhogid_prot_rec_query, address_rhogs_folder +"HOG_B"+ str(rhogid_num).zfill(7)+".fa", "fasta")
            # rhogids_prot_records_oma = []
            # for hog_elements in oma_db.member_of_fam(rhogid_num):   # this gets the member of roothog 2 (HOG:000002)
            #    prot_hog_element = ProteinEntry(oma_db, hog_elements)
            #    #print(prot_hog_element.omaid, prot_hog_element.hog_family_nr, len(prot_hog_element.sequence),prot_hog_element.sequence[0])
            #    rhogids_prot_records_oma.append(SeqRecord(Seq(prot_hog_element.sequence), id=prot_hog_element.omaid))
            # rhogids_prot_records_both= rhogids_prot_records_oma +  rhogid_prot_records_query
            # rhogids_prot_records.append(rhogids_prot_records_both)

    # print("all hogs   (>2 <100) has written.",len(rhogids_prot_records_query),len(rhogids_list),
    # len(rhogid_prot_records_query), len(rhogid_prot_records_query[0]))

    current_time = datetime.now().strftime("%H:%M:%S")
    print(current_time, "- Sequences of roothogs are writtend as fasta file in " + address_rhogs_folder)

    return rhogid_num_list, rhogids_prot_records_query
