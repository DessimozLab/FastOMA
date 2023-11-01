

import os.path

from ._utils_subhog import logger_hog
from . import _utils_roothog
from . import _config


"""
proteomes of species as fasta files in /proteome/
omamer's output of  in /hogmap/
hog and HOG are used interchangeably here. 
"""


# if not os.path.exists(_config.in_folder):
#    os.mkdir(_config.in_folder) # s

# in_folder+"omamer_database/oma_path/OmaServer.h5"
# logger_hog.info("rHOG inferece has started. The oma database address is in " + _config.oma_database_address)
# (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(_config.oma_database_address)
# (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species)



def infer_roothogs():
    _config.set_configs()

    # print(_config.in_folder)
    print(_config.logger_level)
    # import sys
    # sys.exit(0)

    #folder = "/scratch/smajidi1/euk_omamer200.dev8_2/test/hogmap/FastOMA-main/testdata/in_folder/"
    species_names, prot_recs_lists,fasta_format_keep = _utils_roothog.parse_proteomes() # optional input folder
    prot_recs_all = _utils_roothog.add_species_name_prot_id(species_names, prot_recs_lists)

    hogmaps = _utils_roothog.parse_hogmap_omamer(species_names,fasta_format_keep) # optional input folder

    splice_files =  os.path.exists("./splice/")
    if splice_files:
        isoform_by_gene_all = _utils_roothog.parse_isoform_file(species_names)
        isoform_selected,  isoform_not_selected = _utils_roothog.find_nonbest_isoform(species_names,isoform_by_gene_all,hogmaps)
        _utils_roothog.write_isoform_selected(isoform_by_gene_all, isoform_selected,prot_recs_lists)
        # for each isoform file, there will be a file ending with _selected_isoforms.tsv
        hogmaps = _utils_roothog.handle_splice(hogmaps,isoform_not_selected)


    rhogs_prots = _utils_roothog.group_prots_roothogs(hogmaps)

    rhogs_prots = _utils_roothog.handle_singleton(rhogs_prots,hogmaps)
    rhogs_prots = _utils_roothog.merge_rhogs(hogmaps, rhogs_prots)
    rhogs_prots = _utils_roothog.roothogs_postprocess(hogmaps, rhogs_prots)
    address_rhogs_folder = "./omamer_rhogs/"
    rhogid_written_list = _utils_roothog.write_rhog(rhogs_prots, prot_recs_all, address_rhogs_folder, min_rhog_size=2)
        #(rhogs_prot_records, address_rhogs_folder, min_rhog_size=2)

    #(query_species_names, query_prot_recs) = _utils_roothog.parse_proteome()
    #print("query_species_names ", len(query_prot_recs))
    #query_prot_recs = _utils_roothog.add_species_name_gene_id(query_prot_recs, query_species_names)
    #hogmap_allspecies_elements = _utils_roothog.parse_hogmap_omamer(query_species_names)
    # (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
    #  prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies,
    #  prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
    #query_prot_recs_filt = _utils_roothog.filter_prot_mapped(query_species_names, query_prot_recs, query_prot_names_species_mapped)
    #logger_hog.info("size of query_prot_recs_filt is " + str(len(query_prot_recs_filt)) + " " + str(len(query_prot_recs_filt[0])))

    # splice_files =  os.path.exists("./splice/")
    # if splice_files:
    #     isoform_by_gene_all = _utils_roothog.parse_isoform_file(query_species_names)
    #     #query_species_names[0], isoform_by_gene_all[0], isoform_by_gene_all[2][:2]
    #     not_selected_isofroms_all = _utils_roothog.find_nonbest_isoform(hogmap_allspecies_elements, isoform_by_gene_all)
    #     prots_hogmap_hogid_allspecies_ ,  query_prot_recs_filt_ = _utils_roothog.handle_splice(prots_hogmap_hogid_allspecies, query_prot_recs_filt, not_selected_isofroms_all, query_prot_names_species_mapped)
    # else:
    #     prots_hogmap_hogid_allspecies_= prots_hogmap_hogid_allspecies
    #     query_prot_recs_filt_ = query_prot_recs_filt


    #rhogids_list, rhogids_prot_records_query = _utils_roothog.group_prots_roothogs(prots_hogmap_hogid_allspecies_, query_species_names, query_prot_recs_filt_)
    # rhogid_list_raw=utils_rhog.write_rhog(rhogids_list,rhogids_prot_records_query, _config.in_folder+"rhogs/". "rhogs_raw",2)
    #rhogids_list_filt, rhogids_prot_records_query_filt = _utils_roothog.filter_rhog(rhogids_list, rhogids_prot_records_query,prots_hogmap_fscore_allspecies,query_species_names, query_prot_names_species_mapped)
    # for pure usage of this python file, you can set the output folder
    # output_folder_rhog = _config.in_folder + "rhogs_all/" // omamer_rhogs
    # using nextflow
    # import sys
    #output_folder_rhog = "./omamer_rhogs/" #  sys.argv[1]  #
    #rhogid_list_filt1 = _utils_roothog.write_rhog(rhogids_list_filt, rhogids_prot_records_query_filt, output_folder_rhog, 2)  # min_rhog_size, max_rhog_size
    # the list of rhogid_list_filt1 only contain represetative roothog for merged rhogs

