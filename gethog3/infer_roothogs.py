
from ._utils_subhog import logger_hog
from . import _utils_roothog
from . import _config

"""
proteomes of species as fasta files in /proteome/
omamer's output of  in /hogmap/
hog and HOG are used interchangeably here. 
"""


# if not os.path.exists(_config.working_folder):
#    os.mkdir(_config.working_folder) # s

# working_folder+"omamer_database/oma_path/OmaServer.h5"
# logger_hog.info("rHOG inferece has started. The oma database address is in " + _config.oma_database_address)
# (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(_config.oma_database_address)
# (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species)



def infer_roothogs():
    _config.set_configs()

    # print(_config.working_folder)
    print(_config.logger_level)
    # import sys
    # sys.exit(0)

    (query_species_names, query_prot_recs) = _utils_roothog.parse_proteome()

    print("query_species_names ", len(query_prot_recs))

    query_prot_recs = _utils_roothog.add_species_name_gene_id(query_prot_recs, query_species_names)
    hogmap_allspecies_elements = _utils_roothog.parse_hogmap_omamer(query_species_names)

    (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
     prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies,
     prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

    query_prot_recs_filt = _utils_roothog.filter_prot_mapped(query_species_names, query_prot_recs, query_prot_names_species_mapped)

    logger_hog.info("size of query_prot_recs_filt is " + str(len(query_prot_recs_filt)) + " " + str(
        len(query_prot_recs_filt[0])))

    rhogids_list, rhogids_prot_records_query = _utils_roothog.group_prots_roothogs(prots_hogmap_hogid_allspecies,
                                                                                query_species_names,
                                                                                query_prot_recs_filt)

    # rhogid_num_list_raw=utils_rhog.write_rhog(rhogids_list,rhogids_prot_records_query, _config.working_folder+"rhogs/". "rhogs_raw",2)

    rhogids_list_filt, rhogids_prot_records_query_filt = _utils_roothog.filter_rhog(rhogids_list,
                                                                                 rhogids_prot_records_query,
                                                                                 prots_hogmap_fscore_allspecies,
                                                                                 query_species_names,
                                                                                 query_prot_names_species_mapped)
    # for pure usage of this python file, you can set the output folder
    # output_folder_rhog = _config.working_folder + "rhogs_all/"
    # using nextflow



    # import sys
    output_folder_rhog = "./" #  sys.argv[1]  #
    rhogid_num_list_filt1 = _utils_roothog.write_rhog(rhogids_list_filt, rhogids_prot_records_query_filt, output_folder_rhog, 2)  # min_rhog_size, max_rhog_size
