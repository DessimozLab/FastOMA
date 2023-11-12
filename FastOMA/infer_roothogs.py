

import os.path
from shutil import which

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

    hogmaps, unmapped = _utils_roothog.parse_hogmap_omamer(species_names,fasta_format_keep) # optional input folder

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
    rhogs_prots = _utils_roothog.filter_big_roothogs(hogmaps, rhogs_prots)


    address_rhogs_folder = "./temp_omamer_rhogs/"
    min_rhog_size =2
    rhogid_written_list = _utils_roothog.write_rhog(rhogs_prots, prot_recs_all, address_rhogs_folder, min_rhog_size)
    linclust_available=which("mmseqs") #  True #
    # if memseqs is not installed the output will be empty / None
    if linclust_available:
        num_unmapped_singleton = _utils_roothog.collect_unmapped_singleton(rhogs_prots, unmapped, prot_recs_all,  "singleton_unmapped.fa")
        if num_unmapped_singleton:
            result_linclust = _utils_roothog.run_linclust(fasta_to_cluster="singleton_unmapped.fa")
            logger_hog.debug(" linclust is done "+ result_linclust)
            num_clusters = _utils_roothog.write_clusters(address_rhogs_folder, min_rhog_size)
            logger_hog.debug("we wrote "+str(num_clusters)+"  new clusters with linclust ")
