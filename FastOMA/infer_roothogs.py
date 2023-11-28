import os.path
from shutil import which

from . import _utils_roothog
from . import _config
logger = _config.logger_hog

"""
proteomes of species as fasta files in /proteome/
omamer's output of  in /hogmap/
hog and HOG are used interchangeably here. 
"""


def fastoma_infer_roothogs():
    import argparse
    parser = argparse.ArgumentParser(description="checking parameters for FastOMA")
    parser.add_argument("--proteomes", required=True, help="Path to the folder containing the input proteomes")
    parser.add_argument("--splice", help="Path to the folder containing the splice information files")
    parser.add_argument("--hogmap", help="Path to the folder containing the hogmap files")
    parser.add_argument("--out-rhog-folder", required=True, help="Folder where the roothog fasta files are written")
    parser.add_argument('-v', action="count", default=0, help="Increase verbosity to info/debug")
    conf = parser.parse_args()
    logger.setLevel(level=30 - 10 * min(conf.v, 2))
    logger.debug("Arguments: %s", conf)

    species_names, prot_recs_lists, fasta_format_keep = _utils_roothog.parse_proteomes(conf.proteomes)  # optional input folder
    prot_recs_all = _utils_roothog.add_species_name_prot_id(species_names, prot_recs_lists)

    hogmaps, unmapped = _utils_roothog.parse_hogmap_omamer(species_names, fasta_format_keep, folder=conf.hogmap)  # optional input folder

    splice_files = conf.splice is not None and os.path.exists(conf.splice)
    if splice_files:
        isoform_by_gene_all = _utils_roothog.parse_isoform_file(species_names, folder=conf.splice)
        isoform_selected,  isoform_not_selected = _utils_roothog.find_nonbest_isoform(
            species_names, isoform_by_gene_all, hogmaps
        )
        _utils_roothog.write_isoform_selected(isoform_by_gene_all, isoform_selected, prot_recs_lists)
        # for each isoform file, there will be a file ending with _selected_isoforms.tsv
        hogmaps = _utils_roothog.handle_splice(hogmaps, isoform_not_selected)

    rhogs_prots = _utils_roothog.group_prots_roothogs(hogmaps)
    rhogs_prots = _utils_roothog.handle_singleton(rhogs_prots, hogmaps)
    rhogs_prots = _utils_roothog.merge_rhogs(hogmaps, rhogs_prots)
    rhogs_prots = _utils_roothog.filter_big_roothogs(hogmaps, rhogs_prots)


    min_rhog_size = 2
    rhogid_written_list = _utils_roothog.write_rhog(rhogs_prots, prot_recs_all, conf.out_rhog_folder, min_rhog_size)
    linclust_available=which("mmseqs")  # True #
    # if memseqs is not installed the output will be empty / None
    if linclust_available:
        num_unmapped_singleton = _utils_roothog.collect_unmapped_singleton(rhogs_prots, unmapped, prot_recs_all,  "singleton_unmapped.fa")
        if num_unmapped_singleton:
            result_linclust = _utils_roothog.run_linclust(fasta_to_cluster="singleton_unmapped.fa")
            logger.debug(" linclust is done %s", result_linclust)
            num_clusters = _utils_roothog.write_clusters(conf.out_rhog_folder, min_rhog_size)
            logger.debug("we wrote %d new clusters with linclust ", num_clusters)
