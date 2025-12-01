import os.path
import subprocess
from shutil import which

from . import _utils_roothog, logger
from . import __version__ as fastoma_version
from .logging_setup import setup_logging



"""

fastoma-infer-roothogs --proteomes proteome --hogmap hogmap --out-rhog-folder omamer_rhogs -vv

"""


def fastoma_infer_roothogs():
    import argparse
    parser = argparse.ArgumentParser(description="checking parameters for FastOMA")
    parser.add_argument("--version", action="version", version="FastOMA v"+fastoma_version)
    parser.add_argument("--proteomes", required=True, help="Path to the folder containing the input proteomes")
    parser.add_argument("--splice", help="Path to the folder containing the splice information files")
    parser.add_argument("--hogmap", help="Path to the folder containing the hogmap files")
    parser.add_argument("--out-rhog-folder", required=True, help="Folder where the roothog fasta files are written") #out_rhog_folder
    parser.add_argument('-v', action="count", default=0, help="Increase verbosity to info/debug")
    parser.add_argument('--min-sequence-length', required=False, default=50, type=int,
                        help="minimum sequence length. Shorter sequences will be ignored. (Default=50)")

    parser.add_argument("--mergHOG-ratioMax-thresh", required=False, type=float, default=0.8, help="For merging rootHOGs, threshold of ratioMax ") # mergHOG_ratioMax_thresh
    parser.add_argument("--mergHOG-ratioMin-thresh", required=False, type=float, default=0.9, help="For merging rootHOGs, threshold of ratioMin ") # mergHOG_ratioMin_thresh
    parser.add_argument("--mergHOG-shared-thresh", required=False, type=float, default=10, help="For merging rootHOGs, threshold of number shared proteins ") # mergHOG_shared_thresh
    parser.add_argument("--mergHOG-fscore-thresh", required=False, type=float, default=70, help="For merging rootHOGs, threshold of famlut score shared proteins ") # mergHOG_fscore_thresh
    parser.add_argument("--big-rhog-size", required=False, type=int, default=50*1000, help= "For big rootHOGs, we have different heuristics") # big_rhog_size
    parser.add_argument("--big-fscore-thresh", required=False, type=int, default=95, help="For huge rootHOGs, we have different heuristics, like filtering low family score protiens") # big_fscore_thresh

    conf = parser.parse_args()
    setup_logging(conf.v)
    logger.debug("Arguments: %s", conf)

    # Step 1: Parse input data
    logger.info("Parsing proteomes...")
    species_names, prot_recs_lists, fasta_format_keep = _utils_roothog.parse_proteomes(
        conf.proteomes, conf.min_sequence_length
    )
    logger.info("Processing protein records...")
    prot_recs_all = _utils_roothog.add_species_name_prot_id(prot_recs_lists)

    logger.info("Parsing HOGMAP files...")
    hogmaps, unmapped = _utils_roothog.parse_hogmap_omamer(
        prot_recs_lists, fasta_format_keep, folder=conf.hogmap
    )

    # Step 2: Handle splice isoforms
    isoform_data = None
    if conf.splice is not None and os.path.exists(conf.splice): 
        logger.info("Processing splice variants...")
        isoform_data = _utils_roothog.handle_splice_variants(species_names, hogmaps, conf.splice)
        _utils_roothog.write_selected_isoforms(isoform_data, prot_recs_lists)
        hogmaps = isoform_data['filtered_hogmaps']

    # Step 3: Group proteins into root HOGs
    logger.info("Grouping proteins into root HOGs...")
    # rhogs_prots = _utils_roothog.group_prots_roothogs(hogmaps)
    # rhogs_prots = _utils_roothog.handle_singleton(rhogs_prots, hogmaps, conf)
    # rhogs_prots = _utils_roothog.merge_rhogs2(hogmaps, rhogs_prots, conf)
    # rhogs_prots = _utils_roothog.filter_big_roothogs(hogmaps, rhogs_prots, conf)

    rhogs_prots = _utils_roothog.create_root_hogs(hogmaps, conf)


    # Step 4: Save results
    logger.info("Saving results...")
    _utils_roothog.save_gene_id_mapping(prot_recs_all, isoform_data)
    # _utils_roothog.write_root_hogs(rhogs_prots, prot_recs_all, conf.out_rhog_folder)
    
    logger.info(f"Successfully created {len(rhogs_prots)} root HOGs")


    min_rhog_size = 2
    rhogid_written_list = _utils_roothog.write_rhog(rhogs_prots, prot_recs_all, conf.out_rhog_folder, min_rhog_size)
    linclust_available = which("mmseqs")  # True #
    # if memseqs is not installed the output will be empty / None
    if linclust_available:
        try:
            res = subprocess.run([linclust_available, "-h"], check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            logger.error("mmseqs linclust is not working: %s", e)
            logger.error("Skipping clustering of unmapped singletons")
            linclust_available = False
    if linclust_available:
        singleton_unmapped_path = "singleton_unmapped.fa"
        cnt_unmapped_singleton = _utils_roothog.collect_unmapped_singleton(rhogs_prots, unmapped, prot_recs_all, unmapped_singleton_fasta=singleton_unmapped_path)
        if cnt_unmapped_singleton:
            cluster_file = _utils_roothog.run_linclust(fasta_to_cluster="singleton_unmapped.fa")
            num_clusters = _utils_roothog.write_clusters(cluster_file, conf.out_rhog_folder, min_rhog_size)
            logger.debug("we wrote %d new clusters with linclust", num_clusters)
    else:
        logger.info("mmseqs linclust / easy-cluster not available, skipping clustering of unmapped and singletons")


if __name__ == "__main__":
    fastoma_infer_roothogs()
