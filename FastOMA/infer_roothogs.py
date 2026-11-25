import os.path
import sys

from . import _utils_roothog, logger
from .zoo.wrappers import WrapperError
from .zoo.wrappers.clusterers import cluster_sequences, find_mmseqs
from . import _protein_trace
from . import __version__ as fastoma_version
from .logging_setup import setup_logging



"""

fastoma-infer-roothogs --proteomes proteome --hogmap hogmap --out-rhog-folder omamer_rhogs -vv

without OMAmer, e.g. for simulated proteomes, all proteins are clustered de novo with mmseqs:

fastoma-infer-roothogs --proteomes proteome --denovo --out-rhog-folder omamer_rhogs -vv

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

    parser.add_argument("--denovo", action="store_true",
                        help="infer the rootHOGs by de-novo clustering of all proteins with mmseqs instead of "
                             "OMAmer placements (e.g. for simulated data); --hogmap is ignored")
    mmseqs = parser.add_argument_group(
        "clustering", "the unmapped proteins and the singletons (with --denovo: all proteins) are clustered with "
                      "mmseqs easy-cluster; options that are not given keep the mmseqs defaults")
    mmseqs.add_argument("--mmseqs-threads", type=int, help="number of threads of mmseqs (default: all cores)")
    mmseqs.add_argument("--mmseqs-coverage", type=float, help="minimum alignment coverage (mmseqs -c)")
    mmseqs.add_argument("--mmseqs-cov-mode", type=int, help="how the coverage is defined (mmseqs --cov-mode)")
    mmseqs.add_argument("--mmseqs-sensitivity", type=float, help="sensitivity (mmseqs -s)")
    mmseqs.add_argument("--mmseqs-min-seq-id", type=float, help="minimum sequence identity (mmseqs --min-seq-id)")
    parser.add_argument("--protein-trace", action="store_true",
                        help="write protein_events_*.events files, used by fastoma-protein-trace to trace every input protein")

    conf = parser.parse_args()
    setup_logging(conf.v)
    logger.debug("Arguments: %s", conf)
    if conf.protein_trace:
        _protein_trace.enable()

    # Step 1: Parse input data
    logger.info("Parsing proteomes...")
    species_names, prot_recs_lists, fasta_format_keep = _utils_roothog.parse_proteomes(
        conf.proteomes, conf.min_sequence_length
    )
    logger.info("Processing protein records...")
    prot_recs_all = _utils_roothog.add_species_name_prot_id(prot_recs_lists)

    linclust_available = find_mmseqs()
    if conf.denovo:
        if not linclust_available:
            logger.error("De-novo clustering (--denovo) needs mmseqs, which is not on the PATH.")
            sys.exit(1)
        logger.info("De-novo clustering: all proteins are clustered with mmseqs, no OMAmer placements")
        hogmaps, unmapped = _utils_roothog.denovo_unmapped(prot_recs_lists)
    else:
        logger.info("Parsing HOGMAP files...")
        hogmaps, unmapped = _utils_roothog.parse_hogmap_omamer(
            prot_recs_lists, fasta_format_keep, folder=conf.hogmap
        )

    # Step 2: Handle splice isoforms
    isoform_data = None
    if conf.splice is not None and os.path.exists(conf.splice): 
        logger.info("Processing splice variants...")
        isoform_data = _utils_roothog.handle_splice_variants(species_names, hogmaps, conf.splice,
                                                             {sp: set(prots) for sp, prots in prot_recs_all.items()})
        _utils_roothog.write_selected_isoforms(isoform_data, prot_recs_lists)
        hogmaps = isoform_data['filtered_hogmaps']
        unmapped = _utils_roothog.filter_unmapped_by_isoforms(unmapped, isoform_data['non_selected'])
    _utils_roothog.check_all_proteins_accounted(prot_recs_all, hogmaps, unmapped, isoform_data)

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
    cnt_unmapped_singleton = 0
    if linclust_available:
        to_cluster = _utils_roothog.collect_unmapped_singleton(rhogs_prots, unmapped, prot_recs_all,
                                                               unmapped_singleton_fasta="singleton_unmapped.fa")
        cnt_unmapped_singleton = len(to_cluster)
        try:
            clusters = cluster_sequences(to_cluster, executable=linclust_available, threads=conf.mmseqs_threads,
                                         coverage=conf.mmseqs_coverage, cov_mode=conf.mmseqs_cov_mode,
                                         sensitivity=conf.mmseqs_sensitivity, min_seq_id=conf.mmseqs_min_seq_id)
        except WrapperError as e:
            # e.g. killed because it ran out of memory: pass the exit code on, so Nextflow can retry the task
            logger.error("%s", e)
            sys.exit(e.exit_code)
        num_clusters = _utils_roothog.write_clusters(clusters, conf.out_rhog_folder, min_rhog_size)
        logger.info("%d of %d mmseqs clusters written as rootHOGs", num_clusters, len(clusters))
    else:
        logger.warning("mmseqs linclust / easy-cluster not available, skipping clustering of unmapped and singletons")

    # every protein is in exactly one of these (checked above), so the big-rootHOG filter is the rest
    n_input = sum(len(prots) for prots in prot_recs_all.values())
    n_in_rhogs = sum(len(rhogs_prots[rhogid]) for rhogid in rhogid_written_list)
    n_singletons = sum(len(prots) for prots in unmapped.values()) + sum(len(p) == 1 for p in rhogs_prots.values())
    n_isoforms = sum(len(v) for v in _utils_roothog.dropped_isoforms(prot_recs_all, hogmaps, unmapped,
                                                                     isoform_data).values())
    logger.info("%d proteins (after the length filter): %d in %d rootHOGs, %d unmapped or singletons %s, "
                "%d isoforms not selected, %d removed from big rootHOGs", n_input, n_in_rhogs,
                len(rhogid_written_list), n_singletons,
                "sent to clustering" if linclust_available else "not clustered (no mmseqs)",
                n_isoforms, n_input - n_in_rhogs - n_singletons - n_isoforms)


if __name__ == "__main__":
    fastoma_infer_roothogs()
