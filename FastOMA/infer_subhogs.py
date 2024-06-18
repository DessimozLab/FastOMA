import os

from . import _utils_subhog
from . import _infer_subhog
from ._wrappers import logger
from . import __version__ as fastoma_version

"""

fastoma-infer-subhogs  --input-rhog-folder rhogs_rest/0  --output-pickles "pickle_hogs"  \
    --species-tree  species_tree_checked.nwk -vv --parallel # --msa-write --gene-trees-write
    
"""

def fastoma_infer_subhogs():

    import argparse
    parser = argparse.ArgumentParser(description="checking parameters for FastOMA",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", action="version", version="FastOMA v"+fastoma_version)
    parser.add_argument("--input-rhog-folder", required=True, help="Path to the input rootHOG folder.")
    parser.add_argument("--parallel", action='store_true', help="use concurrent parallel per rootHOG")
    parser.add_argument("--species-tree", required=True,
                        help="Path to the input species tree file in newick format")
    parser.add_argument("--output-pickles", required=False, default="pickle_hogs",
                        help="Path to the output folder")

    parser.add_argument("--threshold-dubious-sd", required=False, type=float, default=1/10,
                        help="Threshold to remove proteins in a gene tree due to low species overlap score, not enough evidence for duplication event.") # threshold_dubious_sd
    parser.add_argument("--number-of-samples-per-hog", type=int, default=5,
                        help="Number of representatives (sequences) per HOG. Defaults to ")
    parser.add_argument("--overlap-fragments", required=False, type=float, default=0.15,
                        help="Threshold overlap between two sequences (rows) in MSA to decide whether they are fragments of a gene.")  # overlap_fragments
    parser.add_argument("--gene-rooting-method", required=False, default="midpoint", # gene_rooting_method
                        help="The method used for rooting of gene tree :    midpoint    mad     Nevers_rooting .")
    parser.add_argument("--gene-trees-write", action='store_true',
                        help="writing the all gene trees .")  # the order seems to be nwk_SD_labeled.nwk, dubious_sd0.nwk_SD_labeled.nwk, dubious_sd1.nwk_SD_labeled.nwk
    parser.add_argument("--msa-write", action='store_true',
                        help="writing the raw MSAs (might have more genes that the final gene tree).")
    parser.add_argument("--msa-filter-method",
                        choices=("col-row-threshold", "col-elbow-row-threshold", "trimal"),
                        default="col-row-threshold",
                        help="The method used for filtering MSAs.")
    parser.add_argument("--gap-ratio-row", required=False, type=float, default=0.3,
                        help="For trimming the MSA, the threshold of ratio of gaps for each row.")
    parser.add_argument("--gap-ratio-col", required=False, type=float, default=0.5,
                        help="For trimming the MSA, the threshold of ratio of gaps for each column.")
    parser.add_argument("--min-col-trim", required=False, type=int, default=50,  # todo min rows trim
                        help="min no. columns in msa to consider for filtering")
    parser.add_argument('-v', action="count", default=0, help="Increase verbosity to info/debug")
    conf_infer_subhhogs = parser.parse_args()
    logger.setLevel(level=30 - 10 * min(conf_infer_subhhogs.v, 2))
    logger.debug("Arguments: %s", conf_infer_subhhogs)

    address_rhogs_folder = conf_infer_subhhogs.input_rhog_folder
    # address_rhogs_folder = "./"  # _config.input_rhog_folder
    inferhog_concurrent_on = conf_infer_subhhogs.parallel
    if inferhog_concurrent_on:
        print("parallelization for subhog inference is on.")

    if not os.path.exists(conf_infer_subhhogs.output_pickles):
        os.makedirs(conf_infer_subhhogs.output_pickles)

    pickles_subhog_folder_all = "./" # pickle per taxonomic level

    list_rhog_fastas_files = _utils_subhog.list_rhog_fastas(address_rhogs_folder)
    print("there are ", len(list_rhog_fastas_files), "rhogs in the input folder")

    rhogs_fa_folder = address_rhogs_folder

    list_rhog_fastas_files_rem = _utils_subhog.list_rhog_fastas(address_rhogs_folder)
    print("there are ", len(list_rhog_fastas_files_rem), "rhogs remained in the input folder", list_rhog_fastas_files_rem[:5] )

    hogs_rhog_xml_batch = _infer_subhog.read_infer_xml_rhogs_batch(list_rhog_fastas_files_rem, inferhog_concurrent_on, conf_infer_subhhogs.output_pickles, pickles_subhog_folder_all, rhogs_fa_folder, conf_infer_subhhogs)

    print("finsihed ", address_rhogs_folder)

    threshold_dubious_sd= 0.1


if __name__ == "__main__":
    fastoma_infer_subhogs()