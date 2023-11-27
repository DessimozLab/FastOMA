

import argparse
import sys
import logging
logger_level = "DEBUG"            # DEBUG INFO
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
logger_hog = logging.getLogger("hog")

if logger_level == "INFO":
    logger_hog.setLevel(logging.INFO)
if logger_level == "DEBUG":
    logger_hog.setLevel(logging.DEBUG)  # TRACE  DEBUG INFO  WARN  ERROR  FATAL


input_rhog_folder = "./"
species_tree_address = "species_tree.nwk"
species_tree_checked = "species_tree_checked.nwk"

protein_format_qfo_dataset_before2022 = True

fragment_detection_msa_merge = True  # if this is false and fragment_detection_msa -> we'll remove both fragments in orthology analsys at parent level but report it in orthoxml Dubiousfragment

threshold_dubious_sd = 1/10
overlap_fragments = 0.15

mergHOG_ratioMax_thresh = 0.8
mergHOG_ratioMin_thresh = 0.9
mergHOG_shared_thresh = 70
threshod_f_score_merging = 50
mergHOG_mean_thresh=5000*1000


add_outgroup = False
rooting_method = "midpoint" #"Nevers_rooting" #"midpoint"#"midpoint"  # "midpoint" "mad"


gene_trees_write = False
msa_write = False
gene_trees_write_all = False
msa_write_all = False
keep_subhog_each_pickle = False

big_rhog_size = 50 * 1000
omamer_family_threshold = 0

subsampling_hogclass = True
hogclass_max_num_seq = 20  # 40 subsampling in msa # ver very 2
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 50
hogclass_tresh_ratio_gap_col = 0.6  # 0.8 for very very big

filter_nonchild_rootHOG= False

num_msas_merge_mafft = 100
if hogclass_max_num_seq < 10:
    num_msas_merge_mafft= 10   # small value means always use simple mafft (not mafft merge)

automated_trimAL = False

label_SD_internal = "species_overlap"  # "reconcilation" "species_overlap"
tree_tool = "fasttree"  # "fasttree"  "iqtree"  # for  gene tree with two, we use
rooting_mad_executable_path = "mad"  # /work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/installers/mad/
mmseqs_executable_path ="mmseqs"


inferhog_tresh_ratio_gap_row =0.3 # 0.6   # to have more proteins in the ortho groups 0.1
inferhog_tresh_ratio_gap_col =0.5  # 0.6   # ver very 0.8
inferhog_min_cols_msa_to_filter = 50 #300 #50  # used for msa before gene tree inference and  saving msa in hog class

inferhog_filter_all_msas_row = True
inferhog_resume_rhog = True  # main.py False
# The intermediate files, internal node  pickle files is not working with nextflow. the reason is that the pickles_subhog_folder_all is relative and stored in nextflow_work folder.  this folder can not be used for  the re-submitting
# todo: resume use files pcikels file of a level from previous run
inferhog_resume_subhog = True  # read pickle_subhog  # _infer_subhog.py

inferhog_max_workers_num = 6 # for big rootHOG, we use parallelisation on taxanomic levels of species tree.

inferhog_min_hog_size_xml = 2     # by setting this as 1, pyham won't work on xml output.

big_rhog_filesize_thresh = 400 * 1000
sum_list_rhogs_filesize_thresh = 2 * 1e6


orthoxml_v03 = True


def set_configs():
    parser = argparse.ArgumentParser(description="This is FastOMA ")
    parser.add_argument('--logger-level', default="DEBUG")
    parser.add_argument("--version", action="version", help="Show version and exit.", version="0.1.4",)  # version=__version__
    parser.add_argument('--species-tree-checked', default="species_tree_checked.nwk")
    parser.add_argument('--input-rhog-folder')     # , default="./rhog"
    parser.add_argument('--parallel', action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument('--fragment-detection', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--low-so-detection', action=argparse.BooleanOptionalAction, default=True)

    config_parser = parser.parse_args()
    # Namespace(logger_level=None, in_folder=None)
    setattr(sys.modules[__name__], 'logger_level', config_parser.logger_level)
    setattr(sys.modules[__name__], 'input_rhog_folder', config_parser.input_rhog_folder)
    setattr(sys.modules[__name__], 'parallel', config_parser.parallel)
    setattr(sys.modules[__name__], 'species_tree_checked', config_parser.species_tree_checked) # todo _checked.nwk
    setattr(sys.modules[__name__], 'fragment_detection', config_parser.fragment_detection)
    setattr(sys.modules[__name__], 'low_so_detection', config_parser.low_so_detection)
    print("config_parser 3 ", config_parser)

