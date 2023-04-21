
# global path
import argparse
import sys


# setattr(sys.modules[__name__], 'in_folder', config_parser.in_folder)



# mkdb_parser.add_argument(
#     "--min_fam_size", default=6, help="Only root-HOGs with a protein count passing this threshold are used.",
#     type=int
# )

#  to have more proteins in the ortho groups
#  omamer_fscore_treshold_big_rhog = 0.05
#  inferhog_tresh_ratio_gap_row = 0.1


import logging


logger_level = "DEBUG"            # DEBUG INFO



logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger_hog = logging.getLogger("hog")

if logger_level == "INFO":
    logger_hog.setLevel(logging.INFO)  # DEBUG WARN  INFO
if logger_level == "DEBUG":
    logger_hog.setLevel(logging.DEBUG)  # DEBUG WARN  INFO
# TRACE  DEBUG INFO  WARN  ERROR  FATAL


#
# logging.basicConfig(
#     format='%(asctime)s %(levelname)-8s %(message)s',
#     level=logging.INFO,
#     datefmt='%Y-%m-%d %H:%M:%S')
# logger_hog = logging.getLogger("hog")
#
# if _config.logger_level == "INFO":
#     logger_hog.setLevel(logging.INFO)  # DEBUG WARN  INFO
# if _config.logger_level == "DEBUG":
#     logger_hog.setLevel(logging.DEBUG)  # DEBUG WARN  INFO
# # TRACE  DEBUG INFO  WARN  ERROR  FATAL



input_rhog_folder = "./"

# in_folder = "./in_folder"+ "/"
species_tree_address = "species_tree.nwk"
# no space or special charcter in internal node.
# protein_format_qfo_dataset = True

fragment_detection = False

## output writing files
gene_trees_write = False
msa_write = False
gene_trees_write_all = False
msa_write_all = False
keep_subhog_each_pickle = False


# filtering omamer
omamer_fscore_treshold_big_rhog = 0.5  #  to have more proteins in the ortho groups 0.05  considering for big rhogs
treshold_big_rhog_szie = 1000

## hogclass configs
hogclass_max_num_seq = 20  # subsampling in msa
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 300
hogclass_tresh_ratio_gap_col = 0.2

automated_trimAL = False

label_SD_internal = "species_overlap"  # "reconcilation" "species_overlap"
threshold_dubious_sd = 1/10+0.01
#threshold_sd_suspicious_fragment_ratio = 1/3
tree_tool = "fasttree"  # "fasttree"  "iqtree"  # for  gene tree with two, we use

rooting_method = "midpoint"  # "midpoint" "mad"
rooting_mad_executable_path = "mad" # /work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/installers/mad/

##inferhog
inferhog_tresh_ratio_gap_row = 0.1   # to have more proteins in the ortho groups 0.1
inferhog_tresh_ratio_gap_col = 0.4
inferhog_min_cols_msa_to_filter = 100  # used for msa before gene tree inference and  saving msa in hog class

inferhog_filter_all_msas_row = True


inferhog_resume_rhog = False  # main.py False
# The intermediate files, internal node  pickle files is not working with nextflow
# the reason is that the pickles_subhog_folder_all is relative and stored in nextflow_work folder
# this folder can not be used for  the re-submitting
inferhog_resume_subhog = True  # read pickle_subhog  # _infer_subhog.py

# inferhog_concurrent_on = True now as an argument
inferhog_max_workers_num = 8

## xml
# write_all_prots_in_header = False  # if false writes only those in the hog group
inferhog_min_hog_size_xml = 2     # by setting this as 1, pyham won't work on xml output.

# batch_roothogs
big_rhog_filesize_thresh = 600 * 1000  # 600 would be better
sum_list_rhogs_filesize_thresh = 2 * 1e6

# big_rhog_filesize_thresh = 1.6 * 1000  # 600 would be better
# sum_list_rhogs_filesize_thresh = 5 * 1e3


def set_configs():
    parser = argparse.ArgumentParser(description="This is GETHOG3 ")
    # parser.add_argument('--working-folder', help="in_folder")
    parser.add_argument('--logger-level', default="DEBUG")
    #  $rhogs_big_i - -parrallel

    parser.add_argument("--version", action="version", help="Show version and exit.",
        version="0.0.6",) # version=__version__

    parser.add_argument('--species-tree-address', default="species_tree_test.nwk")
    parser.add_argument('--input-rhog-folder') # , default="./rhog"
    parser.add_argument('--parrallel', default=False)

    config_parser = parser.parse_args()

    # Namespace(logger_level=None, in_folder=None)
    setattr(sys.modules[__name__], 'logger_level', config_parser.logger_level)
    setattr(sys.modules[__name__], 'input_rhog_folder', config_parser.input_rhog_folder)
    setattr(sys.modules[__name__], 'parrallel', config_parser.parrallel)
    setattr(sys.modules[__name__], 'species_tree_address', config_parser.species_tree_address)

    print("config_parser 3 ", config_parser)

'''

TODO
 default doesnt work
'''