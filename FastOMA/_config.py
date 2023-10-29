



import argparse
import sys
import logging
logger_level = "DEBUG"            # DEBUG INFO
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger_hog = logging.getLogger("hog")
if logger_level == "INFO":
    logger_hog.setLevel(logging.INFO)
if logger_level == "DEBUG":
    logger_hog.setLevel(logging.DEBUG)  # TRACE  DEBUG INFO  WARN  ERROR  FATAL
input_rhog_folder = "./"
species_tree_address = "species_tree.nwk"  # no space or special charcter in internal node.
# protein_format_qfo_dataset = True

# the following is not controlled by
# argument in nextflow run 
fragment_detection = True  # this can be changed through adding arguments to commond line run in nextflow run  --low-so-detection --fragment-detection
fragment_detection_msa_merge = True  # if this is false and fragment_detection_msa -> we'll remove both fragments in orthology analyss at parent level but report it in orthoxml Dubiousfragment
low_so_detection = True
# for fragment detection is better to subsampling_hogclass= False but make FastOMA slow
threshold_dubious_sd = 1/10
overlap_fragments = 0.15

mergHOG_ratioMax_thresh = 0.8
mergHOG_ratioMin_thresh = 0.9
mergHOG_shared_thresh = 50

# VARIABLE_threshold_dubious_sd
# threshold_dubious_sd = float(os.getenv(VARIABLE_threshold_dubious_sd, 0.1))

## output writing files
gene_trees_write = True
msa_write = True
gene_trees_write_all = True
msa_write_all = True
keep_subhog_each_pickle = True

big_rhog_size = 60 * 1000
omamer_family_threshold = 90
#
# omamer_fscore_treshold_big_rhog = 0.04 # 0.5 # means no thresold #0.2 #0.5  #  to have more proteins in the ortho groups 0.05  considering for big rhogs
# omamer_treshold_big_rhog_szie = 100 #9000 #100
#
# # for very big rhog, we need to be more stringent
# omamer_treshold_big_rhog_szie2 = 50*1000
# omamer_fscore_treshold_big_rhog2 = 0.6 #0.9

hogclass_max_num_seq = 40  # subsampling in msa # ver very 2
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 50
hogclass_tresh_ratio_gap_col = 0.6  # 0.8 for very very big
# old code after samplign if there are 2 seq sampled, then at least one nongap
subsampling_hogclass = True


num_msas_merge_mafft = 100
if hogclass_max_num_seq < 10:
    num_msas_merge_mafft= 10 # small value means always use simple mafft (not mafft merge)


automated_trimAL = False

label_SD_internal = "species_overlap"  # "reconcilation" "species_overlap"
# threshold_sd_suspicious_fragment_ratio = 1/3
tree_tool = "fasttree"  # "fasttree"  "iqtree"  # for  gene tree with two, we use
rooting_method = "midpoint"  # "midpoint" "mad"
rooting_mad_executable_path = "mad"  # /work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/installers/mad/

##inferhog
inferhog_tresh_ratio_gap_row =0.1 # 0.6   # to have more proteins in the ortho groups 0.1
inferhog_tresh_ratio_gap_col =0.5  # 0.6   # ver very 0.8
inferhog_min_cols_msa_to_filter = 50 #300 #50  # used for msa before gene tree inference and  saving msa in hog class

inferhog_filter_all_msas_row = True
inferhog_resume_rhog = True  # main.py False
# The intermediate files, internal node  pickle files is not working with nextflow
# the reason is that the pickles_subhog_folder_all is relative and stored in nextflow_work folder
# this folder can not be used for  the re-submitting
inferhog_resume_subhog = True  # read pickle_subhog  # _infer_subhog.py

# inferhog_concurrent_on = True now as an argument
inferhog_max_workers_num = 6

## xml
# write_all_prots_in_header = False  # if false writes only those in the hog group
inferhog_min_hog_size_xml = 2     # by setting this as 1, pyham won't work on xml output.

# batch_roothogs
big_rhog_filesize_thresh = 400 * 1000
sum_list_rhogs_filesize_thresh = 2 * 1e6


# big_rhog_filesize_thresh = 1.6 * 1000  # 600 would be better
# sum_list_rhogs_filesize_thresh = 5 * 1e3
orthoxml_v03 = True



def set_configs():
    parser = argparse.ArgumentParser(description="This is GETHOG3 ")     # parser.add_argument('--working-folder', help="in_folder")
    parser.add_argument('--logger-level', default="DEBUG")
    parser.add_argument("--version", action="version", help="Show version and exit.", version="0.0.6",)  # version=__version__
    parser.add_argument('--species-tree-address', default="species_tree_test.nwk")
    parser.add_argument('--input-rhog-folder')     # , default="./rhog"
    parser.add_argument('--parallel', action=argparse.BooleanOptionalAction)
    parser.add_argument('--fragment-detection', action=argparse.BooleanOptionalAction)
    parser.add_argument('--low-so-detection', action=argparse.BooleanOptionalAction)

    config_parser = parser.parse_args()
    # Namespace(logger_level=None, in_folder=None)
    setattr(sys.modules[__name__], 'logger_level', config_parser.logger_level)
    setattr(sys.modules[__name__], 'input_rhog_folder', config_parser.input_rhog_folder)
    setattr(sys.modules[__name__], 'parallel', config_parser.parallel)
    setattr(sys.modules[__name__], 'species_tree_address', config_parser.species_tree_address)
    setattr(sys.modules[__name__], 'fragment_detection', config_parser.fragment_detection)
    setattr(sys.modules[__name__], 'low_so_detection', config_parser.low_so_detection)
    print("config_parser 3 ", config_parser)



'''

TODO
 default doesnt work
'''
