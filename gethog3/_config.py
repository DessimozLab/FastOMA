
# global path
import argparse
import sys

working_folder = "./working_folder"+ "/"

species_tree_address = working_folder + "species_tree.nwk"
                       # no space or special charcter in internal node.
protein_format_qfo_dataset = True

## output writing files
gene_trees_write = False  # this also goes for writing msas
keep_subhog_each_pickle = False

# filtering omamer
omamer_fscore_treshold_big_rhog = 0.5  # 0.2
treshold_big_rhog_szie = 10

## hogclass configs
hogclass_max_num_seq = 5  # subsampling in msa
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 500
hogclass_tresh_ratio_gap_col = 0.2

automated_trimAL = False
lable_SD_internal = "species_overlap"   #  "reconcilation" "species_overlap"
tree_tool = "fasttree" #  "fasttree"  "iqtree"  # for  gene tree with two, we use

rooting_method = "midpoint"  # "midpoint" "mad"
rooting_mad_executable_path = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/installers/mad/mad"

##inferhog
inferhog_tresh_ratio_gap_row = 0.4
inferhog_tresh_ratio_gap_col = 0.4
inferhog_min_cols_msa_to_filter = 100  # used for msa before gene tree inference and  saving msa in hog class

inferhog_filter_all_msas_row = True

inferhog_resume_rhog  = True      # main.py False
inferhog_resume_subhog = True      # read pickle_subhog  # _inferhog.py

# inferhog_concurrent_on = True now as an argument
inferhog_max_workers_num = 8

## xml
# write_all_prots_in_header = False  # if false writes only those in the hog group
inferhog_min_hog_size_xml = 2  # by setting this as 1, pyham won't work on xml output.

logger_level = "DEBUG"  # DEBUG INFO


def set_configs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--working-folder', help="working_folder")
    parser.add_argument('--logger-level')
    config = parser.parse_args()
    print(config)
    setattr(sys.modules[__name__], 'logger_level', config.logger_level)
    setattr(sys.modules[__name__], 'working_folder', config.working_folder)

