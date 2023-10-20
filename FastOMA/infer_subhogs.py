# from xml.dom import minidom
# import xml.etree.ElementTree as ET
from . import _utils_subhog
from . import _infer_subhog
from . import _config
# import os
# from os import listdir
# import os
import sys
# from os import listdir
# import _config

from ._config import logger_hog
import os

# from ._utils import logger_hog
# import ._utils_rhog

def infer_subhogs():
    _config.set_configs()
    logger_hog.debug("logger_level is "+str(_config.logger_level))
    logger_hog.debug("input_rhog_folder is " + str(_config.input_rhog_folder))
    logger_hog.debug("parallel is " + str(_config.parallel))
    logger_hog.debug("species_tree_address is " + str(_config.species_tree_address))
    logger_hog.debug("fragment_detection is " + str(_config.fragment_detection))
    logger_hog.debug("low_so_detection is " + str(_config.low_so_detection))
    logger_hog.debug("inferhog_max_workers_num is " + str(_config.inferhog_max_workers_num))
    logger_hog.debug("inferhog_tresh_ratio_gap_row is " + str(_config.inferhog_tresh_ratio_gap_row))
    logger_hog.debug("inferhog_tresh_ratio_gap_col is " + str(_config.inferhog_tresh_ratio_gap_col))
    logger_hog.debug("inferhog_min_cols_msa_to_filter is " + str(_config.inferhog_min_cols_msa_to_filter))
    logger_hog.debug("big_rhog_size is " + str(_config.big_rhog_size))
    logger_hog.debug("omamer_family_threshold is " + str(_config.omamer_family_threshold))
    logger_hog.debug("hogclass_max_num_seq is " + str(_config.hogclass_max_num_seq))
    logger_hog.debug("hogclass_min_cols_msa_to_filter is " + str(_config.hogclass_min_cols_msa_to_filter))
    logger_hog.debug("inferhog_resume_rhog is " + str(_config.inferhog_resume_rhog))
    logger_hog.debug("inferhog_resume_subhog is " + str(_config.inferhog_resume_subhog))
    logger_hog.debug("fragment_detection is " + str(_config.fragment_detection))
    logger_hog.debug("fragment_detection_msa_merge is " + str(_config.fragment_detection_msa_merge))
    logger_hog.debug("low_so_detection is " + str(_config.low_so_detection))
    logger_hog.debug("threshold_dubious_sd is " + str(_config.threshold_dubious_sd))
    logger_hog.debug("overlap_fragments is " + str(_config.overlap_fragments))
    logger_hog.debug("gene_trees_write is " + str(_config.gene_trees_write))
    logger_hog.debug("msa_write is " + str(_config.msa_write))
    logger_hog.debug("gene_trees_write_all is " + str(_config.gene_trees_write_all))
    logger_hog.debug("msa_write_all is " + str(_config.msa_write_all))
    logger_hog.debug("keep_subhog_each_pickle is " + str(_config.keep_subhog_each_pickle))


    address_rhogs_folder = _config.input_rhog_folder
    inferhog_concurrent_on = False
    inferhog_concurrent_on = _config.parallel #  sys.argv[2]   # "False"  # "False"  #

    if inferhog_concurrent_on:
        print("parallelization for subhog inference is on.")

    pickles_rhog_folder = "./" # pickles_temp/ pickle_rhogs
    if not os.path.exists(pickles_rhog_folder):
        os.makedirs(pickles_rhog_folder)

    # pickles_rhog_folder = "./"
    pickles_subhog_folder_all = "./"

    # inferhog_concurrent_on = inferhog_concurrent_on_string == "True"

    print("input is", address_rhogs_folder)

    list_rhog_fastas_files = _utils_subhog.list_rhog_fastas(address_rhogs_folder)
    print("there are ", len(list_rhog_fastas_files), "rhogs in the input folder")

    rhogs_fa_folder = address_rhogs_folder

    list_rhog_fastas_files_rem = _utils_subhog.list_rhog_fastas(address_rhogs_folder)
    print("there are ", len(list_rhog_fastas_files_rem), "rhogs remained in the input folder", list_rhog_fastas_files_rem[:5] )

    hogs_rhog_xml_batch = _infer_subhog.read_infer_xml_rhogs_batch(list_rhog_fastas_files_rem, inferhog_concurrent_on, pickles_rhog_folder, pickles_subhog_folder_all, rhogs_fa_folder)

    print("finsihed ", address_rhogs_folder)

# TODO all below

"""
  
- use same seed for hog sampling and fasttree/mafft if they have 
- add python code for validate the input make sure the file is there and decent
- check input tree 
_utils_subhog line 92 add check for spaces or chars in 

- print as debug all the variables 
# [_config.oma_database_address, _config.working_folder_root , _config.species_tree_address , _config.working_id ,
 _config.protein_format_qfo_dataset, _config.in_folder, _config.omamer_fscore_treshold_big_rhog,
  _config.treshold_big_rhog_szie, _config.gene_trees_write, _config.keep_subhog_each_pickle, 
  _config.hogclass_max_num_seq, _config.hogclass_min_cols_msa_to_filter, _config.hogclass_tresh_ratio_gap_col, _config.automated_trimAL, _config.lable_SD_internal ,
   _config.rooting_method, _config.rooting_mad_executable_path , _config.inferhog_tresh_ratio_gap_row , _config.inferhog_tresh_ratio_gap_col , _config.inferhog_min_cols_msa_to_filter
    , _config.inferhog_filter_all_msas_row , _config.inferhog_resume_rhog  , _config.inferhog_resume_subhog , _config.inferhog_max_workers_num , _config.inferhog_min_hog_size_xml, _config.logger_level]

FileExistsError: [Errno 17] File exists: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf//pickles_subhog/'

when concurrent has a issue it doesnt stop
Out[6]: {<Future at 0x7f1b48d9afa0 state=finished raised TypeError>: 'KORCO_'}
add eception to show , whn this happens for which taxnomic level and rhog

proteins were associated with known HOGs, and warn if less than ~80 % of proteins are associated with known HOGs.

double check 
     prots_to_remove
     it seems that it always an input argument! 
     

"""


# if address_rhogs_folder.endswith("/"):
#     batch_folder=address_rhogs_folder.split("/")[-2]
# elif "/" in address_rhogs_folder:
#     batch_folder = address_rhogs_folder.split("/")[-1]
# else:
#     batch_folder = address_rhogs_folder
#
# print("rhogs in the input folder", batch_folder)
#
# rhogs_fa_folder = _config.in_folder + rhogs_fa_folder_pure + "/" + batch_folder + "/"
# list_done_rhogid = []
# if _config.inferhog_resume_rhog:
#     list_done_raw = listdir(pickles_rhog_folder)
#     for file in list_done_raw:
#        numr = int(file.split(".")[0].split("_")[1])
#        list_done_rhogid.append(numr)
# #list_done_rhogid = []
# list_rhog_fastas_files_rem = [i for i in list_rhog_fastas_files if i not in list_done_rhogid
#list_rhog_fastas_files_rem = [573914]  # [622140] #[605975] #[594043]

