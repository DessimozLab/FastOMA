# from xml.dom import minidom
# import xml.etree.ElementTree as ET
import _utils
import _inferhog
import os
# from os import listdir
# import os
import sys
# from os import listdir
# import _config

# from ._utils import logger_hog
# import ._utils_rhog


address_rhogs_folder = sys.argv[1]

inferhog_concurrent_on_string = sys.argv[2]   # "False"  # "False"  #
pickles_rhog_folder = "./"
pickles_subhog_folder_all = "./"

inferhog_concurrent_on = inferhog_concurrent_on_string == "True"

print("input is", address_rhogs_folder)

list_rhog_fastas_files = _utils.list_rhog_fastas(address_rhogs_folder)
print("there are ", len(list_rhog_fastas_files), "rhogs in the input folder")

rhogs_fa_folder = address_rhogs_folder


list_rhog_fastas_files_rem = _utils.list_rhog_fastas(address_rhogs_folder)
print("there are ", len(list_rhog_fastas_files_rem), "rhogs remained in the input folder", list_rhog_fastas_files_rem[:5] )

hogs_rhog_xml_batch = _inferhog.read_infer_xml_rhogs_batch(list_rhog_fastas_files_rem, inferhog_concurrent_on, pickles_rhog_folder, pickles_subhog_folder_all, rhogs_fa_folder)

print("finsihed ", address_rhogs_folder)


"""
TODO

print as debug
# [_config.oma_database_address, _config.working_folder_root , _config.species_tree_address , _config.working_id ,
 _config.protein_format_qfo_dataset, _config.working_folder, _config.omamer_fscore_treshold_big_rhog, _config.treshold_big_rhog_szie, _config.gene_trees_write, _config.keep_subhog_each_pickle, _config.hogclass_max_num_seq, _config.hogclass_min_cols_msa_to_filter, _config.hogclass_tresh_ratio_gap_col, _config.automated_trimAL, _config.lable_SD_internal , _config.rooting_method, _config.rooting_mad_executable_path , _config.inferhog_tresh_ratio_gap_row , _config.inferhog_tresh_ratio_gap_col , _config.inferhog_min_cols_msa_to_filter , _config.inferhog_filter_all_msas_row , _config.inferhog_resume_rhog  , _config.inferhog_resume_subhog , _config.inferhog_max_workers_num , _config.inferhog_min_hog_size_xml, _config.logger_level]

species_tree_address = working_folder_root + "tree_fastaname.nwk" # no space or special charcter in internal node, 
arise error or solve it

print("there shouldnt be any space in the tree name internal node name as well")
  '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf//pickles_subhog/rhog_833762/delta/epsilon subdivisions.pickle'

  FileExistsError: [Errno 17] File exists: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf//pickles_subhog/'

when concurrent has a issue it doesnt stop
Out[6]: {<Future at 0x7f1b48d9afa0 state=finished raised TypeError>: 'KORCO_'}
add eception to show , whn this happens for which taxnomic level and rhog

    if len(msa) <= 2:
        wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
        wrapper_tree.options.options['-fastest'].active = True
        
        we don't need tree for msa of 2 !

precuaitions
genetrees     is not with prefix
there shouldnt be any space in the tree name internal node name as well"


- internal node _0 ,  the root name disapears , beacause of ete3 behaviour 

"""


# # address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/working_nf/rhogs_rest/122/"  # bb1
# inferhog_concurrent_on_string = sys.argv[2]   # "False"  # "False"  #
# # mkdir pi_big_rhog pi_big_subhog pi_rest_rhog  pi_rest_subhog genetrees
# prefix_pickle = sys.argv[3]  # "pi_big"
# rhogs_fa_folder_pure = sys.argv[4]   # "rhogs_rest"  # or  "rhogs_big" # sys.argv[4] #
# pickles_rhog_folder = _config.working_folder + "/" + prefix_pickle + "_rhog/"
# pickles_subhog_folder_all = _config.working_folder + "/" + prefix_pickle + "_subhog/"
# pickles_rhog_folder  should be created in nextflow, or previous step of parralell infer subhog,
# workers may conflict of creating folders
# # pickles_subhog_folder = pickles_subhog_folder_all+"/rhog_" + str(rhogid_num) + "/"
# inferhog_concurrent_on = inferhog_concurrent_on_string == "True"  # sys.argv[2]
# # if not os.path.exists(pickles_rhog_folder):
# #     os.mkdir(pickles_rhog_folder)
# #rhogid_num = int(rhog_file.split("/")[-1].split(".")[0].split("_")[1][1:])
# #rhogid_batch = [rhogid_num]
# list_rhog_fastas_files = _utils.list_rhog_fastas(address_rhogs_folder)
# print("there are ", len(list_rhog_fastas_files), "rhogs in the input folder")
# if address_rhogs_folder.endswith("/"):
#     batch_folder=address_rhogs_folder.split("/")[-2]
# elif "/" in address_rhogs_folder:
#     batch_folder = address_rhogs_folder.split("/")[-1]
# else:
#     batch_folder = address_rhogs_folder
#
# print("rhogs in the input folder", batch_folder)
#
# rhogs_fa_folder = _config.working_folder + rhogs_fa_folder_pure + "/" + batch_folder + "/"
#
# list_done_rhogid = []
# if _config.inferhog_resume_rhog:
#     list_done_raw = listdir(pickles_rhog_folder)
#     for file in list_done_raw:
#        numr = int(file.split(".")[0].split("_")[1])
#        list_done_rhogid.append(numr)
# #list_done_rhogid = []
# list_rhog_fastas_files_rem = [i for i in list_rhog_fastas_files if i not in list_done_rhogid
#list_rhog_fastas_files_rem = [573914]  # [622140] #[605975] #[594043]

