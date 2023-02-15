

working_folder_root = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/testgethog/" # "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/"  # bird_hog




#oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/archive/OmaServer.h5"
# oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/archive/OmaServer.h5"

# bird
# working_folder_root = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/"  # bird_hog
# species_tree_address = working_folder_root + "birds370_iqtree_treefile_95bootstrap_internal_name_6let_16Nov_.nwk"
# working_id = "hog3_nov25f/"
# protein_format_qfo_dataset = False

# qfo




#species_tree_address = working_folder_root + "lineage_tree_qfo_2.phyloxml" #phyloxml"
species_tree_address = working_folder_root + "lineage_tree_qfo_2.phyloxml"
                       # "tree_fastaname.nwk"  # no space or special charcter in internal node,

working_id =  "./"  # "working_nf/"
protein_format_qfo_dataset = True

working_folder = working_folder_root + working_id


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


"""
hardcoded folder names
omamer_search in folder working_folder_root
gene_id_pickle_file = working_folder + "gene_id_dic_xml.pickle"
address_rhogs_folder_raw = working_folder + "rhogs_raw/"
pickles_subhog_folder_all = in the nextflow file
gene_trees_folder =  working_folder+"genetrees"
"""
