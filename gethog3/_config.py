


oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/omafast/archive/OmaServer.h5"
working_folder_root = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/"  # bird_hog

species_tree_address = working_folder_root + "lineage_tree_qfo_.nwk" #phyloxml"
working_id = "hog3_nov25/"
working_folder = working_folder_root + working_id
omamer_fscore_treshold_big_rhog = 0.5  # 0.2
treshold_big_rhog_szie = 3000

gene_trees_write = True  # this also goes for writing msas
keep_subhog_each_pickle = True #False
protein_format_qfo_dataset = True





###  dask configs
dask_level = 0   # 1:one level (rhog), 2:both levels (rhog+taxonomic)  3:only taxonomic level  0: no dask
dask_n_core = 1
dask_machine = "slurm"  # "local"  "slurm"
dask_memory_slurm = "40GB"
dask_time_slurm = "05:00:00"
dask_n_jobs = 180

## hogclass configs
hogclass_max_num_seq = 5  # subsampling in msa
hogclass_min_cols_msa_to_filter = hogclass_max_num_seq * 500
hogclass_tresh_ratio_gap_col = 0.2

##inferhog
inferhog_dask_2nd_rhogsize = 200
inferhog_tresh_ratio_gap_row = 0.4
inferhog_tresh_ratio_gap_col = 0.2
inferhog_min_cols_msa_to_filter = 3000  # used for msa before gene tree inference and  saving msa in hog class


## xml
write_all_prots_in_header = False  # if false writes only those in the hog group
inferhog_min_hog_size_xml = 2



"""
hardcoded folder names
omamer_search in folder working_folder_root

gene_id_pickle_file = working_folder + "gene_id_dic_xml.pickle"
address_rhogs_folder_raw = working_folder + "rhogs_raw/"

pickles_subhog_folder_all = _config.working_folder + "/pickles_subhog/"
pickles_rhog_folder = = _config.working_folder + "/pickles_rhog/"

gene_trees_folder =  working_folder+"genetrees"

#file_folders = (address_rhogs_folder_filt, gene_trees_folder, pickle_folder, species_tree_address)
# gene_id_pickle_file = working_folder + "gene_id_dic_xml.pickle"
    name = str(omamer_fscore_treshold_big_rhog)+"_"+str(treshold_big_rhog_szie    
    address_rhogs_folder_filt = working_folder + "rhogs_" + name + "/"
    pickle_folder = working_folder + "pickle_b_"+name+"/"
    hogs_children_level_pickle_folder = working_folder+"/pickle_hog_children3/"
    print(pickle_folder)
    
    output_xml_name = "out_xml_b_"+name+"_.xml"


"""
