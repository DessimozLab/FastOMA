from xml.dom import minidom
import xml.etree.ElementTree as ET
import _utils
import _inferhog
from _utils import logger_hog
import _utils_rhog

if __name__ == '__main__':
    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gene_trees_folder = working_folder + "/gene_trees_/"
    address_rhogs_folder = working_folder + "/rhogs_g10_s100/"  # "  old3/rhog_all/  /rhog_size_g2_s500/" sample_rootHOG
    species_tree_address = working_folder + "/archive/lineage_tree_qfo.phyloxml"
    pickle_folder = working_folder + "/pickle_folder/"
    # format_prot_name = 1  # 0:bird(TYTALB_R04643)  1:qfo(tr|E3JPS4|E3JPS4_PUCGT)
    file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)

    step = "hog"
    print("we are here line16")
    if step == "rhog":
        """
        Structure of folders:
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address= the address to the oma databases
        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is a HOG.
        """
        import pyoma.browser.db as db
        oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
        # working_folder+"omamer_database/oma_path/OmaServer.h5"
        print("rHOG inferece has started. The oma database address is in ", oma_database_address)
        (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
        (query_species_names, query_prot_records_species) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
        query_prot_records_species = _utils_rhog.add_species_name_gene_id(query_prot_records_species, query_species_names, working_folder)
        hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)

        (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
        prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

        query_prot_records_species_filtered = _utils_rhog.filter_prot_mapped(query_species_names,
                                                                              query_prot_records_species,
                                                                              query_prot_names_species_mapped)

        print(len(query_prot_records_species_filtered), len(query_prot_records_species_filtered[0]))
        (rhogid_num_list, rhogids_prot_records_query) = _utils_rhog.group_prots_roothogs(prots_hogmap_hogid_allspecies,
                                                                                         address_rhogs_folder,
                                                                                         query_species_names,
                                                                                         query_prot_records_species_filtered)
        # step = "hog"

    if step == "hog":
        print("we are here line43")
        rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
        logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list)) + ".")

        rhogid_num_list = rhogid_num_list[:15]
        dask_level = 0  # 1:one level (rhog), 3:both levels (rhog+taxonomic)

        print(rhogid_num_list)
        number_roothog = len(rhogid_num_list)
        num_per_parralel = 4
        parralel_num = int(number_roothog/num_per_parralel)
        if number_roothog != parralel_num*num_per_parralel: parralel_num += 1
        rhogid_batch_list = []
        for list_idx in range(parralel_num):
            if list_idx == parralel_num:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
            else:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
            rhogid_batch_list.append(rhogid_num_list_portion)

        if dask_level % 2 == 1:
            from _dask_env import client_dask
            dask_out_list = []
        else:
            hogs_a_rhog_xml_all_list =[]

        for rhogid_batch_idx in range(len(rhogid_batch_list)):
            rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
            logger_hog.info("\n *==* \nNumber of working root hog in the batchid:"+str(rhogid_batch_idx)+" is " + str(len(rhogid_batch)) + ".")

            if dask_level % 2 == 1:
                # len_tresh = 1000
                # for rhogid_num_i in range(len(rhogid_num_list_input)):
                #    rhogid_num = rhogid_num_list_input[rhogid_num_i]
                #    rhogid_len = rhogid_len_list[rhogid_num_i]
                #    if rhogid_len < len_tresh:
                # vars_input_future = client_dask.scatter(vars_input)
                dask_out = client_dask.submit(_inferhog.read_infer_xml_rhogs, rhogid_batch, file_folders, dask_level)
                dask_out_list.append(dask_out)
            else:
                hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhogs(rhogid_batch, file_folders, dask_level)
                hogs_a_rhog_xml_all_list.append(hogs_a_rhog_xml_all)
                print(hogs_a_rhog_xml_all)

        if dask_level % 2 == 1:
            for dask_out in dask_out_list:
                hogs_a_rhog_xml_all = dask_out.result()
                print(hogs_a_rhog_xml_all)
            print("dask out gathered")



    output_xml_name = "out12b.xml"
    _inferhog.collect_write_xml(working_folder, pickle_folder, output_xml_name)

    print("main py is finished.")



"""
to do :
    - get rid of gene_id_name, write in the rhog file 
    - dobule check function merge_subhogs
    
"""
