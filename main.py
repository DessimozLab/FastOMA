from xml.dom import minidom
import xml.etree.ElementTree as ET
import _utils
import _inferhog
from _utils import logger_hog

if __name__ == '__main__':
    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gene_trees_folder = working_folder + "/gene_trees_test/"
    address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # "  old3/rhog_all/  /rhog_size_g2_s500/" sample_rootHOG
    species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"
    pickle_address = working_folder + "/pickle_folder/"
    format_prot_name = 1  # 0:bird(TYTALB_R04643)  1:qfo(tr|E3JPS4|E3JPS4_PUCGT)

    step = "hog"
    print("we are here line16")
    if step == "roothog":
        """
        Structure of folders:
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address= the address to the oma databases
        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is a HOG.
        """
        # import pyoma.browser.db as db
        # oma_database_address = address_working_folder+"omamer_database/oma_path/OmaServer.h5"
        # print("program has started. The oma database address is in ",oma_database_address)
        # (oma_db, list_oma_species) = parse_oma_db(oma_database_address)
        # (query_species_names, query_prot_records_species) = parse_proteome(list_oma_species)
        # query_prot_records_species = add_species_name(query_prot_records_species,query_species_names)
        # hogmap_allspecies_elements = parse_hogmap_omamer(query_species_names)
        # (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
        # prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
        # query_prot_records_species_filtered =  filter_prot_mapped(query_species_names, query_prot_records_species,
        # query_prot_names_species_mapped)
        # print(len(query_prot_records_species_filtered),len(query_prot_records_species_filtered[0]))
        # (rhogid_num_list, rhogids_prot_records_query) = group_prots_roothogs(prots_hogmap_hogid_allspecies,
        # address_rhogs_folder)
        # step == "hog"

    if step == "hog":
        print("we are here line43")
        rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
        logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list)) + ".")

        rhogid_num_list = rhogid_num_list[:98]
        dask_future = True
        dask_future_taxon = False

        print(rhogid_num_list)
        number_roothog = len(rhogid_num_list)
        num_per_parralel = 10
        parralel_num = int(number_roothog/num_per_parralel)
        if number_roothog != parralel_num*num_per_parralel: parralel_num += 1
        rhogid_batch_list = []
        for list_idx in range(parralel_num):
            if list_idx == parralel_num:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
            else:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
            rhogid_batch_list.append(rhogid_num_list_portion)

        if dask_future:
            from _dask_env import client_dask

        dask_out_list = []
        for rhogid_batch_idx in range(len(rhogid_batch_list)):
            rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
            logger_hog.info("\n *==* \nNumber of working root hog in the batchid:"+str(rhogid_batch_idx)+" is " + str(len(rhogid_batch)) + ".")
            gene_id_name = _utils.gene_num_convertor(rhogid_batch, address_rhogs_folder, format_prot_name, rhogid_batch_idx)
            print("Number of genes in the batch is ", len(gene_id_name))
            vars_input = (gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder,
                          pickle_address, dask_future, dask_future_taxon, format_prot_name)

            if dask_future:
                # len_tresh = 1000
                # for rhogid_num_i in range(len(rhogid_num_list_input)):
                #    rhogid_num = rhogid_num_list_input[rhogid_num_i]
                #    rhogid_len = rhogid_len_list[rhogid_num_i]
                #    if rhogid_len < len_tresh:
                # vars_input_future = client_dask.scatter(vars_input)
                dask_out = client_dask.submit(_inferhog.read_infer_xml_rhogs, rhogid_batch, vars_input)
                dask_out_list.append(dask_out)
            else:
                hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhogs(rhogid_batch, vars_input)
                print(hogs_a_rhog_xml_all)

        if dask_future:
            for dask_out in dask_out_list:
                hogs_a_rhog_xml_all = dask_out.result()
                print(hogs_a_rhog_xml_all)
            print("dask out gathered")

    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA", "originVersion": "Nov 2021", "version": "0.3"})  #
    groups_xml = ET.SubElement(orthoxml_file, "groups")
    for hog_xml in hogs_a_rhog_xml_all:
        groups_xml.append(hog_xml)
    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    print(xml_str)
    print("main py is finished.")



"""
to do :
    - get rid of gene_id_name, write in the rhog file 
    - dobule check function merge_subhogs
    
"""
