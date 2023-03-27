

from _utils import logger_hog
import _utils_rhog
import _utils
# from distributed import get_client
# from dask.distributed import rejoin, secede
from Bio import SeqIO


if __name__ == '__main__':
    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gene_trees_folder = "" # in_folder + "/gene_trees_/"
    # check gene_trees_folder exist otherwise mkdir this
    pickle_folder = ""
    species_tree_address = ""


    address_rhogs_folder = working_folder + "/rhog_g10k_fil/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
    address_rhogs_folder_fil = working_folder + "/rhog_g10k_fil3/"
    file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)


    oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/FastOMA/archive/OmaServer.h5"
    # in_folder+"omamer_database/oma_path/OmaServer.h5"
    print("rHOG inferece has started. The oma database address is in ", oma_database_address)
    (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
    (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
    hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)

    (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_fscore_allspecies,
    prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

    print("start finding filtered proteins")
    prot_name_filt = [ ]
    for species_i, query_prot_names in enumerate(query_prot_names_species_mapped):
        print(species_i)
        prots_hogmap_seqlen = prots_hogmap_seqlen_allspecies[species_i]
        prots_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen_allspecies[species_i]
        prots_hogmap_fscore = prots_hogmap_fscore_allspecies[species_i]
        for prot_i, prot_name in enumerate(query_prot_names):
            prot_hogmap_seqlen = prots_hogmap_seqlen[prot_i]
            prot_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen[prot_i]
            prot_hogmap_fscore = prots_hogmap_fscore[prot_i]
            if prot_hogmap_seqlen != "na" and prot_hogmap_subfmedseqlen != "na":
                prot_hogmap_subfmedseqlen = int(prot_hogmap_subfmedseqlen)
                prot_hogmap_seqlen = int(prot_hogmap_seqlen)
                prot_hogmap_fscore= float(prot_hogmap_fscore)
                if prot_hogmap_subfmedseqlen * 0.8 < prot_hogmap_seqlen < prot_hogmap_subfmedseqlen * 1.8 and prot_hogmap_fscore > 0.8:
                    prot_name_filt.append(prot_name)
    print("finished finding filtered proteins")
    rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)

    for rhogid_num in [605945]:# rhogid_num_list
        print("working on rhog", rhogid_num )
        prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        rhog_i_filtered = []
        for rhog in rhog_i:
            if rhog.name.split("||")[0] in prot_name_filt:
                rhog_i_filtered.append(rhog)
        print("writing rhog", rhogid_num)
        SeqIO.write(rhog_i_filtered, address_rhogs_folder_fil + "HOG_B" + str(rhogid_num).zfill(7) + ".fa", "fasta")
    print("finished writing filtered rhogs")


