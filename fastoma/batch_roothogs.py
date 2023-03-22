

# from Bio import SeqIO
import os
from os import makedirs
import shutil
from os import listdir
from . import _config


def list_rhog_fastas(address_rhogs_folder):
    """
     create orthoxml_to_newick.py list of rootHOG IDs  stored in the folder of rHOG .
     input: folder address
     output: list of rhog Id (integer)
    """
    rhog_files = listdir(address_rhogs_folder)
    rhogid_num_list= []
    for rhog_file in rhog_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)

    return rhogid_num_list



def folder_1h_rhog(address_rhogs_folder, output_folder_big, output_folder_rest):

    # in_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/working_nfp/"

    #work_fldor_out = in_folder

    rhogid_num_list = list_rhog_fastas(address_rhogs_folder)
    len(rhogid_num_list)
    dic_rhognum_size = {}
    for rhogid_num in rhogid_num_list[2:]:
        rhog_i_prot_address = address_rhogs_folder + "/HOG_" + str(rhogid_num).zfill(7) + ".fa"
        # rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
        rhog_i_size = os.path.getsize(rhog_i_prot_address)
        dic_rhognum_size[rhogid_num] = rhog_i_size
    len(dic_rhognum_size)


    #list_list_rest_rhog = [[]]
    list_list_rest_size = [[]]
    #list_list_big = []

    # to have at least one big and one rest rhog, this makes nextflow pipline easier,
    # otherwise we need if condition  collect_subhogs won't satisfy if one of them is empty

    list_list_rest_rhog = [[rhogid_num_list[0]],[]]  # each insid list should take 1h to
    list_list_rest_size = [[],[]]
    list_list_big = [rhogid_num_list[1]]



    for rhognum, size in dic_rhognum_size.items():
        # print(rhognum, size)
        if size > _config.big_rhog_filesize_thresh:
            list_list_big.append(rhognum)
        else:
            if sum(list_list_rest_size[-1]) < _config.sum_list_rhogs_filesize_thresh:
                list_list_rest_rhog[-1].append(rhognum)
                list_list_rest_size[-1].append(size)
            else:
                list_list_rest_rhog.append([rhognum])
                list_list_rest_size.append([size])

    # makedirs(work_fldor_out+"rhogs_rest")
    for folder_id, list_rhog in enumerate(list_list_rest_rhog):
        # print(folder_id)
        # output_folder_rest = work_fldor_out + "/rhogs_rest/"
        makedirs( output_folder_rest + str(folder_id))
        for rhogid_num in list_rhog:
            name = "HOG_" + str(rhogid_num).zfill(7) + ".fa"
            folder_rest =output_folder_rest + str(folder_id) + "/"
            shutil.copy(address_rhogs_folder + name, folder_rest + name)

    # makedirs(work_fldor_out+"rhogs_big")
    #  output_folder_big = work_fldor_out + "rhogs_big/"
    for folder_id, rhogid_num in enumerate(list_list_big):
        name = "HOG_" + str(rhogid_num).zfill(7) + ".fa"
        folder_big = output_folder_big + "b" + str(folder_id) + "/"
        makedirs(folder_big)
        shutil.copy(address_rhogs_folder + name, folder_big + name)

    return 1



def batch_roothogs():

    input_rhog = "./rhogs_all/"
    output_folder_big = "./rhogs_big/"
    output_folder_rest = "./rhogs_rest/"
    folder_1h_rhog(input_rhog, output_folder_big, output_folder_rest)

