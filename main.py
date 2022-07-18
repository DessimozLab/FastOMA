

import logging
from datetime import datetime
import os
from os import listdir
from random import sample

from ete3 import Phyloxml
from ete3 import Tree

import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import xml.etree.ElementTree as ET
import itertools
from xml.dom import minidom
#import concurrent.futures

import dask
from dask.distributed import LocalCluster
from dask.distributed import Client
from dask_jobqueue import SLURMCluster

import gc



def list_rhog_fastas(address_rhogs_folder):
    """
     create a list of rootHOG IDs  stored in the folder of rHOG .
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




if __name__ == '__main__':

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    step = "hog"

    logging.basicConfig()
    logger_hog = logging.getLogger("hog")
    logger_hog.setLevel(logging.INFO)  # WARN

    address_rhogs_folder = working_folder + "/sample_rootHOG/" #"/rhog_size_g2_s500/"
    species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"


    print("we are here ")
    if step == "roothog":
        pass
    if step == "hog":
        rhogid_num_list = list_rhog_fastas(address_rhogs_folder)

        print(len(rhogid_num_list),rhogid_num_list[:2])


