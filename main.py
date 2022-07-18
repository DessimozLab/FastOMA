


"""

from datetime import datetime
import os

from random import sample

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
"""


import utils
from utils import logger_hog


if __name__ == '__main__':

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"

    address_rhogs_folder = working_folder + "/sample_rootHOG/" #"/rhog_size_g2_s500/"
    species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"




    step = "hog"
    print("we are here ")
    if step == "roothog":
        pass
    if step == "hog":
        rhogid_num_list = utils.list_rhog_fastas(address_rhogs_folder)

        logger_hog.info("Number of root hog is "+str(len(rhogid_num_list))+".")
        print(rhogid_num_list[:2])
        rhogid_num_list_temp = rhogid_num_list[:299]

        (species_tree) = utils.read_species_tree(species_tree_address)



        print("last")