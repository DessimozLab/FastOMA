# from xml.dom import minidom
# import xml.etree.ElementTree as ET
import _utils
import _inferhog

# from _utils import logger_hog
# import _utils_rhog
#from os import listdir
#import os

import sys





address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/rhogs/1/"  # sys.argv[1]
inferhog_concurrent_on = False  # sys.argv[2]

print("input is", address_rhogs_folder)

#rhogid_num = int(rhog_file.split("/")[-1].split(".")[0].split("_")[1][1:])
#rhogid_batch = [rhogid_num]


list_rhog_fastas_files = _utils.list_rhog_fastas(address_rhogs_folder)
print("there are ",len(list_rhog_fastas_files),"rhogs in the input folder")
folder=address_rhogs_folder # .split("/")[-2]
print("rhogs in the input folder",folder)

from os import listdir
#pickle_folder ="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/gethog3_nf_/pickles_rhog"
#list_done_raw = listdir(pickle_folder)
list_done_rhogid = []
#for file in list_done_raw:
#    numr = int(file.split(".")[0].split("_")[1])
#    list_done_rhogid.append(numr)



list_rhog_fastas_files_rem = [i for i in list_rhog_fastas_files if i not in list_done_rhogid]

print("there are ",len(list_rhog_fastas_files_rem),"rhogs remained in the input folder")

hogs_rhog_xml_batch = _inferhog.read_infer_xml_rhogs_batch(list_rhog_fastas_files_rem,inferhog_concurrent_on, folder)


print("finsihed ",address_rhogs_folder)



# to do   FileNotFoundError: [Errno 2] No such file or directory: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3//gethog3_nf_//pickles_rhog//file_69696.pickle'



