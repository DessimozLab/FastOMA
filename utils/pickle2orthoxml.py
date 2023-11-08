

from xml.dom import minidom
import xml.etree.ElementTree as ET
import pickle

import sys

input_pickle= sys.argv[1]   # "file_D0680685.pickle"
handle=open(input_pickle,'rb')
orthoxml_file = pickle.load(handle)

print(len(orthoxml_file))
xml_str = minidom.parseString(ET.tostring(orthoxml_file[0])).toprettyxml(indent="   ")

with open(input_pickle+"_noheader.orthoxml","w") as out_file:
    out_file.write(xml_str)



