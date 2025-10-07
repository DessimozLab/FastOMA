
import sys
import logging
import pyham


"""
To run
python orthoxml2rootHOGidLevel.py fastoma_out
fastoma_out: where species_tree_checked.nwk and FastOMA_HOGs.orthoxml are located.

output: 
a two-column tsv file, HOGid taxid

"""


fastoma_out=sys.argv[1] 

treeFile=fastoma_out+'/species_tree_checked.nwk'
orthoxmlFile=fastoma_out+'/FastOMA_HOGs.orthoxml'


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger("hog")
logger.setLevel(logging.INFO)

ham_analysis = pyham.Ham(treeFile, orthoxmlFile, use_internal_name=True)
rootHOG_dic=ham_analysis.get_dict_top_level_hogs()
HOGid_level_list=[]
for HOGid, HOG in rootHOG_dic.items():
    level=HOG.__dict__['_properties']['TaxRange']
    HOGid_level_list.append((HOGid,level))
print('We found ' +str(len(HOGid_level_list))+ ' rootHOGs.')

output_f = open('FastOMA_HOGs.level.tsv','w')
output_f.write('##rootHOGid'+'\t'+'TaxRange'+'\n')
for HOGid,level  in HOGid_level_list :
    output_f.write(HOGid+'\t'+level+'\n')
output_f.close()

print('RootHOGs level info is written in FastOMA_HOGs.level.tsv')
