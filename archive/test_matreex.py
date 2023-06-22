from matreex import matreex
import json
import ete3
def make_matreex_html(hog_ids, ham_object, species_tree, outpath, ancestral_level ="Euteleostomi" ):
    sp_tree = ete3.Tree(species_tree, format=1,quoted_node_names=True)
    for n in sp_tree.traverse():
            n.add_features(
                description='', color='',
                mcol='')
    tmp_gene_trees = []
    for hog_id in hog_ids:
        root_hog = ham_object.get_hog_by_id(hog_id)
        merged_gt = matreex.ham2gt(root_hog,hog_id)
        tmp_gene_trees += matreex.get_hog_gene_trees(ancestral_level, sp_tree, merged_gt, hog_id)
    for gt in tmp_gene_trees:
            matreex._add_lost_subtrees(gt, sp_tree)
    gts = matreex.merge_gene_trees(tmp_gene_trees, sp_tree)
    gt_json = matreex.gt2json(gts)
    st_json = matreex.st2json(sp_tree & gts.taxon)

    matreex.write_html(json.dumps(gt_json), json.dumps(st_json),outpath )

    return 1




import pyham        # pyham package for ham analysis
import pandas as pd # pandas for dataframes

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

pickle_folder= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run4_May_omamer/salamin/"
idx= 657554
pickle_files_adress= ["file_"+str(idx)+".pickle"]


output_xml_name= pickle_folder + pickle_files_adress[0]+"_.orthoxml"

print(output_xml_name)
orthoxml_path = output_xml_name
nwk_path= pickle_folder + "../in_folder/species_tree.nwk"
tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
tree_str[:20]



ham_analysis = pyham.Ham(tree_str, output_xml_name, use_internal_name=True)

list_hog_top = list(ham_analysis.top_level_hogs.keys())
list_hog_top

#root_hog_id = "HOG_0013607_sub10165" # "HOG_0024531_sub10729"
root_hog_id= list_hog_top[0]
print(root_hog_id)
# select an HOG
hog =  ham_analysis.get_hog_by_id(root_hog_id)
print(hog)
# create the iHam for it and save it into an html file

output_name = output_xml_name+"_00.html"# "iham_HOG{}.html".format(hog.hog_id)
ham_analysis.create_iHam(hog=hog, outfile=output_name)
output_name


output_name


ham_analysis = pyham.Ham(tree_str, output_xml_name, use_internal_name=True)

list_hog_top = list(ham_analysis.top_level_hogs.keys())
list_hog_top[:10]


hog_ids = list_hog_top  # ['HOG_0613829_sub16986','HOG_0613832_sub19203', 'HOG_0613832_sub19204', 'HOG_0613832_sub19150']
outpath=pickle_folder+ "/matreex_output/"+hog_ids[0]+"_and_rest_matreex.html"


make_matreex_html(hog_ids, ham_analysis, nwk_path, outpath, ancestral_level ="Euteleostomi" )

