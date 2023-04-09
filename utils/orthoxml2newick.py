


import sys
import os
from FastOMA.zoo.hog.convert import orthoxml_to_newick

"""
how to run
    python orthoxml2newick.py my_hogs.orthoxml
"""

input_orthoxml = sys.argv[1]
output_folder = "output_folder_trees"

os.mkdir(output_folder)

trees = orthoxml_to_newick(input_orthoxml)

print("We extracted "+str(len(trees))+" trees from the input HOG orthoxml"+input_orthoxml)

# write them as files
for treeid_hog, tree in trees.items():
    tree_file_i = output_folder+"/tree_"+str(treeid_hog)+".nwk"
    with open(tree_file_i,'w') as handle:
        handle.write(tree)
    handle.close()
    # tree_i.write(format=1, format_root_node=True, outfile=tree_file_i)
print("We wrote "+str(len(trees))+" trees  in nhx format from the input HOG orthoxml"+input_orthoxml+"in "+output_folder)
print("You can visualise each tree using https://beta.phylo.io/viewer/ as extendeed newick format.")
