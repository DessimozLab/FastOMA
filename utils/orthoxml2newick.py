




# copy to zoo folder
#


import sys
sys.path.insert(1, '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/software/gits/zoo/zoo/hog/')

from convert import orthoxml_to_newick


folder= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run28jan/working_nf/"
file= folder+"hog_euk_28jan_all.orthoxml"
trees = orthoxml_to_newick(file)
print(len(trees))

# write them as files
