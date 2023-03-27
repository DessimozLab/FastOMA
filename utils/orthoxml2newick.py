



from FastOMA.zoo.hog.convert import orthoxml_to_newick

input_orthoxml = "folder/my_hogs.orthoxml"
output_folder = "output_folder"


trees = orthoxml_to_newick(input_orthoxml)

print(len(trees))


# write them as files

for tree_idx, tree_i in enumerate(trees):

    tree_file_i = output_folder+"tree_"+str(tree_idx)+".nwk"
    tree_i.write(format=1, outfile=tree_file_i)



