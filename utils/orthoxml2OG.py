
# orthologous per group
# for speciefic taxonomic level,  top level
# option whether include its substes,  groups that emerge after
# for root combiaiton is not alowed
# (1) folder for fasta file
# (2) text file, per group, gene id
#  (3) extract gene markers, single copy ortho roots ,  for species tree reconstruction, , number of gene markers, or min
# this code is for converting an OrthoXML file to a set of Fasta files as Ortholougous groups

from ete3 import Tree
import sys
import os
from FastOMA.zoo.hog.convert import orthoxml_to_newick
from Bio import SeqIO




def max_og_tree(tree):
    for node in tree.traverse("preorder"):
        # for node in xml_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "attriremoved") and n.attriremoved==True):
        if not node.is_leaf() and hasattr(node,"Ev") and node.Ev == 'duplication':       # node.name[:3] == "dup"
            dup_node = node
            children = dup_node.get_children()
            list_num_species = []
            for child in children:
                child_name_leaves = child.get_leaves()
                species_list = []
                for leaf in child_name_leaves:
                    name = leaf.name
                    if name[:3] == "no_":
                        name = leaf.name.split("_")[-1]
                    if name in species_dic:
                        species_name = species_dic[name]
                        species_list.append(species_name)
                    else:
                        print("species not in the dic ",name)
                species_set = set(species_list)
                list_num_species.append(len(species_set))
            index_max_species = list_num_species.index(max(list_num_species))
            # if there are few children with identical number of species, the case would be not a polytomi but two children with one species
            # num_occurence = [1 for i in list_num_species if i == max(list_num_species)]
            # if len(num_occurence) > 1:
            #    print("please check this case with the developer the tool. The tree has polytomy.")
            child_max_species = children[index_max_species]
            children_to_remove = [i for i in children if i != child_max_species]
            for child_to_remove in children_to_remove:
                for i in child_to_remove.get_leaves():
                    i.in_og = "no"


    og_prot_list = []
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            if hasattr(node,"in_og") and node.in_og == "no":
                pass # print(node.name)
            else:
                og_prot_list.append(node.name)

    return og_prot_list




input_orthoxml=sys.argv[1] # "out_folder/output_hog_.orthoxml" 
rhog_all_folder = sys.argv[2] # "out_folder/rhogs_all/" 
fasta_format = ".fa" # of the rhogs_all


output_file = "maximal_og_prot.tsv"


trees, species_dic = orthoxml_to_newick(input_orthoxml, return_gene_to_species=True) # encode_levels_as_nhx=False,  xref_tag="protId",
print("We extracted "+str(len(trees))+" trees  in NHX format from the input HOG orthoxml"+input_orthoxml)


OGs = {}
for hog_id, tree_string in trees.items():

    tree = Tree(tree_string,format=1)
    og_prot_list = max_og_tree(tree)
    OGs[hog_id] = og_prot_list


print("done")


with open(output_file, 'w') as handle:
    for hog_id, og_prot_list in OGs.items():
        line_text = str(hog_id)+"\t"+str(og_prot_list)+"\n"
        handle.write(line_text)
handle.close()

print("We wrote the protein families information in the file "+output_file)


out_folder_ogs = "OGs_maximal/"
os.makedirs(out_folder_ogs)

print("writing "+str(len(OGs))+" OGs as fasta files in folder " +out_folder_ogs )
for hog_id, og_prot_list in OGs.items(): #hog_id="HOG_0667494_sub10524"
    rhog_id = "_".join(hog_id.split("_")[:2]) 

    rhogs_all_address = rhog_all_folder + rhog_id + "."+fasta_format
    rhogs_all_prots = list(SeqIO.parse(rhogs_all_address, "fasta"))

    og_prots = []
    og_prot_list = OGs[hog_id]
    for rhogs_prot in rhogs_all_prots:
        if rhogs_prot.id.split("||")[0] in og_prot_list:
            og_prots.append(rhogs_prot)

    og_id =  "OG_" + hog_id  # one OG per rootHOG      # "/HOG_"+ str(rhogid_num).zfill(7)
    SeqIO.write(og_prots, out_folder_ogs+og_id+".fa", "fasta")   
print("writing done")

