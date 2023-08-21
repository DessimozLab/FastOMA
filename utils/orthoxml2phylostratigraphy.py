
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")


working_folder="./"

nwk_path= working_folder+"in_folder/species_tree.nwk" # species tree should be pruned (no extra leaves)

tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
print(tree_str[:10])

orthoxml_path=folder+"out_folder/output_hog_.orthoxml"
ham_analysis = pyham.Ham(tree_str, orthoxml_path, use_internal_name=True)
print("Ham analysis done") # for a big orthoxml file it can take ~30mins

#phylostratigraphy

#create tree profile, classify all genomes by extant or ancestral, and get % of dup, lost, retained, and gained
treeprofile = ham_analysis.create_tree_profile(outfile= folder+"/out_folder/phylostratigraphy.html")
treemap = treeprofile.compute_tree_profile_full()

