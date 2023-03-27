


# orthologous per group

# for speciefic taxonomic level,  top level
# option whether include its substes,  groups that emerge after
# for root combiaiton is not alowed

# (1) folder for fasta file
# (2) text file, per group, gene id


#  (3) extract gene markers, single copy ortho roots ,  for species tree reconstruction, , number of gene markers, or min



# this code is for converting an OrthoXML file to a set of Fasta files as Ortholougous groups


def xmlTree2OGgenes(xmltree, prot_species_proteome):
    for node in xml_tree.traverse("preorder"):
        # for node in xml_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "attriremoved") and n.attriremoved==True):

        if not node.is_leaf() and node.name[:3] == "dup":
            dup_node = node
            children = dup_node.get_children()
            list_num_species = []
            for child in children:
                child_name_leaves = child.get_leaves()
                species_ls = []
                for leaf in child_name_leaves:
                    name = leaf.name
                    if name[:3] == "no_":
                        name = leaf.name.split("_")[-1]

                    if name in prot_species_proteome:
                        species_name = prot_species_proteome[name]
                    else:
                        print(name)

                    species_ls.append(species_name)
                species_set = set(species_ls)

                list_num_species.append(len(species_set))
            index_max_species = list_num_species.index(max(list_num_species))
            child_max_species = children[index_max_species]

            remov_children = [i for i in children if i != child_max_species]
            for remov_child in remov_children:
                for i in remov_child.get_leaves():
                    i.name = "no_" + i.name

                # remov_child.attriremoved = True

    prot_name_list_OG = []
    for node in xml_tree.traverse("preorder"):
        if node.is_leaf() and node.name[:3] != "no_":
            prot_name_list_OG.append(node.name)

    # print("number of prots in the tree was", len(xml_tree),"but OG list contains",len(prot_name_list_OG))
    return prot_name_list_OG
