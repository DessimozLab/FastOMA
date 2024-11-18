

from xml.dom import minidom
import xml.etree.ElementTree as ET
import pickle
from FastOMA._utils_subhog import read_species_tree
from FastOMA.collect_subhogs import convert_speciestree_to_orthoxml_taxonomy
import sys
from FastOMA.transformer import header_transformer

from FastOMA.collect_subhogs import iter_hogs
from FastOMA.collect_subhogs import update_hogids
from pathlib import Path

```
python pickle2orthoxml.py  "no_header" "file_D0680685.pickle"

python pickle2orthoxml.py "selected_genes"  pickle_folder gene_id_dic_xml.pickle  "species_tree_checked.nwk"     # this will be slow.  gene_id_dic_xml.pickle is in the output of infer_roothogs
```

mode = sys.argv[1]  #"selected_genes" #"no_header" # "selected_genes"  "all_genes"

if mode=="no_header":

    input_pickle= sys.argv[2]   # "file_D0680685.pickle"
    handle=open(input_pickle,'rb')
    orthoxml_file = pickle.load(handle)

    print(len(orthoxml_file))
    xml_str = minidom.parseString(ET.tostring(orthoxml_file[0])).toprettyxml(indent="   ")

    with open(input_pickle+"_noheader.orthoxml","w") as out_file:
        out_file.write(xml_str)

if mode =="selected_genes":

    input_pickle = sys.argv[2] # a folder  of pickles  pickle_folder
    gene_id_pickle_file = sys.argv[3] # generated in infer_roothogs.
    # available in out_folder/temp_output/gene_id_dic_xml.pickle
    # this keeps the gene name and the gene integer ID used in orthoxml.
    species_tree = sys.argv[4] # "species_tree_checked.nwk"

    handle=open(input_pickle,'rb')
    orthoxml_file1 = pickle.load(handle) # todo might have two elements inside?
    gene_int_set = set()
    num_digit = 10  # integer ids  # assumption ?
    for orthoxml_part in orthoxml_file1:
        xml_str = minidom.parseString(ET.tostring(orthoxml_part)).toprettyxml(indent="   ")
        gene_int_set_i = set([int(i[1:num_digit+1]) for i in xml_str.split("geneRef id=")[1:] ])
        gene_int_set.update(gene_int_set_i)

    from datetime import datetime
    fastoma_version= "0"
    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/",
                                                   "origin": "FastOMA " + fastoma_version,
                                                   "originVersion": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                   "version": "0.5"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)  # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
    print("We read the gene_id_name dictionary with %d items", len(gene_id_name))

    speciestree = read_species_tree(species_tree)
    taxonomy, name2taxid = convert_speciestree_to_orthoxml_taxonomy(speciestree)
    print("Now creating the header of orthoxml")

    id_transform_= "noop" #  noop:No transformation, "UniProt":   '>sp|P68250|1433B_BOVIN' --> P68250""")

    id_transformer = header_transformer(id_transform_)

    #  #### create the header of orthoxml ####
    for query_species_name, list_prots in gene_id_name.items():
        first=True
        for (gene_idx_integer, query_prot_name) in list_prots:
            if gene_idx_integer in gene_int_set:
                if first:
                    species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "taxonId": str(name2taxid[query_species_name]), "NCBITaxId": "0"})
                    database_xml = ET.SubElement(species_xml, "database", attrib={"name": "database", "version": "2023"})
                    genes_xml = ET.SubElement(database_xml, "genes")
                    prot_id = id_transformer.transform(query_prot_name)
                    gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": prot_id})
                    first=False
                else:
                    prot_id = id_transformer.transform(query_prot_name)
                    gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": prot_id})



    print("gene_xml is created.")
    # orthoxml_file.append(taxonomy)

    scores = ET.SubElement(orthoxml_file, "scores")
    ET.SubElement(scores, "scoreDef", {"id": "CompletenessScore",
                                       "desc": "Fraction of expected species with genes in the (Sub)HOG"})

    #  #### create the groups of orthoxml   ####
    groups_xml = ET.SubElement(orthoxml_file, "groups")

    with open(input_pickle, 'rb') as handle:
        hogs_a_rhog_xml = pickle.load(handle)
    for idx, hog_a_rhog_xml in enumerate(hogs_a_rhog_xml):
        fam = idx # this could be improved
        groups_xml.append(update_hogids(fam, hog_a_rhog_xml, name2taxid))
    #for fam, hogs_a_rhog_xml in enumerate(iter_hogs(Path(pickle_folder)), start=1):
    #    groups_xml.append(update_hogids(fam, hogs_a_rhog_xml, name2taxid))
    print("converting the xml object to string.")

    output_xml_name= input_pickle+".orthoxml"
    with open(output_xml_name, 'wb') as fh:
        ET.indent(orthoxml_file, space='  ', level=0)
        orthoxml = ET.ElementTree(orthoxml_file)
        orthoxml.write(fh, encoding="utf-8", xml_declaration=True, )
    print("orthoxml is written in %s", output_xml_name)



