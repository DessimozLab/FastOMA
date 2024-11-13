
from FastOMA.infer_roothogs import fastoma_infer_roothogs
from FastOMA._wrappers import logger
from FastOMA.infer_subhogs import fastoma_infer_subhogs


# --low-so-detection --fragment-detection

# --input-rhog-folder ./bb/ --parrallel True  --species-tree species_tree.nwk

#a=2
#fastoma_infer_subhogs()
#  proteome    --hogmap hogmaps   --splice splice  --out-rhog-folder out
import sys
logger.debug("hello ")
folder="pycharm_projects/fastoma_test/"
sys.argv.extend(['--proteomes', folder+"proteome"])
sys.argv.extend(['--hogmap', folder+"hogmaps"])
sys.argv.extend(['--splice', folder+"splice"])
sys.argv.extend(['--out-rhog-folder', folder+"out"])
sys.argv.extend(['-vv'])
fastoma_infer_roothogs()

a=2 # a
#
# from FastOMA.zoo.hog import transform
#
# #from zoo.tree_utils import collapse, gene_species, transform, HOG_coverages
#
# import io
# import lxml.etree
# orthoxml_file = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/benchmark-webservice3/orthoxml/euk_omamer200.dev8_13oct.orthoxml"
#
#
# orthxml_str = []
# with open(orthoxml_file, "r") as f:
#     for i in f:
#         orthxml_str.append(i)
# print(len(orthxml_str))
# dic_gene_integer={}
# for line in orthxml_str:
#     if "gene id" in line:
#         found=False
#         gene_int= line.split("\"")[1]
#         gene_name = line.split("\"")[3]
#         dic_gene_integer[gene_int] = gene_name
#
#
#
# orthoxml_etree=lxml.etree.parse(orthoxml_file)
#
# pw_orthologs_integer = sorted(list(transform.iter_pairwise_relations(orthoxml_etree)))
# # iter_pairwise_relations(obj, rel_type=None    (def:'ortholog' , but possible to use 'paralog')
# print(len(pw_orthologs_integer))
# print(pw_orthologs_integer[:2])
# pw_orthologs_gene =[]
# for pair in pw_orthologs_integer:
#     pw_orthologs_gene.append((dic_gene_integer[pair[0]],dic_gene_integer[pair[1]]))
#
#
#
# print(len(pw_orthologs_gene))
#
# output_file = open(orthoxml_file+"_pairs.tsv","w")
# for  pair in pw_orthologs_gene:
#     output_file.write(pair[0]+"\t"+pair[1]+"\n")
#
# output_file.close()


#
#
# # orthoxml_handle= open(orthoxml_file,"r")
# # orthoxml =""
# # for line in orthoxml_handle:
# #     orthoxml+=line
#
#
# from xml.etree.ElementTree import XMLParser
#
# parser = XMLParser()
# with open(orthoxml_file, 'rb') as xml:
#     for chunk in xml:
#         parser.feed(chunk)
# parser.close()
#
#
# lxml.etree.parse(oxml)
#
# orthoxm= lxml.etree.parse(orthoxml)
#
# # expected = [("1", "2"), ("1", "3"), ("1", "4"), ("1", "5"), ("1", "6"),
# #             ("2", "5"), ("2", "6"), ("3", "4"), ("3", "5"), ("3", "6"),
# #             ("4", "5"), ("4", "6"), ("5", "6")]
# #    self.assertEqual(expected, pw_orthologs)
#
# from xml.etree import ElementTree
# tree = ElementTree.parse(orthoxml_file)
# root = tree.getroot()
