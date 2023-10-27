
from FastOMA.zoo.hog import transform

#from zoo.tree_utils import collapse, gene_species, transform, HOG_coverages

import io
import lxml.etree
import sys
orthoxml_file = sys.argv[1]
#"/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_qfo/benchmark-webservice3/orthoxml/euk_omamer200.dev8_13oct.orthoxml"


orthxml_str = []
with open(orthoxml_file, "r") as f:
    for i in f:
        orthxml_str.append(i)
print(len(orthxml_str))
dic_gene_integer={}
for line in orthxml_str:
    if "gene id" in line:
        found=False
        gene_int= line.split("\"")[1]
        gene_name = line.split("\"")[3]
        dic_gene_integer[gene_int] = gene_name



orthoxml_etree=lxml.etree.parse(orthoxml_file)

pw_orthologs_integer = sorted(list(transform.iter_pairwise_relations(orthoxml_etree)))
# iter_pairwise_relations(obj, rel_type=None    (def:'ortholog' , but possible to use 'paralog')
print(len(pw_orthologs_integer))
print(pw_orthologs_integer[:2])
pw_orthologs_gene =[]
for pair in pw_orthologs_integer:
    pw_orthologs_gene.append((dic_gene_integer[pair[0]],dic_gene_integer[pair[1]]))



print(len(pw_orthologs_gene))
print(pw_orthologs_gene[:2])


output_file = open(orthoxml_file+"_pairs.tsv","w")
for  pair in pw_orthologs_gene:
    output_file.write(pair[0]+"\t"+pair[1]+"\n")

output_file.close()
