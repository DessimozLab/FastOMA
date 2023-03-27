

# to convert orthoxml to rootHOG (gene familes)
#
from Bio import SeqIO

from fastoma.zoo.hog.convert import extract_flat_groups_at_level



input_orthoxml = "folder/my_hogs.orthoxml"
output_folder = "output_folder"


#oxml = #os.path.join(file, "coverage_test_files", simpleEx.orthoxml")
toplevel_groups = []
for grp in extract_flat_groups_at_level(input_orthoxml):
    toplevel_groups.append(set(g.xref for g in grp))

print(len(toplevel_groups))

for toplevel_group_idx , toplevel_group in enumerate(toplevel_groups):
    SeqIO.write(toplevel_group, output_folder+"family_"+str(toplevel_group_idx)+".fasta", "fasta")



#  unversial genes
# prot_name_universal = []
# for group in toplevel_groups:
#     if len(group) > 0.9 * 2181:
#         species = [prot_specis[prot] for prot in group]
#         species_unq = set(species)
#         if len(species_unq) > 0.9 * 2181:
#             prot_name_universal.append(group)
#
# len(prot_name_universal)