



import sys

from FastOMA.zoo.hog import extract_flat_groups_at_level


"""
how to run
    python orthoxml2family.py my_hogs.orthoxml
    
- to convert orthoxml to rootHOG (protein families)
"""

input_orthoxml = sys.argv[1]
output_file = "families_prot.tsv"

toplevel_groups = []
for grp in extract_flat_groups_at_level(input_orthoxml):
    toplevel_groups.append(set(g.xref for g in grp))

# toplevel_groups is a list of sets

print("We extracted "+str(len(toplevel_groups))+" protein families from the input HOG orthoxml"+input_orthoxml)
print("The first one contain "+str(len(toplevel_groups[0]))+" proteins.")

with open(output_file, 'w') as handle:
    for toplevel_group_idx, toplevel_group in enumerate(toplevel_groups):
        line_text = str(toplevel_group_idx)+"\t"+str(toplevel_group)
        handle.write(line_text)
handle.close()

print("We wrote the protein families information in the file "+output_file)


# we need to know the species name of each prot,  as prot_specis dic
# prot_name_universal = []
# for group in toplevel_groups:
#     if len(group) > 0.9 * 2181:
#         species = [prot_specis[prot] for prot in group]
#         species_unq = set(species)
#         if len(species_unq) > 0.9 * 2181:
#             prot_name_universal.append(group)
#
# len(prot_name_universal)