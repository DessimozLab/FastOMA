

# to convert orthoxml to rootHOG (gene familes)
#


from zoo.hog import extract_flat_groups_at_level
folder= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog3_eukaryota/run28jan/working_nf/"
file= folder+"hog_euk_28jan_all.orthoxml"
file

#oxml = #os.path.join(file, "coverage_test_files", simpleEx.orthoxml")
toplevel_groups = []
for grp in extract_flat_groups_at_level(file):
    toplevel_groups.append(set(g.xref for g in grp))

print(len(toplevel_groups))
toplevel_groups[0]


#

write them as fasta files


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