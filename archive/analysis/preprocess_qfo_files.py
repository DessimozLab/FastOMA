

# Proteins in each file belong to the same species.

# change the name of each file based on the species name inside each prot id


from os import listdir
from Bio import SeqIO
import os

working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
prot_folder = working_folder + "/omamer_search_old/proteome/"
project_files = listdir(prot_folder)
query_species_names_old = []
query_species_names_new = []
for file in project_files:
    if file.split(".")[-1] == "fa":
        file_name_split = file.split(".")[:-1]
        query_species_name_old = '.'.join(file_name_split)
        prot_address = prot_folder + query_species_name_old + ".fa"
        prots_record = list(SeqIO.parse(prot_address, "fasta"))
        prot_record = prots_record[0]
        prot_name = prot_record.name  # 'tr|E3JPS4|E3JPS4_PUCGT
        query_species_name_new = prot_name.split("|")[-1].split("_")[-1].strip()
        # if query_species_name_new == 'RAT': query_species_name_new = "RATNO"
        query_species_names_old.append(query_species_name_old)
        query_species_names_new.append(query_species_name_new)

os.mkdir(working_folder+"/omamer_search")
os.mkdir(working_folder+"/omamer_search/proteome/")
os.mkdir(working_folder+"/omamer_search/hogmap")


for idx, query_species_name_old in enumerate(query_species_names_old):
    query_species_name_new = query_species_names_new[idx]

    prot_address_old = working_folder + "omamer_search_old/proteome/" + query_species_name_old + ".fa"
    prot_address_new = working_folder + "omamer_search/proteome/" + query_species_name_new + "_.fa"
    os.system('cp ' + prot_address_old + ' ' + prot_address_new)

    hogmap_address_old = working_folder + "omamer_search_old/hogmap/" + query_species_name_old + ".hogmap"
    hogmap_address_new = working_folder + "omamer_search/hogmap/" + query_species_name_new + "_.hogmap"
    os.system('cp ' + hogmap_address_old + ' ' + hogmap_address_new)


# 13:54:16 - the species DANRE  already exists in the oma database, remove them first



print("done")