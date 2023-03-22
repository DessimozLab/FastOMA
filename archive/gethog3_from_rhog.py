from Bio import SeqIO
import pickle
from os import listdir
import os


# use this python code to convert gethog2 rhogs to fastoma
# check the following folder address
# put the followings in  in in_folder
# 1- gene_id_dic_xml.pickle
# 2- species tree nwk
# 3- rhogs_all


# then run nextfolow pypilei
# nextflow /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/fastoma/archive/gethog3_rhog.nf --input_folder in_folder   --output_folder out_folder -c  /work/FAC/FBM/DBC/cdessim2/default/smajidi1/pycharm_projects/fastoma/nextflow_slurm.config



# folder containing proteomes  species
addr= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog_analysis/DB/"

# output folder
addr2 = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog_analysis/"

# folder contining roothogs
addr_rhogs = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog_analysis/hogs/"

# please creat a folder named   rhogs_all
# output folder
addr3 = addr2 + "rhogs_all/"
#addr3 = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/gethog_analysis/rhogs/"







project_files = listdir(addr)
query_species_names = []
for file in project_files:
    fasta_format = file.split(".")[-1]
    if fasta_format == "fa" or fasta_format == "fasta":
        file_name_split = file.split(".")[:-1]
        query_species_names.append('.'.join(file_name_split))
        fasta_format = file.split(".")[-1]
print(len(query_species_names),query_species_names[-3:])

fasta_format="fa"
query_prot_recs = []
for query_species_names_idx, query_species_name in enumerate(query_species_names):
    prot_address = addr + query_species_name + "."+fasta_format
    prots_record = list(SeqIO.parse(prot_address, "fasta"))
    query_prot_recs.append(prots_record)

query_species_num = len(query_species_names)
print("The are "+str(query_species_num)+" species in the proteome folder.")



def add_species_name_gene_id(addr2,query_prot_recs, query_species_names):
    """
    adding the name of species to each protein record
        - based on file name
    adding gene id number, integer imposed by xml format
    output: updated version of input
    """
    #  _config.in_folder +
    gene_id_pickle_file = addr2+"/gene_id_dic_xml.pickle"
    max_num_prot = int(1e9)
    max_num_prot_per_sp = int(1e6)
    gene_id_name = {}
    for query_species_idx, query_species_name in enumerate(query_species_names):
        query_prot_records = query_prot_recs[query_species_idx]
        gene_counter = max_num_prot + query_species_idx * max_num_prot_per_sp
        gene_id_name[query_species_name] = []
        for query_prot_idx, query_prot_record in enumerate(query_prot_records):
            gene_idx_integer = gene_counter + query_prot_idx
            query_prot_name = query_prot_record.id
            if len(query_prot_name) > 230:
                #logger_hog.info("We are truncating the prot name as it may be problamatic for mafft, " + str(query_prot_name))
                query_prot_name = query_prot_name[:230]
            query_prot_record.id = query_prot_name + "||"+query_species_name+"||"+str(gene_idx_integer)
            gene_id_name[query_species_name].append((gene_idx_integer, query_prot_name))
    # this is used to create the first part of xml file.
    with open(gene_id_pickle_file, 'wb') as handle:
        pickle.dump(gene_id_name, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return query_prot_recs


query_prot_recs2 = add_species_name_gene_id(addr2,query_prot_recs, query_species_names)


max_num_prot = int(1e9)
max_num_prot_per_sp = int(1e6)
gene_id_name_old_new = {}

for query_species_idx, query_species_name in enumerate(query_species_names):
    query_prot_records = query_prot_recs[query_species_idx]
    gene_counter = max_num_prot + query_species_idx * max_num_prot_per_sp
    #gene_id_name[query_species_name] = []
    for query_prot_idx, query_prot_record in enumerate(query_prot_records):
        gene_idx_integer = gene_counter + query_prot_idx
        query_prot_name = query_prot_record.id
        if len(query_prot_name) > 230:
            #   logger_hog.info("We are truncating the prot name as it may be problamatic for mafft, " + str(query_prot_name))
            query_prot_name = query_prot_name[:230]
        query_prot_record_id = query_prot_name + "||"+query_species_name+"||"+str(gene_idx_integer)
        gene_id_name_old_new[query_prot_name] = query_prot_record_id

print("len(gene_id_name_old_new)",len(gene_id_name_old_new))


project_files = listdir(addr_rhogs)
hog_names = []
hog_records = []
for file in project_files:
    fasta_format = file.split(".")[-1]
    if fasta_format == "fa" or fasta_format == "fasta":
        file_name_split = file.split(".")[:-1]
        hog_names.append('.'.join(file_name_split))
        # fasta_format = file.split(".")[-1]

        prot_address = addr + file
        prots_record = list(SeqIO.parse(prot_address, "fasta"))
        hog_records.append(prots_record)

print(len(hog_records), len(hog_names), hog_names[:3], hog_records[-1][:5])


for hog_idx, hog_record in enumerate(hog_records):
    rhogid_num = int(hog_names[hog_idx][3:])
    for prot in hog_record:
        if "||" not in prot.id:
            prot_id_new = gene_id_name_old_new[prot.id]
            prot.id = prot_id_new

    SeqIO.write(hog_record, addr3 + "/HOG_" + str(rhogid_num).zfill(7) + ".fa", "fasta")

