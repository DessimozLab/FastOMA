
import os
import sys

in_folder=sys.argv[1] 
fastoma_out_folder=sys.argv[2] 

"""
To run
python HOGcomposition.py in_folder output
in_folder : where proteome folder is
output: where RootHOGs.tsv is

output: 
species composition table (HOGs in rows, species in columns and each cell to be the protein names of that HOG for each species)
"""
project_files = os.listdir(in_folder+"/proteome/")


fasta_format_keep = ""
species_names = []  # query/input species name based on the file name
for file in project_files:
    species_name, ext = file.rsplit('.', 1)
    if ext in ("fa", "fasta"):
        species_names.append(species_name)
        fasta_format_keep = ext  # last one is stored either fa or fasta
print("number of species:", len(species_names))



prot2species = {} 
for species_name in species_names:
    prot_address = os.path.join(in_folder+"/proteome/", species_name + "." + fasta_format_keep)

    file_prot = open(prot_address,'r')
    for line in file_prot:
        if line.startswith(">"):
            prot_name=line.strip().split("\t")[0].split(" ")[0][1:]
            prot2species[prot_name] = species_name
print("total number of proteins in the fasta files",len(prot2species))

roothog_dic={}
roothog_file= open(fastoma_out_folder+"/RootHOGs.tsv",'r')
for line in roothog_file:
    if line.startswith("RootHOG"):
        continue
    roothog,protein, omamerroothog= line.strip().split("\t")
    if roothog in roothog_dic:
        roothog_dic[roothog].append(protein)
    else:
        roothog_dic[roothog]=[protein]
print("number of HOGs", len(roothog_dic))


file_out=open("HOGcomposition.tsv",'w')
file_out.write('RootHOG'+'\t'+'\t'.join(species_names)+'\n')
for roothog, prots in roothog_dic.items():

    prot_species_dic={}
    for species_name in species_names:
        prot_species_dic[species_name]=[]
    
    for prot in prots:
        species_name=prot2species[prot]
        prot_species_dic[species_name].append(prot)

    prot_species=[]
    for species_name in species_names:
        prots_raw= prot_species_dic[species_name]
        if prots_raw:
            prot_species.append(str(prots_raw)[1:-1])
        else:
            prot_species.append('NA')
        
            
    #prot_species_str=','.join(prot_species)
    file_out.write(roothog+'\t'+'\t'.join(prot_species)+'\n')
    
file_out.close()

print("output is written as file ", "HOGcomposition.tsv")
