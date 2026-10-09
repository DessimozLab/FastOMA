import os
import sys
from Bio import SeqIO



"""
python generate_splice_file_Ensembl.py input_folder

"""



def parse_proteomes(folder=None, min_sequence_length=0):  # list_oma_species
    """
    
    parsing fasta files of proteins located in /proteome/
    using Bio.SeqIO.parse
    Each fasta file is for one species.  The file name is the species name.
    All proteomes should be with the same extension.
    output: prot_recs_lists: a dic with key as species name and  its value list of Biopython record of species.
    """
    # Adapted from https://github.com/DessimozLab/FastOMA/blob/main/FastOMA/_utils_roothog.py
    extension = ".fa"

    project_files = os.listdir(folder)
    print(folder, " as proteome folder, found ", len(project_files), " files.")
    species_names = []  # query/input species name based on the file name
    for file in project_files:
        species_name = file.split(extension)[0]
        #print("%s: %s | %s", file, species_name, ext)
        #if ext in ("fa", "fasta"):
        species_names.append(species_name)
            #fasta_format_keep = ext  # last one is stored either fa or fasta

    # todo accept all fasta formats in the input prtoeome folder, fasta, fa, fna, ..
    prot_recs_lists = {} # key: species name, value is a dic of query protein Biopython records. # 'MYCGE': [SeqRecord(seq=Seq('MDFDK
    #smallprot_recs_lists ={}
    for species_name in species_names:
        prot_address = os.path.join(folder, species_name +  extension)
        prots_record = list(rec for rec in SeqIO.parse(prot_address, "fasta") if len(rec) >= min_sequence_length)
        #prots_record_small = list(rec for rec in SeqIO.parse(prot_address, "fasta") if len(rec) < min_sequence_length)
        # logger.debug(prots_record)
        #print(f"{species_name} contains {len(prots_record)} that are at least {min_sequence_length} long.")
        prot_recs_lists[species_name] = prots_record
        #smallprot_recs_lists[species_name]=prots_record_small

    print("There are ", len(species_names), "species in the proteome folder.")
    return species_names, prot_recs_lists#, fasta_format_keep #, smallprot_recs_lists

def extract_spliced_prots(gene2prot_transcript):
    genes_withmore1=[]
    diff_prot_transcript=[]
    prots_with_splice={}
    for species_name, gene2prot_transcript_sp in gene2prot_transcript.items():
        
        for gene in gene2prot_transcript_sp.keys():
            prot_transcript_list=gene2prot_transcript_sp[gene]
            if len(prot_transcript_list)>1:
                genes_withmore1.append((species_name,gene))
                prots= [prot for  (prot,transcript) in prot_transcript_list]
                transcripts= [transcript for  (prot,transcript) in prot_transcript_list]
                if len(set(prots))!=len(set(transcripts)):
                    print(species_name,gene)
                if set(prots)!=set(transcripts):
                    diff_prot_transcript.append((species_name,gene))
                if species_name in prots_with_splice:
                    prots_with_splice[species_name].append(prots)
                else:
                    prots_with_splice[species_name]=[prots]
    print('There are ',len(genes_withmore1),' genes with more than one transcripts.')
    print('For these genes, there are ',len(diff_prot_transcript),' genes that the protein id is different than transcript id.')
    print('We found ',len(prots_with_splice),' species that have some genes with more than one transcripts.')
    return prots_with_splice




def extract_transcripts(species_names, prot_recs_lists):
    total_splices=0
    gene2prot_transcript={}
    for species_name, prot_recs_list in prot_recs_lists.items():
        gene2prot_transcript[species_name]={} # 
        
        for prot_rec in prot_recs_list:
            prot_name= prot_rec.id # >PSK40689.1 pep primary_assembly:GCA003013735v1:PYFQ01000001:776887:779901:1 gene:C7M61_000337 transcript:C7M61_000337_t1 gene_biotype:protein_coding transcript_biotype:protein_coding description:RIC1 domain-containing protein [Source:UniProtKB/TrEMBL;Acc:A0A2P7YXL8]
            transcriptid = prot_rec.description.split("transcript:")[1].split(" ")[0]
            geneid = prot_rec.description.split("gene:")[1].split(" ")[0]
            protid= prot_rec.id
            if geneid in gene2prot_transcript[species_name]:
                total_splices+=1
                gene2prot_transcript[species_name][geneid].append((protid,transcriptid))
                
            else:
                gene2prot_transcript[species_name][geneid] = [(protid,transcriptid)]
    
                
    print('Total number transcripts that are more than one transcript per gene is', total_splices)
    print('Total number of species',len(gene2prot_transcript),'. For species ',species_name,', the number of genes is ',len(gene2prot_transcript[species_name]))

    return gene2prot_transcript

def write_splice_files(proteome_folder,prots_with_splice,output_folder):

    print('Writing ', len(prots_with_splice), 'splice files')
    os.mkdir(output_folder) 
    for species_name, prots_list in prots_with_splice.items():
        file_splice= open(output_folder+species_name+".splice",'w')
        for prots in prots_list:
            line= ";".join(prots)
            file_splice.write(line+"\n")
        file_splice.close()
    print('Splice files per species are ready in ', output_folder)
    return 1
    


if __name__ == "__main__":


    proteome_folder=  sys.argv[1] 
    
    print("Parsing the proteome folder", proteome_folder)
    species_names, prot_recs_lists = parse_proteomes(proteome_folder)  
    
    gene2prot_transcript = extract_transcripts(species_names, prot_recs_lists)
    
    prots_with_splice = extract_spliced_prots(gene2prot_transcript)
    output_folder= output_folder= proteome_folder+"/../splice_files/"
    write_splice_files(proteome_folder,prots_with_splice, output_folder)
    print("done!")


