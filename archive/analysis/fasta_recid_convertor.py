
from Bio import SeqIO
import sys

hogfile=sys.argv[1] #"HOG0041971.fa" 
records = []
for record in SeqIO.parse(hogfile, "fasta"):
    # id='TS022Sopim_PP_06T000597.1', name='TS022Sopim_PP_06T000597.1', description='TS022Sopim_PP_06T000597.1 TS022Sopim_PP_06T000597.1||TS022||1038013718 
    record.id = record.description.split(" ")[1]
    records.append(record)
SeqIO.write(records, "HOG_"+hogfile[3:], "fasta")
print("the updated fasta file is written" +"HOG_"+hogfile[3:])
