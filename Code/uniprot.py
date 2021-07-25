# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 00:30:59 2021

@author: Abdel-Rahman
"""
from Bio import SeqIO
UniSequence = list(SeqIO.parse("uniprot_sprot.fasta", "fasta"))
print("The sequence of the Unknown fasta file is ",UniSequence[0])
print("The sequence of the Unknown fasta file is ",UniSequence[0].seq)
print("The ID sequence of the Unknown fasta file is ",UniSequence[0].id)
def extract_id(header):
    return header.split('|')[1]        
        
# importing pandas as pd   
import pandas as pd   
       
# list of id, sequence  
id = []
seq=[]
for seq_record in UniSequence:
    id.append(extract_id(seq_record.id))
for seq_record in UniSequence:
    seq.append(seq_record.seq)
# dictionary of lists   
dict = {'id': id,'sequence':seq}   
       
df = pd.DataFrame(dict)  
df.head()
# saving the dataframe  
df.to_csv('protein_dataset.csv')