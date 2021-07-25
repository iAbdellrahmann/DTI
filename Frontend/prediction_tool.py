# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 01:56:10 2021

@author: Abdel-Rahman
"""


import pandas as pd
from Bio import SeqIO
from Bio import Align
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.DataStructs
import networkx as nx
import random
from gensim.models import Word2Vec
import pickle
import warnings
import streamlit as st
import os
#from PIL import Image
#import awesome_streamlit as ast
def write():
    def New_Drug_DDS(sdf):
      supplier = rdkit.Chem.SDMolSupplier('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/structures.sdf')
      molecules = [mol for mol in supplier if mol is not None]
      fingerprints = dict()
      for mol in molecules:
          drugbank_id = mol.GetProp('DATABASE_ID')
          fingerprint = rdkit.Chem.AllChem.GetMorganFingerprint(mol, 2)
          fingerprints[drugbank_id] = fingerprint
      newDrug = rdkit.Chem.SDMolSupplier(sdf)
      newMolecule = [mol for mol in newDrug if mol is not None]
      newFingerprint = dict()
      for mol in newMolecule:
          new_drugbank_id = mol.GetProp('DATABASE_ID')
          new_fingerprint = rdkit.Chem.AllChem.GetMorganFingerprint(mol, 2)
          newFingerprint[new_drugbank_id] = new_fingerprint
      similarity_rows = list()
      for i in fingerprints.keys():
        if (i==new_drugbank_id):
            continue
        else:
            similarity = rdkit.DataStructs.DiceSimilarity(fingerprints[i], newFingerprint[new_drugbank_id])
            similarity = round(similarity, 4)
            if([new_drugbank_id,i, similarity] not in similarity_rows):
              similarity_rows.append([new_drugbank_id,i, similarity,f'https://go.drugbank.com/structures/{new_drugbank_id}/image.svg',f'https://go.drugbank.com/structures/{i}/image.svg'])
      if(len(similarity_rows)!=0):
          similarity_rows = pd.DataFrame(similarity_rows, columns=['compound_1', 'compound_2', 'label','lin_1','lin_2'])
          indexNames = similarity_rows[similarity_rows['label'] <0.5 ].index
          similarity_rows.drop(indexNames , inplace=True)
      return similarity_rows
    
    
    
    def get_randomwalk(node, path_length,G):
        
        random_walk = [node]
        
        for i in range(path_length-1):
            temp = list(G.neighbors(node))
            temp = list(set(temp) - set(random_walk))    
            if len(temp) == 0:
                break
    
            random_node = random.choice(temp)
            random_walk.append(random_node)
            node = random_node
            
        return random_walk
    
    
    def New_Drug_int(sdf):
      ddsf=New_Drug_DDS(sdf)
      all_edge_list_df = pd.read_csv("D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Graph EdgeList/all_edge_list_df.csv")
      all_edge_list_df=all_edge_list_df.iloc[: , 1:]
      newDF=all_edge_list_df.append(ddsf)
      G = nx.Graph()
      G=nx.from_pandas_edgelist(newDF,"compound_1","compound_2","label",create_using=nx.Graph())
      all_nodes = list(G.nodes())
      random_walks = []
      for n in all_nodes:
          for i in range(5):
              random_walks.append(get_randomwalk(n,10,G))
      warnings.filterwarnings('ignore')
      model = Word2Vec(window = 4, sg = 1, hs = 1,vector_size=64, negative = 0, # for negative sampling
                      alpha=0.03, min_alpha=0.0007,
                      seed = 14)
      model.build_vocab(random_walks, progress_per=2)
      model.train(random_walks, total_examples = model.corpus_count, epochs=20, report_delay=1)
      d={}
      for n in G.nodes():
          d[n]=model.wv[n].tolist()
      emb=pd.DataFrame.from_dict(d)
      emb_transposed = emb.T
      new_drug_emb=emb_transposed.loc[list(ddsf['compound_1'])[0]]
      interaction_dataset=pd.read_csv('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Graph EdgeList/uniprot links.csv')
      uniID=list()
      for index,row in interaction_dataset.iterrows():
        if (row['DrugBank ID'] in list(ddsf["compound_2"])):
          uniID.append(row["UniProt ID"])
      poss_interaction=list(set(uniID))
      cols=["F"+ str(i) for i in range(1,65)]
      cols.append('Drug ID')
      cols.append('Target ID')
      df= pd.DataFrame(data=None)
      for i in poss_interaction:
        possible_interaction_df = new_drug_emb.add(emb_transposed.loc[i])
        possible_interaction_df['Drug ID'] = list(ddsf['compound_1'])[0]
        possible_interaction_df['Target ID'] = i
        df= df.append(possible_interaction_df,ignore_index=True)
      df.columns = cols
      prediction_data = df.iloc[:,:-2].values
      load_model = pickle.load(open('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Frontend/svm_model.pkl', 'rb'))
      prediction = load_model.predict(prediction_data)
      df['Label']=prediction
      proteins=list()
      for i in range(len(prediction)):
        if prediction[i]==1:
          pr_id=list(df['Target ID'])[i]
          proteins.append([pr_id,f'https://alphafold.ebi.ac.uk/entry/{pr_id}'])
      true_int=[]
      if(list(ddsf['compound_1'])[0] in list(interaction_dataset['DrugBank ID'])):
          for index, row in interaction_dataset.iterrows():
              if row["DrugBank ID"]==list(ddsf['compound_1'])[0]:
                  true_int.append(row["UniProt ID"])
      if(len(true_int)!=0):
          true_int=pd.DataFrame(true_int)    
          true_int.columns=['Known Interactions']
          true_int     
      return proteins
    
    
    
        
    
    def new_protien_sim(file):
      UniSequence = list(SeqIO.parse(file, "fasta"))
      id=UniSequence[0].id
      id=id.split('|')[1]
      protien=UniSequence[0].seq
      aligner = Align.PairwiseAligner()
      protein_dataset=pd.read_csv('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Inner_join/protein_inner_join.csv')
      protien_similer=[]
      for i in range(len(protein_dataset)):
        if(id==protein_dataset['id'][i]):
              continue
        Score = aligner.score(protien,protein_dataset['sequence'][i])
        sim=Score/max(len(protein_dataset['sequence'][i]),len(protien))
        if (sim>0.4):
          pr_id=protein_dataset['id'][i]
          protien_similer.append([id,pr_id,sim,f'https://alphafold.ebi.ac.uk/entry/{id}',f'https://alphafold.ebi.ac.uk/entry/{pr_id}'])
      if(len(protien_similer)!=0):
          protien_similer=pd.DataFrame(protien_similer)
          protien_similer.columns=['compound_1', 'compound_2', 'label','Link 1','Link 2']
      return protien_similer
    
    
    
    
    def New_protien_int(file):
        protien_sim=new_protien_sim(file)
        interaction_dataset=pd.read_csv('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Graph EdgeList/uniprot links.csv')
        drugID=list()
        for index,row in interaction_dataset.iterrows():
          if (row['UniProt ID'] in list(protien_sim["compound_2"])):
            drugID.append(row["DrugBank ID"])
        all_edge_list_df = pd.read_csv("D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Graph EdgeList/all_edge_list_df.csv")
        all_edge_list_df=all_edge_list_df.iloc[: , 1:]
        newDF=all_edge_list_df.append(protien_sim)
        G = nx.Graph()
        G=nx.from_pandas_edgelist(newDF,"compound_1","compound_2","label",create_using=nx.Graph())
        all_nodes = list(G.nodes())
        random_walks = []
        for n in all_nodes:
            for i in range(5):
                random_walks.append(get_randomwalk(n,10,G))
        warnings.filterwarnings('ignore')
        model = Word2Vec(window = 4, sg = 1, hs = 1,vector_size=64, negative = 0, # for negative sampling
                        alpha=0.03, min_alpha=0.0007,
                        seed = 14)
        model.build_vocab(random_walks, progress_per=2)
        model.train(random_walks, total_examples = model.corpus_count, epochs=20, report_delay=1)
        d={}
        for n in G.nodes():
            d[n]=model.wv[n].tolist()
        emb=pd.DataFrame.from_dict(d)
        emb_transposed = emb.T
        new_drug_emb=emb_transposed.loc[list(protien_sim['compound_1'])[0]]
        poss_interaction=list(set(drugID))
        cols=["F"+ str(i) for i in range(1,65)]
        cols.append('Drug ID')
        cols.append('Target ID')
        df= pd.DataFrame(data=None)
        for i in poss_interaction:
          possible_interaction_df = new_drug_emb.add(emb_transposed.loc[i])
          possible_interaction_df['Drug ID'] = i
          possible_interaction_df['Target ID'] = list(protien_sim['compound_1'])[0]
          df= df.append(possible_interaction_df,ignore_index=True)
        df.columns = cols
        prediction_data = df.iloc[:,:-2].values
        load_model = pickle.load(open('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project/Frontend/svm_model.pkl', 'rb'))
        prediction = load_model.predict(prediction_data)
        df['Label']=prediction
        drugs=list()
        for i in range(len(prediction)):
          if prediction[i]==1:
            dr_id=list(df['Drug ID'])[i]
            drugs.append([dr_id,f'https://go.drugbank.com/structures/{dr_id}/image.svg'])
        true_int=[]
        if(list(protien_sim['compound_1'])[0] in list(interaction_dataset['UniProt ID'])):
            for index, row in interaction_dataset.iterrows():
                if row["UniProt ID"]==list(protien_sim['compound_1'])[0]:
                    true_int.append(row["DrugBank ID"])
        if(len(true_int)!=0):
            true_int=pd.DataFrame(true_int)    
            true_int.columns=['Known Interactions']
            true_int     
        return drugs
    
    
    #image = Image.open('D:/1911/Bioinformatics/Level 4/Second semester/Graduation Project ours/streamlit/logo.jpg')
    #st.image(image, width = 500)

    
    
    def file_selector(folder_path='.'):
        filenames=[_ for _ in os.listdir(folder_path) if _.endswith('.sdf') or _.endswith('.fasta')]
        selected_filename = st.sidebar.selectbox('Select your Protein or Drug', filenames)
        return os.path.join(folder_path, selected_filename)
    
    
    st.title("Methodology")
    st.markdown("""
                * Based on Drug-Drug Similarities, Protein-Protein Similarities and Known Drug-Target Interaction.
                Using Deepwalk Algorithm which is working on top of random walk and Word2Vec Model
                on these Similarities and interaction we will extract features which then passed to SVM
                Model for Prediction.
                """)
    
    st.sidebar.header('Prediction Panel')
    filename = file_selector()
    st.sidebar.write('You selected `%s`' % filename)
    st.sidebar.markdown("""
    [Example SDF file](https://go.drugbank.com/structures/small_molecule_drugs/DB11331.sdf)
    """)
    st.sidebar.markdown("""
    [Example Protein file](https://www.uniprot.org/uniprot/P00734.fasta)
    """)
    if st.sidebar.button('New Interaction For Protein'):
        n_p=New_protien_int(filename)
        df=pd.DataFrame(n_p)
        df.columns=['Drugs',"3D structure"]
        st.write(""" 
                 # Possible Drugs
                 """)
        st.dataframe(df)
        
    if st.sidebar.button('Protein Protein Similarity'):
        p=new_protien_sim(filename)
        p=pd.DataFrame(p)
        p.columns=['Your Protein','Similar Protein','Similarity',"Yor drug's 3D structure","Similar drug 3D structure"]
        st.write(""" 
                 # Protein Protein Similarity
                 """)
        st.dataframe(p)
    
    if st.sidebar.button('New Interaction For Drug'):
        n_d=New_Drug_int(filename)
        df=pd.DataFrame(n_d)
        st.write(""" 
                 # Possible Proteins
                 """)
        df.columns=['Possible Interactions',"3D structure"]
        st.dataframe(df)
    if st.sidebar.button('Drug Drug Similarity'):
        d=New_Drug_DDS(filename)
        d.columns=['Your Drug','Similar Drug','Similarity',"Yor protein's 3D structure","Similar protein 3D structure"]
        st.write(""" 
                 # Drug Drug Similarity
                 """)
        st.dataframe(d)
    