#!/usr/bin/env python
# coding: utf-8

# In[52]:

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import sys


# # Parameters

# In[53]:


max_N=0.05
max_per_sp=2


# In[54]:


ref = sys.argv[1]
# ref = 'NCBI_18s'


# In[55]:


if ref == 'NCBI_18s':
    gb_file = 'NCBI_18s.gb'
    gene = ['rrn18','18S rRNA','18S ribosomal RNA']
    acc_type = ['gene','rRNA']
    min_len=1400; max_len=2400
elif ref == 'NCBI_28s':
    gb_file = 'NCBI_28s.gb'
    gene = ['rrn28','28S rRNA','28S ribosomal RNA']
    acc_type = ['gene','rRNA']
    min_len=3000; max_len=3800
elif ref == 'NCBI_16s':
    gb_file = 'NCBI_16s.gb'
    gene = ['rrn16','16S rRNA','16S ribosomal RNA']
    acc_type = ['gene','rRNA']
    min_len=1200; max_len=1800
elif ref == 'NCBI_23s':
    gb_file = 'NCBI_23s.gb'
    gene = ['rrn23','23S rRNA','23S ribosomal RNA']
    acc_type = ['gene']
    min_len=2500; max_len=2900
elif ref == 'NCBI_rbcL':
    gb_file = 'NCBI_rbcL.gb'
    gene = ['rbcL','rbcl','ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit']
    acc_type = ['CDS']
    min_len=1100; max_len=1500
elif ref == 'NCBI_trnL':
    gb_file = 'NCBI_trnL.gb'
    gene = ['trnL','tRNA-Leu']
    acc_type = ['tRNA']
    min_len=30; max_len=80
elif ref == 'NCBI_ITS1':
    gb_file = 'NCBI_ITS1.gb'
    gene = ['ITS','ITS1','internal transcribed spacer 1']
    acc_type = ['misc_RNA']
    min_len=180; max_len=280
elif ref == 'NCBI_ITS2':
    gb_file = 'NCBI_ITS2.gb'
    gene = ['ITS2','internal transcribed spacer 2']
    acc_type = ['misc_RNA']
    min_len=170; max_len=280
elif ref == 'NCBI_rpl2':
    gb_file = 'NCBI_rpl2.gb'
    gene = ['rpl2','ribosomal protein L2']
    acc_type = ['cds']
    min_len=500; max_len=1500
elif ref == 'NCBI_ndhf':
    gb_file = 'NCBI_ndhf.gb'
    gene = ['ndhf','ndhF']
    acc_type = ['cds']
    min_len=1000; max_len=2500


# # Main

# In[56]:


print(gb_file,gene,acc_type,min_len,max_len)


# In[57]:


def get_qualifier(feature, attribute):
    try:
        return feature.qualifiers[attribute][0]
    except:
        return None 


# In[58]:


# %%time
count=0
for line in open(gb_file): 
    if 'LOCUS' in line:
        count += 1
print(count)


# In[59]:


# %%time
print('reading genbank_file',end='...')
rec_ls = []; rec_rm=[]; rec_count=0
for record in SeqIO.parse(gb_file, "genbank"):
    rec_count += 1
    if record.features:
        for feature in record.features:
            if (feature.type in acc_type):
                seq_dic={}
                seq_dic['Locus'] = record.id
                seq_dic['type'] = feature.type
                if (feature.type in ['gene','CDS']):
                    seq_dic['gene'] = get_qualifier(feature, 'gene')
                elif (feature.type in ['rRNA','tRNA','misc_RNA']):
                    seq_dic['gene'] = get_qualifier(feature, 'product')
                if seq_dic['gene'] in gene:
                    seq_dic['Seq'] = str(feature.location.extract(record).seq)
                    seq_dic['Len'] = len(seq_dic['Seq'])
                    seq_dic['Nn'] = seq_dic['Seq'].count('N')
                    for feature in record.features:
                        if (feature.type == "source"):
                            seq_dic['sci_name'] = get_qualifier(feature, 'organism')
                            seq_dic['mol_type'] = get_qualifier(feature, 'mol_type')
                            seq_dic['TaxID'] = get_qualifier(feature, 'db_xref').replace('taxon:','')
                    rec_ls.append(seq_dic)      
                else:
                    rec_rm.append(seq_dic)
print('read',rec_count,'accessions')


# In[60]:


rec_df = pd.DataFrame(rec_ls)
print(rec_df.shape[0],'entries for',rec_df.sci_name.nunique(),'species')
print(rec_df.groupby('type').size().sort_values(ascending=False).to_dict())
print(rec_df.groupby('gene').size().sort_values(ascending=False).to_dict())


# In[61]:


scut=min_len;
print(rec_df.Len.quantile([.01,.05,.1,0.5,.9,.95,.99]).to_dict())
print(rec_df[rec_df.Len>scut].Len.quantile([.01,.05,.1,0.5,.9,.95,.99]).to_dict())
print(rec_df[rec_df.Len>scut].Len.median()+(rec_df[rec_df.Len>scut].Len.std()*2))
print(rec_df[rec_df.Len>scut].Len.median()-(rec_df[rec_df.Len>scut].Len.std()*2))
rec_df.Len.hist(bins=100);


# In[62]:


rec_df['rN'] = rec_df.Nn/rec_df.Len
print('Removing',rec_df[rec_df.rN>=max_N].shape[0],'accessions with too many Ns')
print(rec_df.shape[0],end=' > ')
rec_df = rec_df[rec_df.rN<max_N]
print(rec_df.shape[0])
print('Removing',rec_df[rec_df.Len<min_len].shape[0],'accessions too small')
print(rec_df.shape[0],end=' > ')
rec_df = rec_df[rec_df.Len>=min_len]
print(rec_df.shape[0])
print('Removing',rec_df[rec_df.Len>max_len].shape[0],'accessions too long')
print(rec_df.shape[0],end=' > ')
rec_df = rec_df[rec_df.Len<=max_len]
print(rec_df.shape[0])


# In[63]:


rm_char='[]()×'
for char in rm_char:
    rec_df['sci_name'] = rec_df['sci_name'].str.replace(char,'')


# In[64]:


print('sending',rec_df.sci_name.nunique(),'species names to WCVP_taxo')
rec_df.groupby('sci_name').head(1).sci_name.to_csv(gb_file.replace('.gb','_NCBI.csv'),index=False)


# In[65]:


print('running wcvp_taxo',end='...')
# print(os.system('python ../../PAFTOL_DB/wcvp_taxo.py ../../PAFTOL_DB/wcvp_v5_jun_2021.txt ' + \
#           gb_file.replace('.gb','_NCBI.csv') + ' -g -s similarity_genus -d divert_genusOK'))
print(os.system('python wcvp_taxo.py wcvp_v5_jun_2021.txt ' +           gb_file.replace('.gb','_NCBI.csv') + ' -g -s similarity_genus -d divert_genusOK'))
wcvp = pd.read_csv(gb_file.replace('.gb','_NCBI_wcvp.csv'))
wcvp = wcvp[wcvp.sci_name.notnull()]
print('found',wcvp.sci_name.nunique(),'species in WCVP')
print(rec_df.shape[0],end=' > ')
rec_df = pd.merge(rec_df.rename(columns={'sci_name':'Ini_sci_name'}),wcvp,how='inner',on='Ini_sci_name')
print(rec_df.shape[0])


# In[66]:


for char in rm_char:
    rec_df['sci_name'] = rec_df['sci_name'].str.replace(char,'')


# In[67]:


sp_count = rec_df.groupby('sci_name').size().to_frame()
print('reducing dataset to max',max_per_sp,'accessions per species, ',(sp_count[0]>2).sum())
print(rec_df.shape[0],end=' > ')
rec_df = rec_df.sort_values('Len',ascending=False).groupby('Locus').head(1).groupby('sci_name').head(max_per_sp)
print(rec_df.shape[0])
rec_df = rec_df.sort_values(['family','genus','sci_name']).reset_index(drop=True)
print('f:',rec_df.family.nunique(),'g:',rec_df.genus.nunique(),'s:',rec_df.sci_name.nunique())


# In[68]:


# rec_df = rec_df[rec_df['type']=='gene']


# In[69]:


types = list(rec_df.type.unique())
print(rec_df.groupby('type').size().to_dict())
rec_fasta=[]
for idx, row in rec_df.iterrows():
    record = SeqRecord(Seq(row.Seq))
    record.id = row.Locus
    record.description = ';gene=' + row.gene + ',type=' + row.type     + ',f=' + row.family + ',g=' + row.genus + ',s=' + row.sci_name + ',ini_s=' + row.Ini_sci_name + ';'
    rec_fasta.append(record)
SeqIO.write(rec_fasta,gb_file.replace('.gb','.fasta'),format='fasta')
rec_df[['Locus','gene','mol_type', 'Len',
          'sci_name', 'kew_id','family', 'genus', 'species', 'infraspecies', 'Duplicates',
          'Ini_sci_name', 'TaxID']].to_csv(gb_file.replace('.gb','_TAXO.csv'),index=False)


# In[70]:


print(rec_df.Len.quantile([.01,.05,.1,0.5,.9,.95,.99]).to_dict())
print(rec_df.Len.median()+(rec_df.Len.std()*2))
print(rec_df.Len.median()-(rec_df.Len.std()*2))
rec_df.Len.hist(bins=50);


# In[71]:


rec_rm_df = pd.DataFrame(rec_rm)
print(rec_rm_df.groupby('type').size().to_dict())

