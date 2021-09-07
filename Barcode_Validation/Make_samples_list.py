#!/usr/bin/env python
# coding: utf-8

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

# In[1]:


import pandas as pd
import argparse; import os


# In[2]:


parser = argparse.ArgumentParser()
parser.add_argument("--db"); parser.add_argument("--DataSource")
opts = parser.parse_args()
db_export_file = opts.db;  DataSource = opts.DataSource


# In[14]:


# Notebook only
# db_export_file = '../PAFTOL_DB/2021-07-05_paftol_export.csv'
# DataSource = 'OneKP'


# In[21]:


# Load data and keep entries for DataSource
db = pd.read_csv(db_export_file)
db.DataSource.replace({'Annotated genome':'AG','Unannotated genome':'UG'},inplace=True)
db = db[db.DataSource==DataSource].astype({'idSequencing':'int','idPaftol':'int'})
print(db.shape[0],DataSource,'samples')
if db.R1FastqFile.isna().sum()>0:
    print(db.R1FastqFile.isna().sum(),'samples have no R1FastqFile and are removed from further analysis')
    db = db[db.R1FastqFile.notnull()] # Will not run samples that have no fastq files. Note: some are single end


# In[23]:


# Make Sample column
if DataSource in ['OneKP','SRA','UG','AG']:
    db['Sample'] = db['ExternalSequenceID']
elif DataSource in ['PAFTOL']:
    db['Sample'] = db['idSequencing'].astype('str').apply(lambda x: 'PAFTOL_' + x.zfill(6))
elif DataSource in ['GAP']:
    db['Sample'] = db['idSequencing'].astype('str').apply(lambda x: 'GAP_' + x.zfill(6))
else:
    print('could not find Datasource',DataSource)


# In[34]:


# Write samples table
print('\nWrite samples table',DataSource + '/' + DataSource + '_samples.csv')
samples_df = db[['Sample','idSequencing', 'idPaftol', 'Family', 'Genus', 'Species']]    .rename(columns={'Family':'family','Genus':'genus','Species':'species'})
samples_df.to_csv(DataSource + '/' + DataSource + '_samples.csv',index=False)
print('First line:\n',samples_df[:1].to_string(index=False))


# In[48]:


# List existing validation cards and output list of samples to blast
samples_done = [filename.replace('BV_','').replace('.csv','') for filename in os.listdir(DataSource + '/Barcode_Validation/')]
print('\nfound',len(samples_done),'validation cards')
samples_todo = db[db.Sample.isin(samples_done)==False]
print(samples_todo.shape[0],'samples to blast, ',db[db.Sample.isin(samples_done)].shape[0],'samples done')


# In[46]:


samples_todo.Sample.to_csv(DataSource + '/Samples_to_barcode.txt',index=False,header=False)

