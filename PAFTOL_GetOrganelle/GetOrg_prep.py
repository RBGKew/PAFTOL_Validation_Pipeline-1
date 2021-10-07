#!/usr/bin/env python
# coding: utf-8

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

# In[1]:


import pandas as pd
import os; import sys
import argparse


# In[2]:


## Parameters
parser = argparse.ArgumentParser(
    description='Blast sample sequences on Barcode database and process the results')
parser.add_argument("--db", type=str, help="latest paftol_export")
parser.add_argument("--DataSource", type=str, help="DataSource (e.g. PAFTOL, SRA)")
parser.add_argument("--rem_search", type=str, help="List completed samples if fasta or log exist")
args = parser.parse_args()

export_file = args.db
DataSource = args.DataSource
rem_search = args.rem_search


# In[3]:


# # For notebook only
# export_file = '../PAFTOL_DB/2021-07-27_paftol_export.csv'
# DataSource = 'GAP'
# rem_search = 'log'


# In[4]:


# Load export for datasource
db = pd.read_csv(export_file)
db = db[(db.DataSource==DataSource) & (db.R1FastqFile.notnull())]
if DataSource == 'PAFTOL':
    fastq_path = '/science/projects/paftol/AllData_symlinks/'
    db['Sample_Name'] = 'PAFTOL_' + db['idSequencing'].astype(int).astype('str').str.zfill(6)
    db['R1_path'] = fastq_path + db.Sample_Name + '_R1.fastq.gz'
    db['R2_path'] = fastq_path + db.Sample_Name + '_R2.fastq.gz'
elif DataSource == 'GAP':
    fastq_path = '/science/projects/paftol/AllData_symlinks/'
    db['Sample_Name'] = 'GAP_' + db['idSequencing'].astype(int).astype('str').str.zfill(6)
    db['R1_path'] = fastq_path + db.Sample_Name + '_R1.fastq.gz'
    db['R2_path'] = fastq_path + db.Sample_Name + '_R2.fastq.gz'
elif DataSource == 'SRA':
    fastq_path = '/data/projects/paftol/SRA_Data/'
    db['Sample_Name'] = db.ExternalSequenceID
    db['R1_path'] = fastq_path + db.R1FastqFile
    db['R2_path'] = fastq_path + db.R2FastqFile
else:
    print('unknown action for',DataSource)
    sys.exit()

print(db.shape[0],'samples in total')


# In[5]:


# List fasta_pt and fasta_nr
db['fasta_pt']=False; db['fasta_nr']=False;
fasta_pt = pd.DataFrame(os.listdir(DataSource + '/fasta_pt/'),columns=['file'])
if fasta_pt.shape[0]>0:
    fasta_pt['Sample_Name'] = fasta_pt.file.str.split('_pt',expand=True)[0]
    db['fasta_pt']=db.Sample_Name.isin(fasta_pt.Sample_Name)
fasta_nr = pd.DataFrame(os.listdir(DataSource + '/fasta_nr/'),columns=['file'])
if fasta_nr.shape[0]>0:
    fasta_nr['Sample_Name'] = fasta_nr.file.str.split('_nr',expand=True)[0]
    db['fasta_nr']=db.Sample_Name.isin(fasta_nr.Sample_Name)
print(db.fasta_pt.sum(),'/',db.shape[0],'pt recovered')
print(db.fasta_nr.sum(),'/',db.shape[0],'nr recovered')


# In[6]:


# Check logs
db['log_pt']=False; db['log_nr']=False;
logs_df = pd.DataFrame(os.listdir(DataSource + '/logs/'),columns=['file'])
if logs_df.shape[0]>0:
    logs_df['Sample_Name'] = logs_df.file.str.split('log_',expand=True)[1]
    logs_df['Type'] = logs_df.Sample_Name.str.split('.',expand=True)[1]
    logs_df['Sample_Name'] = logs_df.Sample_Name.str.split('.',expand=True)[0]
    logs_df['Organelle'] = logs_df.Sample_Name.str.split('_').str[-1]
    logs_df['Sample_Name'] = logs_df.Sample_Name.str.replace('_nr','').str.replace('_pt','')
    logs_df['filesize']=logs_df.file.apply(lambda x: os.stat(DataSource + '/logs/' + x).st_size)
    db['log_pt']=db.Sample_Name.isin(logs_df[(logs_df.Type=='log') & (logs_df.Organelle=='pt')]['Sample_Name'])
    db['log_nr']=db.Sample_Name.isin(logs_df[(logs_df.Type=='log') & (logs_df.Organelle=='nr')]['Sample_Name'])
print(db.log_pt.sum(),'/',db.shape[0],'pt processed')
print(db.log_nr.sum(),'/',db.shape[0],'nr processed')
if logs_df.shape[0]>0:
    db['error_pt']=pd.merge(db, logs_df[(logs_df.Type=='err') & (logs_df.Organelle=='pt')]).filesize > 0
    db['error_nr']=pd.merge(db, logs_df[(logs_df.Type=='err') & (logs_df.Organelle=='nr')]).filesize > 0
    print(db.error_pt.sum(),'/',db.shape[0],'error during pt recovery')
    print(db.error_nr.sum(),'/',db.shape[0],'error during nr recovery')


# In[7]:


if rem_search == 'fasta':
    todo_pt = db[(db.fasta_pt==False)][['Sample_Name','R1_path','R2_path']]
    todo_nr = db[(db.fasta_nr==False)][['Sample_Name','R1_path','R2_path']]
elif rem_search == 'log':
    todo_pt = db[(db.log_pt==False)][['Sample_Name','R1_path','R2_path']]
    todo_nr = db[(db.log_nr==False)][['Sample_Name','R1_path','R2_path']]
if todo_pt.shape[0]>0:
    print('\n',todo_pt.shape[0],DataSource,'samples listed for pt recovery')
if todo_nr.shape[0]>0:
    print('\n',todo_nr.shape[0],DataSource,'samples listed for nr recovery')


# In[ ]:


if todo_pt.shape[0]>0:
    for idx, row in todo_pt.iterrows():
#         print(row['R1_path'],end=':'); print(os.path.exists(row['R1_path']))
        todo_pt.loc[idx,'R1_exist'] = os.path.exists(row['R1_path'])
#         print(row['R2_path'],end=':'); print(os.path.exists(row['R2_path']))
        todo_pt.loc[idx,'R2_exist'] = os.path.exists(row['R2_path'])
todo_pt = todo_pt[(todo_pt.R1_exist) & (todo_pt.R2_exist)]
if todo_pt.shape[0]>0:
    print(todo_pt.shape[0],'paired-end fastq files found')
    todo_pt[['Sample_Name','R1_path','R2_path']].to_csv(DataSource + '/remaining_pt.txt',index=False,header=None)
else:
    print('no fastq file found or no sample to process, remaining list not written')

if todo_nr.shape[0]>0:
    for idx, row in todo_nr.iterrows():
#         print(row['R1_path'],end=':'); print(os.path.exists(row['R1_path']))
        todo_nr.loc[idx,'R1_exist'] = os.path.exists(row['R1_path'])
#         print(row['R2_path'],end=':'); print(os.path.exists(row['R2_path']))
        todo_nr.loc[idx,'R2_exist'] = os.path.exists(row['R2_path'])
todo_nr = todo_nr[(todo_nr.R1_exist) & (todo_nr.R2_exist)]
if todo_nr.shape[0]>0:
    print(todo_nr.shape[0],'paired-end fastq files found')
    todo_nr[['Sample_Name','R1_path','R2_path']].to_csv(DataSource + '/remaining_nr.txt',index=False,header=None)
else:
    print('no fastq file found or no sample to process, remaining list not written')


# In[ ]:


# if todo_pt.shape[0]>0:
#     todo_pt['R1_exist'] = todo_pt.apply(lambda row: os.path.exists(row['R1_path']),axis=1)
#     todo_pt['R2_exist'] = todo_pt.apply(lambda row: os.path.exists(row['R2_path']),axis=1)
# if todo_nr.shape[0]>0:
#     todo_nr['R1_exist'] = todo_nr.apply(lambda row: os.path.exists(row['R1_path']),axis=1)
#     todo_nr['R2_exist'] = todo_nr.apply(lambda row: os.path.exists(row['R2_path']),axis=1)

# if todo_pt.shape[0]>0:
#     todo_pt = todo_pt[(todo_pt.R1_exist) & (todo_pt.R2_exist)]
#     print(todo_pt.shape[0],'paired-end fastq files found')
#     todo_pt[['Sample_Name','R1_path','R2_path']].to_csv(DataSource + '/remaining_pt.txt',index=False,header=None)

# if todo_nr.shape[0]>0:
#     print('\n',todo_nr.shape[0],DataSource,'samples listed for nr recovery')
#     todo_nr = todo_nr[(todo_nr.R1_exist) & (todo_nr.R2_exist)]
#     print(todo_nr.shape[0],'paired-end fastq files found')
#     todo_nr[['Sample_Name','R1_path','R2_path']].to_csv(DataSource + '/remaining_nr.txt',index=False,header=None)

