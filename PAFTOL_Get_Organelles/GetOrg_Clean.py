#!/usr/bin/env python
# coding: utf-8


# # GetOrganelles Cleaning
##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

# In[1]:


import os
import shutil
import sys
import argparse
from Bio import SeqIO


# In[ ]:


parser = argparse.ArgumentParser(
    description='Process GetOrganelles results and clean tmp folders')
parser.add_argument("--path", type=str, help="path of folder to process")
args = parser.parse_args()
path = args.path
print('path to folder:',path)
Sample = path.split('/')[-2].replace('_pt','').replace('_nr','')
print('Sample:',Sample)


# In[32]:


org = path.split('/')[-2].split('_')[-1]
if org in ['pt','nr']:
    print('organelle:',org)
else:
    print('wrong parsing of organelle:',org)
    sys.exit()


# ## Copy best organelle fasta

# In[14]:


fasta_files = [ifile for ifile in os.listdir(path) if ifile.endswith('.fasta')]
if len(fasta_files)==1:
    print('1 fasta file:',fasta_files[0])
    shutil.copyfile(path + fasta_files[0], 'fasta_' + org + '/' + Sample + '_' + org + '.fasta')
elif len(fasta_files)>1:
    print('found',len(fasta_files),'fasta files')
    best_fasta=''
    best_len=0
    for ifasta in fasta_files:
        sum_len=0
        for record in SeqIO.parse(path + ifasta, "fasta"):
            sum_len += len(record.seq)
        if sum_len>best_len:
            best_len=sum_len
            best_fasta=ifasta
    shutil.copyfile(path + best_fasta, 'fasta_' + org + '/' + Sample + '_' + org + '.fasta')
else:
    print('either no fasta or error, exiting.')
    sys.exit()


# ## Remove temp files

# In[7]:


rm_dirs=['filtered_spades','seed']
rm_files=['filtered_1_paired.fq.tar.gz','filtered_2_paired.fq.tar.gz',
          'filtered_1_unpaired.fq.tar.gz','filtered_2_unpaired.fq.tar.gz']
for idir in rm_dirs:
    rm_dir = path + idir
    print(rm_dir)
    if os.path.isdir(rm_dir):
        shutil.rmtree(rm_dir)
for ifile in rm_files:
    rm_file = path + ifile
    print(rm_file)
    if os.path.isfile(rm_file):
        os.remove(rm_file)


# ## Compress and remove folder

# In[ ]:


zip_path='Archives/' + Sample + '_' + org + '.gz'
os.system('tar -zcvf ' + zip_path + ' ' + path)
if os.path.isfile(zip_path):
    print('compressed succesfully to',zip_path,', removing folder',path)
    shutil.rmtree(path) 

