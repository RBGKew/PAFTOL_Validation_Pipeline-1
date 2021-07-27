#!/usr/bin/env python
# coding: utf-8

# # GetOrganelles Cleaning


import os
import shutil
import sys
path=os.getcwd()
file_path=sys.argv[1]
project_name = file_path.split('/')[0]
file_name = file_path.split('/')[1]
# path = path + '\\' + project_name
path = path + '/' + project_name
os.chdir(path)
print('current path:',path)


#Aquiring stats
stats_pt_file = 'stats_pt_' + project_name + '.txt'
stats_nr_file = 'stats_nr_' + project_name + '.txt'
#stats_du_file = 'stats_du_' + project_name + '.txt'
os.system('module load seqkit')
os.system('seqkit stats -j 4 *_pt/embplant_pt*.fasta > ' + stats_pt_file)
os.system('seqkit stats -j 4 *_nr/embplant_nr*.fasta > ' + stats_nr_file)
#os.system('du -a > ' + stats_du_file)



import pandas as pd
import numpy as np


# ## 1. Check logs

getorg_df = pd.read_csv(file_name)
if 'ExternalSequenceID' in getorg_df.columns:
    Dsource='SRA'
    getorg_df=getorg_df.rename(columns={'ExternalSequenceID':'Sample'})
else:
    Dsource='PAFTOL'

print(getorg_df.shape[0],'samples in .csv file')


# Check logs
logs_df = pd.DataFrame(os.listdir('logs/'),columns=['file'])
logs_df['Sample'] = logs_df.file.str.replace('log_','')
logs_df['Type'] = logs_df.Sample.str.split('.',expand=True)[1]
logs_df['Sample'] = logs_df.Sample.str.split('.',expand=True)[0]
if Dsource=='SRA':
    logs_df['Organelle'] = logs_df.Sample.str.split('_',expand=True)[1]
else:
    logs_df['Organelle'] = logs_df.Sample.str.split('_',expand=True)[2]

logs_df['Sample'] = logs_df.Sample.str.replace('_nr','').str.replace('_pt','')
logs_df = logs_df[logs_df.Organelle.isna()==False]
for idx, row in logs_df.iterrows():
    logs_df.loc[idx,'filesize']=os.stat('logs/' + logs_df.loc[idx,'file']).st_size

getorg_df['log']=getorg_df.Sample.isin(logs_df[logs_df.Type=='log']['Sample'])
print(getorg_df.log.sum(),'/',getorg_df.shape[0],'samples processed')
print('continuing with processed samples only')
getorg_df = getorg_df[getorg_df.log]

getorg_df['error_pt']=pd.merge(getorg_df,logs_df[(logs_df.Type=='err') & (logs_df.Organelle=='pt')]).filesize > 0
getorg_df['error_nr']=pd.merge(getorg_df,logs_df[(logs_df.Type=='err') & (logs_df.Organelle=='nr')]).filesize > 0
print(getorg_df.error_pt.sum(),'/',getorg_df.shape[0],'error during plastome recovery')
print(getorg_df.error_nr.sum(),'/',getorg_df.shape[0],'error during ribosomeDNA recovery')


# ## 2. Recovery stats

print(' ')
print('Stats on recovered organelles:')

def proc_seqkit_stats(stats_file):
    seq_df=pd.read_table(stats_file, sep=' +')
    seq_df['Sample']=seq_df.file.str.split('/',expand=True)[0]
    seq_df['sum_len'] = seq_df['sum_len'].str.replace(',','').astype(int)
    seq_df['min_len'] = seq_df['min_len'].str.replace(',','').astype(int)
    seq_df['avg_len'] = seq_df['avg_len'].str.replace(',','').astype(float)
    seq_df['max_len'] = seq_df['max_len'].str.replace(',','').astype(int)
    seq_df['source'] = seq_df.file.str.split('/embplant_',expand=True)[1].str.split('.',expand=True)[0]
    seq_df['K'] = seq_df.file.str.split('.K',expand=True)[1].str.split('.',expand=True)[0]#.astype(int)
    seq_df['assembly'] = seq_df.file.str.split('.graph1',expand=True)[0].str.split('.',expand=True)[2]
    seq_df['repeat_pattern'] = seq_df.file.str.contains(pat = 'repeat_pattern') 
    seq_df['Sample']=seq_df.Sample.str.replace('_' + seq_df.loc[0,'source'],'')
#     seq_df = seq_df.drop(columns=['format','type']).rename(columns={'file':'fasta_file_' + org, 'num_seqs':'num_seqs_' + org,
#          'sum_len':'sum_len_' + org,'min_len':'min_len_' + org,'avg_len':'avg_len_' + org,'max_len':'max_len_' + org})
    return seq_df


print('getting pt stats:')
seq_pt = proc_seqkit_stats(stats_pt_file)
print('found pt fasta for',seq_pt.Sample.nunique(),'Sample')
print((seq_pt.groupby('Sample').size()>1).sum(),'samples have more than 1 file')
print((seq_pt.groupby('Sample').size()>2).sum(),'samples have more than 2 files')
print('Getting the longest sum of contig length')
seq_pt_samp=seq_pt.sort_values('sum_len',ascending=False).groupby('Sample').first().sort_values('Sample').reset_index()
print('SumContigLength:',end='\n\n')
print(seq_pt_samp.sum_len.describe().astype(int),end='\n\n')

print(' ')
print('getting nr stats:')
seq_nr = proc_seqkit_stats(stats_nr_file)
print('found nr fasta for',seq_nr.Sample.nunique(),'Sample')
print((seq_nr.groupby('Sample').size()>1).sum(),'samples have more than 1 file')
print((seq_nr.groupby('Sample').size()>2).sum(),'samples have more than 2 files')
print('Getting the longest sum of contig length')
seq_nr_samp=seq_nr.sort_values('sum_len',ascending=False).groupby('Sample').first().sort_values('Sample').reset_index()
print('SumContigLength:',end='\n\n')
print(seq_nr_samp.sum_len.describe().astype(int),end='\n\n')


print(' ')
getorg_df['pt_exist'] = getorg_df.Sample.isin(seq_pt_samp.Sample)
getorg_df['nr_exist'] = getorg_df.Sample.isin(seq_nr_samp.Sample)
print(getorg_df.pt_exist.sum(),'/',getorg_df.shape[0],'plastomes recovered')
print(getorg_df.nr_exist.sum(),'/',getorg_df.shape[0],'ribosomeDNA recovered')


# ## 3. Clean folders

# to_clean
def clean_org(to_clean, iorg, rm_dirs, rm_files,verbose=True):
	for idx, row in to_clean.iterrows():
			for idir in rm_dirs:
				rm_dir = row.Sample + iorg + '/' + idir
				if os.path.isdir(rm_dir):
					print('removing:', rm_dir)
					shutil.rmtree(rm_dir)
			for ifile in rm_files:
				rm_file = row.Sample + iorg + '/' + ifile
				if os.path.isfile(rm_file):
					print('removing:', rm_file)
					os.remove(rm_file)

rm_dirs=['filtered_spades','seed']
rm_files=['filtered_1_paired.fq.tar.gz','filtered_2_paired.fq.tar.gz',
          'filtered_1_unpaired.fq.tar.gz','filtered_2_unpaired.fq.tar.gz',
          'filtered_1.fq.tar.gz']
		
print(' ')
to_clean_pt = getorg_df[getorg_df.pt_exist]
print('cleaning pt folders for',to_clean_pt.shape[0],'samples')
clean_org(to_clean_pt, '_pt', rm_dirs, rm_files)
to_clean_nr = getorg_df[getorg_df.nr_exist]
print('\ncleaning nr folders for',to_clean_nr.shape[0],'samples')
clean_org(to_clean_nr, '_nr', rm_dirs, rm_files)

# Remove fastq files if pt and nr recovered
to_clean_fastq = getorg_df[(getorg_df.pt_exist) & (getorg_df.nr_exist)]
print('\nRemoving fastq for',to_clean_fastq.shape[0],'samples')
for idx, row in to_clean_fastq.iterrows():
	rm_files=['_1.fastq','_2.fastq','_R1.fastq','_R2.fastq','.fastq']
	for ifile in rm_files:
		rm_file = 'fastq_data/' + row.Sample + ifile
		if os.path.isfile(rm_file):
			print('removing:', rm_file)
			os.remove(rm_file)

# ## 4. Copy files
print(' ')
print('copying files and creating tables')
if os.path.isdir('fasta_pt')==False:
	os.mkdir('fasta_pt')
if os.path.isdir('fasta_nr')==False:
	os.mkdir('fasta_nr')

for idx, row in seq_pt_samp.iterrows():
	shutil.copyfile(row.file, 'fasta_pt/' + row.Sample + '_pt.fasta')
for idx, row in seq_nr_samp.iterrows():
	shutil.copyfile(row.file, 'fasta_nr/' + row.Sample + '_nr.fasta')

seq_pt_samp.to_csv(project_name + '_GetOrg_pt_Results.csv',index=False)
seq_nr_samp.to_csv(project_name + '_GetOrg_nr_Results.csv',index=False)
