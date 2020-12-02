#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from tqdm import tqdm
import numpy as np
import os
import argparse


# ## Genes

# In[2]:


genes_df = pd.DataFrame(list(zip(
    ['BOLD_pln_rbcL','BOLD_pln_rbcLa','BOLD_pln_matK',
     'NCBI_pln_trnL', 'NCBI_pln_trnH',
     'NCBI_pln_18s','NCBI_pln_plastome'], 
    [400,400,400,200,100,1000,1000],
    [90,90,90,90,90,90,0])),
    columns =['gene', 'min_len','min_cov'])
genes_df


# ## Parameters

# In[ ]:


parser = argparse.ArgumentParser(
    description='Script used to concatenate, filter and merge blast output with taxo')
parser.add_argument("project_folder", type=str, help="Name of project, folder")
parser.add_argument("samples_df_filename", type=str, help="path to sample list in the project folder")
args = parser.parse_args()

project_folder = args.project_folder
samples_df_filename = args.samples_df_filename
col_taxo_db=['SeqID','species','genus','family']


# In[3]:


# project_folder='OKP'
# samples_df_filename='OKP_list.txt'
# col_taxo_db=['SeqID','species','genus','family']


# ## Functions

# In[4]:


def clean_sseqid(col_sseqid):
    try:
        col_sseqid2 = col_sseqid.str.split('|',expand=True)[1] # format is e.g. gb|KY652173.1|
        col_sseqid2.loc[col_sseqid2.isna()==True]=col_sseqid[col_sseqid2.isna()==True]
    except:
#         print('no sseqid format change needed')
        col_sseqid2 = col_sseqid
    return col_sseqid2
# tmp=pd.DataFrame(['gb|KY652173.1|','test.1','gb|KY652173.1','KY|652.1'],columns=['test'])
# print(tmp)
# print(clean_sseqid(col_sseqid=tmp.test))


# In[5]:


# Load blast output
def load_blast_file(Sample,blastpath,gene):
    filepath=blastpath + '/' + str(Sample) + '_' + gene + '.out'
    if os.path.exists(filepath): 
        if os.stat(filepath).st_size > 0:
            blast_df=pd.read_csv(filepath, header=None, sep='\t')
            blast_df.columns=['qseqid', 'sseqid', 'pident', 'length', 'slen', 'qlen', 'mismatch', 'gapopen', 'qstart', 
                              'qend', 'sstart', 'send', 'evalue', 'bitscore']
            
            blast_df['Sample']=Sample
            blast_df['scov']=round(blast_df.length/blast_df.slen*100,0)
            blast_df['qcov']=round(blast_df.length/blast_df.qlen*100,0)
            blast_df = blast_df[['Sample','qseqid', 'sseqid',
                        'length', 'slen', 'qlen','scov','qcov',
                        'pident','evalue', 'bitscore']]      
        else:
            return None
    if os.path.exists(filepath)==False:
        return None
    return blast_df


# In[6]:


# Load samples
def load_smpl(smpl_path):
    samples_df=pd.read_csv(smpl_path)
    if ('Sample' not in samples_df.columns):
        print('Could not find sample column. Will read it again as without header and name               first column Sample')
        samples_df=pd.read_csv(smpl_path,header=None).rename(columns={0:'Sample'})
    return samples_df
# samples_df = load_smpl(smpl_path = project_folder + '/' + samples_df_filename)


# In[7]:


# Filter blast output
def filter_blast(blast_filt_df, filter_dict):
    blast_filt_df = blast_filt_df[blast_filt_df.pident>=filter_dict['min_pident']]
    blast_filt_df = blast_filt_df[blast_filt_df.length>=filter_dict['min_length']]
    blast_filt_df = blast_filt_df[blast_filt_df.scov>=filter_dict['min_scov']]
    
    return blast_filt_df


# ## Main

# In[8]:


if __name__ == "__main__":
    print('\n\nProcessing blast output for',genes_df.shape[0],'barcode genes')
    ## Load data
    samples_df = load_smpl(smpl_path = project_folder + '/' + samples_df_filename)
    print(samples_df.shape[0],'samples')
    
    for gene_idx, gene_row in genes_df.iterrows():
        print('\n\nConcatenating blast output of',gene_row.gene)
        ## Load db_taxo
        taxo_db=pd.read_csv('Barcode_DB/' + gene_row.gene + '_TAXO.csv')
        taxo_db = taxo_db[col_taxo_db]
        taxo_db['SeqID']= taxo_db['SeqID'].astype('str')
        print(taxo_db.shape[0],'entries in taxo_db')
        
        filter_dict={'min_pident':95,'min_length':gene_row.min_len,'min_scov':gene_row.min_cov}
        print(filter_dict)
    
        ## Process blast outputs
        blast_taxo=pd.DataFrame(columns=['Sample','qseqid','sseqid','length','slen',
                                         'qlen','scov','qcov','pident','evalue','bitscore',
                                         'SeqID','species','genus','family'])
        blast_taxo.to_csv(project_folder + '/OutBlast_' + gene_row.gene + '.csv',
                                    index=False)
        for idx, row in tqdm(samples_df.iterrows(), total=samples_df.shape[0]):
            # Get blast output
            raw_blast=load_blast_file(Sample = row.Sample, blastpath = project_folder + '/out_blast',
                           gene = gene_row.gene)
            

            # If data, filter blast
            if isinstance(raw_blast,type(None))==False:
                blast_filt_df=filter_blast(raw_blast, filter_dict)
                blast_filt_df['sseqid']= blast_filt_df['sseqid'].astype('str')

                # If matches passed filtering, get taxonomy
                if isinstance(raw_blast,type(None))==False:
                    # Clean sseqid by recovering the element within |
                    blast_filt_df['sseqid']=clean_sseqid(col_sseqid=blast_filt_df['sseqid'])
                    # Make subset to speed up merging
                    taxo_db_subset = taxo_db[taxo_db.SeqID.isin(blast_filt_df.sseqid)]
                    blast_taxo = pd.merge(left=blast_filt_df,right=taxo_db_subset,
                                          how='left',left_on='sseqid',right_on='SeqID')
                    #Reset index
                    blast_taxo = blast_taxo.reset_index().drop(columns=['index'])
                    #Reducing size of qseqid
                    blast_taxo['qseqid']=blast_taxo['qseqid'].str.slice(0,12)
                    #Write output
                    blast_taxo.to_csv(project_folder + '/OutBlast_' + gene_row.gene + '.csv', 
                                mode='a', index=False, header=False)
                    if blast_taxo.SeqID.isna().sum()>0:
                        print('WARNING:',blast_taxo.SeqID.isna().sum(),
                              'have no match between sseqid and SeqID')
                        print(blast_taxo[blast_taxo.SeqID.isna()]['sseqid'])
                    #Append to DF
#                     blast_df = blast_df.append(blast_taxo, True)
        # Info on output
#         print('Found',blast_df.shape[0],'entries for',blast_df.Sample.nunique(),'Samples')
        # Write 
        print('wrote output as ' + project_folder + '/OutBlast_' + gene_row.gene + '.csv')


# In[21]:


df=pd.DataFrame([70,100,90,50,50,20,100, 95],columns=['score'])
df['rank']=df.score.rank(ascending=False,method='min')
df.sort_values('rank')

