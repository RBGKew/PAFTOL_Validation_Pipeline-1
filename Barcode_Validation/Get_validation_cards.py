#!/usr/bin/env python
# coding: utf-8

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

# In[65]:


import pandas as pd
from tqdm import tqdm
import numpy as np
import os
import argparse
import warnings
warnings.filterwarnings('ignore')


# ## Parameters

# In[66]:


## Parameters
parser = argparse.ArgumentParser(
    description='Blast sample sequences on Barcode database and process the results')
parser.add_argument("--samples_file", type=str, help="spreadsheet of samples with their taxonomy")
parser.add_argument("--sample", type=str, help="sample for which a barcode validation will be produced")
parser.add_argument("--barcodes_table", type=str, help="spreadsheet of barcode tests with parameters")
args = parser.parse_args()

samples_file = args.samples_file
sample = args.sample
barcode_tests_file = args.barcodes_table


# In[67]:


# ## Test in Notebook
# # os.chdir('PAFTOL/')
# samples_file = 'PAFTOL_samples.csv'
# sample = 'PAFTOL_007573'
# barcode_tests_file = '../Barcode_DB/Barcode_Tests.csv'


# In[68]:


barcode_DB_dir = os.path.split(barcode_tests_file)[0] +'/'
col_taxo_db=['Locus','species','genus','family']
taxo_ranks=['genus','family']
val_col_order = ['Test','tax_level','taxo','taxo_in_db','Blast','Nmatch','match','rank_pid','rank_bsc','pid','len','scov','qcov',
                 'best','best_pid','best_score','best_scov','best_qcov','NseqID']


# ## Functions

# In[69]:


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


# In[70]:


# Load blast output
def load_blast_file(blastpath):
    if os.path.exists(blastpath): 
        if os.stat(blastpath).st_size > 0:
            blast_df=pd.read_csv(blastpath, header=None, sep='\t')
            blast_df.columns=['qseqid', 'sseqid', 'pident', 'length', 'slen', 'qlen', 'mismatch', 'gapopen', 'qstart', 
                              'qend', 'sstart', 'send', 'evalue', 'bitscore']
            blast_df.pident = blast_df.pident.round(2)
            blast_df['scov']=round( (abs(blast_df.sstart-blast_df.send) +1) / blast_df.slen*100,1)
            blast_df['qcov']=round( (abs(blast_df.qstart-blast_df.qend) +1) / blast_df.qlen*100,1)
            blast_df = blast_df[['qseqid', 'sseqid','length', 'slen', 'qlen','scov','qcov','pident','evalue', 'bitscore']]      
        else:
            return None
    if os.path.exists(blastpath)==False:
        return None
    return blast_df


# In[71]:


# Filter blast output
def filter_blast(blast_filt_df, filter_dict):
    blast_filt_df = blast_filt_df[blast_filt_df.pident>=filter_dict['min_pident']]
    blast_filt_df = blast_filt_df[blast_filt_df.length>=filter_dict['min_length']]
    blast_filt_df = blast_filt_df[blast_filt_df.scov>=filter_dict['min_scov']]
    return blast_filt_df


# In[72]:


def get_blast_results(validic, blast_sample):
    # Calculate ranks
    blast_sample['rank_pid']=blast_sample.pident.rank(method='min',ascending=False)
    blast_sample['rank_bsc']=blast_sample.bitscore.rank(method='min',ascending=False)

    # Validation data
    validic['match'] = validic['taxo'] in list(blast_sample[validic['tax_level']])

    # if there is a match, collect rank, pc identity, matching length and coverage
    if validic['match']==True: 
        idx_tax = blast_sample[blast_sample[validic['tax_level']] == validic['taxo']].index
        validic['rank_pid'] = blast_sample.loc[idx_tax[0],'rank_pid']
        validic['rank_bsc'] = blast_sample.loc[idx_tax[0],'rank_bsc']
        validic['pid'] = blast_sample.loc[idx_tax[0],'pident']
        validic['len'] = blast_sample.loc[idx_tax[0],'length']
        validic['scov'] = blast_sample.loc[idx_tax[0],'scov']
        validic['qcov'] = blast_sample.loc[idx_tax[0],'qcov']

    # Collect best hit and general stats
    validic['best'] = blast_sample.loc[0,validic['tax_level']]
    validic['best_pid'] = blast_sample.loc[0,'pident']
    validic['best_score'] = blast_sample.loc[0,'bitscore']
    validic['best_scov'] = blast_sample.loc[0,'scov']
    validic['best_qcov'] = blast_sample.loc[0,'qcov']
    validic['Nmatch'] = blast_sample.shape[0]
    validic['NseqID'] = blast_sample.sseqid.nunique()
    return validic


# In[73]:


def load_taxo_db(genes_df):
    all_taxo_db = pd.DataFrame()
    for gene_idx, gene_row in genes_df.iterrows():
        ## Load db_taxo
        taxo_db = pd.read_csv(barcode_DB_dir + gene_row.Barcode + '_TAXO.csv')
        taxo_db = taxo_db[col_taxo_db]; taxo_db['Locus']= taxo_db['Locus'].astype('str')
        taxo_db['Barcode'] = gene_row.Barcode
        all_taxo_db = pd.concat([all_taxo_db,taxo_db],ignore_index=True)
        del taxo_db
    return all_taxo_db


# ## Main

# In[74]:


samples_df = pd.read_csv(samples_file)
sample_dic = samples_df[samples_df.Sample==sample].to_dict(orient='records')[0]
sample_dic['genus']


# In[75]:


if __name__ == "__main__":
    ## Load data
    genes_df = pd.read_csv(barcode_tests_file)
    samples_df = pd.read_csv(samples_file)
    # Convert sample info (taxonomy) as a dictionary
    sample_dic = samples_df[samples_df.Sample==sample].to_dict(orient='records')[0]
    print('\n\nProcessing blast output for',genes_df.shape[0],'barcode tests')
    print(genes_df)

    # Don't run sample's validation if it already exists
    validation_file = 'Barcode_Validation/BV_' + sample + '.csv'
    results_blast_df = pd.DataFrame()
    ## Load all db_taxo
    all_taxo_db = load_taxo_db(genes_df)
    # For each gene,
    for gene_idx, gene_row in genes_df.iterrows():
        # barcode db_taxo
        taxo_db = all_taxo_db[all_taxo_db.Barcode==gene_row.Barcode].reset_index(drop=True)
        # For each taxonomic level (genus, family)
        for itax in taxo_ranks:
            validic = {'Test':gene_row.Barcode, 'tax_level': itax, 'taxo': sample_dic[itax]}
            
            ## If taxo in DB
            if validic['taxo'] in list(taxo_db[itax]):
                validic['taxo_in_db'] = True
                
                # Get blast output
                raw_blast=load_blast_file(blastpath = 'out_blast/' + sample + '-' + gene_row.Barcode + '.out')
                
                ## If blast output, filter it
                if isinstance(raw_blast,type(None))==False:
                    validic['Blast'] = True
                    filter_dict={'min_pident':gene_row.blast_pid,'min_length':gene_row.min_len,'min_scov':gene_row.min_cov}
                    blast_filt_df=filter_blast(raw_blast, filter_dict)
                    blast_filt_df['sseqid']= blast_filt_df['sseqid'].astype('str')
                    
                    ## If match after filtering, get validation results
                    if blast_filt_df.shape[0]>0:
                        validic['Nmatch'] = blast_filt_df.shape[0]
                        
                        # Clean sseqid by recovering the element within |
                        blast_filt_df['sseqid']=clean_sseqid(col_sseqid=blast_filt_df['sseqid'])
                        # Make subset to speed up merging
                        taxo_db_subset = taxo_db[taxo_db.Locus.isin(blast_filt_df.sseqid)]
                        blast_taxo = pd.merge(left=blast_filt_df,right=taxo_db_subset,
                                              how='left',left_on='sseqid',right_on='Locus')
                        #Sort values and Reset index
                        blast_taxo = blast_taxo.sort_values('pident',ascending=False).reset_index(drop=True)
                        #Reducing length of qseqid
                        blast_taxo['qseqid']=blast_taxo['qseqid'].str.slice(0,12)
                        if blast_taxo.Locus.isna().sum()>0:
                            print('WARNING:',blast_taxo.Locus.isna().sum(),
                                  'matches are missing their taxonomy', gene_row.Barcode, itax)
                            print(blast_taxo[blast_taxo.Locus.isna()]['sseqid'].to_dict())
                            
                        validic = get_blast_results(validic = validic, blast_sample = blast_taxo)
                        
                    ## If NO match after filtering
                    else:
                        validic['Nmatch'] = 0
                        
                ## If NO blast output
                else:
                    validic['Blast'] = False      
                    
            ## If taxo NOT in DB        
            elif validic['taxo'] not in list(taxo_db[itax]):
                validic['taxo_in_db'] = False
                
            results_blast_df = pd.concat([results_blast_df, pd.DataFrame.from_dict(validic,orient='index').transpose()])

    # Write Output if at least one blast file was found
    if results_blast_df.Blast.sum()>0:
        tmp_val_col = [icol for icol in val_col_order if icol in results_blast_df.columns]
        results_blast_df[tmp_val_col].to_csv(validation_file,index=False)
    else:
        print('No Blast files found for',sample)

