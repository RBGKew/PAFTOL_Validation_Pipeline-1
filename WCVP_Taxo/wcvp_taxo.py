#!/usr/bin/env python
# coding: utf-8

# 
# # wcvp_taxo
# wcvp_taxo is a python3 script for matching and resolving scientific names against the WCVP database (https://wcvp.science.kew.org/)
# 
# ## Input
# ### A. Input files
# The script requires two input tables: The WCVP database and a file with species names to match on WCVP
# 1. **WCVP database**: must be downloaded from http://sftp.kew.org/pub/data-repositories/WCVP/. It will be filtered and save by the script in pickle format. If you whish to update the WCVP database, deleted the .pkl file.
# 2. **Sample file**: This spreadsheet must be in **.csv** format and contain at least one column with the scientific names you wish to match in WCVP. By default the script will look for a column named **scientific_name**. Otherwise it will look for a column called **Species**. If the species name is spread in two columns **(Genus, Species)**, the script with recognize it automatically.
# 
# ### B. Parameters
# These parameters are optional and can be accessed with python wcvp_taxo.py -h
# - **-g, --resolve_genus**: Find taxa for scientific names written in genus sp. format
# - **-s, --similar_tax_method**: Find most similar taxa for misspelled taxa. <br>
# Possibles values are: 
# 	- **similarity_genus**: Search for similar scientific name in WCVP assuming genus is correct (fast)
# 	- **similarity**: Search for similar scientific name in WCVP (slow)
# 	- **request_kewmatch**: Search for similar scientific name using kewmatch (online) (ok if less than <200 queries)
# - **-d, --duplicate_action**. Action to take when multiple wcvp entries match the provided scientific_name. <br>
# Possibles values are: 
# 	- **rank**: reduce duplicates by prioritizing accepted > unplaced > synonym > homotypic_synonym  taxonomic status (keep first entry). 
# 	- **divert**: divert duplicates to _duplicates.csv
# 	- **divert_taxonOK**: divert duplicates to _duplicates.csv, unless all matching entries have the same taxon name in WCVP (keep first entry)
# 	- **divert_speciesOK**: divert duplicates to _duplicates.csv, unless all matching entries have the same species name in WCVP (keep first entry)
# 	- **divert_genusOK**: divert duplicates to _duplicates.csv, unless all matching entries have the same genus name in WCVP (keep first entry and rename as genus sp.)
# - **-oc, --only_changes**: Output file only contains IDs that have a different taxonomy than provided (species, genus or family if provided)
# - **-os, --simple_output**: Output file is simplified to 4 columns: ID, kew-id, Ini_sci_name, sci_name
# - **-v, --verbose**: verbose output in console
# 
# 
# ## Example
# ```console
# python wcvp_taxo.py wcvp_export.txt sample_file.csv -g -s similarity_genus -d divert_taxonOK
# python wcvp_taxo.py wcvp_export.txt sample_file.csv
# python wcvp_taxo.py wcvp_export.txt sample_file.csv -oc -os -s similarity --verbose -d divert
# python wcvp_taxo.py wcvp_export.txt sample_file.csv -g -s similarity -d rank --verbose
# ```
# 
# ## Output
# For the example above, the script will output the following tables:
# * **sample_file_wcvp.csv**: Samples for which the scientific name are resolved.
# * **sample_file_duplicates.csv**: Samples for which the scientific name matched multiple WCVP entries.
# * **sample_file_unresolved.csv**: Samples for which the scientific name did not match any WCVP entries.
# 
# 
# ## Pipeline
# ### Pre-processing
# * Load wcvp database. If only text file exist, saving as .pkl.
# * Find column containing scientific names. scientific_name or sci_name (default), Species or Genus + Species otherwise.
# * Search for column with unique IDs. First column in table will be selected. Creates column with unique IDs if it doesn't exist. Will not pick sci_name or Species as ID.
# 
# ### Initial checks
# * Check if Ini_scinames are written as Genus sp.
# * Check if Ini_scinames exist in WCVP
# * Optional. Find similar names if not in WCVP
# * Check if Ini_scinames have duplicate entries
# * Proceed to matching for valid scientific names
# 
# ### Matching & Resolving
# 1. Find accepted and unplaced matches.
# 2. Resolves synonyms and homotypic synonyms. 
# 3. Resolve duplicates.
# 4. Output tables
# 
# ## Dependencies
# pandas, tqdm<br>
# for similarity: difflib, requests, ast<br>
# numpy, os, argparse, sys

# In[1]:


import pandas as pd
from tqdm import tqdm
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'
import os
import argparse
import sys


# ## Parameters

# In[21]:


parser = argparse.ArgumentParser(
    description='Script used to match species names with wcvp. Requires at least the paths to the wcvp\
                file and to a .csv file containing scientific names (species or genus + species)')
parser.add_argument("wcvp_path", type=str, 
                    help="path to wcvp_export.txt, \
                    download from http://sftp.kew.org/pub/data-repositories/WCVP/")
parser.add_argument("df_path", type=str, 
                    help="path to spreadsheet in .csv format. Note output will be in the same folder")
parser.add_argument("-g", "--resolve_genus", 
                    help="Optional. find taxa for scientific names written in genus sp. format", 
                    action="store_true", default=False)
parser.add_argument("-s",'--similar_tax_method', 
        help="Optional. Find most similar taxa for misspelled taxa. possibles values are: \
        similarity_genus, similarity, request_kew", action="store", default=None)
parser.add_argument("-d",'--duplicate_action', 
        help="Optional. Action to take when multiple wcvp taxon match to a sci_name. possibles values are: \
        rank, divert, divert_taxonOK, divert_speciesOK, divert_genusOK.\
        \n\n rank: reduce duplicates by prioritizing accepted > unplaced > synonym > homotypic synonym \
        taxonomic status. \n\n divert: flag duplicates, remove them from _wcvp.csv output and write them to _duplicates.csv", 
                    action="store", default='rank')
parser.add_argument("-oc", "--only_changes", 
                    help="Optional. Output file only contains IDs that have a different taxonomy than provided", 
                    action="store_true", default=False)
parser.add_argument("-od", "--output_duplicates", 
                    help="Optional. Output a separate file for duplicates as _duplicates.csv", 
                    action="store_true", default=False)
parser.add_argument("-os", "--simple_output", 
                    help="Optional. Output file is simplified to 4 columns: ID, kew-id, Ini_sci_name, sci_name", 
                    action="store_true", default=False)
parser.add_argument("-v", "--verbose", 
                    help="Optional. verbose output in console", 
                    action="store_true", default=False)
args = parser.parse_args()

wcvp_path = args.wcvp_path
df_path = args.df_path
resolve_genus=args.resolve_genus
find_most_similar=args.find_most_similar
dupl_action=args.duplicate_action
only_changes=args.only_changes
simple_output=args.simple_output
verbose=args.verbose

status_keep=['Accepted','Unplaced']


# In[13]:


# ## Jupyter Notebook 
# wcvp_path='wcvp_v3_nov_2020.txt'
# # df_path='sample_file.csv'
# # df_path='../paftol_samples/Release_DEC_DF.csv'
# # df_path='../phylotree/Release_1/all_samples_r1.csv'
# # df_path='../db_plants/2020-12-18_paftol_export.csv'
# # df_path='../Barcoding/Barcode_DB_OLDv2/NCBI_pln_trnH_taxo.csv'
# df_path='../Public_Data/public_sciname_list.csv'
# # df_path='../paftol_samples/ENA_submission/paftol_notaxid.csv'
# resolve_genus=False
# find_most_similar='similarity'
# dupl_action='rank'
# verbose=False
# only_changes=False
# simple_output=False
# status_keep=['Accepted','Unplaced']


# ## Functions

# ### Data processing functions

# In[3]:


# Load wcvp file and save as pickle for faster loading
def load_wcvp(wcvp_path):
    print('Loading WCVP...',end='')
    # Load pickel
    if os.path.exists(wcvp_path.replace('.txt','.pkl')):
        print('found .pkl...',end='')
        wcvp = pd.read_pickle(wcvp_path.replace('.txt','.pkl'))
    elif os.path.exists(wcvp_path):
        wcvp = pd.read_table(wcvp_path,sep='|',encoding='utf-8')
        print('found .txt, ',end='')
        # Remove extra columns
        wcvp = wcvp.drop(columns=['parent_kew_id','parent_name','parent_authors'])
        print('saving to .pkl...',end='')
        wcvp.to_pickle(wcvp_path.replace('.txt','.pkl'))
    else:
        print('could not find',wcvp_path)
        sys.exit()
    print(wcvp.shape[0],'entries')
    
    return wcvp

def load_df(df_path):
    print('Loading dataset...',end='')
    try:
        smpl_df = pd.read_csv(df_path,encoding='utf-8')
        print(smpl_df.shape[0],'entries')
        return smpl_df
    except:
        print('could not find',df_path)
        sys.exit()


# In[4]:


#Define ID column
def GetIDcol(df):
    #Check for columns with all unique values
    col_unique=(df.nunique()==df.shape[0]).to_frame().reset_index().rename(columns={0:'unique','index':'column'})
    col_unique = col_unique[col_unique.unique==True]
    col_unique = col_unique[~col_unique['column'].isin(['Ini_sci_name','Ini_Genus','Ini_Species'])]
    if col_unique.shape[0]>0:
        print('found',col_unique.shape[0],'ID column:',end='')
        colsID=list(col_unique['column'])
        colID=colsID[0]
        print(colID)
    else:
        print('No ID column, create ID from index')
        colID='ID'
    #Return new col with ID 
    return colID


# In[5]:


#Find which column contains the scientific name to match
def define_sci_name(smpl_df, verbose=False):
    col_taxo=list(smpl_df.columns[smpl_df.columns.str.contains(
                        'family|genus|species|infraspecies|sci_name|scientific_name',case=False)])
    for itaxo in col_taxo:
        smpl_df = smpl_df.rename(columns = {itaxo:'Ini_' + itaxo.capitalize()})
        if verbose:
            print('renaming ' + itaxo + ' to Ini_' + itaxo, end=', ')
    # Use sci_name if provided
    if 'Ini_Sci_name' in smpl_df.columns:
        print('\nScientific Name is sci_name') 
        smpl_df = smpl_df.rename(columns = {'Ini_Sci_name':'Ini_sci_name'})
    elif 'Ini_Scientific_name' in smpl_df.columns:
        print('\nScientific Name is scientific_name') 
        smpl_df = smpl_df.rename(columns = {'Ini_Scientific_name':'Ini_sci_name'})
    else:
        # Identify is scientific name is in 1 or two columns
        try:
            avg_word_sp=smpl_df['Ini_Species'].str.split().str.len().mean()
            print('avg words in Ini_Species:',round(avg_word_sp,1))
            if round(avg_word_sp)==1:
                print('Scientific Name (Ini_sci_name) is Ini_Genus + Ini_Species')
                smpl_df['Ini_sci_name'] = smpl_df['Ini_Genus'] + ' ' + smpl_df['Ini_Species']
            elif round(avg_word_sp)>=2:
                print('Scientific Name (Ini_sci_name) is Ini_Species')
                smpl_df['Ini_sci_name'] = smpl_df['Ini_Species']
        except:
            print('ERROR: Could not identify species column')
            sys.exit()
    
    return smpl_df


# ### WCVP related functions

# In[6]:


def get_by_taxon_name(df, wcvp):
    tmp_wcvp=wcvp[wcvp.taxon_name.isin(df.sci_name)]
    match = pd.merge(df, tmp_wcvp, how='inner', left_on='sci_name', right_on='taxon_name')
    return match


# In[7]:


def get_by_kew_id(df, wcvp):
    tmp_wcvp=wcvp[wcvp.kew_id.isin(df.kew_id)]
    match = pd.merge(df, tmp_wcvp, how='inner', on='kew_id')
    return match


# In[8]:


#Find closely matching scientific name using difflib.get_close_matches if scientific name was not found
def find_sim(sci_name, wcvp, only_from_genus=True):
    if sci_name==sci_name:
        # Search for similar sci_name with same genus
        smpl_genus=sci_name.split(' ')[0]
        wcvp_gen=wcvp[wcvp.genus.isin([smpl_genus])]
        if wcvp_gen.shape[0]>0:
            sim_tax = difflib.get_close_matches(sci_name, 
                                                wcvp_gen.taxon_name.astype(str), n=1, cutoff=.9)
            if len(sim_tax)>0:
                return sim_tax[0]
         # If didn't work, search for similar sci_name
        else:
            if only_from_genus==False:
                sim_tax = difflib.get_close_matches(sci_name, 
                                            wcvp.taxon_name.astype(str), n=1, cutoff=.9)
                if len(sim_tax)>0:
                    return sim_tax[0]
                else:
                    return None
    else:
        return None
# print(find_sim('Combretum mussaendiflora', wcvp))
# print(find_sim('Scaveola humilis', wcvp, only_from_genus=False))


# In[9]:


#Find closely matching scientific name using kew namematching system
def kew_namematch(sci_name, verbose=False):
    url = "http://namematch.science.kew.org/api/v2/powo/csv"
    payload = '{\"column\": 0,\"headers\": false,\"outputAllColumns\": true,\"currentChunk\": 0,\"data\": [[\"' +        sci_name + '\"]]}'
    headers = {
      'Content-Type': 'application/json'
    }
    try:
        response = requests.request("POST", url, headers=headers, data = payload)
        content=str(response.text).replace("b'",'').replace("'",'')
        if verbose:
            print(content)
        stats_dict=ast.literal_eval(content)['stats']
        if stats_dict['matched']>0:
            try:
                record_dict=ast.literal_eval(content)['records']
                rec_df = pd.DataFrame(record_dict[1], index =record_dict[0]).T
                return rec_df.loc[0,'Scientific Name']
            except:
                print('ERROR: Failed to convert to dict -',sci_name)
                return None
    except:
        print('ERROR: No valid response -',sci_name)
        return None
# print(kew_namematch('Combretum mussaendiflora',verbose=True))


# In[10]:


#Find closely matching scientific name
def get_sim(df, wcvp, find_most_similar, verbose=False):
    print('\nLooking for most similar names')
    df['Similar_sci_name']=np.nan
    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        if find_most_similar=='similarity_genus': 
            df.loc[idx,'Similar_sci_name']=find_sim(row.sci_name, wcvp, only_from_genus=True)
        elif find_most_similar=='similarity': 
            df.loc[idx,'Similar_sci_name']=find_sim(row.sci_name, wcvp, only_from_genus=False)
        elif find_most_similar=='kewmatch': 
            df.loc[idx,'Similar_sci_name']=kew_namematch(row.sci_name)
        if verbose:
            print(idx,row.sci_name,':',df.loc[idx,'Similar_sci_name'])
    
    df['InWCVP']=(df.Similar_sci_name.isna()==False)
    df['Similar_match']=(df.Similar_sci_name.isna()==False)
    
    # Recover accepted and synonym taxa for found similar names
    df = df.drop(columns=['sci_name']).rename(columns={'Similar_sci_name':'sci_name'})
        
    return df


# In[11]:


def get_duplicates_type(df):
    df = df.sort_values('ID').reset_index().drop(columns='index')
    df['Duplicate_type']='Different_taxa'
    ID_ls=df.ID.unique()
    for iID in ID_ls:
        tmp_df=df[df.ID==iID]
        # Check if all entries have the same genus
        if tmp_df.genus.nunique()==1:
            df.loc[tmp_df.index,'Duplicate_type']='Same_Genus'
        # Check if all entries have the same species
        if len(set(return_dupl.genus + ' ' + return_dupl.species))==1:
            df.loc[tmp_df.index,'Duplicate_type']='Same_Species'
        # Check if all entries have the same taxon name
        if tmp_df.taxon_name.nunique()==1:
            df.loc[tmp_df.index,'Duplicate_type']='Same_Taxon'
            
    return df


# ## Main

# In[14]:


if __name__ == "__main__":
    print('\n\n##### wcvp_taxo v0.4 ##### \nAuthor:   Kevin Leempoel \nLast update: 2021-02-18\n')
    
    print(wcvp_path, df_path, 'g:', resolve_genus, ' s:', find_most_similar, ' d:', dupl_action,
      ' oc:', only_changes, ' os:', simple_output, ' v:', verbose)
    
    #Load libraries depending on the similarity method
    if find_most_similar in ['similarity_genus','similarity']:
        import difflib
    if find_most_similar=='request_kew':
        import requests
        import ast
        
    ## Loading and preparing data
    print('\n\nLoading and preparing data')
    wcvp = load_wcvp(wcvp_path)
    smpl_df = load_df(df_path)
    # Find scientific names
    smpl_df = define_sci_name(smpl_df)
    # Select or make ID column
    colID=GetIDcol(smpl_df)
    if colID=='ID':
        smpl_df['ID']=smpl_df.index
    else:
        smpl_df['ID']=smpl_df[colID]
    
    
    ## Initial checks
    smpl_df['sci_name']=smpl_df['Ini_sci_name']
    print('\n\nInitial checks')
    
    # Check if Ini_scinames are written as Genus sp.
    smpl_df['Genus_sp'] = smpl_df['sci_name'].str.split(' ',expand=True)[1].isin(['sp.','sp'])
    if smpl_df.Genus_sp.sum()>0:
        smpl_df.loc[smpl_df.Genus_sp,'sci_name'] = smpl_df.loc[smpl_df.Genus_sp,'sci_name'].str.split(' ',expand=True)[0]
    if resolve_genus:
        print('Genus sp:',(smpl_df.Genus_sp==True).sum(),'IDs with Genus sp. as sci_name')
    
    # Check if Ini_scinames exist in WCVP
    smpl_df['InWCVP']=smpl_df.sci_name.isin(wcvp.taxon_name)
    print('Missing taxa:',(smpl_df.InWCVP==False).sum(),'IDs not in WCVP')
    
    # Optional. Find similar names if not in WCVP
    if find_most_similar in ['similarity_genus','similarity','request_kew']:
        resolved_sim = get_sim(smpl_df[smpl_df.InWCVP==False],wcvp=wcvp,find_most_similar=find_most_similar,verbose=verbose)
        smpl_df = pd.concat([smpl_df[~smpl_df.ID.isin(resolved_sim.ID)], resolved_sim])
        print('find_most_similar: found',smpl_df.Similar_match.sum(),'IDs by similarity')
        
    # Check if Ini_scinames have duplicate entries
    dupl_taxon_names = wcvp[wcvp.taxon_name.isin(smpl_df.sci_name)]['taxon_name'].to_frame().groupby('taxon_name').size()                    .to_frame().reset_index().rename(columns={0:'count'})
    dupl_taxon_names = dupl_taxon_names[dupl_taxon_names['count']>1].taxon_name
    smpl_df['Duplicates']=smpl_df.sci_name.isin(dupl_taxon_names)
    print('Duplicates:',(smpl_df.Duplicates==True).sum(),'IDs matching multiple entries in WCVP')
    
    # Simpler dataframe 
    if smpl_df[~smpl_df.InWCVP].shape[0]>0:
        print('No match for',smpl_df[~smpl_df.InWCVP].shape[0],'IDs')
    smpl_dfs = smpl_df[smpl_df.InWCVP][['ID','sci_name','Duplicates']]
    
    
    
    ## get WCVP taxons
    # Recover accepted and unplaced taxa
    print('\n\nMatching & Resolving')
    match = get_by_taxon_name(smpl_dfs[(smpl_dfs.Duplicates==False)], wcvp)
    return_df = match[match.taxonomic_status.isin(status_keep)]
    print('After direct matching: found match for',return_df.shape[0],'IDs')
    
    # Resolving synonyms
    synonyms = match[match.taxonomic_status.isin(['Synonym','Homotypic_Synonym'])]                .rename(columns={'kew_id':'Ini_kew_id','accepted_kew_id':'kew_id','taxonomic_status':'Ini_taxonomic_status'})
    cols_syn=list(smpl_dfs.columns) + ['Ini_kew_id','Ini_taxonomic_status']
    return_syn = get_by_kew_id(df = synonyms[cols_syn + ['kew_id']], wcvp = wcvp)
    return_df=pd.concat([return_df,return_syn]).reset_index().drop(columns='index')
    print('After resolving synonyms: found match for',return_df.shape[0],'IDs')
    
    
    
    ## Resolving duplicates
    match=get_by_taxon_name(smpl_dfs[(smpl_dfs.Duplicates==True)],wcvp)
    return_dupl = match[match.taxonomic_status.isin(status_keep)]
    
    # Resolving synonyms in duplicates
    synonyms = match[match.taxonomic_status.isin(['Synonym','Homotypic_Synonym'])]                .rename(columns={'kew_id':'Ini_kew_id','accepted_kew_id':'kew_id','taxonomic_status':'Ini_taxonomic_status'})
    return_syn = get_by_kew_id(df = synonyms[cols_syn + ['kew_id']], wcvp = wcvp)
    return_dupl=pd.concat([return_dupl,return_syn]).reset_index().drop(columns='index')
    
    # Action on duplicates
    dupl_df = get_duplicates_type(return_dupl)
    print('Duplicates: found', dupl_df.shape[0],'entries, for',dupl_df.ID.nunique(),
          'duplicated ID, corresponding to', dupl_df.kew_id.nunique(),'kew_id')
    if dupl_action=='rank':
        # Sort matches by ranked taxonomic status and ID
        dupl_df['taxo_rank']=dupl_df.taxonomic_status
        dupl_df.loc[dupl_df['Ini_taxonomic_status'].isna()==False,'taxo_rank']=                    dupl_df.loc[dupl_df['Ini_taxonomic_status'].isna()==False,'Ini_taxonomic_status']
        dupl_df['taxo_rank'] = dupl_df['taxo_rank']            .replace({'Accepted':1, 'Unplaced':2, 'Synonym':3,'Homotypic_Synonym':4,'Artificial Hybrid':5})
        dupl_df = dupl_df.sort_values('taxo_rank').drop_duplicates('ID').reset_index().drop(columns=['taxo_rank','index'])
        return_df = pd.concat([return_df,dupl_df])
    else:
        if dupl_action in ['divert_taxonOK','divert_speciesOK','divert_genusOK']:
            if dupl_action=='divert_taxonOK':
                keep_dup=['Same_Taxon']
            elif dupl_action=='divert_speciesOK':
                keep_dup=['Same_Taxon','Same_Species']
            elif dupl_action=='divert_genusOK':
                keep_dup=['Same_Taxon','Same_Species','Same_Genus']
                
            return_dupl = dupl_df[dupl_df.Duplicate_type.isin(keep_dup)].drop_duplicates('ID')
            return_df = pd.concat([return_df,return_dupl])
            dupl_df = dupl_df[~dupl_df.Duplicate_type.isin(keep_dup)]
        else:
            return_df = return_df[~return_df.ID.isin(dupl_df.ID)]
            return_df['Duplicate_type']=np.nan
        
        # Output remaining duplicates
        if dupl_df.shape[0]>0:  
            print('Unresolved duplicates: writing separate spreadsheet with', dupl_df.shape[0],'entries, for',
                  dupl_df.ID.nunique(),'duplicated ID, corresponding to', dupl_df.kew_id.nunique(),'kew_id')
            dupl_df2 = pd.merge(smpl_df.drop(columns=['InWCVP']),
                                dupl_df.drop(columns=['Duplicates','sci_name']),how='inner',on='ID')
            dupl_df2 = dupl_df2.sort_values('ID').reset_index().drop(columns=['index','accepted_kew_id',
                     'accepted_name','accepted_authors','reviewed']).rename(columns={'sci_name':'sci_name_query'})
            dupl_df2.to_csv(df_path.replace('.csv','_duplicates.csv'),index=False,encoding='utf-8')
    print('After resolving duplicates: found match for',return_df.shape[0],'IDs')
    
    
    ## Cleaning and merging DF
    smpl_df = pd.merge(smpl_df.drop(columns=['InWCVP']),return_df.drop(columns=['Duplicates','sci_name']),how='left',on='ID')
    # Filter duplicates and unresolved taxa
    print(smpl_df.shape)
    if dupl_action in dupl_action in ['divert_taxonOK','divert_speciesOK','divert_genusOK']:
        smpl_df = smpl_df[~((smpl_df.Duplicates==True) & (smpl_df.Duplicate_type.isin(keep_dup)==False))]
        print(smpl_df.shape)
    if find_most_similar in ['similarity_genus','similarity','request_kew']:
        unresolved = smpl_df[smpl_df.Similar_match==False]
        print(unresolved.shape)
        unresolved.to_csv(df_path.replace('.csv','_unresolved.csv'),index=False,encoding='utf-8')
        smpl_df = smpl_df[~smpl_df.ID.isin(unresolved.ID)]
        print(smpl_df.shape)
    # Modify taxonomy for genus sp.
    if dupl_action=='rank':
        mod_gensp=smpl_df[(smpl_df.Genus_sp)].index
    else:
        mod_gensp=smpl_df[(smpl_df.Genus_sp) | (smpl_df.Duplicate_type=='Same_Genus')].index
    smpl_df.loc[mod_gensp,'taxon_name']=smpl_df.loc[mod_gensp,'genus'] + ' sp.'
    smpl_df.loc[mod_gensp,'species']=np.nan
    #Sort table by ID
    out_df = smpl_df.sort_values('ID').reset_index()                .drop(columns=['index','Genus_sp','accepted_kew_id','accepted_name','accepted_authors','reviewed'])                .rename(columns={'sci_name':'sci_name_query','taxon_name':'sci_name'})
    
    
    
    ## Output options
    # Output All
    if only_changes==False:
        if simple_output:
            out_df = out_df[[colID,'kew_id','Ini_sci_name','sci_name']]
        else:
            out_df = out_df.drop(columns='ID')
        out_df.to_csv(df_path.replace('.csv','_wcvp.csv'),index=False,encoding='utf-8')
    # Output changes only    
    elif only_changes:
        out_df['Same_sci_name']=(out_df.Ini_sci_name==out_df.sci_name)
        if 'Ini_Family' in out_df:
            out_df['Same_family']=(out_df.Ini_Family==out_df.family)
            print(out_df.groupby('Same_family').size())
            out_df = out_df[(out_df.Same_family==False) | (out_df.Same_sci_name==False)]                    .drop(columns=['Same_sci_name','Same_family'])
        else:
            out_df = out_df[(out_df.Same_sci_name==False)].drop(columns='Same_sci_name')
            
        if simple_output:
            out_df = out_df[[colID,'kew_id','Ini_sci_name','sci_name']]
        else:
            out_df = out_df.drop(columns='ID')

        out_df = out_df[(out_df.kew_id.notnull())]
        print('Only_changes:',out_df.shape[0],'IDs have changed taxonomy')
        out_df.to_csv(df_path.replace('.csv','_wcvp_changes.csv'),index=False,encoding='utf-8')
        
    print('Done!')


# In[33]:


smpl_df


# In[34]:


smpl_df.loc[smpl_df.Genus_sp,'sci_name']


# In[35]:


smpl_df.Genus_sp


# In[36]:


smpl_df.loc[smpl_df.Genus_sp,'sci_name'].str.split(' ',expand=True)


# In[37]:


out_df


# In[38]:


out_df.isna().sum()


# In[39]:


smpl_df

