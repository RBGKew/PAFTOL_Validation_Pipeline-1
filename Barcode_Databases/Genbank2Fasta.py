#!/usr/bin/env python
# coding: utf-8

# In[179]:

try:
    print('loading dependencies',end='...')
    from Bio import SeqIO
    import pandas as pd
    import sys
    import numpy as np
    print('done')
except:
    print('libraries missing')



# In[180]:


input_file = sys.argv[1]
name = input_file.split('.')[0]
output_fasta = (name + '_raw.fasta').replace('gb_files','tmp')
output_table = (name + '_raw_SeqID.csv').replace('gb_files','tmp')

print('input .gb is',input_file)
print('output .fasta is ',output_fasta)
print('converting file',end='...')
# In[181]:


def get_qualifier(record, attribute):
    try:
        tmp_feat=seq_record.features[0]
        tmp_qual=tmp_feat.qualifiers
        return tmp_qual[attribute][0]
    except:
        return None 


# In[183]:


# https://gist.github.com/tammoippen/4474e838e969bf177155231ebba52386
def crappyhist(a, bins=50, width=140):
    h, b = np.histogram(a, bins)

    for i in range (0, bins):
        print('{:12.5f}  | {:{width}s} {}'.format(
            b[i], 
            '#'*int(width*h[i]/np.amax(h)), 
            h[i], 
            width=width))
    print('{:12.5f}  |'.format(b[bins]))


# In[184]:
output_handle = open(output_fasta, "w")

file_dic={}
for seq_record in SeqIO.parse(input_file, "genbank"):
    seq_dic={}
    seq_dic['species'] = get_qualifier(seq_record, 'organism')
    seq_dic['TaxID']= get_qualifier(seq_record, 'db_xref').replace('taxon:','')
    seq_dic['mol_type'] = get_qualifier(seq_record, 'mol_type')
    seq_dic['length'] = len(seq_record.seq)
    
    try:
        #seq_dic['TaxID'] = int(seq_dic['TaxID'])
        #new_seqid = correct_accessions(seq_record.id)
        #new_seqid = seq_record.id
        #new_description = str(seq_dic).replace('{','').replace('}','').replace(': ','=').replace("'",'')

        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id, seq_record.description, seq_record.seq))

        file_dic[seq_record.id]=seq_dic
        #del new_seqid, new_description
		
    except:
        print('skipping ' + seq_record.id)
    
    
    del seq_dic

output_handle.close()
print('done')

# In[185]:


df = pd.DataFrame.from_dict(file_dic,orient='index').reset_index().rename(columns = {'index':'SeqID'})
df.to_csv(output_table,index=False)


# In[186]:


print('found',df.shape[0],'records')
print('found',df['SeqID'].nunique(),'unique accessions')
print('found',df['species'].nunique(),'unique species names')
print('found',df['TaxID'].nunique(),'unique TaxIDs')
print(df.groupby('mol_type').size().sort_values(ascending=False))


# In[187]:


crappyhist(df['length'], bins=20, width=100)

