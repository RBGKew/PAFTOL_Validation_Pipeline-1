# DNA barcode validation
## Creation or recovery of reference databases
- Create tools to download barcodes from NCBI
- Script for ecoPCR db creation
- Annotated taxonomy of database
- Summary statistics from barcode dbs

## Barcode recovery from samples
- Get Organelles script
- Cutadapt

## Multi-locus barcoding
- Script to blast samples against multiple databases
- Collect results for each blast

## Sample Validation
- Comparison of sample against collected results

# Dependencies
```python
- entrez-direct 
- seqkit
- cutadapt
- blast
- python3
	- Pandas
	- Numpy
	- Bio
	- TQDM
	- Seaborn
```

# Databases
## .A NCBI to Barcode DB
Barcode databases are created from NCBI in 3 steps
1. Downloading the accessions as a .fasta using a query
2. Trimming sequences using barcode primers
3. Create a .fasta database with taxonomic information




# Barcode recovery from samples
## Recover plastom and ribosomal DNA from raw reads with GetOrganelles
seqkit stats -j 4 out_fasta_*/*.fasta > PAFTOL_PILOT404_Organelles.stats


## Multi-locus barcoding
```console
ls okp_data_fasta > okp_data_fasta_list.txt
sbatch blast_barcodes_array_v2.sh /mnt/shared/scratch/kleempoe/paftol/org_barcode/OKP/ okp_data_fasta ../Barcode_DB/db_ls.txt
```

# Sample Validation
## 1. Taxonomic Validation
With this function, samples are:
* Validated if the species is within the first N blast matches (N is 1 by default for highest confidence)
* Invalidated if not within the first N blast matches
* Untested if the species is not in the database
* Untested if there are no matches

This validation can be done at any taxonomic level but its goal was first to validate a species.
Providing a minimum identity is recommended

## 2. Multi-rank Validation
With this function, samples are:
* Validated if the lowest taxonomic rank is within the first N blast matches (N is 1 by default)
* Invalidated if not within the first X blast matches. It will not be tested at higher taxonomic ranks.
* If the taxonomic rank is not in the database, the test is repeated at the next taxonomic rank (species, genus and family by default)
* Untested if the highest rank is not in the DB
* Untested if there are no matches

 1 - Best Hit, Rigorous
Sample is valid if it matches the best hit at given taxonomic level, provided it exist in the database at that taxonomic level. If it fails at the lowest taxonomic rank, it is not tested at higher taxonomic ranks. 


## 3. Validation by Minimum Identity
With this function, samples are:
* Validated if the taxa is in the blast matches, regarless of its ranking
* Invalidated if there is no match at all
* Untested if the taxa is not in the database
* Untested if no sequences passed filtering

This validation is useful to detect contaminated samples, for species absent from databases, and to overcome biases of the procedure (see below).
This test is usually performed at genus and family level. 
Providing a minimum identity is recommended.


## Possible reasons for test failure that are not biological or lab errors:
* Many sequences in reference databases are mislabelled, particularly at species level (up to 20% according to studies)
* Longer sequence match has a higher bitscore
* Shorter sequence match has a higher identity
* Heteroplasmy