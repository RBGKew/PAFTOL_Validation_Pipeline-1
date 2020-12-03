# DNA Barcode Validation - Barcode Databases

We queried the NCBI nucleotide database for the following barcoding loci: 18s, trnH-psbA and trnL, and downloaded the resulting file in genbank format using Entrez-direct (Kans 2010). trnL and trnH sequences were trimmed in Cutadapt (Martin 2011) using universal barcode primers with 10% mismatch tolerated. For trnH-psbA, we used primers psbAf and trnH2 (Tate and Simpson 2003; Sang, Crawford, and Stuessy 1997), and for trnL, primers trnL-c and trnL-d (Taberlet et al. 2007). All NCBI references were filtered based on length, with a minimum and maximum length set for each barcode. Finally, we added the most recent release of plastid data as a reference database of whole plastomes (https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/).

## Barcode recovery and specific processing

For NCBI queries, the script `ncbi_query_genbank_download.sh` downloads an NCBI database in genbank format using entrez-direct and converts it to a fasta file and an list of accessions containing Accession ID, organism name and taxonomic ID.

### NCBI queries and trimming for 18s, trnH-psbA and trnL
Queries were performed on 29/10/2020
```shell
# 18s
sbatch ncbi_query_genbank_download.sh $my_dir '"18S ribosomal RNA"[All Fields] AND plants[filter] AND ("1500"[SLEN] : "2000"[SLEN])' NCBI_pln_18s

min_len=1500
max_len=1900
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta


# trnH
sbatch ncbi_query_genbank_download.sh $my_dir 'trnH[All Fields] AND plants[filter] AND ("250"[SLEN] : "650"[SLEN])' NCBI_pln_trnH

PriF=GTTATGCATGAACGTAATGCTC
PriRRev=GGATTGTGAATCCACCATGCGCG
min_len=100
max_len=1000
cutadapt -a $PriF...$PriRRev --revcomp -j $ncpu \
	--minimum-length $min_len --maximum-length $max_len \
	-o tmp/"$barcode_db"_trim.fasta -e 0.15 --discard-untrimmed tmp/"$barcode_db"_raw.fasta


# trnL
sbatch ncbi_query_genbank_download.sh $my_dir 'trnL[All Fields] AND plants[filter]' NCBI_pln_trnL

PriF=CGAAATCGGTAGACGCTACG
PriRRev=GTTCAAGTCCCTCTATCCCC
min_len=200
max_len=800
cutadapt -a $PriF...$PriRRev --revcomp -j $ncpu \
	--minimum-length $min_len --maximum-length $max_len \
	-o tmp/"$barcode_db"_trim.fasta -e 0.15 --discard-untrimmed tmp/"$barcode_db"_raw.fasta
```

### Whole Plastomes NCBI
Whole plastomes were recovered and concatenated in a single file from https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.gbff.gz (accessed on 29/10/2020)
Then further processed as follows:
```shell
min_len=10000
max_len=250000
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta
```

### BOLD databases
The BOLD Database was downloaded from bold (accessed on 29/10/2020), and sequences for rbcLa and matK were extracted using the following python jupyter notebook `Processing_BOLD_v2.ipynb`.
Sequences were then further filtered as follows:
```shell
# BOLD rbcLa
min_len=450
max_len=750
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta

# BOLD matK
min_len=500
max_len=1600
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta
```



## Common processing
### Trimming
As all databases were all trimmed either with primers or by length, the list of accessions was consequently updated using `seqkit` and `update_table_from_list.py`
```shell
seqkit fx2tab -j $ncpu -i -n tmp/"$barcode_db"_trim.fasta > tmp/"$barcode_db"_seqID.txt
wc -l tmp/"$barcode_db"_seqID.txt
python update_table_from_list.py tmp/"$barcode_db"_seqID.txt tmp/"$barcode_db"_raw_SeqID.csv
wc -l tmp/"$barcode_db"_raw_SeqID.csv
```

### Taxonomy for database
Use [WCVP_taxo](../WCVP_Taxo/) script to recover the taxonomy for a given .fasta db, and writes a new fasta file containing only accessions for which taxons could be resolved against WCVP. Scientific names in genus sp. format were resolved (option -g), as well as scientific names for which duplicate entries all mapped to the same genus (option -d divert_genusOK). Sequences with unresolved names or that matched duplicate entries with different genera names were discarded. 
WCVP database is http://sftp.kew.org/pub/data-repositories/WCVP/wcvp_v2_jun_2020.zip
```shell
python wcvp_taxo_v03.py wcvp_export.txt tmp/"$barcode_db"_raw_SeqID_update.csv -g -d divert_genusOK
wc -l tmp/"$barcode_db"_raw_SeqID_update_wcvp.csv
cp tmp/"$barcode_db"_raw_SeqID_update_wcvp.csv "$barcode_db"_TAXO.csv
awk -F"," '{print $1}' "$barcode_db"_TAXO.csv > tmp/"$barcode_db"_seqID.txt
wc -l tmp/"$barcode_db"_seqID.txt
seqtk subseq tmp/"$barcode_db"_trim.fasta tmp/"$barcode_db"_seqID.txt > "$barcode_db".fasta

#Verification
seqkit stats -j $ncpu "$barcode_db".fasta
wc -l "$barcode_db"_TAXO.csv
seqkit stats -j $ncpu tmp/"$barcode_db"*.fasta
```

### Converting to blastDB
```shell
for idb in BOLD_pln_matK  BOLD_pln_rbcLa NCBI_pln_trnL NCBI_pln_trnH NCBI_pln_18s NCBI_pln_plastome
do
	echo $idb
	makeblastdb -in $idb.fasta -dbtype nucl -parse_seqids
done
```