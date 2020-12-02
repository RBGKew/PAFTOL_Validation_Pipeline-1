# DNA Barcode Validation - Barcode Databases

## Barcode recovery and specific processing
For NCBI queries, the script `ncbi_query_genbank_download.sh` downloads an NCBI database in genbank format and converts it to a fasta file and an list of accessions containing Accession ID, organism name and taxonomic ID.

### NCBI queries and trimming for 18s, trnH-psbA and trnL
```console
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
```console
min_len=10000
max_len=250000
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta
```

### BOLD databases
for BOLD Databases, the 
```console
# BOLD rbcL
min_len=500
max_len=1500
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta

# BOLD matK
min_len=500
max_len=1600
seqkit seq -j $ncpu --min-len $min_len --max-len $max_len tmp/"$barcode_db"_raw.fasta > tmp/"$barcode_db"_trim.fasta
```



## Common processing
### Trimming
As all databases were all trimmed either with primers or by length, the list of accessions was consequently updated using `update_table_from_list.py`
```console
seqkit fx2tab -j $ncpu -i -n tmp/"$barcode_db"_trim.fasta > tmp/"$barcode_db"_seqID.txt
wc -l tmp/"$barcode_db"_seqID.txt
python update_table_from_list.py tmp/"$barcode_db"_seqID.txt tmp/"$barcode_db"_raw_SeqID.csv
wc -l tmp/"$barcode_db"_raw_SeqID.csv
```

### Taxonomy for database
Use wcvp_taxo_v03.py script to recover the taxonomy for a given .fasta db. Will write a new fasta file containing only accessions for which a taxonomy was recovered.

```console
python wcvp_taxo_v03.py wcvp_export.txt tmp/"$barcode_db"_raw_SeqID_update.csv -g
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