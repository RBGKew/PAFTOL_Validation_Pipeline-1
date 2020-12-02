#!/bin/bash
#SBATCH --job-name="genbank_download"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=32000

project_dir=$1
cd $project_dir

my_query=$2
filename=$3

# Download NCBI data as fasta using a query
esearch -db nucleotide \
-query "$my_query" \
| efetch -format gb > gb_files/$filename.gb

# Convert to fasta and output TaxID, SeqID, Organism name
source activate py36
python Genbank2Fasta.py gb_files/$filename.gb

# Get stats from fasta file
seqkit stats "$filename"_ini.fasta

# Count accessions
wc "$filename"_SeqID.csv
