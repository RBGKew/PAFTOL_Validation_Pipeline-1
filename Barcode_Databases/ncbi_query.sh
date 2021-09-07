#!/bin/bash
#SBATCH --job-name="genbank_download"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=32000

my_query=$1
filename=$2

# Download NCBI data as fasta using a query
esearch -db nucleotide -query "$my_query" | efetch -format gb > gb_files/$filename.gb