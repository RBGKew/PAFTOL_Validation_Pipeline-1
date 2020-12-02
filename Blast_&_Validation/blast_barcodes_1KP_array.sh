#!/bin/bash
#SBATCH --job-name="blast_barcodes"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=12000
#SBATCH --array=1-769%40


project_dir=$1
sample_file=$2
ncpu=4

cd $project_dir
mkdir out_blast

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sample_file)

echo "$sample"

for idb in NCBI_pln_18s NCBI_pln_trnL BOLD_pln_rbcL BOLD_pln_rbcLa BOLD_pln_matK NCBI_pln_trnH;
do
	echo $idb
	blastn  -query data_fasta/"$sample".fasta -db ../Barcode_DB/"$idb".fasta \
	-perc_identity 95 -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads $ncpu -max_target_seqs 2000 \
	-out out_blast/"$sample"_"$idb".out
done

for idb in NCBI_pln_plastome;
do
	echo $idb
	blastn  -query data_fasta/"$sample".fasta -db ../Barcode_DB/"$idb".fasta \
	-perc_identity 95 -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads $ncpu -max_target_seqs 500 \
	-out out_blast/"$sample"_"$idb".out
done
