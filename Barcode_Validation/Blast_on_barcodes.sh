#!/bin/bash
#SBATCH --job-name="blast_barcodes"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=8000

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################
ncpu=4

project_dir=$1
type=$2
cd $project_dir

barcodes_table=../Barcode_DB/Barcode_Tests.csv
Samples_file=Samples_to_barcode.txt
# sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $Samples_file)
sample=$(sed -n 1p $Samples_file)

if [ $type == contigs ] 
then
	fasta_pt_file=in_fasta/"$sample".fasta
	fasta_nr_file=in_fasta/"$sample".fasta
elif [ $type == pt_nr ] 
then
	fasta_pt_file=fasta_pt/"$sample"_pt.fasta
	fasta_nr_file=fasta_nr/"$sample"_nr.fasta
fi

echo "sample:$sample,pt_fasta:$fasta_pt_file,nr_fasta:$fasta_nr_file"

sed 1d $barcodes_table | while read iline; do
	idb="$(cut -d',' -f1 <<<"$iline")"
	type="$(cut -d',' -f4 <<<"$iline")"
	max_blast="$(cut -d',' -f5 <<<"$iline")"
	blast_pid="$(cut -d',' -f6 <<<"$iline")"
	echo "Barcode Test:$idb,type:$type,blast_max_matches:$max_blast,blast_min_pid:$blast_pid"
	
	if [ $type == nr ] && [ -f $fasta_nr_file ]; then
		blastn  -query $fasta_nr_file -db ../Barcode_DB/"$idb".fasta \
			-perc_identity $blast_pid -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
			-num_threads $ncpu -max_target_seqs $max_blast \
			-out out_blast/"$sample"-"$idb".out
	elif [ $type == pt ] && [ -f $fasta_pt_file ]; then
		blastn  -query $fasta_pt_file -db ../Barcode_DB/"$idb".fasta \
			-perc_identity $blast_pid -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
			-num_threads $ncpu -max_target_seqs $max_blast \
			-out out_blast/"$sample"-"$idb".out
	else
	  echo "ERROR $idb, invalid type or no fasta file"
	fi
done

python ../Get_validation_cards.py --sample $sample --samples_file "$project_dir"_samples.csv --barcodes_table $barcodes_table