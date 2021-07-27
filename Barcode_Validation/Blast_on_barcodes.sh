#!/bin/bash
#SBATCH --job-name="blast_barcodes"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --mem=8000
ncpu=4

project_dir=$1
type=$2
cd $project_dir

Samples_file=Samples_to_barcode.txt
iline=$(sed -n "$SLURM_ARRAY_TASK_ID"p $Samples_file)
echo $iline
sample="$(cut -d',' -f1 <<<"$iline")"


if [ $type == contigs ] 
then
	fasta_pt_file=in_fasta/"$sample".fasta
	fasta_nr_file=in_fasta/"$sample".fasta
elif [ $type == pt_nr ] 
then
	fasta_pt_file=fasta_pt/"$sample"_pt.fasta
	fasta_nr_file=fasta_nr/"$sample"_nr.fasta
fi

echo "$sample"
echo "$fasta_pt_file"
echo "$fasta_nr_file"

     
if [ -f $fasta_nr_file ]; then
	for idb in NCBI_pln_18s;
	do
		echo $idb
		blastn  -query $fasta_nr_file -db ../Barcode_DB/"$idb".fasta \
		-perc_identity 95 -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
		-num_threads $ncpu -max_target_seqs 2000 \
		-out out_blast/"$sample"-"$idb".out
	done

else
   echo "File $fasta_nr_file does not exist."
fi

if [ -f $fasta_pt_file ]; then

for idb in NCBI_pln_trnL NCBI_pln_trnH BOLD_pln_rbcLa BOLD_pln_matK;
do
	echo $idb
	blastn  -query $fasta_pt_file -db ../Barcode_DB/"$idb".fasta \
	-perc_identity 95 -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads $ncpu -max_target_seqs 2000 \
	-out out_blast/"$sample"-"$idb".out
done

for idb in NCBI_pln_plastome;
do
	echo $idb
	blastn  -query $fasta_pt_file -db ../Barcode_DB/"$idb".fasta \
	-perc_identity 95 -outfmt "6 qseqid sseqid pident length slen qlen mismatch gapopen qstart qend sstart send evalue bitscore" \
	-num_threads $ncpu -max_target_seqs 500 \
	-out out_blast/"$sample"-"$idb".out
done

else
   echo "File $fasta_pt_file does not exist."
fi