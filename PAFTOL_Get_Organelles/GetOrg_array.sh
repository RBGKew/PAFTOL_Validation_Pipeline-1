#!/bin/bash

#SBATCH --job-name="GetOrg"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem=80000
#SBATCH --ntasks=1

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

ncpu=4

module load python/3.7.9
module load blast bowtie2 spades
module load getorganelle


sample_file=$1
org=$2

iline=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sample_file)
echo $iline
sample="$(cut -d',' -f1 <<<"$iline")"
echo $sample
file_path_R1="$(cut -d',' -f2 <<<"$iline")"
file_R1=`basename "$file_path_R1"`; file_R1=${file_R1/.gz/}
echo $file_R1
gunzip -c Data/$file_R1.gz > tmp_fastq_"$org"/$file_R1

file_path_R2="$(cut -d',' -f3 <<<"$iline")"
file_R2=`basename "$file_path_R2"`; file_R2=${file_R2/.gz/}
echo $file_R2
gunzip -c Data/$file_R2.gz > tmp_fastq_"$org"/$file_R2


if [ $org == pt ] 
then
	get_organelle_from_reads.py -1 tmp_fastq_"$org"/$file_R1 -2 tmp_fastq_"$org"/$file_R2 -o GetOrg/"$sample"_pt \
	-R 20 -k 21,45,65,85,105 -t $ncpu -F embplant_pt --zip-files > \
	logs/log_${sample}_pt.log 2> logs/log_${sample}_pt.err
elif [ $org == nr ]
then
	 get_organelle_from_reads.py -1 tmp_fastq_"$org"/$file_R1 -2 tmp_fastq_"$org"/$file_R2 -o GetOrg/"$sample"_nr \
	-R 10 -k 35,85,115 -t $ncpu -F embplant_nr --zip-files > \
	logs/log_${sample}_nr.log 2> logs/log_${sample}_nr.err
fi

rm tmp_fastq_"$org"/$file_R1; rm tmp_fastq_"$org"/$file_R2
python ../GetOrg_Clean.py --path GetOrg/"$sample"_"$org"/