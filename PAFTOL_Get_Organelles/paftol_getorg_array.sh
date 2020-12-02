#!/bin/bash

#SBATCH --job-name="getorg_04"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium
#SBATCH --mem=32000
#SBATCH --ntasks=1
#SBATCH --array=1-1673%18

source activate py36

project_dir=$1
sample_file=$2
ncpu=4

mkdir logs
mkdir out_fasta_pt
mkdir out_fasta_nr

ifasta=$(sed -n "$SLURM_ARRAY_TASK_ID"p $sample_file)
echo $ifasta
sample="$(cut -d',' -f1 <<<"$ifasta")"
echo $sample
file_R1="$(cut -d',' -f2 <<<"$ifasta")"
echo $file_R1
file_R2="$(cut -d',' -f3 <<<"$ifasta")"
echo $file_R2

get_organelle_from_reads.py -1 $file_R1 -2 $file_R2 -o "$sample"_pt \
-R 20 -k 21,45,65,85,105 -t $ncpu -F embplant_pt --zip-files > \
logs/log_${sample}_pt.log 2> logs/log_${sample}_pt.err

get_organelle_from_reads.py -1 $file_R1 -2 $file_R2 -o "$sample"_nr \
-R 10 -k 35,85,115 -t $ncpu -F embplant_nr --zip-files > \
logs/log_${sample}_nr.log 2> logs/log_${sample}_nr.err

