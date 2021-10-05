#!/bin/bash

#SBATCH --job-name="GetOrg"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --mem=80000
#SBATCH --ntasks=1
ncpu=4

module load python/3
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


file_path_R2="$(cut -d',' -f3 <<<"$iline")"

if [ ! -z "$file_path_R2" ]; then
	file_R2=`basename "$file_path_R2"`; file_R2=${file_R2/.gz/}
	echo $file_R2
	
	if [ $org == pt ] 
	then
		get_organelle_from_reads.py -1 Data/$file_R1.gz -2 Data/$file_R2.gz -o GetOrg/"$sample"_pt \
		--max-reads 536870912 -R 20 -k 21,45,65,85,105 -t $ncpu -F embplant_pt --zip-files > \
		logs/log_${sample}_pt.log 2> logs/log_${sample}_pt.err
	elif [ $org == nr ]
	then
		 get_organelle_from_reads.py -1 Data/$file_R1.gz -2 Data/$file_R2.gz -o GetOrg/"$sample"_nr \
		--max-reads 536870912 -R 10 -k 35,85,115 -t $ncpu -F embplant_nr --zip-files > \
		logs/log_${sample}_nr.log 2> logs/log_${sample}_nr.err
	fi
	
else
	echo "Single-end Mode"
	
	if [ $org == pt ] 
	then
		get_organelle_from_reads.py -u Data/$file_R1.gz -o GetOrg/"$sample"_pt \
		--max-reads 536870912 -R 20 -k 21,45,65,85,105 -t $ncpu -F embplant_pt --zip-files > \
		logs/log_${sample}_pt.log 2> logs/log_${sample}_pt.err
	elif [ $org == nr ]
	then
		get_organelle_from_reads.py -u Data/$file_R1.gz -o GetOrg/"$sample"_nr \
		--max-reads 536870912 -R 10 -k 35,85,115 -t $ncpu -F embplant_nr --zip-files > \
		logs/log_${sample}_nr.log 2> logs/log_${sample}_nr.err
	fi
fi

python ../GetOrg_Clean.py --path GetOrg/"$sample"_"$org"/



