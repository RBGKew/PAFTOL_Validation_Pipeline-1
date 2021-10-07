#!/bin/bash
module load python/3

# Only arguments needed are DataSource and paftol_export
# e.g: ./GetOrg_Pipeline.sh 2021-07-27_paftol_export.csv "PAFTOL"
DataSource=$2
paftol_export=$1
rem_search="fasta" # log or fasta
slurmThrottle=5


## Make lists of remaining  samples that have no organelles recovered
rm -f $DataSource/remaining_pt.txt
rm -f $DataSource/remaining_nr.txt;
python GetOrg_prep.py --db $paftol_export --DataSource $DataSource --rem_search $rem_search


## Go to dir and create working directories
cd $DataSource
mkdir -p GetOrg; mkdir -p logs; mkdir -p fasta_pt; mkdir -p fasta_nr; mkdir -p Archives;


## Copy .fastq.gz files in currrent Data directory
mkdir -p Data;
while read iline; do
	file_path_R1="$(cut -d',' -f2 <<<"$iline")"
	file_R1=`basename "$file_path_R1"`; file_R1=${file_R1/.gz/}
	echo $file_R1
	cp -n $file_path_R1 Data/$file_R1.gz

	file_path_R2="$(cut -d',' -f3 <<<"$iline")"
	if [ ! -z "$file_path_R2" ]; then
		file_R2=`basename "$file_path_R2"`; file_R2=${file_R2/.gz/}
		echo $file_R2
		cp -n $file_path_R2 Data/$file_R2.gz
	else
		echo "No R2"
	fi
done < remaining_pt.txt

while read iline; do
	file_path_R1="$(cut -d',' -f2 <<<"$iline")"
	file_R1=`basename "$file_path_R1"`; file_R1=${file_R1/.gz/}
	echo $file_R1
	cp -n $file_path_R1 Data/$file_R1.gz

	file_path_R2="$(cut -d',' -f3 <<<"$iline")"
	if [ ! -z "$file_path_R2" ]; then
		file_R2=`basename "$file_path_R2"`; file_R2=${file_R2/.gz/}
		echo $file_R2
		cp -n $file_path_R2 Data/$file_R2.gz
	else
		echo "No R2"
	fi
done < remaining_nr.txt


## Launch remaining pt
a=($(wc 'remaining_pt.txt')); Ns_pt=${a[0]}; echo $Ns_pt
if (( $Ns_pt > 0 )); then
	sbatch --array=1-${Ns_pt}%$slurmThrottle ../GetOrg_array.sh remaining_pt.txt "pt"
fi

## Launch remaining nr
a=($(wc 'remaining_nr.txt')); Ns_nr=${a[0]}; echo $Ns_nr
if (( $Ns_nr > 0 )); then
	sbatch --array=1-${Ns_nr}%$slurmThrottle ../GetOrg_array.sh remaining_nr.txt "nr"
fi
