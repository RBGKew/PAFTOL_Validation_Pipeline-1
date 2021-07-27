#!/bin/bash
module load python/3.7.9
# Only arguments are DataSource and paftol_export
DataSource=$1
paftol_export=$2


# Make lists of remaining of samples that have no organelles recovered
python GetOrg_prep.py --db $paftol_export --DataSource $DataSource

cd $DataSource
mkdir logs; mkdir fasta_pt; mkdir fasta_nr; mkdir fastq_pt; mkdir fastq_nr;

# Launch remaining pt
a=($(wc 'remaining_pt.txt')); Ns_pt=${a[0]}; echo $Ns_pt
sbatch -p all --array=1-${Ns_pt}%$slurmThrottle ../GetOrg_array.sh remaining_pt.txt
# Launch remaining nr
a=($(wc 'remaining_nr.txt')); Ns_nr=${a[0]}; echo $Ns_nr
sbatch -p all --array=1-${Ns_nr}%$slurmThrottle ../GetOrg_array.sh remaining_nr.txt
