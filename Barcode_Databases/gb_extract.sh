#!/bin/bash
#SBATCH --job-name="gb_ext"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --mem=32000

source activate py36 

# python GB_extract.py NCBI_18s
# python GB_extract.py NCBI_26s
# python GB_extract.py NCBI_23s
# python GB_extract.py NCBI_rbcL
# python GB_extract.py NCBI_trnL
# python GB_extract.py NCBI_ITS1
# python GB_extract.py NCBI_ITS2
# python GB_extract.py NCBI_rpl2
python GB_extract.py NCBI_ndhf
