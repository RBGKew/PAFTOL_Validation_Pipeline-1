#!/bin/bash
#SBATCH --job-name="Barcode_Validation"
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --mem=4000

##################################
# Author: Kevin Leempoel

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################

### Command ###
# ./Barcode_Validation.sh 2021-07-27_paftol_export.csv 'OneKP'

source activate py36
slurmThrottle=10

### Inputs:
# latest paftol_export. Needs the following fields : idSequencing, idPaftol, DataSource, Family, Genus, Species
paftol_export=$1
# DataSource: PAFTOL, SRA, GAP, OneKP, AG, UG
DataSource=$2
if [[ $DataSource == OneKP || $DataSource == AG || $DataSource == UG ]]
then
type="contigs"
elif [[ $DataSource == PAFTOL || $DataSource == SRA || $DataSource == GAP ]]
then
type="pt_nr"
else
  echo "ERROR, invalid datasource $DataSource"
  exit 0
fi

### Directories
# All Datasource directories are expected to be in the same working directory as paftol_export and should look like OneKP/in_fasta, OneKP/blast, PAFTOL/fasta_pt, PAFTOL/fasta_nr etc.
mkdir -p $DataSource/out_blast
mkdir -p $DataSource/Barcode_Validation

### List samples to blast by DataSource. 
# Samples with existing validation cards will be omitted
python Make_samples_list.py --db $paftol_export --DataSource $DataSource

Nsamples=($(wc -l $DataSource/'Samples_to_barcode.txt'))
echo "blast $Nsamples samples on barcode databases" 
if (( $Nsamples > 0 )); then
	sbatch -p short --array=1-${Nsamples}%$slurmThrottle Blast_on_barcodes.sh $DataSource $type
fi