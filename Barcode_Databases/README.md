# DNA Barcode Validation - Barcode Databases

Barcode databases were built from BOLD and NCBI repositories.

## BOLD databases

The BOLD Database was downloaded from bold (accessed on 03/08/2021), and sequences for cpDNA rbcL, rbcLa, and matK as well as rDNA ITS2 were extracted using a custom notebook (`Processing_BOLD.ipynb`).

## NCBI databases

NCBI nucleotide database was queried on 30/07/2021 for the following barcoding loci: cpDNA 23s, cpDNA 16s, rDNA 18s. 16s, trnH-psbA and trnL, and downloaded in GenBank format, with the following queries:

```shell
sbatch ncbi_query.sh '"18S ribosomal RNA"[All Fields] OR "rrn18"[All Fields] AND "Spermatophyta"[Organism] AND ("1000"[SLEN] : "300000"[SLEN])' NCBI_18s

sbatch ncbi_query.sh '"23S ribosomal RNA"[All Fields] OR "rrn23"[All Fields] AND "Spermatophyta"[Organism] AND ("0"[SLEN] : "300000"[SLEN]) AND chloroplast[filter]' NCBI_23s

sbatch ncbi_query.sh '"16S ribosomal RNA"[All Fields] OR "rrn16"[All Fields] AND "Spermatophyta"[Organism] AND ("0"[SLEN] : "300000"[SLEN]) AND chloroplast[filter]' NCBI_16s
```

Genbank files were processed in a custom script `GB_extract.py`, in which genes or rRNA were extracted. All NCBI references were filtered based on length, with a minimum and maximum length set for each barcode. See the beginning of the script for filter values.

Finally, we added the most recent release of plastid data as a reference database of whole plastomes (https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/).

## Taxonomy checks against WCVP
The taxonomy of accessions was resolved against WCVP (last accession: [wcvp_v5_jun_2021.zip](http://sftp.kew.org/pub/data-repositories/WCVP/wcvp_v5_jun_2021.zip)) using our custom script [WCVP_taxo](../WCVP_Taxo/), and written as a new .fasta file and list of accessions. Scientific names in genus sp. format were resolved (option -g), as well as scientific names for which duplicate entries all mapped to the same genus (option -d divert_genusOK). Sequences with unresolved names or that matched duplicate entries with different genera names were discarded. WCVP database is http://sftp.kew.org/pub/data-repositories/WCVP/wcvp_v5_jun_2021.zip

Note that a maximum of two accessions per species were kept in each reference database.

The fasta file is accompanied by a list of accessions containing Accession ID, organism name and taxonomic ID (*_TAXO.csv files).