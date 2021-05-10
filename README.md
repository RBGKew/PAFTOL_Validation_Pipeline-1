# PAFTOL Validation Pipeline
- [DNA Barcoding](#dna-barcoding)
  * [Barcode Databases](#barcode-databases)
  * [Samples](#samples)
    + [PAFTOL](#paftol)
    + [1KP](#1kp)
  * [Taxonomic standardization](#taxonomic-standardization)
  * [Sample Validation](#sample-validation)
    + [1. Blast](#1-blast)
    + [2. Test](#2-test)
    + [3. Validation](#3-validation)
  * [Dependencies](#dependencies)
- [Phylogeny Placement](#phylogeny-placement)
- [Validation decisions](#validation-decisions)

We implemented two validation steps, ran in parallel to verify the family of PAFTOL samples. The two steps are: 

* 1) an *in-silico* DNA barcoding validation, which uses plastid and nuclear ribosomal data for DNA-based identification.
* 2) a phylogenetic validation, which checks the placement of each sample at family level in a preliminary tree relative to its expected position. 


![Family_Validation](Family_Validation.jpg)

# DNA Barcoding
## Barcode Databases
To improve the testability and accuracy of the validation procedure—due to partial recovery of organellar DNA and uneven coverage across families for reference barcodes—Six individual barcode reference databases were built from the NCBI nucleotide and BOLD databases ([https://www.ncbi.nlm.nih.gov/nuccore; ](https://www.ncbi.nlm.nih.gov/nuccore)https://www.boldsystems.org/, accessed on 29/10/2020), one for the whole plastome, and the remaining five for specific loci (ribosomal 18S, as well as plastid *rbc*L, *mat*K, *trn*L, and *trn*H-*psb*A).
The creation and curation of barcode databases is described in detail [here](Barcode_Databases/)

## Samples
### PAFTOL
For PAFTOL samples, plastomes and ribosomal DNA were recovered from raw reads using `GetOrganelles` (Jin et al. 2020). In both cases, recommended parameters were used (https://github.com/Kinggerm/GetOrganelle#recipes; i.e. -R 20 -k 21,45,65,85,105 for plastomes, and -R 10 -k 35,85,115 for nuclear ribosomes). Our GetOrganelle script is in [PAFTOL_Get_Organelles](PAFTOL_Get_Organelles/)

### 1KP
As 353 target genes were also recovered from transcriptomes of the One Thousand Plant Transcriptomes Initiative (Leebens-Mack et al. 2019), we also performed validation by barcoding on 766 1KP samples. Here, we used the original transcriptome fasta files. 

## Taxonomic standardization

Species names in PAFTOL, 1KP and barcode databases were all standardized against the World Checklist of Vascular Plants (https://wcvp.science.kew.org/)  using a [custom python script](WCVP_Taxo/). As such all species reference databases and samples were standardized against the same taxonomy in the validation pipeline.

## Sample Validation
### 1. Blast

Sample sequences were queried against barcode databases using `BLASTn` (Camacho et al. 2009) if their family was present in the database. BLAST results were further filtered with a minimum identity >95%, minimum length (based on barcode length) and minimum coverage of reference locus >=90%. BLAST matches were then ranked by identity.

**Source**|**Barcode**|**Min. length**|**Min. coverage**
:-----:|:-----:|:-----:|:-----:
NCBI|18s|1000|90
NCBI|plastomes|1000|0
NCBI|trnH-psbA|200|90
NCBI|trnL|200|90
BOLD|matK|300|90
BOLD|rbcLa|300|90


```shell
sbatch blast_barcodes_PAFTOL_array.sh PAFTOL_Samples/ paftol_ls.txt
python Process_blast.py PAFTOL_Samples paftol_ls.txt

sbatch blast_barcodes_1KP_array.sh OKP/ okp_ls.txt
python Process_blast.py OKP okp_ls.txt
```

### 2. Test

Up to six barcode tests were thus executed per sample. A sample passed an individual test if the first ranked BLAST match (ranked by identity) confirmed its original family identification, and failed otherwise. Note that controls could only be completed if the specimen’s family was present in the barcode databases and if at least one BLAST match remained after filtering. 

### 3. Validation

The final result of the barcode validation following the six individual barcode tests were determined as follows:  

1. Confirmed: One or more barcode test confirm the family identification of a sample.
2. Rejected: More than ½ of the barcode tests confirm the same incorrect family identification (requires at least two barcode tests).
3. Inconclusive otherwise.

The barcode validation scripts can be found [here](Blast_&_Validation/)

## Dependencies

```python
- entrez-direct 
- seqkit
- cutadapt
- blast
- python3
	- Pandas
	- Numpy
	- Bio
	- TQDM
	- Seaborn
```

# Phylogeny Placement
To conduct phylogenetic validation, a preliminary phylogenetic tree was built using the complete, unvalidated release dataset, following the phylogenetic methods described below. We then assessed which nodes best represents each family in the tree. For every node in the tree, two metrics were calculated: firstly, the proportion of samples belonging to the taxon that subtend the node, and secondly, the proportion of samples subtending the node that belong to the taxon. The two metrics were multiplied to produce a combined score, and for each family, the highest scoring node was subsequently considered to best represent the taxon in the tree (allowing the identification of outlying samples). Where a node has a value of 1, the taxon that it represents is monophyletic in the tree.

The phylogenetic validation of family identification of each sample was determined as:

1. Confirmed: if belonging to the family whose best scoring node had a value >0.5 and found under this node, 
2. Rejected: if identified as belonging to the family whose best scoring node had a value >0.5 but not falling under the node 
3. Inconclusive: if belonging to the family whose best scoring node had a value of ≤0.5.

Note that for families represented in the a tree by a single sample, the validation was performed with respect to their orders. If the order was represented by a single sample, the validation result was coded as inconclusive.

The phylogenetic validation script can be found [here](Phylogeny_Placement/)

# Validation Decisions
The outputs of the phylogenetic and barcoding validation were combined to identify specimens for automatic inclusion and exclusion from the final tree, and where a decision on inclusion/exclusion was subject to expert review (see figure above).
