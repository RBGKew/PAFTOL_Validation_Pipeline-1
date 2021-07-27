# Analyse PAFTOL tree and label nodes

## Method

The program has two main tasks: First: to identify the nodes corresponding to higher taxa, such as family or order; Second: to confirm or reject samples at these taxonomic ranks, based on their placement within or outside these nodes, respectively. 

1. For every node in the tree, two metrics were calculated: (i) the proportion of samples belonging to a taxon (e.g. family) that subtended the node, and (ii) the proportion of samples subtending the node that belonged to that taxon. The two metrics were then multiplied to produce a combined score. For each taxon, the highest scoring node was subsequently considered to best represent the taxon in the tree. Each node corresponding to a taxon were categorized as (i) well-resolved when reaching a score of 1, (ii) resolved when <1 but >0.5, (iii) poorly resolved if the score was <0.5. 

2. The phylogenetic validation of family identification of each sample was determined as:

   a.   Confirmed: if identified as belonging to a family whose best scoring node had a value >0.5 and found as a descendant of this node in the tree.

   b.   Rejected: if identified as belonging to a family whose best scoring node had a value >0.5 but not found as a descendant of this node.

   c.   Inconclusive: if belonging to a family whose best scoring node had a value ≤0.5.

Note that for families represented in the tree by a single sample, the validation was performed with respect to their orders. If the order was represented by a single sample, the sample was considered untestable and coded as inconclusive.

## Usage

**Inputs:**

* opt_tree 		Input tree (required)
* opt_dup 		is an (optional) list of PAFTOL entries to disregard
* opt_root 		is the node ID of the true root. Need to find this manually based on content of tree 

```perl
perl pp.pl -tree treefile.nwk [-dup.dup.txt] -good g.txt -alien a.txt -outlier o.txt -specimen s.txt -tree2 new_treefile.nwk [-root 100000] [-order] [-well] [-help] > output.txt
```


**Options:**

* tree=s  		Input (unrooted) tree in Newick format
* dup=s  		Optional list of nodes to ignore
* good=s  		List of all specimens not needed for manual review (if running on a pre-tree)
* bad=s  		List of badly resolving higher taxa
* alien=s  		List of all specimens intruding in well-defined families
* outlier=s  	List of all taxa outlying an ancestral taxon
* specimen=s 	List of score for how well each specimen matches to its family
* tree2=s  		File to which a simplified (rooted) family-level tree will be written
* root=i  		Node to be used to root the new tree
* order  		Specify this to write a simplified order-level tree instead of a family level tree
* well  		Specify this and specimens will only be considered outliers to well-defined parents
* help  		Prints out usage

**Outputs:**

The program produces:

* A list of outliers (i.e. samples not found under the node best describing their family) and aliens (i.e. samples found in a resolved family but not belonging to it). Note that by default the program produces a summary of the tree for families, but it can used for other taxonomic levels (e.g. order, genus)
* Taxonomic clashes (inconsistent families and orders) 
* The proportion of monophyletic nodes at each level of the taxonomy
* The proportion of monophyletic nodes at each level of the taxonomy
* mean %consistency of each taxonomic head node at each maximum taxonomic group size
* family match score for each species
* summary, rooted species and family level trees

