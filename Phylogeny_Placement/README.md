# Analyse PAFTOL tree and label nodes

## Usage
**Inputs:**
* opt_tree 		Input tree (required)
* opt_dup 		is an (optional) list of PAFTOL entries to disregard
* opt_root 		is the node ID of the true root. Need to find this manually based on content of tree 

**Outputs:**
* default output is to produce summary tree for families
* number (proportion) of monophyletic nodes at each level of the taxonomy
* proportion of monophyletic taxonomic nodes at each maximum taxonomic group size
* mean %consistency of each taxonomic head node at each maximum taxonomic group size
* taxonomic clashes (inconsistent families and orders)
* lists of outliers and aliens
* family match score for each species
* summary, rooted species and family level trees


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

```
perl pp.pl -tree treefile.nwk [-dup.dup.txt] -good g.txt -alien a.txt -outlier o.txt -specimen s.txt -tree2 new_treefile.nwk [-root 100000] [-order] [-well] [-help] > output.txt
```
