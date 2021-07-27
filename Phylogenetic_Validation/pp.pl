#!/bin/perl

# analyse PAFTOL tree and label nodes

# some of the outputs generated include:
# number (proportion) of monophyletic nodes at each level of the taxonomy
# proportion of monophyletic taxonomic nodes at each maximum taxonomic group size
# mean %consistency of each taxonomic head node at each maximum taxonomic group size
# taxonomic clashes (inconsistent families and orders)
# lists of outliers and aliens
# family match score for each species
# summary, rooted species and family level trees

# inputs are opt_tree, opt_dup and opt_root
#   opt_dup is an (optional) list of PAFTOL entries to disrgard
#   opt_root is the node ID of the true root: need to find this manually based on content
#     of tree 
# outputs written to various files
# order flag is used to produce summary tree for orders 
#   default is to produce summary tree for families
# well flag is to produce statistics for outliers to well-defined taxa only

use vars qw($opt_tree
			$opt_dup
			$opt_good
			$opt_bad
            $opt_alien
            $opt_outlier
            $opt_specimen 
            $opt_stats
            $opt_bootstrap
            $opt_tree2 
            $opt_order
            $opt_well
            $opt_help);
            
use Bio::TreeIO;
use Data::Dumper;
use Getopt::Long;
use sort 'stable';
use strict;

&GetOptions("tree=s",    
            "dup=s",    
            "good=s",     
            "bad=s",      
            "alien=s",    
            "outlier=s", 
            "specimen=s",
            "stats=s",
            "bootstrap=s",
            "tree2=s", 
            "order",      
            "well",
            "help");      

########################################################
# USAGE
#
my $USAGE =<<USAGE;

     Usage:

         perl pp.pl -tree treefile.nwk [-dup.dup.txt] -good g.txt -alien a.txt -outlier o.txt -specimen s.txt -stats st.txt -tree2 new_treefile.nwk [-bootstrap bs.txt] [-order] [-well] [-help] > output.txt

         where:
                tree=s  Input (unrooted) tree in Newick format
                 dup=s  Optional list of nodes to ignore
                good=s  List of all specimens not needed for manual review (if running on a pre-tree)
                 bad=s  List of badly resolving higher taxa
               alien=s  List of all specimens intruding in well-defined families
             outlier=s  List of all taxa outlying an ancestral taxon
            specimen=s  List of score for how well each specimen matches to its family
               stats=s  List of taxonomic coherence scores for labelled nodes and associated bootstrap values
           bootstrap=s  List of bootstrap values by position in taxonomy
               tree2=s  File to which a simplified (rooted) family-level tree will be written
                 order  Specify this to write a simplified order-level tree instead of a family level tree
                  well  Specify this and specimens will only be considered outliers to well-defined parents
                  help  Prints out this helpful message

USAGE
#
######################################################

if ($opt_help) {

	print "$USAGE\n";
    exit 0;
}

my %barred;

# filter out a list of accessions to be ignored

if ($opt_dup) {
	
	open (IN, $opt_dup) || die "Could not open $opt_dup\n";

	while (<IN>) {
	
		# each row of the file ends with the identifier of a PAFTOL sample that should
		# not have been included in the tree

		chomp;
		my ($id) = /(\d+)$/;
		$barred{$id} = 1;
	}
}

my $treeio = Bio::TreeIO->new(-format => 'newick', 
                              -file => $opt_tree, 
                              -internal_node_id => 'bootstrap');

my @parents;                        
my $id;
my $maxId = 0;
my %term2leaf;
my %nodeId2term;
my %nodeId2node;
my %nodeId2specimen;
my %nodeId2ancestorId;
my %nodeId2descendantId;
my %taxon2ancestor;
my %leafId2specimen;
my %specimen2leafId;
my %taxon2labelledNodeId;

my @ranks = ('ORDER', 'FAMILY', 'GENUS', 'SPECIES');

if (my $tree = $treeio->next_tree) {

 	my @leaves = $tree->get_leaf_nodes;
 	
 	print STDERR "Raw leaf count: ", scalar @leaves, "\n";
 	
 	# parse tree leaves, and record: 
 	# for each node, the number of species under each taxonomic classification 
 	# attached to it (will be 1 in each case in the beginning as we are only looking at 
 	# leaves)
 	# for each higher level taxonomic term, associate a list of species names
 	# and create a list of parents to be looked at next
 	
 	my $duplicateCount = 0;
 	my %duplicateSpecies;
 	my ($a, $b, $c, $d, $e, $f); # variables to count numbers of previously blacklisted
 	                             # samples
 	
 	LEAF: for my $leaf (@leaves) {
  	
  		my $description = $leaf->id;
  		
  		# $description =~ s/_sp\./ sp/g; # can't use field seperator within fields
  		
  		my ($order, $family, $genus, $species, $leafId, $blacklist) = 
  			split '_', $description;
  		
  		if (($species eq '') || ($order eq 'unknown')) {
  		
			die if $order eq 'unknwown'; # should have fixed all these problems in latest tree
		}
		
		if ($barred{$leafId}) {
		
			$leaf -> id($leafId); # update leaf ID so we can ignore these ones later
			next LEAF;
		}
		
		# had asked Paul to add trailing digits to duplicates but have now decided that
		# this is a bad idea - better to discover duplicates myself
		
		$species =~ s/\d+$//;	
			
		# count duplicates
		
  		my $species = $genus . ' ' . $species;
  		
  		if ($specimen2leafId{$species}) {
  		
  			$duplicateCount ++;
  			$duplicateSpecies{$species} = 1;
  		}
  		
  		# create a specimen name to associate with each leaf - will be the same as the
  		# species name except where there are multiple specimens per species
  	
  		my $specimen = $species;
  		
  		while ($specimen2leafId{$specimen}) {
  		
  			$specimen .= 'X';
  		}
  		
  		$specimen .= $blacklist;
  		
  		## stats on previously blacklisted samples
  		
  		if ($leafId =~ /^\d+$/) {
  		
  			# PAFTOL samples
  		
  			if ($blacklist eq '') {
  			
  				$a++;
  				
  			} else {
  			
  				$b++;
  			}
  		
  		} elsif ($leafId =~ /^\D+$/) {
  		
  			# 1KP samples
  			
  			if ($blacklist eq '') {
  			
  				$c++;
  			
  			} else {
  			
  				$d++;
  			}
  		
  		} else {
  		
  			# ENA samples
  		
  			if ($blacklist eq '') {
  			
  				$e++;
  			
  			} else {
  			
  				$f++;
  			}
  		}
  		
  		print STDERR "Duplicate leaf ID:", $leafId, "\n" if $nodeId2node{$leafId};
  		die if $nodeId2node{$leafId};
  		$leafId .= 'X' if $nodeId2node{$leafId}; 		
  		$maxId = $leafId if $maxId < $leafId;
  		$leaf -> id($leafId); # original metadata gets overwritten by numerical ID
  		$nodeId2node{$leafId} = $leaf;
  		${$nodeId2specimen{$leafId}}{$specimen} = 1;
  		$leafId2specimen{$leafId} = $specimen;
  		$specimen2leafId{$specimen} = $leafId;
		${${${$nodeId2term{'SPECIES'}}{$leafId}}{$species}}{$specimen} = 1;
  		${${${$nodeId2term{'GENUS'}}{$leafId}}{$genus}}{$specimen} = 1;
		${${${$nodeId2term{'FAMILY'}}{$leafId}}{$family}}{$specimen} = 1;
		${${${$nodeId2term{'ORDER'}}{$leafId}}{$order}}{$specimen} = 1;
  		push @{${$term2leaf{'SPECIES'}}{$species}}, $leaf;
  		push @{${$term2leaf{'GENUS'}}{$genus}}, $leaf;
  		push @{${$term2leaf{'FAMILY'}}{$family}}, $leaf;
  		push @{${$term2leaf{'ORDER'}}{$order}}, $leaf;
  		my $parent = $leaf -> ancestor;
  		push @parents, $parent;
  		  		
  		# each leaf is the "labelled node" for the specimen
  		
  		${$taxon2labelledNodeId{'SPECIMEN'}}{$specimen} = $leafId;
  		
  		# store taxonomic relationships
  		
  		${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'SPECIES'} = $species;
  		${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'GENUS'} = $genus;
  		${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'FAMILY'} = $family;
  		${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'ORDER'} = $order;
  		${${$taxon2ancestor{'SPECIES'}}{$species}}{'GENUS'} = $genus;
  		${${$taxon2ancestor{'SPECIES'}}{$species}}{'FAMILY'} = $family;
  		${${$taxon2ancestor{'SPECIES'}}{$species}}{'ORDER'} = $order;
  		${${$taxon2ancestor{'GENUS'}}{$genus}}{'FAMILY'} = $family;
  		${${$taxon2ancestor{'GENUS'}}{$genus}}{'ORDER'} = $order;
  		${${$taxon2ancestor{'FAMILY'}}{$family}}{'ORDER'} = $order;  
    }
    
    print STDERR "Accepted node count: ", scalar keys %nodeId2node, "\n";
    print STDERR join "*", $a, $b, $c, $d, $e, $f, "\n";
    
    for my $rank (keys %term2leaf) {
    
    	print STDERR $rank, "\t", scalar keys %{$term2leaf{$rank}}, "\n";
    }
    
    print STDERR "Duplicates: ", 
    			  $duplicateCount, 
    			  "\tDuplicate species: ", 
    			  scalar keys %duplicateSpecies, "\n";
    
    # work through tree from tips upward, mapping internal nodes to data from their 
    # leaves
    
    $id = $maxId + 1;
      	
  	PARENT: while (scalar @parents > 0) {
  	
  		my $parent = shift @parents;		
  		if ($parent -> id ne '') {
  		
  			next PARENT; # node has already been done
  			
  		} else {
  		
  		 my $check;
  		
  			for my $child ($parent -> each_Descendent) {
  		
  		 		if ($child -> id eq '') {
  		 	
  		 			# node has unseen children
  		 			# handle child first then try parent again
  		 	
  		 			unshift(@parents, $child, $parent);
  		 			next PARENT;
  		 		}
  			}
  		
  			my $text;
  			$parent -> id($id); # assign an ID for this node
  			
  			if ($check == 1) {
  			
  				print STDERR $id , "\n";		
  			}
  					
  			$nodeId2node{$id} = $parent;
  			
  			# transfer child's data to parent
  			
  			for my $child ($parent -> each_Descendent) {
  
  				my $childId = $child -> id;
  				
  				for my $rank (@ranks) {
  				
  					for my $term (keys %{${$nodeId2term{$rank}}{$childId}}) {
  				
  						for my $specimen 
  						  (keys %{${${$nodeId2term{$rank}}{$childId}}{$term}}) {
  				
  							# transfer terms associated with child node to parent node
  							# also keep track of which leaves originated these terms
  							
  							${${${$nodeId2term{$rank}}{$id}}{$term}}{$specimen} = 1;
  							
  							# a hash of all leaves mapped to the ID. Use this later to 
                            # measure number of leaves per node
  							
  							${$nodeId2specimen{$id}}{$specimen} = 1;
  						 }
  					}
  				}
  				
  				# keep track of all ancestor-descendent relationships

  				${$nodeId2ancestorId{$childId}}{$id} = 1;
  				${$nodeId2descendantId{$id}}{$childId} = 1;
  				
  				for my $grandchildId (keys %{$nodeId2descendantId{$childId}}) {
  				
  					${$nodeId2ancestorId{$grandchildId}}{$id} = 1;
  					${$nodeId2descendantId{$id}}{$grandchildId} = 1;
  				}
  			}
  			
  			$id ++;
  			
  			if (defined $parent -> ancestor) {
  			
  				unshift @parents, $parent -> ancestor;
  			}
  		}		
  	}
  	 	
  	# print out data per taxonomic term
  	
  	my (%mono, %para, %taxonSize2phylyCount, %nodeCoverageByTaxon, %singleMember);
  	my %labelledNodeId2taxon;
  	print join "\t", 'RANK', 'NAME', 'SIZE', 'M/P', 'HEAD NODE MATCH';
    print "\n";
    my %singletons;
    
	for my $rank (sort keys %term2leaf) {
	
		print STDERR $rank, "\n";
		
		TERM: for my $term (sort keys %{$term2leaf{$rank}}) {
	
			my @leaves =  @{${$term2leaf{$rank}}{$term}};			
			my $taxonSize = scalar @leaves; # list of leaf nodes associated with this term
			
			# statistics are only meaningful where we have more than 1 leaf per node
			
			if ($taxonSize >= 2) { 
								
				print $rank, "\t", $term, "\t", $taxonSize, "\t";
				my $lca = $tree->get_lca(-nodes => \@leaves); # find LCA for this term
				
				print STDERR $term, "\t", $lca;
					
				if ($lca -> id eq '') {
				
					# how does this happen?
					
					print STDERR "Help!\n";
				}
			
				# how many terms does this node map to?  If 1, that term is monophyletic
				# print this and store answer by size of LCA
					
				if (scalar keys %{${$nodeId2term{$rank}}{$lca -> id}} == 1) {
		
					print "M";
					$mono{$rank} ++;
			 		${${$taxonSize2phylyCount{$rank}}{$taxonSize}}{'M'} ++;
		
				} else {
		
					print "P";
					$para{$rank} ++;
			 		${${$taxonSize2phylyCount{$rank}}{$taxonSize}}{'P'} ++;
				}
			
				# what proportion of nodes under the LCA belong to the specified taxon?
				# 1 if monophyletic
				# print this and store answer by size of LCA
								
				my $nodeCoverageByTaxon = $taxonSize / 
										  scalar keys %{$nodeId2specimen{$lca -> id}}; 
				
				print "\t", $nodeCoverageByTaxon;
				push @{${$nodeCoverageByTaxon{$rank}}{$taxonSize}}, $nodeCoverageByTaxon;
				print "\n";
			
			} else {
			
				# later on in the program, we need labelled nodes for all genera,
				# families and orders - not just those represented by > 1 species
				
				# we can directly label their leaf nodes - whereas for those taxa with 
				# > 1 species we will need to select one of their ancestral nodes to
				# label (will do this later) 
				
				my $leaf = ${${$term2leaf{$rank}}{$term}}[0];
				my $leafId = $leaf -> id;
			 	push @{${$labelledNodeId2taxon{$rank}}{$leafId}}, $term;	
			    ${$taxon2labelledNodeId{$rank}}{$term} = $leafId;
			    ${$singletons{$rank}}{$term} = 1;
			}
		}
	}
	
	print STDERR "SINGLETONS\n";
	
	for my $rank (keys %singletons) {
	
		print STDERR $rank, "\t", scalar keys %{$singletons{$rank}};
		print STDERR "\n";
		
		if ($rank eq 'ORDER') {
		
			for my $taxon (keys %{$singletons{$rank}}) {
				
				print STDERR $taxon, "\t";
			}
			
			print STDERR "\n";
		}
	}

	# print out summary statistics re monophyly
	
	print "\n";
	print join "\t", 'RANK', 'M', 'P';
	print "\n";
		
	for my $rank (keys %mono) {

		print join "\t", $rank, $mono{$rank}, $para{$rank};
		print "\n";
	}
	
	# print out summary statistics for taxa grouped by size
	# i.e. what are the proportion of monophyletic taxa, 
	# and the average coverage of LCAs by the taxa they are LCAs of?
			
	print "\n";
	print join "\t", 'RANK', 'SIZE', 'PROPORTION M', 'AV. HEAD NODE COVERAGE';
	print "\n";

	for my $rank (keys %taxonSize2phylyCount) {
	
		# for each size of taxon (sorted in order), calculate cumulative statistics 
		# (for monophyly) for this and all smaller sizes
	
		my %cumulativeTaxonCount; # records cumulative count;
		my @nodeCoverageByTaxonList;
		my %taxonSize2nodeCoverageByTaxon;
	
		my @sizes = sort {$a <=> $b} keys %{$taxonSize2phylyCount{$rank}};
		
		for (my $n = 0; $n < scalar @sizes; $n++) {
		
			if ($n == 0) {
			
				# smallest size
			
				${$cumulativeTaxonCount{$sizes[$n]}}{'M'} = 
					${${$taxonSize2phylyCount{$rank}}{$sizes[$n]}}{'M'};
					
				${$cumulativeTaxonCount{$sizes[$n]}}{'P'} = 
					${${$taxonSize2phylyCount{$rank}}{$sizes[$n]}}{'P'};
			
			} else {
		
				# all other sizes - keep track of total number of taxa covered in all
				# groups of this size or smaller
		
				${$cumulativeTaxonCount{$sizes[$n]}}{'M'} = 
					${$cumulativeTaxonCount{$sizes[$n - 1]}}{'M'} +
			        ${${$taxonSize2phylyCount{$rank}}{$sizes[$n]}}{'M'};	
			                             
				${$cumulativeTaxonCount{$sizes[$n]}}{'P'} =  
				    ${$cumulativeTaxonCount{$sizes[$n - 1]}}{'P'} +
			        ${${$taxonSize2phylyCount{$rank}}{$sizes[$n]}}{'P'};
			}
			
			# add the coverage statistics for all taxa of the latest size to a list 
			# already containing the coverage stats for all smaller taxa, and caclulate
			# running total
			
			push @nodeCoverageByTaxonList, @{${$nodeCoverageByTaxon{$rank}}{$sizes[$n]}};
			my $total;
			
			for my $coverage (@nodeCoverageByTaxonList) {
			
				$total += $coverage;
			}
			
			# record the average (per LCA) over the current size range
			
			$taxonSize2nodeCoverageByTaxon{$sizes[$n]} = 
				$total /
				scalar @nodeCoverageByTaxonList;
		}
		
		for my $size (sort {$a <=> $b} keys %cumulativeTaxonCount) {
		
			# calculate monophyly ratio for the current range
		
			my $monophylyRatio = ${$cumulativeTaxonCount{$size}}{'M'} / 
			                     (${$cumulativeTaxonCount{$size}}{'M'} + 
			                     ${$cumulativeTaxonCount{$size}}{'P'});
			
			# print out monophyly and coverage statistics for the range
			
			print join "\t", $rank, 
							 $size, 
							 $monophylyRatio, 
							 $taxonSize2nodeCoverageByTaxon{$size};
			
			print "\n";
		}
	}
	
	# what nodes best represent each taxa?
	# for all nodes that have >=1 leaf nodes belonging to a taxon, calculate:
	# coverage of node by taxa, and coverage of taxa by node 
	# if a taxon is monophyletic, both numbers will be 1 at the LCA
	# multiply the two statistics, and prefer node with the largest product
	
	my %term2score;
	my %term2nodeCoverageByTerm;
	my %term2termCoverageByNode;
	
	# for each rank
	
	for my $rank (keys %nodeId2term) {
	
		# for each node
	
		for my $nodeId (keys %{$nodeId2term{$rank}}) {
		
			# for each term of the specified rank that has at least 1 leaf coming under
			# that node
			
			TERM: for my $term (keys %{${$nodeId2term{$rank}}{$nodeId}}) {
			
				# if only one child comes under that node, there is no score to calculate
			
				if (scalar @{${$term2leaf{$rank}}{$term}} >=  2) {
				
					my $match = scalar keys %{${${$nodeId2term{$rank}}{$nodeId}}{$term}};
					
					my $nodeCoverageByTerm = 
						$match / (scalar keys %{$nodeId2specimen{$nodeId}});
					
					my $termCoverageByNode = 
						$match / scalar @{${$term2leaf{$rank}}{$term}};
					
					${${$term2score{$rank}}{$term}}{$nodeId} = $nodeCoverageByTerm * 
															   $termCoverageByNode;
					
					${${$term2nodeCoverageByTerm{$rank}}{$term}}{$nodeId} = 
						$nodeCoverageByTerm;
						
					${${$term2termCoverageByNode{$rank}}{$term}}{$nodeId} = 
						$termCoverageByNode;				
				}
			}
		}
	} 
	
	print "\n";
	print join "\t", 'RANK', 'TERM', 'NODE', 'SIZE', 'SCORE';
	print "\n";	
	
	# $opt_report is a report of species that don't fit good higher level taxa
		
	open (ALIEN, ">$opt_alien") || die "Could not open $opt_alien\n";
	print ALIEN join "\t",'TERM', 'ALIEN';
	open (BAD, "> $opt_bad") || die "Could not open $opt_bad\n";
	print ALIEN "\n";
	
	# file for printing taxonomic coherence metric and support values
	
	open (STATS, "> $opt_stats") || die "Could not open $opt_stats"; 
	
	my %misplaced;
	my $nn = 0;
	my %familyScoreCount;
	my %term2warning;
	my (%specimenScoreCount, %wellResolved);
	
	## label nodes and print their statistics, also look for outliers and aliens
	
	RANK: for my $rank (keys %term2score) {
	
		for my $term (keys %{$term2score{$rank}}) {
		
			next RANK if scalar @{${$term2leaf{$rank}}{$term}} < 2;
		
			# for this term, find the best node and print out its statistics
			
			my $nodeId;
		
			# sorting algorithm has oddly stopped working
			# can't figure out why
		
			# my @nodeIds = sort {((${${$term2score{$rank}}{$term}}{$b} <=> ${${$term2score{$rank}}{$term}}{$a}) || ($b cmp $a)} keys %{${$term2score{$rank}}{$term}};
			# my $nodeId = $nodeIds[0];
				
			# manual sort below until this is figured out
			
			my $nodeId;
			my $warning;
			
			for my $node (keys %{${$term2score{$rank}}{$term}}) {
			
				if (${${$term2score{$rank}}{$term}}{$node} > 
				    ${${$term2score{$rank}}{$term}}{$nodeId}) {
					
					$nodeId = $node;
					$warning = 0;
				
				} elsif (${${$term2score{$rank}}{$term}}{$node} ==
				         ${${$term2score{$rank}}{$term}}{$nodeId}) {
				       
			   		if ($node gt $nodeId) {
				   
				   		$nodeId = $node;
				   	} 
				   
				   	$warning = 1;    
				}
			}
			
			# warning if two best nodes score equally well, or if best node scores at
			# less than or equal to half
			
			if (($warning == 1) || (${${$term2score{$rank}}{$term}}{$nodeId} <= 0.5)) {
			
				print BAD join "\t", $rank,
				                     $term, 
					                 ${${$term2score{$rank}}{$term}}{$nodeId};
				
				for my $leaf (@{${$term2leaf{$rank}}{$term}}) {
			
					my $leafId = $leaf -> id;
					
					print BAD "\t" . $leafId . ': ' . $leafId2specimen{$leafId};
				}
				
				print BAD "\n";
				
				${$term2warning{$rank}}{$term} = 1;
			}
			
			# node statistics
			  							  
			print join "\t", $rank, 
			                 $term, 
			                 $nodeId, 
			                 scalar @{${$term2leaf{$rank}}{$term}},
			                 ${${$term2score{$rank}}{$term}}{$nodeId};
			  
			print "\n";
		
			# summary statistics re family coherence
			
			if ($rank eq 'FAMILY') {
				
				my $simplifiedScore = int(${${$term2score{$rank}}{$term}}{$nodeId} * 10);
				$familyScoreCount{$simplifiedScore} ++;
				
				$specimenScoreCount{$simplifiedScore} += 
					scalar @{${$term2leaf{$rank}}{$term}};;
			}
			
			# store linkage of taxa and assigned nodes
			
			push @{${$labelledNodeId2taxon{$rank}}{$nodeId}}, $term;	
			${$taxon2labelledNodeId{$rank}}{$term} = $nodeId;
			
			# if the taxon hasn't mapped onto a leaf node, print out its coherence score
			# and the bootstrap value of the corresponding node
			
			if (! $leafId2specimen{$nodeId}) {
			
				print STATS join "\t", $rank, 
			       				   	   ${${$term2score{$rank}}{$term}}{$nodeId},
			      				   	   $nodeId2node{$nodeId} -> bootstrap;
			      				   
				print STATS "\n";
			}
					
			# define some nodes as well-resolved for the purpose of analysing outliers
			
			if (${${$term2termCoverageByNode{$rank}}{$term}}{$nodeId} >= 0.9) {
			
				${$wellResolved{$rank}}{$term} = 1;
			}
			
			# list of alien specimens intruding in mostly well resolved taxa
			# note: 3 differences in handling of aliens and outliers:
			# 1. different criteria for whether ancestral node is well-resolved
			# 2. aliens are only recorded for well-resolved ancestors
			# 3. only considering aliens at specimen level
			
			if ((${${$term2nodeCoverageByTerm{$rank}}{$term}}{$nodeId} >= 0.9) &&
			    (${${$term2nodeCoverageByTerm{$rank}}{$term}}{$nodeId} < 1)) {
				
				for my $otherTerm (keys %{${$nodeId2term{$rank}}{$nodeId}}) {
				
					# if a member of a non-matching taxon of the same rank is mapped to
					# this node
				
					if ($otherTerm ne $term) {
					
						# print out all species belonging to that taxon that come under
						# this node
					
						for my $specimen 
							(keys %{${${$nodeId2term{$rank}}{$nodeId}}{$otherTerm}}) {
							
							print ALIEN join "\t", 
							                  'A', 
							                  'SPECIMEN', 
							                  $specimen, 
							                  $specimen2leafId{$specimen}, 
							                  'FAMILY (' . $otherTerm . ' *)';
							
						    print ALIEN "\n";
							
						}
					} 
				}
			}					
		}
	}
	
	# print out some summary stats from the previous section
	
	for my $score (sort {$a<=>$b} keys %familyScoreCount) {
	
		print STDERR join "\t", $score, 
		                        $familyScoreCount{$score}, 
		                        $specimenScoreCount{$score};
		                        
		print STDERR "\n";
	
	}
	
	for my $rank (keys %wellResolved) {
	
		print STDERR "Well resolved: ", 
		             $rank, 
		             "\t", 
		             scalar keys %{$wellResolved{$rank}}, 
		             "\n";
	}

	# explore taxa that have placed as outliers
	# integrate this with a report on taxa that are the only members of their families or 
	# orders (which cannot be outliers by definition)

	my %rankScore;
	my %underranked; # list of families nor placed under their orders
					 # use this later when drawing simplified tree
					 
	for my $rank (keys %term2warning) {
	
		my %rankCount;
	
		for my $term (keys %{$term2warning{$rank}}) {
		
			$rankCount{$rank} += scalar @{${$term2leaf{$rank}}{$term}};
		}
	
		print STDERR join "\t", $rank, 
								scalar keys %{$term2warning{$rank}}, 
								$rankCount{$rank};
		
		print STDERR "\n";
	}
		
	open (OUT, ">$opt_outlier") || die "Could not open $opt_outlier\n";
	open (GOOD, "> $opt_good") || die "Could not open $opt_good.txt\n";
	
	my $n = 0;
	
	for my $rank (keys %taxon2labelledNodeId) {
	
		TAXON: for my $taxon (keys %{$taxon2labelledNodeId{$rank}}) {
	
			$n++ if $rank eq 'SPECIMEN';
			my %match;
			
			# labelled node for this taxon
			
			my $node = ${$taxon2labelledNodeId{$rank}}{$taxon}; 
		
			## check againt self and all ancestral nodes for incompatible matches
			
			# start by adding self to a list of its ancestors		
			
			${$nodeId2ancestorId{$node}}{$node} = 1;
			
			my @levels;

			push @levels, 'SPECIES' if $rank eq 'SPECIMEN';
			push @levels, 'GENUS' if (($rank eq 'SPECIMEN') || ($rank eq 'SPECIES'));
			
			push @levels, 'FAMILY' if (($rank eq 'SPECIMEN') ||
			                           ($rank eq 'SPECIES') || 
			                           ($rank eq 'GENUS'));
			                           
			push @levels, 'ORDER' if (($rank eq 'SPECIMEN') || 
			                          ($rank eq 'SPECIES') || 
			                          ($rank eq 'GENUS') || 
			                          ($rank eq 'FAMILY'));
			           
			## gather statistics on which parent ranks the taxa is an outlier of
			# for each ancestral node
			
			for my $testNode (keys %{$nodeId2ancestorId{$node}}) {
	
				# does it correspond to an ancestral taxon at any relevant rank?
				
				for my $level (@levels) {
				
					my $parent = ${${$taxon2ancestor{$rank}}{$taxon}}{$level};
				
					# if $opt_well is specified, only count mistmatches to well-defined
					# parents
										
					if (($opt_well) && (! ${$wellResolved{$level}}{$parent})) {
					
						$match{$level} ++
					
					} elsif (defined ${$labelledNodeId2taxon{$level}}{$testNode}) {
				
						for my $species 
							(@{${$labelledNodeId2taxon{$level}}{$testNode}}) {
					
							if ($species eq $parent) {
						
								$match{$level} ++;
							}
						}
					}
				}
			}
			
			my $rankKey;
			my $target;
			$target = 4 if $rank eq 'SPECIMEN';
			$target = 3 if $rank eq 'SPECIES';
			$target = 2 if $rank eq 'GENUS';
			$target = 1 if $rank eq 'FAMILY';
			$rankKey += 1000 if $match{'SPECIES'} > 0;
			$rankKey += 100 if $match{'GENUS'} > 0;
			$rankKey += 10 if $match{'FAMILY'} > 0;
			$rankKey += 1 if $match{'ORDER'} > 0;
							
			while (length $rankKey < $target) {
				
				 $rankKey = '0' . $rankKey;
			}
			
			# print out a list of all outliers, and what nodes they are outliers of
			
			if ($rank ne 'ORDER') {
			
				my $reportText;
				my $mismatchCount = 0;
				my $singletonCount = 0;
				my @levels;
				push @levels, 'SPECIES' if length $rankKey == 4;
				push @levels, 'GENUS' if  length $rankKey >= 3;
				push @levels, 'FAMILY'  if  length $rankKey >= 2;
				push @levels, 'ORDER'  if  length $rankKey >= 1;
				
				# consider each ancestral taxa of the current taxon
								
				for my $level (@levels) {
				
					my $ancestor = ${${$taxon2ancestor{$rank}}{$taxon}}{$level};
					my $text;
					$text = ' *' if ${$wellResolved{$level}}{$ancestor};
					$text .= ' !' if ${$term2warning{$level}}{$ancestor};
					
					# is the taxon an outlier wrt the ancestral taxon?
					# if ancestral taxon is very poorly resolved, treat all its children
					# as outliers
					
					if ((($level eq 'SPECIES') && ($rankKey =~ /0\d\d\d$/)) ||
					    (($level eq 'GENUS') && ($rankKey =~ /0\d\d$/)) ||
					    (($level eq 'FAMILY') && ($rankKey =~ /0\d$/)) ||
					    (($level eq 'ORDER') && ($rankKey =~ /0$/)) ||
					    (${$term2warning{$level}}{$ancestor}))  {
						
						$reportText .= "\t" . 
					         	       $level . 
					         	       " (" . 
					           	       ${${$taxon2ancestor{$rank}}{$taxon}}{$level} . 
					          	       $text . 
					          	       ')';
					          	       
					    $mismatchCount++;
					    
					    if (($rank eq 'FAMILY') && ($level eq 'ORDER')) {
					    
					    	$underranked{$taxon} == 1;
					    }
					}
					
					# is the taxon the only representative of the ancestral taxon?
					# only consider this at family/order level
					# (most genera and species are only represented by a single specimen)
					
					if ((($level eq 'FAMILY') || ($level eq 'ORDER')) &&
					     (${$singletons{$level}}{$ancestor})) {
					
						$reportText .= "\t" . 
					                   $level . 
					                   " (" . 
					                   ${${$taxon2ancestor{$rank}}{$taxon}}{$level} .
					                   $text . 
					                   ' 1)';  
					                   
					    $singletonCount ++;          			                  
					}
				}
				
				my $problem = '-'; # keep track of cases to review or reject
				
				if ($reportText ne '') {
				
					# next action
					
					if (($rank eq 'SPECIMEN') &&
					    ($reportText =~ /FAMILY/)) {
					    
						my $family = ${${$taxon2ancestor{$rank}}{$taxon}}{'FAMILY'};
						my $order = ${${$taxon2ancestor{$rank}}{$taxon}}{'ORDER'};
					    
					    if (${$singletons{'FAMILY'}}{$family} == 1) {
					        
					        if ($reportText !~ /ORDER/) {
					        
					        	$problem = 'C'; # can't confirm at family level but order is
					                        	# fine
					                        	
					        } elsif ((${$singletons{'ORDER'}}{$order} == 1) ||
					                 (${$term2warning{'ORDER'}}{$order} == 1)) {
					        
					        	$problem = 'U'; # can't resolve at order level either
					        
					        }  else {
					        
					        	$problem = 'R';
					        }
					    
					    } elsif (${$term2warning{'FAMILY'}}{$family}) {
					    
							$problem = 'U'; # low scoring taxa, review all members
					    
					    } else {
					    
					    	$problem = 'R'; # straightforward reject
					    }
					}
				
					my $errorKey;
					$errorKey = $problem;
					$errorKey .= 'O' if $mismatchCount > 0;
					$errorKey .= 'S' if $singletonCount > 0;
					
					if (${$term2warning{$rank}}{$taxon}) {
					
						$taxon .= ' !';
					}
					
					print OUT join "\t", $errorKey, 
										 $rank, 
										 $taxon, 
										 $specimen2leafId{$taxon};
					
					print OUT $reportText, "\n";
				} 
				
				# print out a list of good specimens
				
				if (($rank eq 'SPECIMEN') && (($problem eq '-') || ($problem eq 'C'))) {
				
					# list of all specimens which are not outliers of their families,
					# not the only member of their family
					
					my $family = ${${$taxon2ancestor{$rank}}{$taxon}}{'FAMILY'};
					my $familyText = 'FAMILY (' . $family;
					$familyText .= ' *' if ${$wellResolved{'FAMILY'}}{$family};
					$familyText .= ')';         
					
					print GOOD join "\t", "G",
										  $rank,
										  $taxon,
										  $specimen2leafId{$taxon},
										  $familyText;
					print GOOD "\n";
										  
				}
			}
			
			${$rankScore{$rank}}{$rankKey} ++;
		}
	}
	
	print STDERR $n, "\n";

	# print summary of outlier classes
	
	print STDERR "TAXON SCORES BY RANK\n";
	
	for my $rank (keys %rankScore) {
	
		for my $key (sort keys %{$rankScore{$rank}}) {
		
			print STDERR join "\t", $rank, $key, ${$rankScore{$rank}}{$key}, "\n";
		}
	}
	
	# assign a family match score for each specimen - this is the score of the best
	# scoring ancestral node (wrt to the family) for the taxon
	
	my %specimen2score;
	
	SPECIMEN: for my $specimen (keys %specimen2leafId) {
	
		my $family = ${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'FAMILY'};
		next SPECIMEN if scalar @{${$term2leaf{'FAMILY'}}{$family}} < 2;
		my $leafId = $specimen2leafId{$specimen};
		my $bestScore = 0;
	
		for my $ancestor (keys %{$nodeId2ancestorId{$leafId}}) {
	
			my $score = ${${$term2score{'FAMILY'}}{$family}}{$ancestor};
			$bestScore = $score if $score > $bestScore;
		}
		
		$specimen2score{$specimen} = $bestScore;
	}
	
	# print out a score for how well each mis-fitting species fits to its family
	
	open (SPECIMEN, ">$opt_specimen") || die "Could not open $opt_specimen\n";
	
	for my $specimen (sort {$specimen2score{$a} <=> $specimen2score{$b}} 
					  keys %specimen2score) {
		
		my $family = ${${$taxon2ancestor{'SPECIMEN'}}{$specimen}}{'FAMILY'};			 
		my $familyNode = ${$taxon2labelledNodeId{'FAMILY'}}{$family};
		my $familyScore = ${${$term2score{'FAMILY'}}{$family}}{$familyNode};
		
		if ($specimen2score{$specimen} != $familyScore) {
	
			print SPECIMEN join "\t", $specimen, 
								 	  $family, 
								 	  $specimen2score{$specimen},
								 	  $familyScore;				 
		
			print SPECIMEN "\n";
		
		} 
	}
	
	## print out node sets (for use in figure)
	# classes are sub-generic, generic, sub-familial, familial, sub-ordinal, ordinal, 
	# supra-ordinal
	
	my @ranks = ('GENUS', 'FAMILY', 'ORDER');
	
	if ($opt_bootstrap) {
	
		open (OUT, "> $opt_bootstrap") || die "Could not open $opt_bootstrap";
		my %recorded;
		
		for my $rank (@ranks) {
	
			for my $term (keys %{$taxon2labelledNodeId{$rank}}) {
		
				my $nodeId = ${$taxon2labelledNodeId{$rank}}{$term};
			
				for my $descendantNodeId (keys %{$nodeId2descendantId{$nodeId}}) {
				
					#print STDERR $descendantNodeId , "\n";
			
					if ((! $leafId2specimen{$descendantNodeId}) && 
				    	(! $recorded{$descendantNodeId}) &&
				    	(! $barred{$descendantNodeId})) {
				    	
				    	my $score = $nodeId2node{$descendantNodeId} -> bootstrap;
						print OUT 'SUB-' . $rank, "\t", $score, "\n";
						$recorded{$descendantNodeId} = 1;
					}
				}
			}	
		
			for my $term (keys %{$taxon2labelledNodeId{$rank}}) {
		
				my $nodeId = ${$taxon2labelledNodeId{$rank}}{$term};
			
				if ((! $leafId2specimen{$nodeId}) && 
					(! $recorded{$nodeId}) &&
					(! $barred{$nodeId})) {
			
					my $score = $nodeId2node{$nodeId} -> bootstrap;
					print OUT $rank, "\t", $score, "\n";
					$recorded{$nodeId} = 1;
				}
			}
		}
	
		# supra-ordinal nodes
	
		my $root = $tree -> get_root_node;
		my $rootId = $root -> id;
	
		for my $descendantNodeId (keys %{$nodeId2descendantId{$rootId}}) {
	
			if ((! $leafId2specimen{$descendantNodeId}) && 
				(! $recorded{$descendantNodeId}) &&
				(! $barred{$descendantNodeId})) {
	 	
	 			my $score = $nodeId2node{$descendantNodeId} -> bootstrap;
	 			print OUT "SUPRA-ORDER\t", $score, "\n" if $score ne ''; 		
	 		}
		}
	}
	
	## print out a rooted tree containing only families and orders
	# $opt_order allows families to also be excluded
	
	# first step, find labelled nodes for the included taxa and prepare a description
	# for each to be shown in the tree

	my %reducedTree;
	my %descendentCount;
	my $root;
	my %nodeId2text;
	my $outlierNodeId;
	my ($inconsistent, $consistent) = 0;
	
	if (! $opt_order) {
	
		for my $nodeId (keys %{$labelledNodeId2taxon{'FAMILY'}}) {
	
			my @texts;
		
			for my $term (@{${$labelledNodeId2taxon{'FAMILY'}}{$nodeId}}) {
		
				# for families not found under their orders, use a special separator
				
				my $separator = '_';
		
				if ($underranked{$term}) {
			
					$separator = '!';
					$inconsistent++;
				
				} else {
				
					$consistent ++;
				}
		
				my $text = $term . $separator;
				
				# score for the labelled node (rounded down to the nearest 0.1)
				# will be empty for single member families for which no score will have
				# been calculated)
				
				$text .= ((int (10* ${${$term2score{'FAMILY'}}{$term}}{$nodeId})) / 10) 
			  	  if defined ${$term2score{'FAMILY'}}{$term};		          
			
				push @texts, $text;
			}
		
			$nodeId2text{$nodeId} .= join ',', @texts; 
			
			print STDERR $nodeId, "\t", $nodeId2text{$nodeId}, "\n";
		}
	}
	
	# now do something similar for orders
	
	for my $nodeId (sort keys %{$labelledNodeId2taxon{'ORDER'}}) {
	
		my @texts;
		
		for my $term (@{${$labelledNodeId2taxon{'ORDER'}}{$nodeId}}) {
			
			my $text = uc $term . '_';
			
			$text .= ((int (10* ${${$term2score{'ORDER'}}{$term}}{$nodeId})) / 10) 
			  if defined ${$term2score{'ORDER'}}{$term};   
			        
			push @texts, $text;
			$outlierNodeId = $nodeId if $text eq 'AMBORELLALES_';
		}
		
		$nodeId2text{$nodeId} .= join ',', @texts; 
	}
	
	# create untipped tree (i.e. lose nodes downstream of interesting nodes)
	# by working up from labelled nodes identified in previous section
	# when we hit a node on a previously followed path, boost its count then stop
	# traversing
	# thus, nodes with a count > 1 are junction nodes between paths
	
	my $n = 1;
	
	LABEL: for my $nodeId (sort keys %nodeId2text) {
		
		$n++;
		
		while (! $reducedTree{$nodeId}) {
		
			my $node = $nodeId2node{$nodeId};
		
			if (! defined $node -> ancestor) {

				$root = $nodeId;
				next LABEL;
			}
	
			my $ancestorId = $node->ancestor->id;
			$reducedTree{$nodeId} = $ancestorId;
			
			$nodeId = $ancestorId;
			$descendentCount{$nodeId} ++;
					
		}	
	}
	
	print STDERR ($n -1), "\t", scalar keys %reducedTree, "\n";
	
	# parse unwanted nodes from tree 
	# work up from labelled nodes and remove ancestors that are neither labelled nor 
	# junctions
	
	my %done;
	
	NODE: for my $nodeId (sort keys %nodeId2text) {
		
		while ($reducedTree{$nodeId}) {
			
			next NODE if $done{$nodeId};
		
			if (my $parent = $reducedTree{$nodeId}) {
			
				# parent node is not a junction and only has 1 descendent
		
				if (($descendentCount{$parent} < 2) && 
			    	((! defined ${$labelledNodeId2taxon{'FAMILY'}}{$parent}) || 
			    	 ($opt_order)) && 
			     	(! defined ${$labelledNodeId2taxon{'ORDER'}}{$parent})) {
					
					# splice parent from tree
					
					my $grandparent = $reducedTree{$parent};
					$reducedTree{$nodeId} = $grandparent;
					delete $reducedTree{$parent};
				
				} else {
					
					$done{$nodeId} = 1;
			
					if ($reducedTree{$parent}) {
					
						$nodeId = $parent;
					
					} else {
					
						# parent is root, we're done 
						
						next NODE;
					}
				}
			
			} else {
			
				# in case the root has been labelled
				
				$done{$nodeId} = 1;
				next NODE;
			}
		}
	}

	# parse tree a third time, create reverse links (parent -> child)
	
	my %reducedTreeReversed;
	undef %done;
	
	NODE: for my $nodeId (sort keys %reducedTree) {
		
		push @{$reducedTreeReversed{$reducedTree{$nodeId}}}, $nodeId;
	}
		
	## now create a properly rooted tree
	
	# in a tree with non-angiosperm outliers, the root should be inserted between the 
	# angiosperms and the gymnosperms, and these should only be linked to each other via
	# the root
	
	my %revisedTreeReversed;
	my $trueRootNodeId = 10000;
	my $angiospermNodeId = 18703;
	my $gymnospermNodeId = 18702;
	my @newList;	
	push @{$reducedTreeReversed{$trueRootNodeId}}, $angiospermNodeId, $gymnospermNodeId;
	
	undef $reducedTree{$angiospermNodeId} 
		if $reducedTree{$angiospermNodeId} eq $gymnospermNodeId;
		
	undef $reducedTree{$gymnospermNodeId} 
		if $reducedTree{$gymnospermNodeId} eq $angiospermNodeId;
		
	my @newAngiospermList;
	
	for my $nodeId (@{$reducedTreeReversed{$angiospermNodeId}}) {
	
		push @newAngiospermList, $nodeId if $nodeId ne $gymnospermNodeId;
	}
	
	undef  @{$reducedTreeReversed{$angiospermNodeId}};
	push @{$reducedTreeReversed{$angiospermNodeId}}, @newAngiospermList;
	
	my @newGymnospermList;
	
	for my $nodeId (@{$reducedTreeReversed{$gymnospermNodeId}}) {
	
		push @newGymnospermList, $nodeId if $nodeId ne $angiospermNodeId;
	}
	
	undef @{$reducedTreeReversed{$gymnospermNodeId}};
	push @{$reducedTreeReversed{$gymnospermNodeId}}, @newGymnospermList;
	
	# now work through the tree starting from the new root
	# both upwards and downwards relationships wrt the original root are downwards 
	# relationships wrt the new root
	
	my @tips = ($trueRootNodeId);
	
	while (scalar @tips > 0) {
		
		my $tip = shift @tips;
		
		print STDERR "TIP: ", $tip, "\n";
		
		my @newTips;
		
		if (defined $reducedTreeReversed{$tip}) {
		
			push @newTips, @{$reducedTreeReversed{$tip}};
			
			print STDERR "Adding from reversed tree ", join " ", @{$reducedTreeReversed{$tip}}, "\n";
		}
		
		push @newTips, $reducedTree{$tip} if $reducedTree{$tip};
		
		print STDERR "Adding from tree ", $reducedTree{$tip}, "\n";
		
		for my $newTip (@newTips) {
					
			# because we are looking up and down, need to make sure we don't go backwards
		
			if (! $revisedTreeReversed{$newTip}) {
			
				 push @{$revisedTreeReversed{$tip}}, $newTip;	
				 push @tips, $newTip;
		    }
		}
	}
		
	# print out tree in Newick format
	
	my $string = _processNode($trueRootNodeId, '', \%revisedTreeReversed, \%nodeId2text);
	open (TREE, "> $opt_tree2") || die "Could not open $opt_tree2\n";
	print TREE $string, ";\n";
}

sub _processNode($$$) {

	# simple recursive algorithm for printing out Newick tree

	my ($node, $string, $children, $nodeId2text) = @_;
	my %children = %$children;
	my %nodeId2text = %$nodeId2text;
	
	if ((! defined $children{$node}) || (scalar @{%children{$node}} == 0)) {
	
		# node has no children, or node has no children left
	
		$string .= ',' if (($string ne '') && ($string !~ /\($/));
		
		if ($nodeId2text{$node}) {
		
			$string .= $nodeId2text{$node};
		
		} else {
		
			$string .= $node;
		}
	
	} else {
	
		$string .= ',' if (($string ne '') && ($string !~ /\($/));
		$string .= '(';
		
		# deal with the node's children before returning to the node
		
		while (scalar @{$children{$node}} > 0) {
		
			my $child = shift @{$children{$node}};
			$string = _processNode($child, $string, $children, $nodeId2text);
		}
		
		$string .=  ")";
		
		if ($nodeId2text{$node}) {
		
			$string .= $nodeId2text{$node};
		
		} else {
		
			$string .= $node;
		}
	}
	
	return $string;
}

__END__
