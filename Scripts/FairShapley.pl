#!/usr/bin/perl

#FairShapley.pl
#
#This script is a modified version of the script FairShapley.pl introduced in
#Wicke, K. and Fischer, M. Comparing the rankings obtained from two biodiversity indices: the Fair Proportion Index and the Shapley Value. 
#Journal of Theoretical Biology (2017), Volume 430, 207-214, doi.org/10.1016/j.jtbi.2017.07.010
#
#Compared to the original version, it additionally outputs the average FP index of each species and prints it to an output file specified by the user (lines 226-231).
#Note that the calculation of this average assumes that all input trees contain the same set of species.
#If this is not the case, the average returned might not be correct and needs to be adjusted manually.


#Modules needed for running the script
use strict;
use warnings;
use Try::Tiny;				#Error Handling
use Getopt::Long;			#Command line arguments
use Bio::TreeIO;			#Parser for tree files
use Bio::Tree::Tree; 		#An Implementation of TreeI Interface
use Bio::Tree::Node;		#A simple Tree node



#Subroutines used:
sub calc_fp;
sub calc_sv_rooted;
sub calc_sv_unrooted;
sub print_fp;
sub print_sv;




#Description of how to use the script
my $usage = <<'ENDUSAGE';
FairShapley.pl
	Compute the Fair Proportion Index and/or (modified) Shapley Value 
	for phylogenetic trees in Newick format. 
	Note, that the FP-Index is only defined for rooted trees, while 
	the SV-index can be calculated for both rooted and unrooted trees.

	Please make sure that BioPerl is installed on your machine!

SYNOPSIS:
FairShapley.pl --in=filename (--out=filename) (--values=value)
	
OPTIONS:

	--out=filename
		Generates an output file containing the computed values.

	--values=value
		Choose the diversity index to calculate.
		Options for value:
			0 Calculation of FP
			1 Calculation of SV (modified SV if the tree is rooted, original SV if the tree is unrooted)
			2 Calculation of FP and SV unrooted on rooted trees
			3 Calculation of FP and (modified) SV 
		Per default the FP and the (modified) SV are calculated (option 3).

DESCRIPTION:

	Example: FairShapley.pl --in=myTree out=myResults --values=0

ENDUSAGE



#Variables for the command line arguments
my ($input,$output, $values, $help);
$values=3; 	#per default both FP and SV are calculated


GetOptions('in=s'=>\$input,'out=s'=>\$output,'values=i'=>\$values,'help'=>\$help);

#Print Usage
if ($help){
    print $usage;
    exit(1);
}



#Throw error if no input filename is given by the user
if (!defined($input)){
    print "Missing tree input filename.\n$usage";
    exit(1);
}

#Throw error if an input filename is given, but the file doesn't exist
open(my $in, "<", $input) or die "Could not open file '$input'. $!";



#Check whether a name of an output file is given.
#If so, check further if this file already exists and let the user decide whether to overwrite or append the file or to stop the script
our $terminal = 0; 
if (!defined($output)){
    $terminal = 1;
}

elsif (-e "$output"){
	print "\nThe file $output already exists. \nIf you want to overwrite the file, type OVERWRITE,\nif you want to append the new data to the old file, type APPEND. \nElse the program will die.\n";
	my $answer = <STDIN>;
	chomp($answer);
	if($answer eq "OVERWRITE" ){
		unlink ($output) or die;	
	}
	elsif(not $answer eq "APPEND"){
		die;
	}
}


###############################################################################################################################################################

#Read in the file containing the tree(s)
my $treeio = Bio::TreeIO->new(-format => 'newick',
                            -file   => $input);


our $tree_counter = 0; 	#counts the number of trees in the input file
our %average_FP;



#Loop over all trees and process them one after another
try{
	while(our $tree = $treeio->next_tree ){

		$tree_counter ++;
        print "Tree_counter".$tree_counter."\n";
		
		#Decide whether the tree is rooted or unrooted 
		#(Bio::TreeIO always sets the top node as root node even if the tree is unrooted. 
		#To decide whether this node really is a root in the sense of Phylogenetics, simply count the children of the node.
		#If the number of children equals 2, the tree is rooted, otherwise, it's unrooted.)
		my $candidate_root = $tree->get_root_node();
		my @children_of_root = $candidate_root->each_Descendent();
		my $rooted;
		if (@children_of_root == 2){
			$rooted = 1;
		}
		else{
			$rooted = 0;
		}  	  


		#Store important information about the tree
		our @nodes = $tree->get_nodes(); 		#Array of nodes
		our @leaves = $tree->get_leaf_nodes();	#Array of leaves
		our $tree_size = @leaves; 				#Number of leaves
		our $root = $candidate_root; 			#root

	
		#Calculation of the diversity indices depending on whether the tree is rooted or not and on the options chosen
	
		#Case 1: Calculation of FP
		if($values == 0){
			if ($rooted == 1){
				my $tree_fp = calc_fp();
                my %fp_hash = %{$tree_fp};
                foreach (keys %fp_hash){
                    $average_FP{$_} += $fp_hash{$_};
                }
				print_fp($tree_fp);
			}
			else{
				print "Your input tree is unrooted. The FP-Index can't be calculated for unrooted trees!\n";
			}
		}

		#Case 2: Calculation of SV 
		elsif ($values == 1){
			if($rooted == 1){
				my $tree_fp = calc_fp();	#in the rooted case SV is derived from FP!
				my $tree_sv_rooted = calc_sv_rooted($tree_fp);
				print_sv($tree_sv_rooted);
			}
			else{
				my $tree_sv_unrooted = calc_sv_unrooted();
				print_sv($tree_sv_unrooted);
			}
		}


		#Case 3: Calculation of FP and the unrooted SV on rooted trees
		elsif ($values == 2){
			if($rooted == 1){
		
				my $tree_fp = calc_fp();
				print_fp($tree_fp);

				my $tree_sv_unrooted = calc_sv_unrooted();
				print_sv($tree_sv_unrooted);
			}
			else {
				print "\nYour input tree (tree ".$tree_counter.") is unrooted.\nTherefore neither the FP-Index nor the unrooted SV on rooted trees can be calculated!\n";
			}
		}


		#Case 4: Calculation of FP and SV (default)
		else{
			if($rooted == 1){

				my $tree_fp = calc_fp();
				print_fp($tree_fp);
		
				my $tree_sv_rooted = calc_sv_rooted($tree_fp);
				print_sv($tree_sv_rooted);
			}
			else {
				print "\nYour input tree (tree ".$tree_counter.") is unrooted. The FP-Index can't be calculated for unrooted trees!\n";

				my $tree_sv_unrooted = calc_sv_unrooted();
				print_sv($tree_sv_unrooted);
			}
		}


	}

    
	#Print the average FP index (assuming that all trees in the input contain the same set of species)
     
	open (OUTPUT, ">>$output") or die ("Could not open file $output\n");
		print OUTPUT "\nThe average FP indices for your set of ".$tree_counter. " input trees are : \n";
			foreach (keys %average_FP){
				print OUTPUT $_."\t".$average_FP{$_}/$tree_counter."\n";
			}
	close(OUTPUT);




###############################################################################################################################################################
	#Subroutine for calculation the Fair Proportion Index
	sub calc_fp{
	
		#Loop over all nodes and store their possible contribution to the FP-Indix of a leaf
		my %FP_contribution;	#Hash to store the contribution to FP for every edge
		my %fp_hash;	#Hash to store the FP-Indices for each leaf

		foreach(our @nodes){
			if($_->branch_length()){	#This test is necessary as the root node has no value assigned to it
				my @descendent_nodes = $_->get_all_Descendents(); 	#functions, which returns all descendents of a given node
				my $descendent_leaves = 0;

				#Count how many of the descendents are leaves
				foreach (@descendent_nodes){
					if ($_->is_Leaf()){
						$descendent_leaves++;
					}
				}

				#If the node has no descendent nodes, i.e. it's a leaf, its contribution is its own branch length
				if ($descendent_leaves == 0){
					$FP_contribution{$_} = $_->branch_length();
				}

				#Else the contribution of a node is its branch_length divided by the numer of descendent leaves
				else {
					$FP_contribution{$_} = $_->branch_length()/$descendent_leaves;
				}
			}
		}
			

		#Loop over all leaves to calculate the FP-Index for each leaf
		foreach (our @leaves){

			my $fp_index = 0;
			my @lineage = our $tree->get_lineage_nodes($_);
			splice @lineage,0,1;	#Remove the root, since it doesn't have a branch length assigned to it

			#Loop over all lineage nodes except the first one (root) and sum up their contribution to the FP-Index of the current leaf
			foreach (@lineage){
				$fp_index = $fp_index + $FP_contribution{$_};
			}
			
			$fp_index = $fp_index + $FP_contribution{$_}; 	#add contribution of the leaf itself
			$fp_hash{$_->id()} = $fp_index;
		}
	
		return  \%fp_hash;
	}		

###############################################################################################################################################################
	#Subroutine for calculation of the Shapley Value of rooted trees (derived from the FP-Index: SV(a) = FP(a)-PD(a)/n)
	sub calc_sv_rooted{

		my %hash = %{$_[0]};	#FP-Indices passed to the subroutine

		my $mod_value;			#Modified Shapley Value
		my %hash_sv_rooted;		#Hash to store the modified Shapley Values for all leaves
		my %hash_pd;			#Hash to store the PD for every leaf

		#Loop over all leaves and calculate the distance to the root, i.e. the PD of the leaves
		foreach (our @leaves){
			my $PD = 0;
			my @lineage = our $tree->get_lineage_nodes($_);
			splice @lineage,0,1;	#Remove the root, since it doesn't have a branch length assigned to it

			foreach (@lineage){
				$PD = $PD + $_->branch_length();
			}
			$PD = $PD + $_->branch_length();
			$hash_pd{$_->id()} = $PD;
		}

		
		#Calculate the modified Shaples Values from the Fair Proportion Indices by substraction of PD(Leaf)/n
		while ((my $key, my $value) = each %hash){
			$mod_value = $value - $hash_pd{$key}/our $tree_size;
			$hash_sv_rooted{$key} = $mod_value;
		}
	
		return \%hash_sv_rooted;
	}

###############################################################################################################################################################
	#Subroutine for calculation of the Shapley Value of unrooted trees
	sub calc_sv_unrooted{
		my %hash_sv_unrooted;
		#my @sv_unrooted;	#Array to store the Shapley-Values of the unrooted tree
		my %far_sets;		#For each edge (node) store those leaves, which will be separated from the tree when removing this edge (= far set)
							#$far_sets{$edge} = reference to array containing descendent leaves of this edge
	
		my %branch_lengths; 	#Hash to store the branch_lengths of each edge/node as it's somehow not possible to assess them later
	
		#Loop over all nodes to fill the hash %far_sets
		foreach (our @nodes){
			if ($_->branch_length()){	#Check whether the current node has a value assigned to it
				$branch_lengths{$_} = $_->branch_length(); 	#Store the branch_length
			
				my @descendent_nodes = ();      #Store all descendent nodes of a node
				my @descendent_leaves = (); 	#Store all descendent leaves of a node
			
				#If the current node is a leaf the array @descendent_leaves will only contain this element. 
				#The set of nodes separated from the tree when removing this node/edge, will also only contain this node itself.
				if ($_->is_Leaf()){
					push(@descendent_leaves, $_);
					$far_sets{$_} = \@descendent_leaves;
				}
			
				#If the current node is not a leave, one has to go through all descendent nodes of this node and filter them for leaves.
				#Those leaves will then form the "far set".
				else{
					@descendent_nodes =  $_->get_all_Descendents();
					foreach(@descendent_nodes){
						if ($_->is_Leaf()){
							push(@descendent_leaves, $_);
						}
					}
					$far_sets{$_} = \@descendent_leaves;
				}
			}
		}

	
		#Loop over all leaves to calculate the Shapley-Values. For each leaf loop over all edges to do so.
		foreach (our @leaves){
		
			my $current_leaf = $_;
			my $far_set = 0;	
			my $containing_set = 0; 	
			my $shapley_value = 0;
		
			while (my ($key,$value) = each %far_sets){	#Key:Node, Value:Array reference
				my $branch_length = $branch_lengths{$key};

				#If the current leaf is among the descendents of the current edge $_, the size of the containing_set will simply
				#be the size of the number of descendent leaves of this edge, i.e. @{$value}. The size of the far_set is calculated by subtracting
				#the size of the containg_set from the tree size. 
				if (grep $_ eq $current_leaf, @{$value}){
					$containing_set = @{$value};
					$far_set = our $tree_size - $containing_set;
					$shapley_value = $shapley_value + $branch_length * $far_set / ($tree_size * $containing_set);
				}

				#If the current leaf is not among the descendents of the current edge $_, those descendents will form the far_set. 
				#Thus the size of the far_set i @{$value} and again the size of the containing_set calculates as $tree_size - $far_set.
				else{
					$far_set = @{$value};
					$containing_set = our $tree_size - $far_set;
					$shapley_value = $shapley_value + $branch_length * $far_set / ($tree_size * $containing_set);
				}
			}

			$hash_sv_unrooted{$_->id()} = $shapley_value;
		}

		return \%hash_sv_unrooted;
	}

	###############################################################################################################################################################

	sub print_fp{
		my %hash = %{$_[0]};
		#my @sorted_keys = sort{$hash{$b} <=> $hash{$a}} keys %hash;
		if ($terminal == 1){
			print "\nThe Fair Proportion Indices for tree ".$tree_counter." are : \n";
			foreach (keys %hash){
				printf("%-30s%-20s\n",$_,$hash{$_});
			}
		}
		else {
			open (OUTPUT, ">>$output") or die ("Could not open file $output\n");
			print OUTPUT "\nThe FP indices for tree ".$tree_counter." are : \n";
			foreach (keys %hash){
				#printf OUTPUT "%-30s%-20s\n",$_,$hash{$_};
				print OUTPUT $_."\t".$hash{$_}."\n";
			}
			close(OUTPUT);
		}
	}
	


	sub print_sv{
		my %hash = %{$_[0]};
		my @sorted_keys = sort{$hash{$b} <=> $hash{$a}} keys %hash;
		if ($terminal == 1){
			print "\nThe Shapley Values for tree ".$tree_counter." are : \n";
			foreach (@sorted_keys){
				printf("%-30s%-20s\n",$_,$hash{$_});
			}
		}
		else {
			open (OUTPUT, ">>$output") or die ("Could not open file $output\n");
			print OUTPUT "\nThe Shapley Values for tree ".$tree_counter." are : \n";
			foreach (@sorted_keys){
				#printf OUTPUT "%-30s%-20s\n",$_,$hash{$_};
				print OUTPUT $_."\t".$hash{$_}."\t";
			}
			close(OUTPUT);
		}
	}

}
catch{}

finally{
	if(@_){
		print "An error has occured: @_\nPlease, make sure that your input file is in Newick format!\n";
	}
}		

