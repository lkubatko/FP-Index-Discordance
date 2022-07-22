#RootTrees.jl
#Script to root unrooted trees by specifying an outgroup using the Julia package PhyloNetworks: 
#https://crsl4.github.io/PhyloNetworks.jl/latest/

using PhyloNetworks


#Example usage for the mammal data set: 
#Root 447 gene trees stored in the file "447GeneTreesPhyMLGTR.txt" using taxon "Gal" as the outgroup
#(Original source for gene trees: Supplementary text S2 to 
#"Mark S. Springer, John Gatesy: The gene tree delusion, Molecular Phylogenetics and Evolution, 94(A), 1-33, 2016").


#Read in multiple trees from a file
trees = readMultiTopology("447GeneTreesPhyMLGTR.txt")

#Loop over the gene trees stored in "trees" and root them using "Gal" as the outgroup; save rooted trees in file "447GeneTreesPhyMLGTR_rooted.txt"
for i in 1:447
	rootedtree = rootatnode!(trees[i],"Gal")
	writeTopology(rootedtree, "447GeneTreesPhyMLGTR_rooted.txt"; append=true)
end




