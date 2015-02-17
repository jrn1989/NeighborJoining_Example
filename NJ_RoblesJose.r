# Data Mining Project - "Neighbor Joining"
# José Robles Noriega
# DMKM - Bucharest 2013 

# To check if the required package is installed
if("ape" %in% rownames(installed.packages()) == FALSE) {install.packages("ape")};


# ape library that provides functions for reading, writing, plotting, and manipulating phylogenetic trees
library("ape")

# Attention: Dont forget to set up the working directory source file location
setwd("~/Desktop/InBucharest/DataMining/RoblesJose_NeighborJoining_DataMiningProject")

# First dataset

# Data reader
dnaSequences_1 <- read.dna("fusADNAAlign.fasta.txt", format = "fasta")

# To compute pairwise distance for all the sequences
distances_1 <- dist.dna(dnaSequences_1, model = "K80")

# To compute NJ 
nj_tree_1 <- nj(distances_1)

# Plotting results
png(file="nj_1.png", units="in", width=11, height=8.5, res=300)
plot(ladderize(nj_tree_1),"p",cex=0.4,dir="d",edge.width=0.5)
dev.off()
summary.phylo(nj_tree_1)

# Bootstaping
png(file="nj_boot_1.png", units="in", width=11, height=8.5, res=300)
myboot_1 <- boot.phylo(nj_tree_1, dnaSequences_1, function(xx) nj(dist.dna(xx)), B=100,block=20)
plot.phylo(ladderize(nj_tree_1),"p",cex=0.4,dir="d",edge.width=0.5)   # plot the unrooted phylogenetic tree
nodelabels(myboot_1,cex=0.5)    # plot the bootstrap values
nj_tree_1$node.label <- myboot_1
dev.off()

# Second dataset

dnaSequences_2 <- read.dna("woodmouseDNAAlign.fasta.txt", format = "fasta")

# To compute pairwise distance for all the sequences
distances_2 <- dist.dna(dnaSequences_2, model = "K80")
# To compute NJ 
nj_tree_2 <- nj(distances_2)

# Plotting results
png(file="nj_2.png", width=11, height=8.5, units="in",res=300)
plot(ladderize(nj_tree_2),"p",edge.width=0.5)
dev.off()
summary.phylo(nj_tree_2)


# Bootstaping
png(file="nj_boot_2.png", units="in", width=11, height=8.5, res=300)
myboot_2 <- boot.phylo(nj_tree_2, dnaSequences_2, function(xx) nj(dist.dna(xx)), B=100,block=20)
plot.phylo(ladderize(nj_tree_2),"p",edge.width=0.5)   # plot the unrooted phylogenetic tree
nodelabels(myboot_2)    # plot the bootstrap values
nj_tree_2$node.label <- myboot_2
dev.off()
