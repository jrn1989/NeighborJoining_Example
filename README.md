# README 
by: José Robles - December 2013

In this small project we introduce the Neighbor Joining (NJ) Algorithm, along with a couple of examples in R.

## Files:

* fusADNAAlign.fasta: first data set - see **Datasets** section for more details

* woodmouseDNAAlign.fasta: second data set - see **Datasets** section for more details

* nj_1.png: image of the phylogenetic tree for the 1st dataset

* nj_2.png: image of the phylogenetic tree for the 2nd dataset

* nj_boot_1.png: image of the phylogenetic tree with bootstrap analysis for the 1st dataset

* nj_boot_2.png: image of the phylogenetic tree with bootstrap analysis for the 2nd dataset

* NJ_RoblesJose.r: r script to calculate the phylogenetic trees with Neighbor Joining algorithm and bootstrap analysis for the given datasets.


## Context:

In biology, one of the subjects of study is phylogenetics. More precisely, this is the study of all the evolutionary relationships that exist among all the living organisms. In the past, scientists used to rely on comparing physical characteristics such as color, size, number of legs etc. in order to build those relationships. Nowadays, sequences of ADN and proteins are used to make this comparisons, always under the hypothesis that similar species are genetically closer. The result of phylogenetic studies is a hypothesis about the evolutionary history of taxonomic groups and one of the ways to represent it is through a phylogenetic tree:

* Leaf nodes represent the existing species

* Internal nodes represent the ancestors

*  The root represent the oldest ancestor from the point of view of evolution
The main practical application of the Neighbor Joining algorithm is precisely to build those trees.

## NJ:

NJ is a bottom-up clustering method created by Naruya Saitou and Masatoshi Nei in 1987 [1]. The input data is a distance matrix, that is being modified during each iteration. The separation between each pair of nodes is changed according the average divergence from all other nodes and the tree is built by adding a link to the pair of nodes that are least distant in the modified matrix. When 2 nodes are linked, their common ancestor is added to the tree and the terminal nodes with their respective branches are removed from the tree. This pruning process converts the newly added common ancestor into a terminal node on a tree of reduced size. So at each iteration, two terminal nodes are replaced by a new node. We calculate a new distance matrix by taking the distances from the new node and the other terminal nodes and the process is repeated. The algorithm stops when two nodes remain, separated by a single branch.


## Known Results and Issues:

The main advantage of NJ is that it is relatively fast, with a polynomial complexity O(n3) making it suitable for large datasets and bootstraping, being widely used by molecular evolutionist. It allows the rapid inclusion of all homologous sequences of the dataset in a single tree [7]. A good example of this is the Ribosomal Database Project[5] that provides to the research community annotated aRNA sequences and some tools to analyse them, being NJ one of the algorithms to build phylogenetic trees. Also, com- pared to other algorithms, NJ has the advantage that it does not assume that lineages evolve at the same time

To build a distance matrix it is necessary to use some alignment technique in the sequences and this can cause some lost of information. Another drawback is that it can assign negative lengths to some branches.

## Datasets:

The first dataset identifies the bacterial lineage corresponding to the gene fusA. This data, with 149 ADN sequences, comes from the paper authored by Scott Santos and Howard Ochman (University of Arizona, Tucson), available from the Ribosomal Database Project (additional genes and sequences developed since publication can be found there [4].) 

Example of DNA:

```
>Acinetobacter_ADP1_fusA
------------------------------------------------------------ATGGCCCGCCAAACCCCAAT
TAGTAATTACCGTAACATCGGTATTTCTGCTCACATCGATGCAGGTAAAACAACTACAACAGAACGTATTTTGTTCTACA
CAGGTGTATCTCACAAAATTGGTGAAGTACACGAT---------GGTGCAGCAACAATGGACTGGATGGAACAAGAGCAA
GAACGCGGTATTACAATTACCTCGGCTGCTACAACTTGTTTCTGGTCTGGTATGGGTAACCAATTTGCA-----------
```

The second dataset is a set of 15 sequences of the mitochondrial gene cytochrome b of the woodmouse (Apodemus sylvaticus) which is a subset of the data analysed by [2]. An example of ADN sequence for this dataset looks like this:

```
>No305
nttcgaaaaa cacacccact actaaaantt atcagtcact ccttcatcga cttaccagct
ccatctaaca tttcatcatg atgaaacttc ggctcattac taggagtctg cctaataatc
caaatcctca caggcttatt cctagcaata cactacacat cagacacaat aacagcattc
tcttcagtaa cccatatttg tcgagacgta aattatggct gactaattcg atacatacat
```

## Results and Evaluation:

Applying the NJ algorithm we obtain the phylogenetic tree shown in figure 1, along with the following statistics:

```
Phylogenetic tree: nj_tree
  Number of tips: 149
  Number of nodes: 147
  Branch lengths:
    mean: 0.06383912
    variance: 0.005324346
    distribution summary:
      Min.    1st Qu.     Median    3rd Qu.       Max.
-6.641e-17  5.546e-03  3.767e-02  1.070e-01  3.841e-01
  No root edge.
  First ten tip labels: Acinetobacter_ADP1_fusA
                        Agro_tumefaciens_C58_fusA
                        Aquifex_aeolicus_fusA
                        Bac_anthracis_0581_fusA
                        Bac_anthracis_AmesfusA
                        Bac_anthracis_Sterne_fusA
                        Bac_cereus_10987_fusA
                        Bac_cereus_14579_fusA
                        Bac_halodurans_fusA
                        Bac_subtilis_fusA
  No node labels.
```  


![](https://github.com/jrn1989/NeighborJoining_Example/blob/master/nj_1.png)

To evaluate the resulting tree we can perform a resampling method such as bootstrap, see next figure. Here, the numbers in blue boxes are bootstrap values for the nodes in the tree.

![](https://github.com/jrn1989/NeighborJoining_Example/blob/master/nj_boot_1.png)


This value is and indicator of how much confidence we have in the group defined by that node in the tree. With a high value, near to 100%, we are very confident that the group defined by that node is correct; on the other hand, with a very low value we are not so confident.
In this case, the bootstrap values were calculated with 100 random resamples, each one of them with 20 randomly sampled columns from the alignment. Each resample forms a sort of fake alignment of its own and a phylogenetic can be based upon this. These 100 trees are known as the ”bootstrap trees”. Then, for each group that we see in the original phylogenetic tree, we count in how many of the 100 bootstrap trees it appears: this is the bootstrap value for the group in the original tree.

And for the second dataset, we obtain:

```  
Phylogenetic tree: nj_tree_2
  Number of tips: 15
  Number of nodes: 13
  Branch lengths:
    mean: 0.002430697
    variance: 5.135314e-06
    distribution summary:
      Min.    1st Qu.     Median    3rd Qu.       Max.
-3.161e-05  7.977e-04  1.715e-03  3.087e-03  9.407e-03
  No root edge.
  First ten tip labels: No305
                        No304
                        No306
                        No0906S
                        No0908S
                        No0909S
                        No0910S
                        No0912S
                        No0913S
                        No1103S
```  

![](https://github.com/jrn1989/NeighborJoining_Example/blob/master/nj_2.png)

![](https://github.com/jrn1989/NeighborJoining_Example/blob/master/nj_boot_2.png)

## Conclusions:

Neighbor Joining is a bottom-up clustering mainly used in the reconstruction of phylogenetic trees. It is relatively simple and fast which makes it suitable for bootstrapping and for analyzing large datasets. Unfortunately, it relies on the intermediate step of building distance measures which is not always very ac- curate. Even though there are methods that can build a phylogenetic the from the ADN sequences based on maximum likelyhood (without calculating distance measures) it stills being used and improved[8].


## References

[1] Saitou, Naruya, and Masatoshi Nei. ”The neighbor-joining method: a new method for reconstructing phylogenetic trees.” Molecular biology and evolution 4.4 (1987): 406-425.

[2] Michaux, J. R., Magnanou, E., Paradis, E., Nieberding, C. and Libois, R. (2003) Mitochondrial phylogeography of the Woodmouse (Apodemus sylvaticus) in the Western Palearctic region. Molecular Ecology, 12, 685–697.

[3] Santos, Scott R., and Howard Ochman. ”Identification and phylogenetic sorting of bacterial lineages with universally conserved genes and proteins.” Environmental Microbiology 6.7 (2004): 754-759.

[4] http://rdp.cme.msu.edu/misc/resources.jsp

[5] http://rdp.cme.msu.edu/index.jsp

[6] http://a-little-book-of-r-for-bioinformatics.readthedocs.org/en/latest/

[7] http://www.icp.ucl.ac.be/~opperd/private/neighbor.html

[8] Pearson, William R., Gabriel Robins, and Tongtong Zhang. ”Generalized neighbor-joining: more reliable phylogenetic tree reconstruction.” Molecular biology and evolution 16.6 (1999): 806-816.
