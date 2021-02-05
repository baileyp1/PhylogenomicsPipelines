# PhylogenomicsPipelines

This repository contains two pipelines to perform phylogenomic analysis. One pipeline recovers genes from sample Illumina read data and the second pipeline performs phylogenetic analysis on the recovered genes to obtain a species tree. They will run on Linux (and with the [Slurm](https://slurm.schedmd.com/) job manager, if installed) and MacOS.

This software has been used in the following work to construct and analyse the [Kew Tree of Life](https://treeoflife.kew.org/) ([PAFTOL](https://www.kew.org/science/our-science/projects/plant-and-fungal-trees-of-life) project):
```
Baker et al (2021) A comprehensive phylogenomic platform for exploring the angiosperm tree of life (in preparation)
``` 

For gene recovery use
```
recover_genes_from_all_samples.sh
```
For phylogenetic analysis use
```
make_species_trees_pipeline.sh
```
On the command line type the name of the program for brief instructions on its use. For the phylogenetic analysis pipeline it is best to start in a fresh directory.
 
##  Additional software required
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum). 
<!--
For the Java programs (Trimmomatic, ASTRAL) it is best to set up an alias, named as described below -->

For gene recovery:
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), exactly version 0.39 (or alter the version number in script recover_genes_from_one_sample.sh, line ~54)
* [Paftools](https://github.com/RBGKew/pypaftol) or [HybPiper](https://github.com/mossmatters/HybPiper)
* If using HybPiper, [seqtk](https://github.com/lh3/seqtk)

For phylogenetic analysis:
* Python 2.7
* [seqtk](https://github.com/lh3/seqtk)
* bc, a Linux command line utility
* [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
* [MAFFT ](https://mafft.cbrc.jp/alignment/software/) or [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) (UPP requires SEPP and PASTA as explained in the link)
* [FastTree](http://www.microbesonline.org/fasttree/), [RAxML-NG](https://github.com/amkozlov/raxml-ng) or [IQ-TREE version 2](http://www.iqtree.org)
* [RAxML](https://github.com/stamatak/standard-RAxML) (used for building a species tree using a concatenated alignment)
* [Newick Utilities](http://cegg.unige.ch/newick_utils)
* [ASTRAL](https://github.com/smirarab/ASTRAL), exactly version 5.7.4 (or alter the version number in script make_species_trees.sh, line ~235)
* [AMAS.py](https://github.com/marekborowiec/AMAS) (check it has the 'trim' option) and/or [trimAl](http://trimal.cgenomics.org/) (for trimming if those options are selected)
* [TreeShrink](http://trimal.cgenomics.org/)
* [R](https://www.r-project.org/) v3.4.x or higher  (used by TreeShrink and trimAl)
 
## How to use, outputs and further details  
## recover_genes_from_all_samples.sh
### Example
A typical example is shown at the bottom of the command line help. 

## make_species_trees_pipeline.sh 
### Example
A basic example is shown at the bottom of the command line help. A more extensive analysis is presented below (options for this pipeline in brackets) with jobs set to run via Slurm job manager on a High Performance Computing (HPC) Linux cluster. 

Build genes trees from sample fasta files, formatted as described for option -a, by aligning the input DNA sequence for each gene with UPP (option -A),  filtering out genes with low sequence coverage across the alignment (option -F), removing very rare insertions (option -K), building each gene tree with IQTREE-2 using the Ultrafast bootstrap option (option -q), using TreeShrink to identify usually long branches in the gene trees (option -T), then collapsing nodes with bootstrap values < 30 % (option -L) before building a species tree with ASTRAL. Finally, FASTTREE and RAxML are then used to reconstruct a supermatrix tree built from a concatenated set of the UPP gene alignments:
```
make_species_trees_pipeline.sh \
-a \
-D 'dna' \
-A upp \
-t <path_to>/taxon_info_for_tree_labels.csv \
-g <path_to>/<file_with_target_geneIds_ONLY.txt> \
-F '60 0' \
-K 0.003 \
-q iqtree2-B1000-nm1000 \
-T \
-L 30 \
-c 26 \
-C 4 \
-Q long \
-Y long \
-R 5500 \
-U 12000 \
-V 0 \
-W 0 \
-H 1 \
-X 5921,5596:6:20000 \
*.fasta \
> make_species_trees_pipeline.log 2>&1 &
```
Note that in the above command there must not be a space character after the back slash, there must be a space before the back slash and any option values that contain spaces need to be quoted e.g. option -F. 

