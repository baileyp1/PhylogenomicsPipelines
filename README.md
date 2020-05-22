# PhylogenomicsPipelines

This repository contains two pipelines to perform phylogenomics analysis. One pipeline recovers genes from sample Illumina read data and the second pipeline performs phylogenetic analysis on the recovered genes to obtain a species tree. They run on Linux (with the Slurm job manager, if installed) and MacOS.

For gene recovery use
```
recover_genes_from_all_samples.sh
```
For phylogenetic analysis use
```
make_species_trees_pipeline.sh
```
On the command line type the name of the program for brief instructions on its use.
 
## Required software dependencies
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum)

For gene recovery:
* Trimmomatic, exactly version 0.39
* Paftools or HybPiper
* If using HybPiper, seqtk

For phylogenetic analysis:
* seqtk
* blast
* bc, a Linux command line utility
* exonerate
* MAFFT
* fasttree or raxml-ng
* RAxML
* Newick Utilities
* Astral, exactly version 5.7.3
* AMAS.py

<!-- Next section: Further Details on how to use and outputs   
Show commandline and typical usage with explantion
e.g. throttle set to 1

Further advice on software dependencies
e.g. adding to $PATH

Explanation of how pipelines work
The aim here is to describe the main steps in detail
   -->





