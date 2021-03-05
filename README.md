
#  PhylogenomicsPipelines

This repository contains two pipelines to perform phylogenomic analysis. One pipeline recovers genes from sample Illumina read data and the second pipeline performs phylogenetic analysis on the recovered genes to obtain a species tree. They will run on Linux (and with the [Slurm](https://slurm.schedmd.com/) job manager, if installed) and MacOS.

This software has been used in the following work to construct and analyse the [Kew Tree of Life](https://treeoflife.kew.org/) ([PAFTOL](https://www.kew.org/science/our-science/projects/plant-and-fungal-trees-of-life) project):

Baker et al (2021) A comprehensive phylogenomic platform for exploring the angiosperm tree of life (submitted to [bioRxiv](https://doi.org/10.1101/2021.02.22.431589) and Systematic Biology)

For gene recovery use
```
recover_genes_from_all_samples.sh
```
For phylogenetic analysis use
```
make_species_trees_pipeline.sh
```


To run the scripts correctly, they need to be made available from any directory on the command line. Do this by adding the path of the PhylogenomicsPipelines directory (wherever you want that to be) to the .bash_profile file or equivalent file in your home directory - e.g.:
```bash
PATH=$PATH:<full_path_to_your_home_directory>/ProgramFiles/PhylogenomicsPipelines
```
Log back in to pick up the new addition to the path, then on the command line type the name of the programs shown above to see brief instructions on their use.
 
##  Additional software required
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum). For the Java programs (e.g. Trimmomatic, ASTRAL) it is necessary to set up a global variable (add to your .bash_profile or equivalent file) as follows:
```bash
export TRIMMOMATIC=<path_to_executable>/Trimmomatic-0.39/trimmomatic-0.39.jar
export ASTRAL=<path_to_executable>/astral/Astral/astral.5.7.4.jar
```

For gene recovery:
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Paftools](https://github.com/RBGKew/pypaftol) or [HybPiper](https://github.com/mossmatters/HybPiper)
* If using HybPiper, [seqtk](https://github.com/lh3/seqtk)
<!--
For generating basic statistics following gene recovery e.g. read depth, percent of reads on-target:
-->  

For phylogenetic analysis:
* Python 2.7
* [seqtk](https://github.com/lh3/seqtk)
* bc, a Linux command line utility
* [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
* [MAFFT ](https://mafft.cbrc.jp/alignment/software/) or [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) (UPP requires SEPP and PASTA as explained in the link)
* [PAL2NAL](http://www.bork.embl.de/pal2nal/) (if option '-D codon' is used)
* [FastTree](http://www.microbesonline.org/fasttree/), [RAxML-NG](https://github.com/amkozlov/raxml-ng) or [IQ-TREE version 2](http://www.iqtree.org)
* [RAxML](https://github.com/stamatak/standard-RAxML) (used for building a species tree using a concatenated alignment)
* [Newick Utilities](http://cegg.unige.ch/newick_utils)
* [ASTRAL](https://github.com/smirarab/ASTRAL)
* [AMAS.py](https://github.com/marekborowiec/AMAS) (ensure the version has the 'trim' option) and/or [trimAl](http://trimal.cgenomics.org/) (for trimming if those options are used)
* [TreeShrink](https://github.com/uym2/TreeShrink) (if option is used)
* [R](https://www.r-project.org/) v3.4.x or higher  (if TreeShrink or trimAl are being used)
 
## How to use, outputs and further details  
## recover_genes_from_all_samples.sh
### Example
A typical example is shown at the bottom of the command line help. 

## make_species_trees_pipeline.sh 
### Example
A basic example is shown at the bottom of the command line help. A more extensive analysis is presented below (options for this pipeline in brackets) with jobs set to run via Slurm job manager on a High Performance Computing (HPC) Linux cluster. 

Build genes trees from sample fasta files, formatted as described for option -a, by aligning the input DNA sequence for each gene with UPP (option -A),  filtering out genes with low sequence coverage across the alignment (option -F), removing very rare insertions (option -K), building each gene tree with IQTREE-2 using the Ultrafast bootstrap option (option -q), using TreeShrink to identify usually long branches in the gene trees (option -T), then collapsing nodes with bootstrap values < 30 % (option -L) before building a species tree with ASTRAL. Finally, FastTree and RAxML are then used to reconstruct a supermatrix tree built from a concatenated set of the UPP gene alignments:
```bash
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

The Slurm options, -Q, -Y, -V and -W, will need to be altered depending on how Slurm is set up. Slurm memory options, -R and -U, will need to be altered depending on the size of the data set; also refer to option -X which helps with memory efficiency by  allowing the user to increase memory for specific genes that require significantly more memory. The values set in the above example were appropriate for building gene trees then a species tree with ASTRAL containing up to 3,200 samples (Refer to the citation above). 

### Outputs 
The main output files are listed below in output order. Additional files will exist depending on the option(s) used, the option(s) value and the underlying programs triggered by the option(s).

* \<geneName\>.[dna|protein|codon].fasta<br>
  Unaligned sequence file, one gene from (all) samples in the data set<br>
  
* \<geneNameId\>.[dna|protein|codon].aln.for_tree.fasta<br>
  Aligned sequence file, after any filtering and/or trimming. Used to build the species tree

* \<geneNameId\>.[dna|protein|codon].gene_tree_USE_THIS.nwk<br>
  Gene tree file in Newick format
  
* make_gene_trees[-slurm_array_id].log<br>
  Log file containing any program outputs and information 

* Summary files

  * \<prefix>_summary_of_sample_quality.txt
  * \<prefix\>.summary_gene_recovery.txt
  * \<prefix>.summary_stats.txt

* When sequence filtering or TreeShrink options are used, a new directory is created and the retained sequences from those steps are re-aligned automatically. A new directory is created that contains the outputs already described above called:<br>
  after_[reAlnFilterSeqs|treeshrink]\_USE\_THIS_[dna|protein|codon]
 <!-- 
      * run3b_dna_gene_trees_for_coelescence_phylo_bs_less_30_rmed.nwk - present the new file now   -->

* Species tree Newick files
  * \<prefix>.[dna|protein|codon].species\_tree.astral_-t2.nwk<br>
  Tree with all support values reported at each node
  * \<prefix>.[dna|protein|codon].species\_tree.astral_pp1_value.nwk<br>
 Tree with the local posterior probabilities only reported at the nodes
  * \<prefix>.[dna|protein|codon].\<method>.species_tree_USE_THIS.nwk<br>
  Tree with annotated of the tree tip labels if option -t used<br>
  * \<prefix>.fasttree_[dna|protein|codon]_species_tree_USE_THIS.nwk<br>
 Tree using concatenated gene alignments with FastTree
  * \<prefix>.\_RAxML_100\_bs_[dna|protein|codon]_species_tree_USE_THIS.nwk<br>
 Tree using concatenated gene alignments with RAxML
 
* Summary files<br>
  * run3b.amas_summary.txt<br>
    summary of gene trees reported by AMAS.py
  * run3b_summary_tree_stats.txt<br>
    summary tree statistics
---
<span style="font-size:small;">Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew</span>




