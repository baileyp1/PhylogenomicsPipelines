# PhylogenomicsPipelines

This repository contains two pipelines to perform phylogenomic analysis. One pipeline recovers genes from sample Illumina read data and the second pipeline performs phylogenetic analysis on the recovered genes to obtain a species tree. They will run on Linux (and with the [Slurm](https://slurm.schedmd.com/) job manager, if installed) and MacOS.

This software has been used in the following work to construct and analyse the [Kew Tree of Life](https://treeoflife.kew.org/) ([PAFTOL](https://www.kew.org/science/our-science/projects/plant-and-fungal-trees-of-life) project):

Baker et al (2022) A comprehensive phylogenomic platform for exploring the angiosperm tree of life (submitted to [bioRxiv](https://doi.org/10.1101/2021.02.22.431589) and now published in [Systematic Biology](https://doi.org/10.1093/sysbio/syab035))

For gene recovery use \&nbsp &nbsp&nbsp&nbsp
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
Log back in to pick up the new addition to the path or type
```bash
source ~/.bash_profile
```
then on the command line type the name of the programs shown above to see brief instructions on their use. At least for the phylogenetic analysis script it is best to work in a new empty directory. The main reason (amongst a few others) is that the  script uses wild cards in several places to pick up a set of files and may pick up extra files that also match, however most file name sets are now checked for and deleted at the start of the analysis (except in the case of option -b).
 
##  Additional software required
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum). If the Java programs are being used (i.e. Trimmomatic, ASTRAL[-MP], Jalview, Picard) it is necessary to set global variables with exactly these names: TRIMMOMATIC, ASTRAL, ASTRALMP, JALVIEW, PICARD. Add them to your .bash_profile or (equivalent file) as follows - e.g.
```bash
export TRIMMOMATIC=<path_to_executable>/Trimmomatic-<VERSION>/trimmomatic-<VERSION>.jar
export ASTRAL=<path_to_executable>/astral/Astral/astral.<VERSION>.jar
```
Note that for ASTRAL-MP, an extra Java flag is also required to be set:
```bash
export ASTRALMPLIB=-Djava.library.path=<path_to_executable>/ASTRAL-MP/Astral/lib
```

For gene recovery (if known, specific version requirements are shown in brackets):
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Paftools](https://github.com/RBGKew/pypaftol) or [HybPiper](https://github.com/mossmatters/HybPiper) (version 1.3.1 plus patch for issue [41](https://github.com/mossmatters/HybPiper/issues/41))
* If using HybPiper, Perl
* If using HybPiper, [seqtk](https://github.com/lh3/seqtk) (version 1.3)
* If using option -S:
  * [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) (version 2.4)
  * [BWA](http://bio-bwa.sourceforge.net)
  * [Samtools](http://www.htslib.org) (version 0.7.17 or must contain 'coverage' and 'depth' programs)
  * [Picard](https://broadinstitute.github.io/picard)

For phylogenetic analysis (if known, specific version requirements are shown in brackets):
* Python
* [seqtk](https://github.com/lh3/seqtk) (version 1.3)
* [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) (version 2.4)
* [MAFFT ](https://mafft.cbrc.jp/alignment/software/) or [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md). UPP requires SEPP and PASTA as explained in the link
* [Jalview](https://www.jalview.org),  to view gene alignments (optional)
* [PAL2NAL](http://www.bork.embl.de/pal2nal/) if option '-D codon' is used
* [FastTree](http://www.microbesonline.org/fasttree/), [RAxML-NG](https://github.com/amkozlov/raxml-ng) or [IQ-TREE version 2](http://www.iqtree.org)
* [RAxML](https://github.com/stamatak/standard-RAxML) (used for building a species tree using a concatenated alignment)
* [Newick Utilities](http://cegg.unige.ch/newick_utils)
* [ASTRAL](https://github.com/smirarab/ASTRAL). For [ASTRAL-MP](https://github.com/smirarab/ASTRAL/tree/MP), clone and switch to the ASTRAL-MP branch.
* [AMAS.py](https://github.com/marekborowiec/AMAS) (ensure the version has the 'trim' option) and/or [trimAl](http://trimal.cgenomics.org/) (for trimming if those options are used)
* [TreeShrink](https://github.com/uym2/TreeShrink) if option is used (version 1.3.5+)
* [R](https://www.r-project.org/) (version 4.0)  if the TreeShrink option <!--or trimAl are--> is used

One way of setting up the above software requirements is described [here](https://github.com/baileyp1/PhylogenomicsPipelines/blob/master/Species_Tree_Readme.md). 
 
## How to run example:  inputs, options, outputs and log files
## recover_genes_from_all_samples.sh
### Example of running recover_genes_from_all_samples.sh
A typical example is shown at the bottom of the command line help.<br>
### Outputs from recover_genes_from_all_samples.sh
The main output file from the gene recovery pipeline is an unaligned sequence file with all recovered genes from one sample within a directory dedicated to each sample:
* \<SampleName\>.fasta (Paftools)<br>
  * FASTA header format:<br> 
    \>geneNameId (option -a)
* \<SampleName\>_all_genes.fasta (HybPiper)
  * FASTA header format:<br> 
    \>sampleNameId-geneNameId (default)<br>

Gene recovery stats output (option -S):
* \<SampleName\>_gene_recovery_stats.txt<br>

The overall_gene_recovery_stats.sh script creates a table of results from the above file for multiple samples. 
  
## make_species_trees_pipeline.sh 
### Example of running make_species_trees_pipeline.sh 
A basic example is shown at the bottom of the command line help and can be used as a template to add more command options. A more extensive analysis is presented below (options for this pipeline in brackets) with jobs set to run via Slurm job manager on a High Performance Computing (HPC) Linux cluster.

Build genes trees from sample fasta files, formatted as described for option -a, by aligning the input DNA sequence for each gene with UPP (option -A),  filtering out genes with low sequence coverage across the alignment (option -F), removing very rare insertions (option -K), building each gene tree with FastTree (option -q), using TreeShrink to identify unusually long branches in the gene trees (option -T), then collapsing nodes with bootstrap values < 15 % (option -L) before building species trees with the programs specified in option -s, first ASTRAL, then FastTree is used to reconstruct a supermatrix tree built from a concatenated set of the UPP gene alignments:
```bash
make_species_trees_pipeline.sh \
-b
-a \
-D 'dna' \
-A upp \
-M '-M 80 -T 0.15' \
-t <path_to>/taxon_info_for_tree_labels.csv \
-g <path_to>/<file_with_target_geneIds_ONLY.txt> \
-F '60 0' \
-O 0 \
-K 0.003 \
-q fasttree \
-T \
-L 15 \
-s 'astralmp fasttree' \
-o mid-point \
-C 3 \
-c 32 \
-Q long \
-Y himem \
-R 13000 \
-U 1500000 \
-V 0 \
-W 0 \
-p release_3.0_tree \
-X 5921,5596,5333,5614:12:70000 \
*.fasta \
> make_species_trees_pipeline.log 2>&1 &
```
Note that in the above command there must not be a space character after the back slash, there must be a space before the back slash and any option values that contain spaces need to be quoted e.g. option -F. 

The Slurm options, -Q, -Y, -V and -W, will need to be altered depending on how Slurm is set up. Slurm memory options, -R and -U, will need to be altered depending on the size of the data set; also refer to option -X which helps with memory efficiency by allowing the user to increase memory for specific genes that require significantly more memory. Note that for one very large gene alignment (Angiosperms353 gene 5921), the memory for PASTA (used by UPP) had to be increased from the default by setting the --max-mem-mb parameter to 10GB (--max-mem-mb=10000) under the '[pasta]' section in the ~/.sepp/upp.config file.

The values set in the above example were appropriate for building gene trees then a species tree with ASTRAL-MP containing up to 10,709 samples over a 6 day period. Refer to the citation above for an overview of the methods used and use the KewTreeOfLife link above to see version 3.0 of the tree display and any method changes that apply to this version of the tree.

### Outputs from make_species_trees_pipeline.sh  
The main output files from the phylogenetic analysis pipeline are listed below in output order. Additional files will exist depending on the option(s) used, the option(s) value and the underlying programs triggered by the option(s).

* \<geneNameId\>.[dna|protein|codon].fasta<br>
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
## Further details
### Installing, compiling and running the additional software
If using an existing computer set up, many of the additional programs required might already be installed. If so, just ensure that the program version is up to date. The requirement for a specific version is listed in the software section above, if known. Minor programs listed as optional and not installed are ignored and their outputs will not exist e.g. Jalview for gene alignment images.

### How to run further:  inputs and options
### make_species_trees_pipeline.sh
### Inputs
### Preparing input fasta files
It might be necessary to select  only a subset of input files in the same directory. This is possible by replacing this line in the command example above:
```
*.fasta \
```
and instead using a sample names list file, one line per sample, with the following code:
```
`cat <sample_list_file> | awk '{print "<your_path_to_fasta_files>/" $1 ".fasta"}' `
```
This command creates a file path for each sample name in the list. ‘cat' prints the list and awk creates the path to each file. The two back ticks fire off the command resulting in a space-separated file list which in turn is accepted by the main program.<br>
Note: for this command to work the sample identifiers must match at least part of the filename and the awk string might need to be adapted further.

To check that the file paths actually exist, try to list them like so:
```
cat <sample_list_file> | awk '{print "<path_to_files>/" $1 <your_path_to_fasta_files>/" $1 ".fasta"}' | xargs ls -l
```
### Options
**Option -q**  name of phylogeny program for gene trees from DNA sequences.

  Values:<br>
* iqtree2<br>
default IQTREE2 run with model testing (slow for big alignments) and Ufboot bootstrapping<br>
* iqtree2-B1000-nm1000<br>
  single model of evolution (GTR for DNA, JTT for protein) and Ufboot bootstrapping with 1000 tree iterations<br>
* iqtree2-B1000-nm110<br>
  same as for 'iqtree2-B1000-nm1000’ except number of iterations is set to 110. Approximately 10x faster but less thorough tree search
* iqtree2-B1000-nm200 (originally set up for testing run time)    
  same as for 'iqtree2-B1000-nm1000’ except number of iterations is set to 200. Approximately 5x faster but less thorough tree search (originally set up for testing run time)

**Option -s** name of phylogeny program(s) to use for the species tree(s)
* When 'astral' is used, memory is fixed to 12GB and is sufficient to build a species tree with approximately 3,200 samples. For larger data sets it would be better to use ASTRAL-MP for faster analysis but it will need  more memory. Memory for ASTRAL-MP should be set via option -U.

**Option -H** Slurm array throttle for gene trees<br>
* By default, this option is set to 1 only (i.e. runs one gene) which is useful to make sure the set up is running as expected for one gene. However, it means that the Slurm throttle has to be altered by hand to run more genes in parallel e.g. for 50 genes:
```bash
scontrol update arraytaskthrottle=50 job=<slurm_job_id>
```
* If any type of sequence filtering is on (options -F, -I, -T), then upon gene re-alignment, the throttle is set to 50 by default.

* If option -X is used, an additional slurm array is set up for gene trees that need significantly higher memory and/or cpu for the alignment or tree steps relative to other genes in the data set. There are some! The throttle will also need to be increased for this array too. It is not necessary to remove these genes from the main list.

---
<span style="font-size:small;">Copyright © 2023 The Board of Trustees of the Royal Botanic Gardens, Kew</span>



