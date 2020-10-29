# PhylogenomicsPipelines

This repository contains two pipelines to perform phylogenomics analysis. One pipeline recovers genes from sample Illumina read data and the second pipeline performs phylogenetic analysis on the recovered genes to obtain a species tree. They can run on Linux (with the Slurm job manager, if installed) and MacOS.

For gene recovery use
```
recover_genes_from_all_samples.sh
```
For phylogenetic analysis use
```
make_species_trees_pipeline.sh
```
On the command line type the name of the program for brief instructions on its use. For the phylogenetic analysis pipeline you should start in a fresh directory.
 
## Required software dependencies
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum)

For gene recovery:
* Trimmomatic, exactly version 0.39 (or alter the version number in script recover_genes_from_one_sample.sh, line ~50)
* Paftools or HybPiper
* If using HybPiper, seqtk

For phylogenetic analysis:
* Python 2.7
* [seqtk](https://github.com/lh3/seqtk)
* bc, a Linux command line utility
* [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
* [MAFFT ](https://mafft.cbrc.jp/alignment/software/) or [UPP](https://github.com/smirarab/sepp/blob/master/README.UPP.md) (UPP requires SEPP and PASTA as explained in the link)
* [FastTree](http://www.microbesonline.org/fasttree/), [raxml-ng](https://github.com/amkozlov/raxml-ng) or [IQ-TREE version 2](http://www.iqtree.org)
* [RAxML](https://github.com/stamatak/standard-RAxML) (used for building a species tree using a concatenated alignment)
* [Newick Utilities](http://cegg.unige.ch/newick_utils)
* [ASTRAL](https://github.com/smirarab/ASTRAL), exactly version 5.7.4 (or alter the version number in script make_species_trees.sh, line ~151)
* [AMAS.py](https://github.com/marekborowiec/AMAS) and/or [trimAl](http://trimal.cgenomics.org/) (for trimming if those options are selected)
* [TreeShrink](http://trimal.cgenomics.org/)
* R (used only by treeshrink and trimAl)
 
## How to use, outputs and further details  
### recover_genes_from_all_samples.sh
<!--
All options (NB - problem with presenting all options is what happens if they change - have to update here as well -OK?)
```
-s   add sample name and fastq file names via a csv table file (must have a header line);
     format: SampleName,R1FastqName,R2FastqName (required option) ]
```
File containing a list of samples and data files to recover genes from.
```
-f   FULL path to all samples(required); N.B. no filenames, just the full path to them, not a relative path and no wild cards! ]
```
Can be a full path or relative path to a single directory containing all samples to use
```
  [ -t   file name of target genes in fasta format (required option) ]
```
NEED to place an example targets file in the example folder  
```
[ -a   file name of adaptors in fasta format (required option) ]
```
You need to know what adaptors to use for your sequencing data. An example file is provided for the example data provided (Trimmomatic, palidindromic mode). STILL NEED TO ADD THIS.
```
  [ -y   Hyb-Seq program; options are 'paftools' or 'hybpiper' (default=paftools) ]
```
```
  [ -p   directory prefix for each sample (default=Sample) ]
```


```
  [ -c   number of cpu to use (default=4) ]
  [ -m   Slurm memory to use (in MB); Paftools requires >> 20000, HybPiper requires << 20000 (default=20000) ]
  [ -q   Slurm partition (queue) to use (default=fast) ]
  [ -h   prints usage and description ]
  [ -v   program version]
```
Typical example of how to run:
```
<path to>/recover_genes_from_all_samples.sh \
-s <table_file.csv> \
-t <angiosperms353TargetsFile.fasta> \
-f <fastq_files_path> \
-a <illumina_adaptors.fasta> \
-p Sample \
-c 4 \
-m 20000 \
-q main \
> recover_genes_from_all_samples.log 2>&1 &
```
Can either specify the full path to the bash script or set the path in your bash_profile file so that you can just type the name of the program on the commandline. A directory for each sample is created using the prefix specified by -p (default value is 'Sample_'). You will need to alter the Slurm options, especially the -q option to reflect the names used by your Slurm setup.

Show commandline and typical usage with explantion
e.g. throttle set to 1

-->

### make_species_trees_pipeline.sh 
### Example
A basic example is shown at the bottom of the command line help. A more extensive analysis is presented below.

Build genes trees from sample fasta files, formatted as described for option -a, by aligning the input DNA sequence for each gene with UPP (option -A),  filtering out genes with low sequence coverage across the alignment (option -F), removing very rare insertions (option -K), building each gene tree with IQTREE-2 using the Ultrafast bootstrap option (option -q), using TreeShrink to identify usually long branches in the gene trees (option -T), then collapsing nodes with bootstrap values < 10% before building a species tree with ASTRAL (option -L). Finally, FASTTREE and RAxML are then used to reconstruct a supermatrix tree build from a concatenated set of the UPP gene alignments:
```
make_species_trees_pipeline.sh \
-a \
-D 'dna' \
-A upp \
-M '--retree 2' \
-t <path_to>/taxon_info_for_tree_labels.csv \
-g <path_to>/<file_with_target_geneIds_ONLY.txt> -F '60 30' \
-K 0.003 \
-q iqtree2 \
-T \
-L 30 \
-c 8 \
-C 6 \
-Q long \
-R 15000 \
-U 12000 \
-V 0 \
-W 0 \
-H 1 \
*.fasta \
> make_species_trees_pipeline.log 2>&1 &
```
Note that in the above command there must not be a space character after the back slash, there must be a space before the back slash and any option values that contain spaces need to be quoted e.g. option -M. The Slurm options, -Q, -R, -U, -V and -W) will need to be altered depending on 


<!--
### Outputs
The main output files are listed below. Some files will only exist where relevant options have been selected    [ turn list below into a table ]
* \<geneNameId\>.dna.fasta<br>
  Gene-wise raw DNA sequence file, one gene from (all) samples in the sample set
* \<geneNameId\>.protein.fasta<br>
  Gene-wise raw protein sequence file, one gene from (all) samples in the sample set (optional)	
|	|		|-- <geneNameId>.dna.aln.for_tree.fasta					Gene DNA sequence alignment file
|	|		|-- <geneNameId>.protein.aln.for_tree.fasta				Gene protein sequence alignment file
|	|		|-- <geneNameId>.dna.gene_tree_USE_THIS.nwk				Gene tree file in Newick format				
|	|		|-- run<species_tree_runId>.summary_gene_recovery.txt		
|	|		|-- run<species_tree_runId>.summary_stats.txt
|	|		|-- run<species_tree_runId>.dna.<method>.species_tree_USE_THIS.nwk	Species tree
|	|		|-- new directoriues made



## Further info on software dependencies
Disucss it workd on different OS - developed specifically on Linux x 2 and the Darwin on macbook  HighSierra Howver I hope the pipelinces can run on other closely related OS. I coudl spend lotso tim developin and checking other OS but instead but instead it's more effficient to just curl up quietly and disappear into /dev/null. 
e.g. adding to $PATH
NB - using AMAS.py the trim option needs to exist for pipleiene to work - if it is not, you may have an older version so download the latest Github repo 

## Explanation of how pipelines work
The aim here is to describe the main steps in detail and enough to trouble shoot certain errors.
It is possible to re-run the script in the same directory but if there are files that you actually want to remove before re-running you shoudl start in a fresh direcotry. Wildcards are used to pick up the set fo fiel requried for a particualr step e.g. if you realise that you no longer need a fasta file you need to remove the modified.fasta file from the directory otherwise it will be included again. Second example: if you remove one or more genes from the analysis then existing files from the previous run will be picker up again
   -->





