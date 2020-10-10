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
On the command line type the name of the program for brief instructions on its use. For the phylogenetic analysis pipeline you should start with a new directory.
 
## Required software dependencies
The following programs need to be installed and available from the command line by typing the native program name. Some but not all of them are easily available from software installers (e.g. bioconda, brew, apt, yum)

For gene recovery:
* Trimmomatic, exactly version 0.39 (or alter the version number in script recover_genes_from_one_sample.sh, line ~50)
* Paftools or HybPiper
* If using HybPiper, seqtk

For phylogenetic analysis:
* python 2.7
* seqtk
* bc, a Linux command line utility
* exonerate
* MAFFT
* fasttree, raxml-ng or [IQ-TREE version 2](http://www.iqtree.org)
* RAxML
* Newick Utilities
* [ASTRAL](https://github.com/smirarab/ASTRAL), exactly version 5.7.3 (or alter the version number in script make_species_trees.sh, line ~151)
* AMAS.py and/or trimAl (for trimming if those options are selected)
* treeshrink
* R (used only by treeshrink and trimAl if )

## Further Details options, how to use and outputs  
To finish
<!--
### recover_genes_from_all_samples.sh
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



### make_species_trees_pipeline.sh  
Options
```  don't really need this one here
	-h             shows this message
```
```
	-v             program version
```
The version number is taken from the GitHub tag
CHANGE program to pipeline.
	-a             add sample name onto the fasta header from the input fasta file name.
	               Expected gene identifier format in the fasta header: >geneId (no hyphen '-' characters)
	               NB - the files with the modified fasta headers are written to files with name ending \_modified.fasta_
```
	-t <csv file>  add sample name and other info from a comma separated value (csv) table file into the tree leaf labels.
	               Format of table row: sample_name/identifier, followed by any species information as required. 
	               Sample_name/identifier must match the sample fasta file name (minus any [dot] ending suffix e.g. .fasta)
```
Table can have a header line or not as desired.

	-u             add contig length info onto species tree tips (requires option -t)
	-g <file>      file (including path to it) containing list of gene names (required option)
	-i             make gene trees only
	-j             make species trees only
	-f <float>     fraction of [well conserved regions in] the alignment covered by a sample sequence.
	               Minimum to tolerate (default=0.6; 0 would mean no filtering, i.e. include sequence of any length)
	-s <float>     fraction of samples in each gene tree.
	               Minumum to tolerate (default=0.6; 0 would mean no filtering, include all available samples in each gene tree)
	-p <string>    prefix for the output filenames e.g. taxonomic group (default=tree_pipeline)
	-q <string>    name of phylogeny program for gene trees from DNA sequences.
                   	Options are, fastest to slowest: fasttree, raxml-ng (default=fasttree)
	-r <string>    name of phylogeny program for gene trees from protein sequences.
                   	If required, options are, fastest to slowest: fasttree, raxml-ng (no default)
	-c <integer>   number of cpu to use for RAxML in supermatrix method (default=8)


A typical example:
<path to>/make_species_trees_pipeline.sh \\
-g <geneListFile> \\
-t <sampleTreeTipInfoFile> \\
-f 0.6 \\
-s 0.3 \\
-q fasttree \\
-c 8 \\
<path_to_recovered_genes_from_samples>/*.fasta \\
> make_species_trees_pipeline.log 2>&1 &


## Further advice on software dependencies
Disucss it workd on different OS - developed specifically on Linux x 2 and the Darwin on macbook  HighSierra Howver I hope the pipelinces can run on other closely related OS. I coudl spend lotso tim developin and checking other OS but instead but instead it's more effficient to just curl up quietly and disappear into /dev/null. 
e.g. adding to $PATH
NB - using AMAS.py the trim option needs to exist for pipleiene to work - if it is not, you may have an older version so download the latest Github repo 

## Explanation of how pipelines work
The aim here is to describe the main steps in detail and enough to trouble shoot certain errors.
It is possible to re-run the script in the same directory but if there are files that you actually want to remove before re-running you shoudl start in a fresh direcotry. Wildcards are used to pick up the set fo fiel requried for a particualr step e.g. if you realise that you no longer need a fasta file you need to remove the modified.fasta file from the directory otherwise it will be included again. Second example: if you remove one or more genes from the analysis then existing files from the previous run will be picker up again
   -->





