#!/bin/bash

################################
# make_species_trees_pipeline.sh

# Author:   Paul Bailey
################################
shopt -s failglob


####################################################################
# Set up commandline parameters and fetch their values using getopts
####################################################################

# Variables for any command line flags needing a default value:
# INPUT FILE OPTIONS:
useGenewiseFiles=no
addSampleName=no
sampleTableFile=no			# 10.5.2020 - only just added in - check - ensures cmdline parameter is always accupied which is critical
option_u=no

# ALIGNMENT OPTIONS:
seqType=dna
alnProgram=mafft 				# 27.7.2020 - may want to merge this option with mafftAlgorithm so it woudl be quoted liek so: 'mafft --retree 2'
mafftAlgorithm='--retree 2'     # Hidden option - '--maxiterate 1000' #'--retree 1'   '--retree 2' - need to ry and merge with the alnProgram option somehow

# FILTERING AND TRIMMING OPTIONS:
filterSeqs1=no			# My filter sequencing option (option 1)
filterSeqs2=no
fractnAlnCovrg=0.6			# NB - variable is used via filterSeqs1, not used if filterSeqs1 is not used
fractnSamples=0.3			# NB - variable is used via filterSeqs1, not used if filterSeqs1 is not used
trimAln1=no					# Filter alignment columns with optrimAl
trimAln2=no					# Filter alignment columns to remove rarer insertions; maximum limit of percent occupancy to trim at


maxColOccThreshold=90		### 7.9.2020 - now testing lower values e.g. 30 and 15 - still need to check input value; also shoudl it be minColOccLenThreshold
							### Hidden variable - length unit in bases

# PHYLOGENY OPTIONS:
geneTreesOnly=no
speciesTreesOnly=no
phyloProgramDNA=fasttree
phyloProgramPROT=no			# Work around to specify any program so software testing code will not crash! Ensures cmd parameter is always occupied which is critical
treeshrink=no
collapseNodes=no			# option -L

# OTHER OPTIONS:
# Hidden options (i.e. not apparent from the help menu but they always have a value so can be and have to be used in downstream scripts).
# Need to check them if I make them public
fileNamePrefix=tree_pipeline
fractnMaxColOcc=0.7
slurmThrottle=50			# Was 1, set to 50 as now have two rounds of alignment 
cpuGeneTree=1				# still keep this separate from the supermatrix tree; then I can specify loads more for the supermatrix.
							# Maybe hava a comma separated list e.g. 4,30 for gene and species tree.
cpu=8						# number of cpu to use for RAxML in supermatrix method
partitionName=main   		# Values depend on the cluster being used so good to have a flagged option for this
							# NB - make_species_tree.sh uses 'long' queue - need an extra variable for that 
geneTreeSlurmMem=20000		# option -R		 
speciesTreeSlurmMem=50000	# option -U; I think I need 500000 GB mem for the RAxML large tree (Slurm)
geneTreeSlurmTime=1-00:00	# SBATCH -t 0-36:00;  A time limit of zero requests that no time limit be imposed - like the mem option
speciesTreeSlurmTime=0		# SBATCH -t 0-36:00
####-D  upload details to PAFTOL database (internal use only) 	


function usage	{

cat << EOF

Program description: makes species trees from fasta files, either of two types:
                     1. recovered genes (DNA) from multiple samples, one fasta file per sample (default)
                        By Default, the fasta header format MUST BE: >sampleId-geneId but one other option (-a) is available.
                     2. gene-wise (DNA) fasta files, one fasta file per gene (option -G) 

Usage:               make_species_trees_pipeline.sh [OPTIONS] fastafile1 fastafile2 fastafile3 ...

OPTIONS <value>:
    -h               shows this message
    -v               program version
INPUT FILE OPTIONS:
    -G               make gene trees, starting from gene-wise fasta files rather than files containing all genes per sample
                     Gene name/identifier must be identical to the sample fasta file name (minus any [dot] ending suffix e.g. .fasta),
                     else change the gene name list in option -g so that it is.	
                     Fasta header line format MUST BE: >sampleId
    -g <file>        file (including path to it) containing list of gene names only (required option)        
    -a               add sample name onto the fasta header from the input fasta file name.
                     Expected gene identifier format in the input fasta header: >geneId (no hyphen '-' characters allowed)
    -t <csv file>    add sample name and other info from a comma separated value (csv) table file into the tree leaf labels.
                     Format of table row: sample_name/identifier, followed by any species information (include sample_name/identifier again if required) 
    -u               add contig length info onto species tree tips (requires option -t)
                     Sample_name/identifier must be identical to the sample fasta file name (minus any [dot] ending suffix e.g. .fasta)
                     Not available in gene-wise mode (option -G)
    -p <string>      prefix for the output filenames e.g. taxonomic group (default=tree_pipeline)
ALIGNMENT OPTIONS:
    -D <string>      sequence type to use: dna, protein, codon (default=dna). N.B. use with multiple types must be quoted (e.g. 'dna protein')
                     codon is input DNA aligned but guided by a protein alignment
    -A <string>      alignment program to use: mafft, UPP 	[ but for which residue type - make for all for now ]
    -M <string>      if using mafft, specify algorithm to use in quotes i.e. '--retree 1', '--retree 2', '--maxiterate 1000' etc (default='--retree 2')
FILTERING AND TRIMMING OPTIONS:
    -F <2 integers>  filter sequences option 1. Format: '<1> <2>' (no default; N.B. values must be quoted)
                     1. filters sequences in gene alignments by coverage, in this option by the percent of well conserved regions (columns) 
                        in the alignment covered by a sample sequence. A well conserved column is one with > 70 % residue occupancy
                        Minimum percent to tolerate (advised=60; 0 would mean no filtering, i.e. include sequence of any length)
                     2. percent of samples in each gene tree.
                        Minumum percent to tolerate (advised=30; 0 would mean no filtering, include all available samples in each gene tree)
                     Don't use with option -I else option -F only will applied.
    -I <integer>     filter sequences option 2 (no default; advised=60). Filters sequences from the alignment based on their length as a 
                     percent of the median length of all sequences in the data set. Don't use with option -F else option -F only will applied.
    -J               trim alignment option 1. Filter columns with OpTrimAL (Shee el al 2020, https://doi.org/10.3389/fpls.2020.00258)
    -K <fraction>    trim alignment option 2. Filter columns to remove rarer insertions (no default; example: 0.003 will remove columns with 0.3 % occupancy)
PHYLOGENY OPTIONS:
    -i               make gene trees only
    -j               make species trees only. Gene trees must already exist in run directory
    -q <string>      name of phylogeny program for gene trees from DNA sequences.
                   	 Options are, fastest to slowest: fasttree, iqtree2, raxml-ng (default=fasttree)	
    -r <string>      name of phylogeny program for gene trees from protein sequences.
                   	 If required, options are, fastest to slowest: fasttree, iqtree, raxml-ng (no default)
    -S <string>      name of phylogeny program for supermatrix approach (concatenated gene alignments).
                     If required, options are, fastest to slowest: fasttree, raxml-ng (no default)
    -T               use Treeshrink on gene trees (followed by re-alignment)
    -L <integer>     collapse gene tree nodes with bootstrap support less than <integer> percent (no default)
OTHER OPTIONS: 
    -C <integer>     number of cpu to use for genetrees; NB - not convinced >1 cpu works robustly for raxml-ng with small datasets! (default=1)              	-
    -c <integer>     number of cpu to use for RAxML in supermatrix method (default=8)
    -R <integer>     Slurm memory to use (in MB) for gene trees (default=0; means no limit is imposed)
    -U <integer>     Slurm memory to use (in MB) for species trees (default=0; means no limit is imposed)
	-V <string>      Slurm time limit to use for gene trees. Format: <days>-<hours>:<minutes> e.g. 1-0:0 is 1 day (default=0; means no limit is imposed)
	-W <string>      Slurm time limit to use for species trees. Format: <days>-<hours>:<minutes> e.g. 1-0:0 is 1 day (default=0; means no limit is imposed)
    -Q <string>      Slurm partition (queue) to use (default=medium) ]
    -H <integer>	 Slurm throttle - not available to change here (default=50)


A basic example run is described below:
make genes trees from sample fasta files, formatted as described for option -a, by aligning the input
DNA sequence for each gene with the MAFFT --retree 2 algorithm , building each gene tree with FASTTREE, 
then making species trees with ASTRAL, FASTTREE and RAxML (last two program make a supermatrix tree from 
concatenated gene alignments):

make_species_trees_pipeline.sh \\
-a \\
-D 'dna' \\
-A mafft \\
-M '--retree 2' \\
-t <sampleTreeTipInfoFile> \\
-g <geneListFile> \\
-F '60 30' \\
-q fasttree \\
-c 8 \\
-C 1 \\
-Q medium \\
<path_to_recovered_genes_from_samples>/*.fasta \\
> make_species_trees_pipeline.log 2>&1 &

EOF
}


# Hidden parameters which work if you know about them:
# -m  <float> maximum column occupancy (default=0.7; 0.7 means that columns with 70 % residue occupancy will be counted)


# Extensions to options not inlcuded yet:
# option -F 3rd field. Use the well conserved regions in each alignment only to build an additional tree (0=no, 1=yes)
#                        NBNB - NEED TO WORK OUT HOW THIS WORKS FOR DNA AND PROTEIN - LOTS OF FILES NOW   
#                        - working dir: conservedDNARegionsTree, conservedProteinRegionsTree, conservedCodonRegionsTree


#echo User inputs:    ### For testing only 
while getopts "hvat:ug:ijGF:m:p:M:q:r:TC:c:d:Q:A:D:O:L:I:JK:R:U:V:W:"  OPTION; do	# Remaining options - try capital letters!

	#echo -$OPTION $OPTARG	### For testing only - could try to run through options again below 
	 
	case $OPTION in

		h) usage; exit 1 ;;
		v) echo "make_species_trees_pipeline.sh  version 0.0.1"; exit ;;
		a) addSampleName=yes ;;
		t) sampleTableFile=$OPTARG ;;
		u) option_u=yes ;;
#ALIGNMENT OPTIONS
		D) seqType=$OPTARG ;;
		A) alnProgram=$OPTARG ;;

		g) geneListFile=$OPTARG ;;
		i) geneTreesOnly=yes ;;
		j) speciesTreesOnly=yes ;;
		G) useGenewiseFiles=yes ;;
		F) filterSeqs1=$OPTARG ;;
		#f) fractnAlnCovrg=$OPTARG ;;
		m) fractnMaxColOcc=$OPTARG ;;
		#s) fractnSamples=$OPTARG ;;
		p) fileNamePrefix=$OPTARG ;;
		M) mafftAlgorithm="$OPTARG" ;;
		q) phyloProgramDNA=$OPTARG ;;
		r) phyloProgramPROT=$OPTARG ;;
		T) treeshrink=yes ;;
		C) cpuGeneTree=$OPTARG ;;
		c) cpu=$OPTARG ;;
		Q) partitionName=$OPTARG ;;
		O) maxColOccThreshold=$OPTARG ;;
		L) collapseNodes=$OPTARG ;;
		I) filterSeqs2=$OPTARG ;;
		J) trimAln1=yes ;;
		K) trimAln2=$OPTARG ;;
		R) geneTreeSlurmMem=$OPTARG ;;
		U) speciesTreeSlurmMem=$OPTARG ;;
		V) geneTreeSlurmTime=$OPTARG ;;
		W) speciesTreeSlurmTime=$OPTARG ;;
		?)  echo This option is not allowed. Read the usage summary below.
      	    echo
      	    usage; exit 1 ;;
      	:) echo "Invalid option: $OPTARG requires an argument" 1>&2 ;;		# 26.7.2020 - Added - captures error for options that require an argument - useful - actually doesn't work if option requiring a value is followed by a flag instead
    esac
done


###########################################
# Set up access to all scripts for pipeline
###########################################
# This script accesses several more scripts in the same directory.
# Unless users set the path to this directory they will not be acessible.
# So, pre-fixing the extra scripts with the path of this script.
# NB - I think $0 is always the full path. 
fullPathToWrapperScript=`echo $0 `
#echo  $fullPathToWrapperScript
pathToScripts=`dirname $fullPathToWrapperScript `


#############################
# User input checks and setup
#############################

# Prints usage if no parameters are given:
if [ "$#" -lt 1 ]; then usage; exit 1; fi


###################################### 
# Check required software dependancies
######################################
softwareList=(seqtk fastatranslate fastalength $phyloProgramDNA $phyloProgramPROT nw_ed java astral.5.7.3.jar AMAS.py raxmlHPC-PTHREADS-SSE3)	# fasttree raxml-ng now tested in $phyloProgramDNA and $phyloProgramPROT 
### Difficult to test the following softwares in this way in an array - could try to test separately: 'bc --help' 'est2genome --help' 'mafft --help'
###		28.1.2020 - try to double quote the array to keep these cmds with spaces together!
### NB - astral.5.6.3.jar also may not be straight forward to check!!!! Also fasttreeMP not on Macbook - fasttree is the minimum.
### NBNBNB - in recover_genes_from_all_samples.sh, running seqtk like this makes the script stop, but here it doesn't seem to!!!! It's because there a no set options in this script
#echo ${softwareList[@]}
### NB - 29.8.2020 - only need to check the software actually being used e.g. upp, maffft, phylog programs - could check software when checking the  
for software in ${softwareList[@]}; do
	if [[ $software != 'no' ]]; then			# The default state $phyloProgramPROT and some other software - so don't want to test that!
		echo Testing $software is installed...
		$software >/dev/null 2>&1
		if [[ $? == 127 ]]; then
			# Exit code 127 is for "command not found"
			### 16.3.2020 - NB - this doesn't work as script exits in the above conditional if program is not found!
			echo "ERROR: Not available: $software"
			exit
		fi
	fi
done


# Check positional parameters are present i.e. the gene recovery fasta file(s) are present:
# Command line variables summary:
echo '$# == ' $#										# Total number of all parameters (excludes script name, includes flags and their values, excludes free parameters (ones with no flags)) 
echo \$OPTIND == $OPTIND								# Position of the first free parameter after any options - free parameters must come after any optional parameters.
#echo 'Value of first free parameter: ' ${@:$OPTIND:1}	# Lists the value of the first free parameter from the $@ variable; ${@:$OPTIND:2} will access the first two free parameters.
echo 'Values of all free parameters: '${@:$OPTIND:$#}	# Therefore ${@:$OPTIND:$#} will access all the free parameters
echo $(( $# - $OPTIND + 1 ))							# Number of free parameters, in this case the number of samples
numbrSamples=$(( $# - $OPTIND + 1 ))
#echo $numbrSamples


if [ $(( $# - $OPTIND + 1 )) -lt 4 ]; then				### 30.3.2020 AND speciesTreesOnly == no to handle scritp just processing species trees
														###		Would need to think how to specify the input file(s) !!!!!!!!!!! 
														### 4.7.2020 AND now if gene-wise files == no are entered - can have less than one of those right?
														### 12.8.2020 - Just have an if else clause and say: less than 1 files but in gene-wise mode so OK 
    echo
    echo "ERROR: you need to input fasta files containing recovered genes from at least four species!"
    exit 1
fi


### 16.1.2020 - on hold for the moment
# Set up directory for temporary outputs e.g. gene alignments, gene trees 
###tempOutputsDir=temp_outputs_`date | awk '{print $3 "_" $2 "_" $6 "_" $4}' | sed 's/://g' `$RANDOM
# Create unique temporary directory:
###mkdir $tempOutputsDir


#######################
# Check the input files
#######################

# 1. Check the input files are in fasta format
###for file in ${@:$OPTIND:$#}; do 
	###echo "Check fasta file format here - still to do:" $file

	# Still to do:
	# First just check they are all files - checks commandline parameters haven't leaked through
	# In fact it is really important to check also that file exists and has contents because there may have been no genes recovered!!
	#	actually it may be ok that the file is empty 
	# 1.See IT3F_frame.sh script e.g. file needs to start with >
	# 2.Check there are no empty seq lines - see my 46.notes on assembly checks
	# 3.Check equal # header lines and seq lines
	# 4.nw_ed doesn't like identifiers containing [ and ] chars so should change them or exit with error

	### 15.1.2019 - also print all file paths to a file to select later - poor and OK/good quality - BUT - check it will work
###done
	### 27.4.2020 - I think gene names mustn't contain dots (for extracting a gene from all samples - see line ~59 in make_gene_trees.sh)
 	###             Also for HybPiper output files, sample names mustn't have dots! 



### 7.10.2019 - NOT YET STRESS TESTED ALL SITUATIONS
### 16.7.2020 - also does it work for dna-wise files; also go through logic - doesn't make too much sense now.
# 2. First ensure all filenames supplied are unique:
numbrDuplicateNames=`for file in ${@:$OPTIND:$#}; do 
					uniqueSampleId=$(basename $file | awk -F '.' '{print $1}' )
					done | echo $uniqueSampleId | sort | uniq -c | awk '$1 > 1' | wc -l `
#echo numbrDuplicateNames: $numbrDuplicateNames
if [[ $numbrDuplicateNames -ge 1 ]]; then echo "ERROR: there are $numbrDuplicateNames duplicate names in the input fasta files."; exit 1; fi


# 3. Test that the input files and their paths haven't reached ARG_MAX.
#    getconf ARG_MAX = 262,144 (Macbook) and 2,097,152 on Hatta cluster
lengthOfFreeParameters=`echo ${@:$OPTIND:$#} | awk '{sum+=length($0)} END {print sum}' `
echo lengthOfFreeParameters: $lengthOfFreeParameters
systemARG_MAX=`getconf ARG_MAX `
if [ $lengthOfFreeParameters -gt $systemARG_MAX ]; then
	echo
	echo 'ERROR: the length of the input fasta files an their paths has exceeded the system ARG_MAX variable.'
	echo 'Try to reduce the length of the path to the files or increase the system ARG_MAX value somehow.'
	exit
fi


########################################################################
# Prepare recovered genes fasta headers to contain gene Id and sample Id
########################################################################
# Formats considering:
# 1. >sampleId-geneId (default)
# 2. >geneId (option -a)
#    User has the sample files named in a unique way, then the file name is added to the 
#    fasta header as the sample name: >geneId sampleId
# 3. User supplies a table containing the species info for each sample to include in the tree tips (option -t)
#	 This info has to be matched with some unique info in the fasta filename - same as option -a
#	 It would have to always be the same field in the table and be unique
#	 NB - I think that the csv file can contain as many samples as you want but script will only use the samples with fasts
#	 files entered on the command line 


# First check that only one of option -a and -t are being used:
###.10.5.2020 - don't need this check anymore  
# if [[ $addSampleName == 'yes' &&  -s $sampleTableFile ]]; then
# 	echo Both option -a and option -t selected - re-run with only one of these options - exiting. 
# fi


####################
# Code for option -a
####################
if [[ $addSampleName == 'yes' && $useGenewiseFiles  != 'yes' ]]; then

	echo Will add sample name from filename: $addSampleName


	### 15.1.2020 - if you accidentally have the correct format with this option, then it will try to process samples
	### but the geneids will also contiain the sample id so no genes will be recognised
	### Just have a conditional to exit if any gene name contains a hyphen

 	# Prepare fasta files:
	for file in ${@:$OPTIND:$#}; do

    	# Get filename and chop off the file ending if there is one:
    	uniqueSampleId=`basename $file | awk -F '.' '{print $1}' `
    	#echo $uniqueSampleId
   
    	# Add filename i.e. species id to the fasta header lines:
    	cat  $file | awk -v sampleId=$uniqueSampleId '{if($1 ~ /^>/) {print $1 " " sampleId} else {print $0}}' \
    	> ${uniqueSampleId}_modified.fasta
    	### NB - need to change awk code when altering format
		### NB - keep an eye on the max line length allowed w.r.t. seqtk
		### 14.7.2020 - consider to do what I've done just below for ${geneName}_dna.fasta and check  ${uniqueSampleId}_modified.fasta doesn't already exist
		###				Actually not such a good idea - means that people can't over-run previuous run 

		# If an input fasta file name doesn't exist then the 'modified.fasta' created above will not exist.
		### 9.7.2020 - Can't remember why it is important to test - I think it because otherwise you end up with failed gene tree runs!  
		# Testing whether file exists here:
		if [[ ! -s ${uniqueSampleId}_modified.fasta ]]; then
			echo "ERROR: this input fasta file does not exist or is empty: ${uniqueSampleId}_modified.fasta"
			### Consider to add: If it exists and is empty, rename this file so that it is not picked up again e.g. to <prefix>.fasta_not_used"
			exit 1
		fi
    done
### NB - 16.6.2020 - do I need to check the above files like I do for the format below? - I think I do because even though I'm creating it above
### they may be in a different format e.g. gene-wise files!!!!!!!! Just need to check with code below # samples in the data set
    # Concatenate fasta files and put seqs on a single line:
	fileList=`ls *_modified.fasta `
	echo $fileList
	cat $fileList | seqtk seq -l 0 /dev/fd/0 > all_samples_concatenated.fasta
	geneFile=all_samples_concatenated.fasta


elif [[ $useGenewiseFiles == 'yes' && $addSampleName != 'yes' ]]; then

	echo Will use prepared gene-wise files directly: $useGenewiseFiles

	# Copy the fasta files to exactly <geneName>_dna.fasta, then they are ready for the alignment step.
	# The fasta filename (minus the dot suffix) needs to be exactly the name of the gene in the genelist.
	for file in ${@:$OPTIND:$#}; do

		if [[ ! -s $file ]]; then
			echo "ERROR: this input fasta file does not exist or is empty: $file"
			exit 1
		fi

		# Get filename and chop off the file ending if there is one:
    	geneName=`basename $file | awk -F '.' '{print $1}' `
    	#echo $geneName
		cp $file ${geneName}_dna.fasta
		### 11.8.2020 - consider to reanme to e.g.  ${geneName}_dna_for_aln.fasta and move rather than copy (one less file)
		if [[ $? == 1 ]]; then
			### I've kept this in - it avoids _dna.fasta file stamping on itself if these files are submitted by the user 
			echo "ERROR: if running pipeline (possibly again) in gene-wise mode, please run in a another directory
and input gene-wise fasta files with a relative path (probably from a previous run) e.g. ../old_run/*_dna.fasta"
			exit 1
		fi
	done
	geneFile='use_genewise_files'	# Used to change conditional statements

	# Need to recalculate numbrsamples variable here - up till here in script, $numbrSamples contains the number of gene-wise files entered.
	# NB - Fasta file format here is: >sampelId.
    # NB - on MacOS awk inserts a blank line between output lines so removing them with grep -v '^$'.
    numbrSamples=`cat *_dna.fasta | awk '{if($1 ~ /^>/)  {print $1} }' | grep -v  '^$' | sort -u | wc -l `
else
	# Format is already ready to go (>sampleId-geneId).
   
	### Still deciding on whether to change to this format completely throughout...
	### For the moment, altering the ">species-gene" format back to my "gene species" format:
	for file in ${@:$OPTIND:$#}; do 

### 7.10.2019 - Do I need to do this here????? 
	 	# Get filename and chop off the file ending if there is one:
     	uniqueSampleId=`basename $file | awk -F '.' '{print $1}' `
     	#echo Sample Id: $uniqueSampleId

	 	# Convert from fasta format >sampleId-geneId to >geneId sampleId (my internal format).

	 	# First separate out the fasta Id field (anything before the first white space),
	 	# then check the dash format, checking both sides of the dash,
	 	# more than one dash in gene or sample ids may result in non-unique sample/gene ids, if so script will exit.
		### NB - this step might need thinking about more if you wanted to also include paralogs - each paralog should still be unique but not the speciesId!! 
### 7.10.2019 - NOT YET TESTED ALL SITUATIONS
     	# There should only be one sample identifier per file, so checking whether sample identifier is unique:
     	beforeDashCheck=`cat $file | awk '{if($1 ~ /^>/)  {print $1} }' |  awk -F '-' '{print $1}' | sort | uniq -c | awk '$1 == 1' | wc -l`
     	# All gene identifiers should be unique, so checking if there are any duplicates"
     	afterDashCheck=`cat $file | awk '{if($1 ~ /^>/)  {print $1} }' |  awk -F '-' '{print $2}' | sort | uniq -c | awk '$1 > 1' | wc -l `
     	#echo "beforeDashCheck (sampleId): " $beforeDashCheck 
     	#echo "afterDashCheck (geneId): " $afterDashCheck 
     	if [[ $beforeDashCheck -gt 1 ]]; then echo "ERROR: more than one sample identifier detected on fasta header line (there should only be one sample identifier for the default format) for this sample: $file.
     		Also check that fasta header lines have this format: >sampleId-geneId"; exit 1
     	elif [[ $afterDashCheck -gt 1 ]]; then echo "ERROR: gene identfiers should be unique, one or more not unique for this sample: $file \nAlso check that fasta header lines have this format: >sampleId-geneId"; exit 1
		else cat $file | awk '{if($1 ~ /^>/)  {print $1} else {print $0}}' \
			| awk -F '-' '{if($1 ~ /^>/) {{gsub(/>/,"",$1)} {print ">" $2 " " $1}} else {print $0}}' \
			> ${uniqueSampleId}_modified.fasta
			# NB - fasta format now: >geneId sampleId

			# If an input fasta file name doesn't exist then the 'modified.fasta' filename created above will not exist.
			# Testing whether file exists here:
			if [[ ! -s ${uniqueSampleId}_modified.fasta ]]; then
				echo "ERROR: this input fasta file does not exist or is empty: ${uniqueSampleId}_modified.fasta"
				exit 1
			fi
		fi

		### [[To convert from 'gene species' my format to species-gene format:
		###cat CKDK_Ochnaceae_Ochna_serrulata.fasta | awk '{if($1 ~ /^>/)  {{gsub(/>/,"",$1)} {print ">" $2 "-" $1}} else {print $0}}' >species-gene-format.fasta ]]
		### NB - need to change awk code when altering format in next script
		### NB - keep an eye on the max line length allowed w.r.t. seqtk
    done
	
    # Finally, test whether there are the same number of sample identifiers on the fasta header lines as there are samples submitted:
### 7.10.2019 - NOT YET TESTED ALL SITUATIONS
    # NB - Fasta file format is now: >geneId sampelId.
    # NB - on MacOS awk inserts a blank line between output lines so removing them with grep -v '^$'.
    numbrSamplesInFile=`cat *_modified.fasta | awk '{if($1 ~ /^>/)  {print $2} }' | grep -v  '^$' | sort -u | wc -l `
    ### NB - if file can get through here with no sample field occupied (it can't I don't think), then you would end up counting genes, then number of genes could equal number of samples (unlikely).
    if [ $numbrSamples -ne $numbrSamplesInFile ]; then echo "ERROR: number of samples counted according to the fasta header lines is different to the number of input files"; exit 1
    else
		# Concatenate fasta files and put seqs on a single line:
		fileList=`ls *_modified.fasta `
		#echo $fileList
		cat $fileList | seqtk seq -l 0 /dev/fd/0 > all_samples_concatenated.fasta
		geneFile=all_samples_concatenated.fasta
	fi
fi
#exit


####################
# Code for option -t
####################
# First need to determien whether option is even selected.
if [[ $sampleTableFile != 'no' ]]; then
	if [[ -s $sampleTableFile ]]; then
		echo Will use this file containing text for the tree leaves: $sampleTableFile
		# First check that the table is a csv file (plus header - is possible?)
		### I think I can check here the csv file - must start with the table header! - only check needed - maybe check that it has ',' chars (?)
		


		# Using the filename to associate sample with the correct row of the csv file provided.
		# The user needs to create a list of filenames and add the tree tip info to each row.
		# (How else could this be done? Maybe getting user to supply a list of sample names with <sample_names>.fasta file
		# would be an alternative (like gene recovery step), but it's also convenient not to have to prepare another file.)

		# Prepare a mapfile for Newick Utils to switch in sample info onto the tree tips:
		addTreeTipInfoFromTable.py  $sampleTableFile  tree_tip_info_mapfile.txt  ${@:$OPTIND:$#}
		### 29.1.2020 - not now using ${@:$OPTIND:$#}
		# Check that tree_tip_info_mapfile.txt now exists, if not an identifier is not unique in the table:
		if [[ ! -s tree_tip_info_mapfile.txt ]]; then
			echo "The sample table (from option -t) contains a sample Id that is not unique in the table. Exiting..."
			exit
		fi
	else
    	echo File containing text for the tree leaves does not exist or is empty: $sampleTableFile - exiting.
		exit 1
	fi
fi


if [ -s $geneListFile ]; then echo ""
	### NB - 26.4.2020 - just realised that the gene names should not contain any dot chars - see note in make_gene_trees.sh ~ line 59
else echo "ERROR: the gene list file (option -g) does not exist or is empty: $geneListFile"; exit; fi


### 5.10.2020 - no longer need these two checks - being done further done
# NB - Trying to check that value is an integer - how does this behave when adding 'hello'? looks OK; BUT allows through chars A -F - why?
# if (( $(bc <<< "$fractnAlnCovrg < 0 ") || $(bc <<< "$fractnAlnCovrg > 1 ") )); then
# 	echo "ERROR1: -f option needs to be a fraction between 0 and 1 - exiting"
# 	exit
# fi  


# if (( $(bc <<< "$fractnSamples < 0") || $(bc <<< "$fractnSamples > 1") )); then 
# 	echo "ERROR: -s option needs to be a fraction between 0 and 1 - exiting"
# 	exit
# fi


if (( $(bc <<< "$fractnMaxColOcc < 0") || $(bc <<< "$fractnMaxColOcc > 1") )); then 
	echo "ERROR: -m option needs to be a fraction between 0 and 1 - exiting"
	exit
fi


# Check first for whether $cpu is an integer (this is a logical fudge for detecting whether an integer or string) 
### For maxColOcc_pc 
#if [ "$cpu" -eq "$cpu" ] 2>/dev/null
#then echo ""
#else echo "ERROR: option -c should be an integer - exiting "; exit; fi  


# Check first for whether $cpu is an integer (this is a logical fudge for detecting whether an integer or string) 
if [[ "$cpu" -eq "$cpu" || "$cpuGeneTree" -eq "$cpuGeneTree" ]] 2>/dev/null
then echo ""
else echo "ERROR: option -c should contain an integer value - exiting "; exit; fi

if [ -z $fileNamePrefix ]; then 	#echo ""
	echo "ERROR: option -p used but an output file name prefix string was not supplied - exiting"
	exit
	### This doesn't really work because if there is no paramter value the next parameter becomes the value which is a string!
fi

#if [[ -z $phyloProgram || $phyloProgram != 'fasttree' || $phyloProgram != 'raxml-ng' ]]; then echo ""
#	echo "ERROR: option -q used but a phylogeny program name was not supplied or was incorrectly spelt - exiting"
	### 16.3.2020 - this really doesn't work - is triggered no matter what! Now have a solution in make_gene_trees.sh
	### This doesn't really work because if there is no paramter value the next parameter becomes the value which is a string!
#fi


# Input checks for the sequence type option:
dnaSelected=no
proteinSelected=no
codonSelected=no
if [[ `echo $seqType | grep -o 'dna' ` == 'dna' ]];then 
	dnaSelected=yes
fi
if [[ `echo $seqType | grep -o 'protein' ` == 'protein' ]];then 
	proteinSelected=yes
fi
if [[ `echo $seqType | grep -o 'codon' ` == 'codon' ]];then 
	codonSelected=yes
fi
#echo dnaSelected=$dnaSelected proteinSelected=$proteinSelected codonSelected=$codonSelected
if [[ $dnaSelected == 'no' &&  $proteinSelected == 'no' && $codonSelected == 'no' ]]; then
	echo "ERROR: No sequence type (option -D) was entered or recognised."; exit
fi


# Input checks for the alignment type option -A:
if [[ "$alnProgram" != 'mafft' && "$alnProgram" != 'upp' ]];then 
	echo "ERROR: No alignment program (option -A) was entered or recognised."; exit
fi


# Check ilter sequence options (option -F):
if [[ $filterSeqs1 != 'no' ]]; then
	numbrFields=`echo $filterSeqs1 | awk '{print NF}' `
	if [[ $numbrFields -eq 2 ]]; then
		# Check variable has a sensible integer value range, then whether it is an integer value (uses a logical fudge for detecting whether an integer or string, including non-alphanumeric characters)  
		# then whether there is a value (not sure whether last check is necessary now I check for an integer)
		fractnAlnCovrg=`echo $filterSeqs1 | awk '$1 >= 0 && $1 <= 100 {print $1}' `
		if [[ "$fractnAlnCovrg" -eq "$fractnAlnCovrg" ]] 2>/dev/null; then echo ""
		else echo "ERROR: 1st value of option -F should contain an integer value between 1 and 100 - exiting "; exit; fi
		if [[ -z "$fractnAlnCovrg" ]]; then echo "ERROR: option -F should contain an integer value between 1 and 100 - exiting "; exit; fi
		# Eventually could use these percentages directly, but meanwhile changing into a fraction to use with the existing code:	
		fractnAlnCovrg=`echo $filterSeqs1 | awk 'fractn=$1/100 {print fractn}' `

		fractnSamples=`echo $filterSeqs1 | awk '$2 >= 0 && $2 <= 100 {print $2}' `
		if [[ "$fractnSamples" -eq "$fractnSamples" ]] 2>/dev/null; then echo ""
		else echo "ERROR: 2nd value of option -F should contain an integer value - exiting "; exit; fi
		if [[ -z "$fractnSamples" ]]; then echo "ERROR: option -F should contain an integer value between 1 and 100 - exiting "; exit; fi
		fractnSamples=`echo $filterSeqs1 | awk 'fractn=$2/100 {print fractn}' `
	else
		echo "ERROR: Filter sequences option (option -F) was selected but should contain two values."; exit
	fi
fi


# Check option to collapse gene tree nodes (option -L):
if [[ $collapseNodes != 'no' ]]; then
	# First check variable is an integer (uses a logical fudge for detecting whether an integer or string, including non-alphanumeric characters) 
	# then check variable has a sensible integer value range, then whether there is a value (not sure whether last check is necessary now I have the first integer check)
	if [[ "$collapseNodes" -eq "$collapseNodes" ]] 2>/dev/null; then echo ""
	else echo "ERROR: option -L should contain an integer value - exiting "; exit; fi
	collapseNodes=`echo $collapseNodes | awk '$1 >= 0 && $1 <= 100 {print $1}' `
	if [[ -z "$collapseNodes" ]]; then echo "ERROR: option -L should contain an integer value - exiting "; exit; fi
fi


# Check option to filter alignment columns to remove rarer insertions (trimAln2, option -K):
if [[ $trimAln2 != 'no' ]]; then
	# First check variable is an integer (uses a logical fudge for detecting whether an integer or string, including non-alphanumeric characters) 
	# then check variable has a sensible integer value range, then whether there is a value (not sure whether last check is necessary now I have the first integer check)
	###if [[ "$trimAln2" -eq "$trimAln2" ]] 2>/dev/null; then echo ""
	###else echo "ERROR: option -K should contain an integer value - exiting "; exit; fi
	### 5.10.2020 - actually can't perform the above integer check because the value could be less than 1 % - have to trust the user!
	trimAln2=`echo $trimAln2 | awk '$1 > 0 && $1 <= 100 {print $1}' `
	if [[ -z "$trimAln2" ]]; then echo "ERROR: option -K should contain an integer or fraction value - exiting "; exit; fi
	### Could convert a user percent value to a fraction expected by AMAS.py to standardize the inputs, somthing like: 
	### fractnAlnCovrg=`echo $filterSeqs1 | awk 'fractn=$1/100 {print fractn}' `
fi


# Some parameters selected for checking:
echo "################"
echo "Options selected"
echo "################"
echo dnaSelected: $dnaSelected
echo proteinSelected: $proteinSelected

echo 'filter sequence option 1 (option -F): ' $filterSeqs1
echo 'filter sequence option 2 (option -I): ' $filterSeqs2
#echo '	fractnAlnCovrg: ' $fractnAlnCovrg
#echo '	fractnSamples: ' $fractnSamples
echo 'collapseNodes (option -L): ' $collapseNodes
echo 'trimAln1 (option -J): ' $trimAln1
echo 'trimAln2 (option -K): ' $trimAln2
echo 'treeshrink: ' $treeshrink
echo 'filterSeqs1: ' $filterSeqs1
echo 'geneTreeSlurmMem: ' $geneTreeSlurmMem
echo 'speciesTreeSlurmMem: ' $speciesTreeSlurmMem
echo 'geneTreeSlurmTime: ' $geneTreeSlurmTime
echo 'speciesTreeSlurmTime: ' $speciesTreeSlurmTime
echo
#exit


# In case analysis is re-run in the same working directory AND the filtering/thresholds have beeen changed,
# need to delete the existing alignment and tree files first, otherwise files for some genes may 
# be used from a previous run, even if they have been filtered out in the current run. Probably is
# only relevant for small datasets:
#### 29.9.2020 - still get an error message though!!!!!
#### 1.10.2020 - hey why don't I just put all the other files under the umbrella of the first line, there will always be some *.aln.for_tree.fasta files!
### OH NO THAT WON'T WORK!!!!!!
if ls *.aln.for_tree.fasta >/dev/null 2>&1; then rm *.aln.for_tree.fasta *gene_tree_USE_THIS.nwk ; fi
# Also, now I have added filterShortSeqs() function I need to check whether any of the filtered/trimmed fasta files need removing:
if ls *.aln.after_filter1.fasta >/dev/null 2>&1; then rm *.aln.after_filter1.fasta ; fi
if ls *.aln.after_filter2.fasta >/dev/null 2>&1; then rm *.aln.after_filter2.fasta ; fi
if ls *.aln.after_trim1.fasta >/dev/null 2>&1; then rm *.aln.after_trim1.fasta ; fi
if ls *.aln.after_trim2.fasta >/dev/null 2>&1; then rm *.aln.after_trim2.fasta ; fi
# Other filter/trim output files should go here.
# Also (!), if there is no filtering or trimming done then the file going into 
# tree building is *.dna.aln.fasta so should be removed as well:
if ls *.aln.fasta >/dev/null 2>&1; then rm *.aln.fasta ; fi
# Also w.r.t. the realignment step, also need to delete these files (they mustn't linger around, otherwise gene set will be used, if if already filtered out)
if ls after_reAlnFilterSeqs_USE_THIS/*after_filterSeqs_dna.fasta >/dev/null 2>&1; then rm after_reAlnFilterSeqs_USE_THIS/*after_filterSeqs_dna.fasta ; fi
if ls after_reAlnFilterSeqs_USE_THIS/*after_treeshrink_dna.fasta >/dev/null 2>&1; then rm after_reAlnFilterSeqs_USE_THIS/*after_treeshrink_dna.fasta ; fi


##########################
# End of user input checks
##########################


echo 'Detecting operating system: '
os=`uname `
echo $os
if [ $os == 'Darwin' ]; then
    echo 'Running as if on Mac OS (Darwin)...'
elif [ $os == 'Linux' ]; then
	echo 'Running as if on Linux...'
    slurm=`sbatch -V | grep ^slurm | wc -l `
    if [ $slurm -eq 1 ]; then 
        echo 'Slurm scheduler detected, running via slurm...'
        ### 28.6.2019 - could create a variable here that contains call to "lsrun -J <$jobName>"
        ###	remain empty if slurm not used - need to be able to set $jobName for each program.
	fi
else
    echo 'Unrecognised OS that this software can run on, exiting.'
    exit
fi
#echo slurm: $slurm
#exit


### NB - noted previously that as well as $numbrSamples/species input, I need a list of samples ids - see tree script


# Move into $tempOutputsDir and do all the work in there:
###cd $tempOutputsDir

###########################
echo 'Making gene trees...'
###########################
if [[ $os == 'Darwin' && $speciesTreesOnly == 'no' ]]; then

    cat $geneListFile | \
	while read line ; do
		geneId=`echo $line | tr -d '\n' `	# Need to remove the line return!
		exePrefix="/usr/bin/time -l"		# this time command gets the RSS memory, -l flag doesn't work on Cluster
		echo
		echo
		echo ####################
		echo Processing gene $geneId
		echo ####################
		$pathToScripts/make_gene_trees.sh \
		$geneId \
		$geneFile \
		$fractnAlnCovrg \
		$phyloProgramDNA \
		$phyloProgramPROT \
		$fractnMaxColOcc \
		$cpuGeneTree \
		"$mafftAlgorithm" \
		"$exePrefix" \
		"$alnProgram" \
		$dnaSelected \
		$proteinSelected \
		$codonSelected \
		"$filterSeqs1" \
		$pathToScripts \
		$maxColOccThreshold \
		$filterSeqs2 \
		$trimAln1 \
		$trimAln2 \
		"$treeshrink"
	done > make_gene_trees.log 2>&1
	#exit
	if [[ $filterSeqs1 != 'no' ]]; then
		# Need to name the files to use in the assess script, depends on sequence type selected. 
		if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then
			seqType=protein
			###alnFileSuffix=${seqType}.aln.for_tree.fasta		# before AMAS trim - consider to add this or just do stats at end of all filtering+trimming 
			alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta  # after AMAS trim
			## Could add other filenames used in script (?) 
		else
			seqType=dna
			alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta
		fi
		### 6.10.2020 - NBNB - what about codon aln files - I assume they are required as well if used????????
		echo seqType: $seqType
		echo alnFileForTreeSuffix: $alnFileForTreeSuffix
		$pathToScripts/assess_gene_alignments.sh \
		$fractnAlnCovrg \
		$fractnMaxColOcc \
		$fractnSamples \
		$numbrSamples \
		$fileNamePrefix \
		$geneFile \
		$sampleTableFile \
		$option_u \
		$seqType \
		$alnFileForTreeSuffix \
		> assess_gene_alns.log 2>&1
	fi

	if [[ $treeshrink == 'yes' || $filterSeqs1 != 'no' || $filterSeqs2 != 'no' ]]; then
		##########################################
		echo 'Re-aligning gene alignments because TreeShrink option or filter sequences option 1 or 2 is on, then continuing analysis in the "after_treeshrink_USE_THIS" or "after_reAlnFilterSeqs_USE_THIS" directory...'
		##########################################
		run_treeshrink_and_realign.sh \
		"$numbrSamples" \
		"$phyloProgramDNA" \
		"$phyloProgramPROT" \
		"$sampleTableFile" \
		"$geneListFile" \
		"$fractnAlnCovrg" \
		"$fractnMaxColOcc" \
		"$fractnSamples" \
		"$mafftAlgorithm" \
		"$cpuGeneTree" \
		"$partitionName" \
		"$pathToScripts" \
		"$geneTreesOnly" \
		"$dnaSelected" \
		"$proteinSelected" \
		"$codonSelected" \
		"$treeshrink" \
		"$filterSeqs1" \
		"$alnProgram" \
		"$maxColOccThreshold" \
		"$filterSeqs2" \
		"$trimAln1" \
		"$trimAln2" \
		> run_treeshrink_and_realign.log 2>&1
		exit	# Species trees will be made after TreeShrink or re-alignment step(s) in nested call to this script, if requested.
	fi

elif [[ $os == 'Linux' && $speciesTreesOnly == 'no' ]]; then
	###exePrefix="/usr/bin/time -v -o g${gene}_mafft_dna_aln_time_and_mem.log"		# NB - this will not work here - need to pick up gene id in Slurm script instead.
    slurm=`sbatch -V | grep ^slurm | wc -l `
    if [ $slurm -eq 1 ]; then
		# Count the # genes to process and fix that number in the Slurm --array parameter.
    	# It has to be set outside the sbatch script I think so that I can automatically set the array size .
    	# It has to be set each time otherwise.
    	numbrGenes=`tail -n+1 $geneListFile | wc -l `
    	numbrGenes=$(( $numbrGenes - 1 ))		# Slurm array needs to start at zero so need to adjust the max value.

    	# Notes on the --array flag - e.g. SBATCH --array=1-${numbrSamples}%50 
    	# 0-352%10	# NB - if input file has a header line, array should start at 1 (but the gene list file doesn't!);
    	# Maximum array size by default is 1000 - can be increased with MaxArraySize (max supported = 4000001, default set to 1001)
    	# NBNB - it's only an administration value set globally i think  - can't be changed by user
		# Prevent jobs going over 1000, otherwise will get empty sample directories and log files - minor issue!:
    	if [ $numbrGenes -ge 1001 ]; then
      		$numbrGenes=1000
      		echo "WARNING: Maximum default size of Slurm array is 1000, processing only the first 1000 samples"
    	elif [ $numbrGenes -lt $slurmThrottle ]; then 
      		slurmThrottle=$numbrGenes
    	fi

        jobInfo=`sbatch -p $partitionName -c $cpuGeneTree -t $geneTreeSlurmTime --mem $geneTreeSlurmMem --array=0-${numbrGenes}%$slurmThrottle  $pathToScripts/slurm_setup_array_to_make_gene_trees.sh \
		$geneFile \
		$geneListFile \
		$fractnAlnCovrg \
		$pathToScripts \
		$phyloProgramDNA \
		$phyloProgramPROT \
		$fractnMaxColOcc \
		$cpuGeneTree \
		"$mafftAlgorithm" \
		"$exePrefix" \
		"$alnProgram" \
		$dnaSelected \
		$proteinSelected \
		$codonSelected \
		"$filterSeqs1" \
		"$maxColOccThreshold" \
		"$filterSeqs2" \
		"$trimAln1" \
		"$trimAln2" \
		"$treeshrink" `

		echo jobInfo: $jobInfo
		jobId=`echo $jobInfo | cut -d ' ' -f 4 `
		echo \$jobId: $jobId - same id as \$SLURM_ARRAY_JOB_ID - from running slurm_setup_array_to_make_gene_trees.sh
		# NB - If jobId variable is used after each call to Slurm then it doesn't matter if
		# a Slurm step is missed out - the jobId last assigned will alway be used
		if [[ $filterSeqs1 != 'no' ]]; then
			# Need to name the files to use in the assess script, depends on sequence type selected. 
			if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then
				seqType=protein
				###alnFileSuffix=${seqType}.aln.for_tree.fasta		# before AMAS trim - consider to add this or just do stats at end of all filtering+trimming 
				alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta  # after AMAS trim
				### Could add other filenames used in script (?) 
			else
				seqType=dna
				alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta
			fi
			echo seqType: $seqType
			echo alnFileForTreeSuffix: $alnFileForTreeSuffix

			jobInfo1=`sbatch -J assess_gene_alns  --dependency=afterok:$jobId -p $partitionName -c 1 -n 1 -o assess_gene_alns.log -e assess_gene_alns.err  $pathToScripts/assess_gene_alignments.sh \
			$fractnAlnCovrg \
			$fractnMaxColOcc \
			$fractnSamples \
			$numbrSamples \
			$fileNamePrefix \
			$geneFile \
			$sampleTableFile \
			$option_u \
			$seqType \
			$alnFileForTreeSuffix `

    		echo jobInfo1: $jobInfo1
			jobId=`echo $jobInfo1 | cut -d ' ' -f 4 `
			echo \$jobId: $jobId - from running assess_gene_alignments.sh
		fi
 

		echo treeshrink: $treeshrink
		echo filterSeqs1: $filterSeqs1
		if [[ $treeshrink == 'yes' || $filterSeqs1 != 'no' ]]; then
			##########################################
			echo 'Re-aligning gene alignments because TreeShrink option or filterSeqs1 option is on, then continuing analysis in the "after_treeshrink_USE_THIS" or "after_reAlnFilterSeqs_USE_THIS" directory...'
			##########################################
			jobInfo2=`sbatch -J trShrnk_realgn  --dependency=afterok:$jobId -p $partitionName -c 1 -n 1 -o run_treeshrink_and_realign.log  -e run_treeshrink_and_realign.err  $pathToScripts/run_treeshrink_and_realign.sh \
			"$numbrSamples" \
			"$phyloProgramDNA" \
			"$phyloProgramPROT" \
			"$sampleTableFile" \
			"$geneListFile" \
			"$fractnAlnCovrg" \
			"$fractnMaxColOcc" \
			"$fractnSamples" \
			"$mafftAlgorithm" \
			"$cpuGeneTree" \
			"$partitionName" \
			"$pathToScripts" \
			"$geneTreesOnly" \
			"$dnaSelected" \
			"$proteinSelected" \
			"$codonSelected" \
			"$treeshrink" \
			"$filterSeqs1" \
			"$alnProgram" \
			"$maxColOccThreshold" \
			"$filterSeqs2" \
			"$trimAln1" \
			"$trimAln2" `
			exit	# Species trees will be made after TreeShrink or re-alignment step(s) in nested call to this script, if requested.
		fi
	else
		cat $geneListFile | \
		while read line ; do
			geneId=`echo $line | tr -d '\n' `	# Need to remove the line return!
			exePrefix="/usr/bin/time -v"		# this time command gets the RSS memory, -l flag doesn't work on Cluster, GNULinux
			echo ####################
			echo Processing gene $geneId
			echo ####################
			$pathToScripts/make_gene_trees.sh \
			$geneId \
			$geneFile \
			$fractnAlnCovrg \
			$phyloProgramDNA \
			$phyloProgramPROT \
			$fractnMaxColOcc \
			$cpuGeneTree \
			"$mafftAlgorithm" \
			"$exePrefix" \
			"$alnProgram" \
			$dnaSelected \
			$proteinSelected \
			$codonSelected \
			"$filterSeqs1" \
			$pathToScripts \
			$maxColOccThreshold \
			"$filterSeqs2" \
			"$trimAln1" \
			"$trimAln2" \
			"$treeshrink"
		done > make_gene_trees.log 2>&1

		if [[ $filterSeqs1 != 'no' ]]; then
			# Need to name the files to use in the assess script, depends on sequence type selected. 
			if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then
				seqType=protein
				###alnFileSuffix=${seqType}.aln.for_tree.fasta		# before AMAS trim - consider to add this or just do stats at end of all filtering+trimming 
				alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta  # after AMAS trim
				## Could add other filenames used in script (?) 
			else
				seqType=dna
				alnFileForTreeSuffix=${seqType}.aln.for_tree.fasta
			fi
			echo seqType for assess script: $seqType
			echo alnFileForTreeSuffix: $alnFileForTreeSuffix
			$pathToScripts/assess_gene_alignments.sh \
			$fractnAlnCovrg \
			$fractnMaxColOcc \
			$fractnSamples \
			$numbrSamples \
			$fileNamePrefix \
			$geneFile \
			$sampleTableFile \
			$option_u \
			$seqType \
			$alnFileForTreeSuffix \
			> assess_gene_alns.log 2>&1
		fi
		echo treeshrink: $treeshrink
		echo filterSeqs1: $filterSeqs1
		if [[ $treeshrink == 'yes' || $filterSeqs1 != 'no' ]]; then
			##########################################
			echo 'Re-aligning gene alignments because TreeShrink option or filterSeqs1 option is on, then continuing analysis in the "after_treeshrink_USE_THIS" or "after_reAlnFilterSeqs_USE_THIS" directory...'
			##########################################
			run_treeshrink_and_realign.sh \
			"$numbrSamples" \
			"$phyloProgramDNA" \
			"$phyloProgramPROT" \
			"$sampleTableFile" \
			"$geneListFile" \
			"$fractnAlnCovrg" \
			"$fractnMaxColOcc" \
			"$fractnSamples" \
			"$mafftAlgorithm" \
			"$cpuGeneTree" \
			"$partitionName" \
			"$pathToScripts" \
			"$geneTreesOnly" \
			"$dnaSelected" \
			"$proteinSelected" \
			"$codonSelected" \
			"$treeshrink" \
			"$filterSeqs1" \
			"$alnProgram" \
			"$maxColOccThreshold" \
			"$filterSeqs2" \
			"$trimAln1" \
			"$trimAln2" \
			> run_treeshrink_and_realign.log 2>&1
			exit	# Species trees will be made after TreeShrink or re-alignment step(s) in nested call to this script, if requested.
		fi
	fi
fi


### Could toggle off/on this step - if no aln, then make pipeline stop at this point
if [ $geneTreesOnly == 'yes' ]; then exit; fi


################################
echo 'Making species tree(s)...'
################################
if [ $os == 'Darwin' ]; then
	exePrefix="/usr/bin/time -l"		# this time command gets the RSS memory, -l flag doesn't work on Cluster
    $pathToScripts/make_species_trees.sh \
	$fractnAlnCovrg \
	$fractnSamples \
	$numbrSamples \
	$fileNamePrefix \
	$geneFile \
	$cpu \
	$phyloProgramDNA \
	$phyloProgramPROT \
	"$exePrefix" \
	$sampleTableFile \
	$dnaSelected \
	$proteinSelected \
	$codonSelected \
	aln.for_tree.fasta \
	> ${fileNamePrefix}_make_species_trees.log 2>&1
elif [ $os == 'Linux' ]; then
	exePrefix="/usr/bin/time -v"
    if [ $slurm -eq 1 ]; then
    	### NB - not sure where to put the $exePrefix!!!!
    	### One option is to put the "time script" cmd in a wrapper but then I need a log file for this sbatch call and delete it from the script header..
    	### I think this is the only way without changing the main script itself.
		echo \$jobId: $jobId - should match previous Slurm step.
		# NB - previous Slurm jobs have to have an exit code of zero to satisfy Slurm --dependancy afterok:$jobId parameter,
		# otherwise would need to use --dependancy afterany:$jobId if exit code coudl be > 0.  
        sbatch --dependency=afterok:$jobId -p long -c $cpu -t $speciesTreeSlurmTime --mem=$speciesTreeSlurmMem -o ${fileNamePrefix}_make_species_trees.log -e ${fileNamePrefix}_make_species_trees.log  $pathToScripts/make_species_trees.sh \
        $fractnAlnCovrg \
        $fractnSamples \
        $numbrSamples \
        $fileNamePrefix \
        $geneFile \
        $cpu \
        $phyloProgramDNA \
        $phyloProgramPROT \
        "$exePrefix" \
        $sampleTableFile \
        $dnaSelected \
		$proteinSelected \
		$codonSelected \
		aln.for_tree.fasta
	else
		$pathToScripts/make_species_trees.sh \
		$fractnAlnCovrg \
		$fractnSamples \
		$numbrSamples \
		$fileNamePrefix \
		$geneFile \
		$cpu \
		$phyloProgramDNA \
		$phyloProgramPROT \
		"$exePrefix" \
		$sampleTableFile \
		$dnaSelected \
		$proteinSelected \
		$codonSelected \
		aln.for_tree.fasta \
		> ${fileNamePrefix}_make_species_trees.log 2>&1
	fi
fi