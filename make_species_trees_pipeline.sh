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
addSampleName=no
sampleTableFile=no			# 10.5.2020 - only just added in - check - ensures cmdline parameter is always accupied which is critical
option_u=no
geneTreesOnly=no
speciesTreesOnly=no
fractnAlnCovrg=0.6
fractnSamples=0.6
fileNamePrefix=tree_pipeline
phyloProgramDNA=fasttree
phyloProgramPROT=no			# Work around to specify any program so software testing code will not crash! Ensures cmd parameter is always occupied which is critical. 
cpu=8						# number of cpu to use for RAxML in supermatrix method
partitionName=main   		# Values depend on the cluster being used so good to have a flagged option for this
							# NB - make_species_tree.sh uses 'long' queue - need an extra variable for that 

# Hidden options (i.e. not apparent from the help menu but they always have a value so can be and have to be used in downstream scripts).
# Need to check them if I make them public
fractnMaxColOcc=0.7
slurmThrottle=1
mafftAlgorithm='--maxiterate 1000'     #'--maxiterate 1000' #'--retree 1'   '--retree 2' - use maxiterate
cpuGeneTree=1					# still keep this separate from the supermatrix tree; then I can specify loads more for the supermatrix.
								# Maybe hava a comma separated list e.g. 4,30 for gene and species tree.
speciesTreeMem=100000			# I think I need 500000 GB mem for the RAxML large tree (Slurm)
####-D  upload details to PAFTOL database (internal use only)   	


function usage	{

cat << EOF

Program description: makes species trees from fasta files containing recovered genes from multiple samples, one fasta file per sample.
                     By Default, the fasta header format MUST BE: >sampleId-geneId but one other option (-a) is available.

Usage: make_species_trees_pipeline.sh [options] fastafile1 fastafile2 fastafile3 ...

OPTIONS:
	-h             shows this message
	-v             program version  
	-a             add sample name onto the fasta header from the input fasta file name.
	               Expected gene identifier format in the fasta header: >geneId (no hyphen '-' characters)
	-t <csv file>  add sample name and other info from a comma separated value (csv) table file into the tree leaf labels.
	               Format of table row: sample_name/identifier, followed by any species information, including sample_name/identifier, as required. 
	-u             add contig length info onto species tree tips (requires option -t)
                   Sample_name/identifier must be identical to the sample fasta file name (minus any [dot] ending suffix e.g. .fasta)	  
	-g <file>      file (including path to it) containing list of gene names (required option)
	-i             make gene trees only
	-j             make species trees only
	-f <float>     fraction of [well conserved regions in] the alignment covered by a sample sequence.
	               Minimum to tolerate (default=0.6; 0 would mean no filtering, i.e. include sequence of any length)
	-s <float>     fraction of samples in each gene tree.
	               Minumum to tolerate (default=0.6; 0 would mean no filtering, include all available samples in each gene tree)
	-p <string>    prefix for the output filenames e.g. taxonomic group (default=tree_pipeline)
	-q <string>    name of phylogeny program for gene trees from DNA sequences.
                   	Options are, fastest to slowest: fasttree, iqtree2, raxml-ng (default=fasttree)
	-r <string>    name of phylogeny program for gene trees from protein sequences.
                   	If required, options are, fastest to slowest: fasttree2, iqtree, raxml-ng (no default)
	-c <integer>   number of cpu to use for RAxML in supermatrix method (default=8)
	-Q <string>    Slurm partition (queue) to use (default=medium) ]


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

EOF
}


# Hidden parameters which work if you know about them:
# -m  <float> maximum column occupancy (default=0.7; 0.7 means that columns with 70 % residue occupancy will be counted)


#echo User inputs:    ### For testing only 
while getopts "hvat:ug:ijf:m:s:p:q:r:c:d:Q:"  OPTION; do	# Remaining options - try capital letters!

	#echo -$OPTION $OPTARG	### For testing only - could try to run through options again below 
	 
	case $OPTION in

		h) usage; exit 1 ;;
		v) echo "make_species_trees_pipeline.sh  version 0.0.1"; exit ;;
		a) addSampleName=yes ;;
		t) sampleTableFile=$OPTARG ;;
		u) option_u=yes ;;
		g) geneListFile=$OPTARG ;;
		i) geneTreesOnly=yes ;;
		j) speciesTreesOnly=yes ;;
		f) fractnAlnCovrg=$OPTARG ;;
		m) fractnMaxColOcc=$OPTARG ;;
		s) fractnSamples=$OPTARG ;;
		p) fileNamePrefix=$OPTARG ;;
		q) phyloProgramDNA=$OPTARG ;;
		r) phyloProgramPROT=$OPTARG ;;
		c) cpu=$OPTARG ;;
		d) cpuGeneTree=$OPTARG ;;
		Q) partitionName=$OPTARG ;;  
		?)  echo This is not allowed. Read the usage summary below.
      	    echo
      	    usage; exit 1 ;;
    esac
done


###########################################
# Set up access to all scripts for pipeline
###########################################
# This script accesses two more scripts in the same directory.
# Unless users set the path to this directory they will not be acessible.
# So, pre-fixing the extra scripts with the path of this script.
# NB - I think $0 is always the full path. 
fullPathToWrapperScript=`echo $0 `
#echo  $fullPathToWrapperScript
pathToScripts=`dirname $fullPathToWrapperScript `


#############################
# User input checks and setup
#############################

# Prints usage if no parameters are given:ß
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
#echo 'Values of all free parameters: '${@:$OPTIND:$#}	# Therefore ${@:$OPTIND:$#} will access all the free parameters
echo $(( $# - $OPTIND + 1 ))							# Number of free parameters, in this case the number of samples
numbrSamples=$(( $# - $OPTIND + 1 ))
#echo $numbrSamples


if [ $(( $# - $OPTIND + 1 )) -lt 4 ]; then				### 30.3.2020 AND speciesTreesOnly == no to handle scritp just processing species trees
    echo
    usage
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
if [[ $addSampleName == 'yes' ]]; then

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
    done

    # Concatenate fasta files and put seqs on a single line:
	fileList=`ls *_modified.fasta `
	echo $fileList
	cat $fileList | seqtk seq -l 0 /dev/fd/0 > all_samples_concatenated.fasta
	geneFile=all_samples_concatenated.fasta
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
     	if [[ $beforeDashCheck -gt 1 ]]; then echo "ERROR: more than one sample identifier detected on fasta header line (there should only be one) for this sample: $file \nAlso check that fasta header lines have this format: >sampleId-geneId"; exit 1
     	elif [[ $afterDashCheck -gt 1 ]]; then echo "ERROR: gene identfiers should be unique, one or more not unique for this sample: $file \nAlso check that fasta header lines have this format: >sampleId-geneId"; exit 1
		else cat $file | awk '{if($1 ~ /^>/)  {print $1} else {print $0}}' \
			| awk -F '-' '{if($1 ~ /^>/) {{gsub(/>/,"",$1)} {print ">" $2 " " $1}} else {print $0}}' \
			> ${uniqueSampleId}_modified.fasta
			# NB - fasta format now: >geneId sampleId
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
if [[ $sampleTableFile != 'no' ]]; then
	if [[ -s $sampleTableFile ]]; then echo ""
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


### 8.5.2020 - DON"T NEED THIS CODE ANY MORE - DELETE
		# for file in ${@:$OPTIND:$#}; do

		# 	# Get filename and chop off the file ending if there is one:
  #   		uniqueSampleId=`basename $file | awk -F '.' '{print $1}' `
  #   		uniqueSampleId=`expr $uniqueSampleId + 0`
  #   		#echo $uniqueSampleId

		# 	$pathToScripts/add_sample_id_from_table.py  $sampleTableFile  $file  $uniqueSampleId  ${uniqueSampleId}_modified.fasta
		# 	### NB - 31.10.2019 - just realised that this for loop should be done within the python script so that the Dict only 
		# 	### has to be set up once; may improve the speed of this step.
		# 	### NB - I could also get the uniqueSampleId in above script after stripping out any zeros - now doing
		# 	### 13.10.2019 - NB - doing this means that if there are > one files with same identifier but differ by leading zeros,
		# 	### these files will overwrite each other and only one will be produced --> avoids duplicate samples going into the analysis.
		# 	### This is convenient for now but consider to be completely strict about sample names, if filenames are not as in the csv file,
		# 	### then they are not found!! However it would seem very silly to have different ids only by the presence of leading zeros!
		# done
  #   	# Concatenate fasta files and put seqs on a single line:
		# fileList=`ls *_modified.fasta `
		# #echo $fileList
		# ### Not required now - seqs coming out of python are already on a single line (I think)
		# cat $fileList | seqtk seq -l 0 /dev/fd/0 > all_samples_concatenated.fasta
		# geneFile=all_samples_concatenated.fasta
	else
    	echo File containing text for the tree leaves does not exist or is empty: $sampleTableFile - exiting.
		exit 1
	fi
fi


if [ -s $geneListFile ]; then echo ""
	### NB - 26.4.2020 - just realised that the gene names should not contain any dot chars - see note in make_gene_trees.sh ~ line 59
else echo "ERROR: the gene list file (option -g) does not exist or is empty: $geneListFile"; exit; fi


# NB - Trying to check that value is an integer - how does this behave when adding 'hello'? look OK; BUT allows through chars A -F - why?
if (( $(bc <<< "$fractnAlnCovrg < 0 ") || $(bc <<< "$fractnAlnCovrg > 1 ") )); then
	echo "ERROR1: -f option needs to be a fraction between 0 and 1 - exiting"
	exit
fi  


if (( $(bc <<< "$fractnSamples < 0") || $(bc <<< "$fractnSamples > 1") )); then 
	echo "ERROR: -s option needs to be a fraction between 0 and 1 - exiting"
	exit
fi


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

if [ -z $fileNamePrefix ]; then echo ""
	echo "ERROR: option -p used but an output file name prefix string was not supplied - exiting"
	exit
	### This doesn't really work because if there is no paramter value the next parameter becomes the value which is a string!
fi

#if [[ -z $phyloProgram || $phyloProgram != 'fasttree' || $phyloProgram != 'raxml-ng' ]]; then echo ""
#	echo "ERROR: option -q used but a phylogeny program name was not supplied or was incorrectly spelt - exiting"
	### 16.3.2020 - this really doesn't work - is triggered no matter what! Now have a solution in make_gene_trees.sh
	### This doesn't really work because if there is no paramter value the next parameter becomes the value which is a string!
#fi

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
		"$exePrefix"
	done > make_gene_tree.log 2>&1
	$pathToScripts/assess_gene_alignments.sh $fractnAlnCovrg $fractnMaxColOcc $fractnSamples $numbrSamples $fileNamePrefix $geneFile tree_tip_info_mapfile.txt $option_u > assess_gene_alns.log 2>&1
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

        jobInfo=`sbatch -p $partitionName -c $cpuGeneTree -t 0-36:00 --mem 10000 --array=0-${numbrGenes}%$slurmThrottle  $pathToScripts/slurm_setup_array_to_make_gene_trees.sh \
		$geneFile \
		$geneListFile \
		$fractnAlnCovrg \
		$pathToScripts \
		$phyloProgramDNA \
		$phyloProgramPROT \
		$fractnMaxColOcc \
		$cpuGeneTree \
		"$mafftAlgorithm" \
		"$exePrefix" `
		echo jobInfo: $jobInfo
		jobId=`echo $jobInfo | cut -d ' ' -f 4 `
		echo \$jobId: $jobId - same id as \$SLURM_ARRAY_JOB_ID 
		jobInfo1=`sbatch -J assess_gene_alns  --dependency=afterok:$jobId -p $partitionName -c 1 -n 1 -o assess_gene_alns.log -e assess_gene_alns.err  $pathToScripts/assess_gene_alignments.sh $fractnAlnCovrg $fractnMaxColOcc $fractnSamples $numbrSamples $fileNamePrefix $geneFile $option_u `
    	echo jobInfo1: $jobInfo1
		jobId1=`echo $jobInfo1 | cut -d ' ' -f 4 `
		echo \$jobId1: $jobId1 - from running assess_gene_alignments.sh
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
			$cpuGeneTrees \
			"$mafftAlgorithm" \
			"$exePrefix"
		done > make_gene_tree.log 2>&1
		$pathToScripts/assess_gene_alignments.sh $fractnAlnCovrg $fractnMaxColOcc $fractnSamples $numbrSamples $fileNamePrefix $geneFile tree_tip_info_mapfile.txt  $option_u
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
	> make_species_trees.log 2>&1
	### 31.8.2019 - I think I can move all intermediate and gene tree files here
	### Actually, this is a lot of files - just delay this and think if I need to do it!!!!
	### I think a better way is to copy all input files to the tmp dir, then cd into it and run everythign there - then put all the final files in the starting dir.
	# mv *_modified.fasta  $tempOutputsDir
	# mv all_samples_concatenated.fasta  $tempOutputsDir
	# g${gene}_dna.fasta
	# g${gene}_protein.fasta
	# g4471_mafft_dna_aln.fasta
	# g4471_mafft_dna_aln_ovr70pc_aln_covrg.txt
	# g4471_mafft_protein_aln.fastag
	# g4471_mafft_protein_aln_ovr70pc_aln_covrg.fasta
	# g4471_mafft_dna_aln_ovr70pc_aln_covrg_fasttree.nwk
	# g4471_mafft_protein_aln_ovr70pc_aln_covrg_fasttree.nwk
	# Malpighiales_fasttree_trees_for_coelescence_phylo.nwk
	# Malpighiales_fasttree_protein_trees_for_coelescence_phylo.nwk
	# mafft_dna_alns_fasta_file_list.txt
	# Malpighiales_fasttree_trees_for_coelescence_phylo_bs_less15pc_rmed.nwk
	# Malpighiales_fasttree_trees_for_coelescence_phylo_astral.nwk
	# mafft_dna_alns_fasta_samples_in_tree.txt
	# samples_submitted.txt
	# RAxML_bootstrap.Malpighiales__raxmlHPC-PTHREADS-SSE
	# RAxML_bestTree.Malpighiales__raxmlHPC-PTHREADS-SSE
	# RAxML_bipartitions.Malpighiales__raxmlHPC-PTHREADS-SSE
	# RAxML_bipartitionsBranchLabels.Malpighiales__raxmlHPC-PTHREADS-SSE


elif [ $os == 'Linux' ]; then
	exePrefix="/usr/bin/time -v"
    if [ $slurm -eq 1 ]; then
    	### NB - not sure where to put the $exePrefix!!!!
    	### One option is to put the "time script" cmd in a wrapper but then I need a log file for this sbatch call and delete it from the script header..
    	### I think this is the only way without changing the main script itself.
													
        sbatch --dependency=afterok:$jobId1 -p long -c $cpu --mem=$speciesTreeMem  $pathToScripts/make_species_trees.sh $fractnAlnCovrg $fractnSamples $numbrSamples $fileNamePrefix $geneFile $cpu $phyloProgramDNA $phyloProgramPROT "$exePrefix" $sampleTableFile
	else
		$pathToScripts/assess_gene_alignments.sh \
		$fractnAlnCovrg \
		$fractnMaxColOcc \
		$fractnSamples \
		$numbrSamples \
		$fileNamePrefix \
		$geneFile \
		tree_tip_info_mapfile.txt \
		$option_u \
		> assess_gene_alns.log 2>&1
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
		> make_species_trees.log 2>&1
		### 31.8.2019 - I think I can move all intermediate and gene tree files here 
	fi
fi
