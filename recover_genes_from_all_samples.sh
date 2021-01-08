#!/bin/bash

###############################
# recover_genes_from_samples.sh

# Based on 51.PAFTOL_gene_recovery_overlapRecover.sh 

# Author:   Paul Bailey
###############################
set -e
set -u
set -o pipefail
shopt -s failglob


####################################################################
# Set up commandline parameters and fetch their values using getopts
####################################################################

# Variables for any command line flags needing a default value:
hybSeqProgram=no            # Removed the default for now: paftools
targetsFile='no'
adapterFasta='no'
samplePrefix=Sample
cpu=4
slurmTime=0-36:00           # SBATCH -t 0-36:00 - was 24h but increased to 36 for the larger samples, then 3 days necessary for just a few PAFTOL samples 
slurmMemory=0               # Using 80,000 MB for most samples and 130000 MB for very large samples - could use 130000 MB on all samples; only 0.3GB required for HybPiper it seems, use 20GB
partitionName=medium        # Values depend on the cluster being used so good to have a flagged option for this
slurmThrottle=1
stats=no                    # stats flag for reporting read coverage and depth
refFilePathForStats=default # gene recovery fasta file path for coverage stats - used for the BWA reference

# Hidden options (i.e. not apparent from the help menu but they always have a value so can be used in downstream scripts):
usePaftolDb=no              # d) usePaftolDb='--usePaftolDb' ;;
spadesCovCutoff=4           # Default=8, 4 might be OK too (currently set to 4); for reads_first.py --cov_cutoff flag - NOT IMPLEMENTED YET





function usage  {

cat << EOF

Program: recover_genes_from_all_samples.sh

Program description: recovers genes from fastq files of multiple samples. Paftools or HybPiper caqn be used to align the reads to a set of
                     target genes, then assemble the reads for each target gene. If Slurm is available samples will be run in parallel

Usage: recover_genes_from_all_samples.sh

OPTIONS <value>:
  -h   
                 prints usage and description
  -v             
  program version
  -s <csv file>  
                 add sample name and fastq file names via a csv table file (must have a header line);
                 format: SampleName,R1FastqName,R2FastqName (required option)
  -f <string>    
                 FULL path to all sample fastq files; N.B. no filenames, just the full path to them, not a relative path and no wild cards! (required option)
  -t <string>    
                 file name of target genes in fasta format (required option)
  -a <string>    
                 file name of adaptors in fasta format (required option)
  -y <string>    
                 Hyb-Seq program; options are: paftools, hybpiper, hybpiper-bwa (required option)
  -S    
                 calculate statistics for gene recovery from Illumina read data (includes reads on-target, read coverage, read depth).
                 This option can be used separately after the gene recoveries have run but the path to the gene recovery fasta 
                 files has to be specified with option -P, if not running in the same location as the original gene recovery run
  -P <string> 
                 Specify FULL path to the gene recovery fasta files (for option -S), but only the part common to all files.
                 This option looks for one of these two scenarios:
                 1. For gene recovery samples in each of their own directory i.e. /<path>/<SampleDirPrefix>_<SampleName>/<SampleName>.fasta, type: /<path> (directory set up, as used by this pipeline) 
                 2. For gene recovery samples all in the same directory i.e. /<path>/<SampleName>.fasta, type: /<path>
                 Note: 'SampleName' needs to match that provided by the sample list in option -s
  -p <string>    
                 directory prefix for each sample (default=Sample)
  -c <integer>   
                 number of cpu to use (default=4)
  -m <integer>   
                 Slurm memory to use (in MB); Paftools requires >> 20000, HybPiper requires << 20000 (default=0, means no limit is imposed in Slurm default mode)
  -T <string>    
                 Slurm time limit, format <days>-<hours>:<minutes> (default=0-36:00 (36 hours)
  -Q <string>    
                 Slurm partition (queue) to use (default=medium)

  -H <integer>
                 Slurm array throttle (default=10; could set to 1, then increase once happy with run with: scontrol update arraytaskthrottle=<integer> job=<jobId>)
 

A typical example:
<path to>/recover_genes_from_all_samples.sh \\
-s <table_file.csv> \\
-t <angiosperms353TargetsFile.fasta> \\
-f <fastq_files_path> \\
-a <illumina_adaptors.fasta> \\
-p Sample \\
-c 4 \\
-m 20000 \\
-Q main \\
> recover_genes_from_all_samples.log 2>&1 &

EOF
}


#echo User inputs:    ### For testing only 
while getopts "hvs:t:f:a:y:p:c:d:H:m:T:Q:SP:"  OPTION; do
 
  #echo -$OPTION $OPTARG    ### For testing only - could try to run through options below 
   
  case $OPTION in

    h) usage; exit 1 ;;
    v) echo "recover_genes_from_all_samples.sh  version 0.0.1"; exit ;;
    s) sampleList=$OPTARG ;;
    t) targetsFile=$OPTARG ;;
    f) paftolDataSymlinksDir=$OPTARG ;;
    a) adapterFasta=$OPTARG ;;
    y) hybSeqProgram=$OPTARG ;;
    p) samplePrefix=$OPTARG ;;
    c) cpu=$OPTARG ;;
    d) usePaftolDb=$OPTARG ;;
    H) slurmThrottle=$OPTARG ;;
    m) slurmMemory=$OPTARG ;;
    T) slurmTime=$OPTARG ;;
    Q) partitionName=$OPTARG ;; 
    S) stats=yes ;;
    P) refFilePathForStats=$OPTARG ;;
    ?)  echo This option does not exist. Read the usage summary below.
            echo
            usage; exit 1 ;;
    esac
done


#############################
# User input checks and setup
#############################
# Prints usage if no parameters are given:
#echo $#
if [ "$#" -lt 1 ]; then usage; exit 1; fi


###################################### 
# Check required software dependancies
######################################
# Only two programs to check currently (Paftools and HybPiper; can use -h and --help, both give a zero exit code)
### 29.3.2020 - addpaftolfastqmaybe do this after all other checks so the print out doesn't appear if there is an issue with input parameters
### 12.5.2020 - also Trimmomatic is required. Others e.g. fastqc are dependancies of Paftools
if [[ $hybSeqProgram == 'hybpiper'* ]]; then
    echo 'Testing HybPiper is installed, will exit with a 127 error if not found.' 
    reads_first.py -h >/dev/null 2>&1    # Will exit here if not found!
    echo 'hybpiper found.'

    #echo 'Testing seqtk is installed...'
    #seqtk >/dev/null 2>&1                # Will exit here if not found! NBNB - also exits here IF FOUND !!!!!!!! Need another solution!
                                          # Maybe run program via the $hybSeqProgram as this works in the trees script - maybe bash is unaware that it is runnign a program!?
                                          # NO - I thinks it's becuase of the set options - maybe should remove them !!!!
elif [[ $hybSeqProgram == 'paftools' ]]; then
    echo 'Testing paftools is installed, will exit with a 127 error if not found.' 
    $hybSeqProgram -h >/dev/null 2>&1   # Will exit here if not found!
    echo 'paftools found.'
fi
if [[ $? == 127 ]]; then
    # Exit code 127 is for "command not found"
    ### 16.3.2020 - this doesn't work as script exits in the above conditional if program is not found!
    ### Readjusted above conditional instead.
    echo "ERROR: Not available: $hybSeqProgram"
    exit
fi


### Check fasta file format here - see species tree script
### For Hybpiper need to ensure there are no period chars (.) in the sample name/ids - for the paralog code


if [ ! -s $sampleList ]; then usage; echo; echo "ERROR: the samples table file (option -s) does not exist or is empty: $sampleList"; exit; fi

if [[ ! -s "$targetsFile" && $hybSeqProgram != 'no' ]]; then usage; echo; echo "ERROR: the target genes file (option -t) does not exist or is empty: $targetsFile"; exit; fi
### 12.5.2020 - Just realsied that I can determine the full path to fiel here then the user just needs to supply relative path - ditto for paftoldataasymlinkdir
###             NB - I have no test for the fastq dir exisitng! NB - I think this harder than i realised - how do you get the full path from a partial path - 8.1.2020 - I think I know this now - see species tree script 

if [[ ! -s "$adapterFasta" && $hybSeqProgram != 'no' ]]; then usage; echo; echo "ERROR: the adaptor file (option -a) does not exist or is empty: $adapterFasta"; exit; fi

if [ -z $samplePrefix ]; then echo ""
  usage; echo "ERROR: option -p used but an output file name prefix string was not supplied - exiting"
  ### This doesn't really work because if there is no paramter value the next parameter becomes the value which is a string!
fi

# Check first for whether $cpu is an integer (this is a logical fudge for detecting whether an integer or string) 
if [ "$cpu" -eq "$cpu" ] 2>/dev/null
then echo ""
else usage; echo "ERROR: option -c should be an integer - exiting "; exit; fi  


if [[ $usePaftolDb != 'PAFTOL' && $usePaftolDb != 'OneKP_Transcripts' && $usePaftolDb != 'OneKP_Reads' && $usePaftolDb != 'SRA' && $usePaftolDb != 'AG' && $usePaftolDb != 'no' ]]; then 
  usage; echo "ERROR: option -d should contain one of the following data sets: PAFTOL, OneKP_Transcripts, OneKP_Reads, SRA or AG - you added \'$usePaftolDb\'. Exiting"
  exit
fi    

##########################
# End of user input checks
##########################


echo 'Detecting operating system: '
os=`uname `
echo $os
### 25.9.2019 - I think these lines could be merged with essentailly same lines below.
if [ $os == 'Darwin' ]; then
    echo 'Running as if on Mac OS (Darwin)...'
elif [ $os == 'Linux' ]; then
  echo 'Running as if on Linux...'
    slurm=`sbatch -V | grep ^slurm | wc -l `
    if [ $slurm -eq 1 ]; then 
        echo 'Slurm scheduler detected, running via slurm...'
        ### 28.6.2019 - could create a variable here that contains call to "lsrun -J <$jobName>"
        ### remain empty if slurm not used - need to be able to set $jobName for each program.
  fi
else
    echo 'Unrecognised OS that this software can run on, exiting.'
    exit
fi
#echo slurm: $slurm
#exit


###########################################
# Set up access to all scripts for pipeline
###########################################
# This script accesses two more scripts in the same directory.
# Unless users set the path to this directory they will not be acessible.
# So, pre-fixing the extra scripts with the path of this script.
# NB - I think $0 is always the full path. 
fullPathToWrapperScript=`echo $0 `
echo  $fullPathToWrapperScript
pathToScripts=`dirname $fullPathToWrapperScript `



###########################
echo 'Processing samples for gene recovery...'
###########################
if [ $os == 'Darwin' ]; then
	exePrefix="/usr/bin/time -l "		# this time command gets the RSS memory - -l flag doesn't work on Cluster - but I hope the cluster version is common to all GNU/Linux. 
	tail -n+2 $sampleList | \
	while read line; do
  	### Keep an eye on whether the $line variable value can break up with any chars
    echo $line
  	$pathToScripts/recover_genes_from_one_sample_test_with_stats.sh "$line"  $targetsFile  $paftolDataSymlinksDir  $adapterFasta  $samplePrefix  $cpu  "$exePrefix" $hybSeqProgram $usePaftolDb $stats $refFilePathForStats 
  done
elif [ $os == 'Linux' ]; then
  exePrefix="/usr/bin/time -v"	# PYTHONPATH works on the Cluster, but on Macbook, it deosn't need to be set (only really need to alter the flag char!)
  slurm=`sbatch -V | grep ^slurm | wc -l `
  if [ $slurm -eq 1 ]; then

    # Count the # samples to process and fix that number in the Slurm --array parameter.\
    # It has to be set outside the sbatch script I think so that I can automatically set the array size .
    # It has to be set each time otherwise.
    numbrSamples=`tail -n+2 $sampleList | wc -l `

    # Notes on the --array flag - e.g. SBATCH --array=1-${numbrSamples}%50 
    # 0-352%10	# NB - if input file has a header line, array should start at 1;
		# 25 samples with 8cpu (80GB mem) is ~ 50% available cpu, 2TB mem (67% of available mem)
    # Maximum array size by default is 1000 - can be increased with MaxArraySize (max supported = 4000001, default set to 1001)
    # NBNB - it's only an administration value set globally i think  - can't be changed by user
		# Prevent jobs going over 1000, otherwise will get empty sample directories and log files - minor issue!:
    if [ $numbrSamples -ge 1001 ]; then
      numbrSamples=1000
      echo "WARNING: Maximum default size of Slurm array is 1000, processing only the first 1000 samples"
    elif [ $numbrSamples -lt $slurmThrottle ]; then 
      slurmThrottle=$numbrSamples
    fi

    jobInfo=`sbatch -p $partitionName -c $cpu -t $slurmTime  --mem $slurmMemory  --array=1-${numbrSamples}%$slurmThrottle  $pathToScripts/slurm_setup_array_to_recover_genes.sh \
    $sampleList \
    $targetsFile \
    $paftolDataSymlinksDir \
    $adapterFasta \
    $samplePrefix \
    $cpu \
    $pathToScripts \
    "$exePrefix" \
    $hybSeqProgram \
    $usePaftolDb \
    $stats \
    $refFilePathForStats `
   	echo jobInfo: $jobInfo			# NB - Don’t need to remember the jobId - unless want to merge with tree pipeline
    jobId=`echo $jobInfo | cut -d ' ' -f 4 `
    echo \$jobId: $jobId - same id as \$SLURM_ARRAY_JOB_ID 
  else
  	exePrefix="/usr/bin/time -v"
    tail -n+2 $sampleList | \
		while read line; do
      ### Keep an eye on whether this variable value can break up with any chars
  		$pathToScripts/recover_genes_from_one_sample.sh "$line"  $targetsFile  $paftolDataSymlinksDir  $adapterFasta  $samplePrefix  $cpu  "$exePrefix" $hybSeqProgram $usePaftolDb $stats $refFilePathForStats 
  	done
  fi
fi
