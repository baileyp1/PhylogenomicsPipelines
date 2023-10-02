#!/bin/bash

##################################
# recover_genes_from_one_sample.sh

# Author: Paul Bailey

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################
set -e
set -u
set -o pipefail
shopt -s failglob


line="$1"
targetsFile=$2
paftolDataSymlinksDir=$3
adapterFasta=$4
samplePrefix=$5
cpu=$6
exePrefix="$7"
hybSeqProgram=$8
usePaftolDb=$9
stats=${10}
refFilePathForStats=${11} 

#echo exePrefix: "$exePrefix"	# check that variable is not split on space
echo hybSeqProgram: $hybSeqProgram
echo stats: $stats
echo refFilePathForStats: $refFilePathForStats
  

sampleId=`echo $line | cut -d ',' -f 1 `    		# I think this will be the PAFTOL_xxxxxx id - same as in the read file names
R1FastqFile=`echo $line | cut -d ',' -f 2 `
R2FastqFile=`echo $line | cut -d ',' -f 3 `
externalSequenceID=`echo $line | cut -d ',' -f 4 `	# For adding the external sequence Id to the paftol_da db

  
if [[ ! -d ${samplePrefix}_$sampleId ]]; then mkdir ${samplePrefix}_$sampleId; fi
cd ${samplePrefix}_$sampleId
echo sampleId: $sampleId
echo externalSequenceID: $externalSequenceID


if [[ -s  $paftolDataSymlinksDir/$R1FastqFile ]]; then ls $paftolDataSymlinksDir/$R1FastqFile; fi
if [ -z "$R2FastqFile" ]; then 
	echo "No R2FastqFile - read trimming and gene recovery will work in single-end mode"
else
	ls $paftolDataSymlinksDir/$R1FastqFile
	ls $paftolDataSymlinksDir/$R2FastqFile
fi
pwd

#if [ $usePaftolDb != 'usePaftolDb' ]; then	- changed - now introducing data set type
if [[ $usePaftolDb == 'no' ]]; then
  			                                        #--nodelist=kppgenomics01.ad.kew.org  # mem normally set to 80000
	#  sbatch -J ${samplePrefix}_${sampleId}_fastqToGenes -p main -t 1-0:00 -c $cpu --mem=80000 -o ${samplePrefix}_${sampleId}_fastqToGenes.log   -e ${samplePrefix}_${sampleId}_fastqToGenes.log_err   --wrap "
	# RUNTIME: For 8 cpu, up to 2 mins; up to 18 GB mem (for 10 samples)
	#          For 4 cpu, up to 24 mins; up to 12 GB mem.    - Total # samples to date: ~2500 - would take 14h using 176 cpu (1 node) 
	#          For 2 cpu, up to 1h50', up to 10 GB mem (tested 40 samples)
	#          Conclusion: the above stats seems to suggest that running this step separately might be the most efficient
	###pathToTrimmomatic=`which trimmomatic-0.39.jar `		# NB - 'which' requires the file to be executable!
	###$exePrefix  java -jar $pathToTrimmomatic PE \		# 12.2.2021 - Changed the way java programs are called to using a global variable
	if [ -z "$R2FastqFile" ]; then
		### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix  java -jar $TRIMMOMATIC SE - OK now so have put back
		$exePrefix java -jar $TRIMMOMATIC SE \
		-threads $cpu \
		-trimlog ${sampleId}_R1_trimmomatic.log \
		$paftolDataSymlinksDir/$R1FastqFile \
		${sampleId}_R1_trimmomatic.fastq.gz \
		ILLUMINACLIP:${adapterFasta}:2:30:10:2:true \
		LEADING:10 \
		TRAILING:10 \
		SLIDINGWINDOW:4:20 \
		MINLEN:40 > ${sampleId}_trimmomatic.log 2>&1
	else
		### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix  java -jar $TRIMMOMATIC PE - OK now so have put back
		$exePrefix java -jar $TRIMMOMATIC PE \
		-threads $cpu \
		-trimlog ${sampleId}_R1_R2_trimmomatic.log \
		$paftolDataSymlinksDir/$R1FastqFile \
		$paftolDataSymlinksDir/$R2FastqFile \
		${sampleId}_R1_trimmomatic.fastq.gz \
		${sampleId}_R1_trimmomatic_unpaired.fastq.gz \
		${sampleId}_R2_trimmomatic.fastq.gz \
		${sampleId}_R2_trimmomatic_unpaired.fastq.gz \
		ILLUMINACLIP:${adapterFasta}:2:30:10:2:true \
		LEADING:10 \
		TRAILING:10 \
		SLIDINGWINDOW:4:20 \
		MINLEN:40 > ${sampleId}_trimmomatic.log 2>&1
	fi
	# 8.3.2020 - changed to palidromic mode - confirm that the extra two parameters is Ok still for normal mode
  

	# Paftools requires unzipped fastq files (NB - if unzipped files already present, files will not be unzipped again (i.e gzip exits with an error), so I've included the -f flag to force uncompression!):
	#srun -J ${sampleId}_unzip_R1_R2 -n 1  gunzip -f ${sampleId}_R1_trimmomatic.fq.gz  ${sampleId}_R2_trimmomatic.fq.gz
	if [[ -s ${sampleId}_R1_trimmomatic.fastq.gz ]]; then gunzip -f ${sampleId}_R1_trimmomatic.fastq.gz; fi  
	if [[ -s ${sampleId}_R2_trimmomatic.fastq.gz ]]; then gunzip -f ${sampleId}_R2_trimmomatic.fastq.gz; fi
fi


### NB - I think that if I check file exists or is zero byte before running this is equivalent to make affect
### convert script to a makefile


usePaftolDbFlag=''
if [ $hybSeqProgram == 'paftools' ]; then

	echo "Running gene recovery with Paftools..."

	if [ $usePaftolDb != 'no' ]; then 
		usePaftolDbFlag='--usePaftolDb'		### NB - 6.7.2021 - I don't think $usePaftolDbFlag needs to be used in the else part of the conditional - it is a completely different command - then I can remove usePaftolDbFlag='' just above.

  		# Separate out the two fields:
  		datasetOrigin=`echo $usePaftolDb | cut -d ':' -f 1 `
  		recoveryRun=`echo $usePaftolDb | cut -d ':' -f 2 `
  		echo "Dataset type: $datasetOrigin"		# Required for addPaftolFastq program
  		echo "Recovery run: $recoveryRun" 		# Required for recoverSeqs program 

### 14.7.2023 - will not need these steps for option -u upload to new paftol db 
		# Remove the .gz ending from the file name for adding to the pafto_da database
		unzippedR1FastqFile=`basename -s .gz $paftolDataSymlinksDir/$R1FastqFile `

		# Need to uncompress raw fastq files to get FastQCStats and upload into the paftol_da database:
		gunzip -f -c $paftolDataSymlinksDir/$R1FastqFile > $unzippedR1FastqFile
		if [ -z "$R2FastqFile" ]; then
			echo "No R2FastqFile"
		else
			unzippedR2FastqFile=`basename -s .gz $paftolDataSymlinksDir/$R2FastqFile `
			gunzip -f -c $paftolDataSymlinksDir/$R2FastqFile > $unzippedR2FastqFile
		fi
		
		
 		# Add sampleId flag and value to paftools command, depending on the data type, PAFTOL or external data.
 		# NB - PAFTOL data doesn't use this flag value but the flag is still required! Paftools checks whether 
 		#      an externalSequenceID exists and fails if one doesn't for a non-paftol sample.
		if [[ -n $externalSequenceID ]]; then
			externalSequenceID="--sampleId $externalSequenceID"
		fi
		# NB --sampleId flag is not now mandatory and paftools checks whether an externalSequenceID exists
		#    and fails if one doesn't for a non-paftol sample.
		export PYTHONPATH=$HOME/lib/python
		if [ -z "$R2FastqFile" ]; then
			paftools --loglevel INFO addPaftolFastq $externalSequenceID \
			--fastqPath $paftolDataSymlinksDir \
			--dataOrigin $datasetOrigin \
			$unzippedR1FastqFile
		else
			paftools --loglevel INFO addPaftolFastq $externalSequenceID \
			--fastqPath $paftolDataSymlinksDir \
			--dataOrigin $datasetOrigin \
			$unzippedR1FastqFile  $unzippedR2FastqFile
		fi
		# NB - the file name must be of this format e.g. PAFTOL_005853_R1.fastq BUT now only for PAFTOL data
		# The fastqPath entered in the database consists of the path and filename e.g. $paftolDataSymlinksDir/$unzippedR1FastqFile
		# --sampleId=$sampleId - must be included for all data types (there must be a value), however in the case of paftol data it is ignored but the flag and a value is still required!
		### 7.6.2021 - NB - I don't think this is true anymore! I think --sampleId and therefore $externalSequenceID can be blank
		# --dataOrigin - dataset origin added via the -d flag of this script/pipeline e.g. -d PAFTOL, -d SRA

		echo "Exit status of paftools addpaftolFastq:" $?
		### NB - paftools recoverSeqs also runs FastQC stats against the raw reads but this is unnessary.
		###     Try to delete this step then I can just remove the raw fastqs here.
### 14.7.2023 - end of requirement for option -d

		# Paftools --usePaftolDb requires the targets filename WITHOUT the path
		# i.e. it insists that targets file is in the pwd (but I could change that):
		cp $targetsFile . 
		targetsFile=`basename $targetsFile `

		# NB - Use of --usePaftolDb flag requires the use of the Paftools trimmomatic flags so have to have a separate paftools recoverSeqs command here.
		#      The Trimmomatic program name needs to be 
		export PYTHONPATH=$HOME/lib/python 			# I had to add this for the cluster ONLY - need to. Check it is OK on Macbook, it should be.
		if [ -z "$R2FastqFile" ]; then
			### 13.10.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix paftools recoverSeqs - OK now so have put back
			$exePrefix paftools recoverSeqs \
			$targetsFile \
			${sampleId}.fasta \
			-f $unzippedR1FastqFile \
			--trimmer trimmomatic \
			--trimmomaticLeadingQuality 10 --trimmomaticTrailingQuality 10 \
			--trimmomaticMinLength 40 \
			--trimmomaticSlidingWindowSize 4 --trimmomaticSlidingWindowQuality 20 \
			--trimmomaticAdapterFname $adapterFasta  \
			--mapper tblastn \
			--assembler overlapSerial \
			--blastNumThreads $cpu \
			--allowInvalidBases \
			--windowSizeReference 50 \
			--relIdentityThresholdReference 0.7 \
			--windowSizeReadOverlap 30 \
			--relIdentityThresholdReadOverlap 0.9 \
			--summaryCsv ${sampleId}_summary.csv \
			$usePaftolDbFlag $recoveryRun \
			> ${sampleId}_overlapSerial.log 2>&1
		else
			### 13.10.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix  paftools recoverSeqs - OK now so have put back
			$exePrefix paftools recoverSeqs \
			$targetsFile \
			${sampleId}.fasta \
			-f $unzippedR1FastqFile \
			-r $unzippedR2FastqFile \
			--trimmer trimmomatic \
			--trimmomaticLeadingQuality 10 --trimmomaticTrailingQuality 10 \
			--trimmomaticMinLength 40 \
			--trimmomaticSlidingWindowSize 4 --trimmomaticSlidingWindowQuality 20 \
			--trimmomaticAdapterFname $adapterFasta  \
			--mapper tblastn \
			--assembler overlapSerial \
			--blastNumThreads $cpu \
			--allowInvalidBases \
			--windowSizeReference 50 \
			--relIdentityThresholdReference 0.7 \
			--windowSizeReadOverlap 30 \
			--relIdentityThresholdReadOverlap 0.9 \
			--summaryCsv ${sampleId}_summary.csv \
			$usePaftolDbFlag $recoveryRun \
			> ${sampleId}_overlapSerial.log 2>&1
		fi

### 14.7.2023 - write to the db here - e.g.:
### if option - u is set : 
### note new flag -u for script here - Python $ENV var instead of creating a file
### make this into a function so i can call it with HybPipr as well
#time paftol_db_add_sample_recoveries_batch_upload.py \
#-c ~/.paftol/paftol_merge_v4_for_recovery_upload_tests.cfg \
#-p Sample \
#-a . \
#-b paftol/AllData_symlinks \
#-f . \
#-g paftol/PAFTOL_gene_recovery_and_trees/run1_paftools \
#-r 1 \
#-u sampleId \
#-C 	#### May want to have this as optional for this parent script as well
#> paftol_da_add_samples_batch_upload.log 2>&1 &


		rm $targetsFile	# If write to database fails, this file doesn't get deleted (c.f. set cmds active), so presence of file is a useful 'marker' for failing to write to db

		# Remove the large fastq files::
		if [[ -s $unzippedR1FastqFile ]]; then rm $unzippedR1FastqFile; fi
		if [[ -s $unzippedR2FastqFile ]]; then rm $unzippedR2FastqFile; fi
	else
		################
		# overlapRecover (OLC approach)
		################
		# RUNTIME: for 4 cpu, up to 16h21' + up to >>23GB mem (average mem = 14GB assessed over 142 samples; a few samples use up to 125GB mem in 8-11h); set to 80GB mem (or lower but will get more fails)
		# Was using:
		#srun -J ${sampleId}_overlapRecover -n 1  -o ${sampleId}_overlapRecover.log  -e ${sampleId}_overlapRecover.log_err \	- issue with srun
		export PYTHONPATH=$HOME/lib/python 			# I had to add this for the cluster ONLY - need to. check it is OK on Macbook, it should be.
		if [ -z "$R2FastqFile" ]; then
			### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix  paftools recoverSeqs - OK now so have put back
			$exePrefix paftools recoverSeqs \
			$targetsFile \
			${sampleId}.fasta \
			-f ${sampleId}_R1_trimmomatic.fastq \
			--mapper tblastn \
			--assembler overlapSerial \
			--blastNumThreads $cpu \
			--allowInvalidBases \
			--windowSizeReference 50 \
			--relIdentityThresholdReference 0.7 \
			--windowSizeReadOverlap 30 \
			--relIdentityThresholdReadOverlap 0.9 \
			--summaryCsv ${sampleId}_summary.csv \
			$usePaftolDbFlag \
			> ${sampleId}_overlapSerial.log 2>&1
		else
			### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix  paftools recoverSeqs - OK now so have put back
			$exePrefix paftools recoverSeqs \
			$targetsFile \
			${sampleId}.fasta \
			-f ${sampleId}_R1_trimmomatic.fastq \
			-r ${sampleId}_R2_trimmomatic.fastq \
			--mapper tblastn \
			--assembler overlapSerial \
			--blastNumThreads $cpu \
			--allowInvalidBases \
			--windowSizeReference 50 \
			--relIdentityThresholdReference 0.7 \
			--windowSizeReadOverlap 30 \
			--relIdentityThresholdReadOverlap 0.9 \
			--summaryCsv ${sampleId}_summary.csv \
			$usePaftolDbFlag \
			> ${sampleId}_overlapSerial.log 2>&1
		fi
		
		# Note:
		# Paftools overlapRecover and recoverSeqs programs (4.9.2019) both output fasta header geneId, not yet sampleId-geneId 
		# --summaryCsv ${sampleId}_summary.csv - requires that the fastq file ends in .fasta not e.g. .fq
		# --tgz ${sampleId}_overlapSerial.tgz  - removed to save disk space
		# recoverSeqs program also requires the --mapper and --assembler options - overlapRecover does not.
		# I think windowSizeReference should be similar to the average size of the trimmed read length (or the N50)
		# NB - 28.2.2020 - I think it's best to keep the output file names very simple and just use the identifier ONLY as the prefix!


		# ### 7.9.2019 - testing single end mode: 
		# gunzip -f ${sampleId}_R1_trimmomatic_unpaired.fastq.gz  ${sampleId}_R2_trimmomatic_unpaired.fastq.gz
		# ### 	NEED TO CHECK IF READ DUPLICSTE NAMES MATTER - I think they do - single end not possible yet
		# cat ${sampleId}_R1_trimmomatic.fastq  ${sampleId}_R2_trimmomatic.fastq \
		# ${sampleId}_R1_trimmomatic_unpaired.fastq  ${sampleId}_R2_trimmomatic_unpaired.fastq \
		#  | sed  's/ \([12]\):[YN]:[0-9]:.\{1,\}/:\1/' \
		# > ${sampleId}_all_trimmed.fastq

		if [[ $stats == 'no' ]]; then
			# Remove the large fastq files:
			if [[ -s ${sampleId}_R1_trimmomatic.fastq ]]; then rm ${sampleId}_R1_trimmomatic.fastq; fi
			if [[ -s ${sampleId}_R1_trimmomatic_unpaired.fastq.gz ]]; then rm ${sampleId}_R1_trimmomatic_unpaired.fastq.gz; fi
			if [[ -s ${sampleId}_R2_trimmomatic.fastq ]]; then rm ${sampleId}_R2_trimmomatic.fastq ${sampleId}_R2_trimmomatic_unpaired.fastq.gz; fi
			if [[ -s ${sampleId}_R1_R2_trimmomatic.log ]];then rm ${sampleId}_R1_R2_trimmomatic.log; fi
			if [[ -s ${sampleId}_R1_trimmomatic.log ]]; then rm ${sampleId}_R1_trimmomatic.log; fi
		fi
	fi
elif [[ $hybSeqProgram == 'hybpiper'* ]]; then

	echo "Running gene recovery with HybPiper..."
	bwa=''
	targetFileFlag='--targetfile_aa' # only used by HybPiper2
	if [[ $hybSeqProgram == 'hybpiper-bwa' || $hybSeqProgram == 'hybpiper2-bwa' ]];then  
		bwa='--bwa'
		targetFileFlag='--targetfile_dna'
		echo Using HybPiper with the --bwa option...
	fi
	# else $bwa remains blank and the default option is used

	# First combine unpaired reads (both single end reads should have unique ids) - but won't I have the same problem as above?!
	unpairedFastqFile=''
	if [[ -n "$R2FastqFile" ]]; then
		gunzip -fc ${sampleId}_R1_trimmomatic_unpaired.fastq.gz ${sampleId}_R2_trimmomatic_unpaired.fastq.gz \
		> ${sampleId}_R1_R2_trimmomatic_unpaired.fastq
		unpairedFastqFile="--unpaired ${sampleId}_R1_R2_trimmomatic_unpaired.fastq"
	fi

	if [[ $hybSeqProgram == 'hybpiper' ]];then 
		
		echo Using HybPiper version 1.3 ...
		### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix reads_first.py --cpu $cpu $bwa - OK now so have put back
		$exePrefix reads_first.py --cpu $cpu $bwa \
		-b $targetsFile \
		-r ${sampleId}_R*_trimmomatic.fastq  \
		--cov_cutoff 4 \
		--prefix ${sampleId} \
		$unpairedFastqFile \
		> ${sampleId}_hybpiper.log 2>&1
		### 11.4.2022 - now using the --unpaired ${sampleId}_R1_R2_trimmomatic_unpaired.fastq \
		# Output: sampleId/geneId/sampleId/sequences/FNA/geneId.FNA; fasta header line: >sampleId
		# NBNB - From what I can make out, --cov_cutoff does seem to operate with spades, even though it says 
		#        flag is for velvetg - set to 4 otherwise default=8
		# NB - 11.4.2020 - trialing --length_pct set to 50

		if [ -s ${sampleId}/genes_with_seqs.txt ]; then
			# Alter the genewise files to be compatible with the make_species_trees_pipeline.sh script,
			# i.e. create sample fasta file containing all genes with faster header format >sampleId-geneId:
			for filePath in ${sampleId}/*/*/sequences/FNA/*.FNA; do
				geneName=`basename $filePath .FNA`
				cat $filePath \
				| awk -v gene=$geneName '{if($1 ~ /^>/) {print $1 "-" gene} else {print $0}}'
			done > ${sampleId}_all_genes.fasta

			### 2.11.2022 - temporary fix - time command doesn't exist after OS updates on the KewHPC nodes: $exePrefix intronerate.py  - OK now so have put back
			$exePrefix intronerate.py --prefix ${sampleId} --addN > ${sampleId}_intronerate.log 2>&1
			# Outputs e.g.:
			# geneId_supercontig.fasta; fasta header line: >sampleId-geneID
			# geneId_introns.fasta; fasta header line: >sampleId-geneID

			# Now concatenate the intronerate outputs:
			cat ${sampleId}/*/*/sequences/intron/*_supercontig.fasta \
			> ${sampleId}_all_genes_supercontig.fasta

			cat ${sampleId}/*/*/sequences/intron/*_introns.fasta \
			> ${sampleId}_all_genes_introns.fasta

			# Create a text file with the sampleId name in it for analysing the outputs for paralogs:
			echo $sampleId > ${sampleId}_name.txt
			# Also need a list of gene names:
			geneList=`cat ${sampleId}/genes_with_seqs.txt | awk '{printf $1 " "}' `
			# Output any paralogs found for each gene
			paralog_investigator.py $sampleId > ${sampleId}_paralog_investigator.log 2>&1
			# Now get all single copy genes AND paralogs into one file and
			# add any paralog suffixes to the geneId as well e.g. .main, .1.
			# Genes from geneId.FNA; fasta header line: >geneId
			# Genes from paralog_retriever.py output: >SampleId.main-geneId NODE_2_length_2887_cov_48.492391,Artocarpus-gene006,0,584,94.83,(-),2346,386
			### NBNB - the paralogs can come out as .0 or .1 etc even though there is only one paralog so they will be treated as different species in the gene trees.
			### Need to ensure that these have the same suffix; still need to check how they cluster in the gene trees though for RAxML.
			### Now need to look out for multiple paralogs...
			# First NOT including paralog suffixes in the geneId:
			for geneName in $geneList; do
				export geneName
				paralog_retriever.py  ${sampleId}_name.txt  $geneName \
				| perl -e '
	        	$paralogCountr = 1;
	        	while ($line = <>) {
	        		chomp $line;
	        		# Records with paralogs identified by detecting .main, .0, .1 etc:
	            	if($line =~ /^>/ && $line =~ /\.(main|\d+)/)    {
	            		($sampleParalog, $restOfLine) = split / /, $line;
	                	($sampleId, $paralogId) = split /\./, $sampleParalog;
	                	if($paralogId eq "main")	{
	                		print $sampleId, "_", $paralogId, "-", $ENV{geneName}, " $restOfLine\n"
	                		}
	                	else	{
	                		print $sampleId, "_p", $paralogCountr, "-", $ENV{geneName}, " $restOfLine\n";
	                		$paralogCountr++;
	                	}
	            	}
	            	elsif($line =~ /^>/)    {
	                	print $line, "-", $ENV{geneName}, "\n"
	            	}
	            	else {
	            		print $line, "\n"
	            	}
	        	}' 
	        done > ${sampleId}_all_single_copy_AND_paralogous_genes_NO_paralog_suffix.fasta  2> ${sampleId}_paralog_number.txt
			# Output fasta header formats in this file:
			# >sampleId_main-geneId		- main paralog chosen by HybPiper
			# >sampleId_p1-geneId		- first additional paralog
			# NB - both paralogs will be used in the make_species_trees_pipeline.sh script (NB - this may break the one gene, one locus assumption for ASTRAL)

			# Now creating geneIds with the paralog suffix, if present, so that no paralogs
			# will be used in the make_species_trees_pipeline.sh script:
			### NB - searching for .main or .\d is quite dangerous if sample names/ids also have dots in them! Unless you were to grab the text after the last dot or inform the user.
			### 17.6.2021 - Instead, you could split into an array, then get the last element (after last dot) and add remainder to $sampleId - easy! Also above file as well
			for geneName in $geneList; do
				export geneName
				paralog_retriever.py  ${sampleId}_name.txt  $geneName \
				| perl -e '
				$paralogCountr = 1;
	        	while ($line = <>) {
	        		chomp $line;
	        		# Records with paralogs identified by detecting .main, .0, .1 etc:
	            	if($line =~ /^>/ && $line =~ /\.(main|\d+)/)    {
	            		($sampleParalog, $restOfLine) = split / /, $line;
	                	($sampleId, $paralogId) = split /\./, $sampleParalog;
	                	if($paralogId eq "main")	{
	                		print $sampleId, "_", $paralogId, "-", $ENV{geneName}, "_", $paralogId, " $restOfLine\n"
	                	}
	                	else	{
	                		print $sampleId, "_p", $paralogCountr, "-", $ENV{geneName}, "_p", $paralogCountr, " $restOfLine\n";
	                		$paralogCountr++;
	                	}
	            	}
	            	elsif($line =~ /^>/)    {
	                	print $line, "-", $ENV{geneName}, "\n"
	            	}
	            	else {
	            		print $line, "\n"
	            	}
	        	}' 
	        done > ${sampleId}_all_single_copy_AND_paralogous_genes_PLUS_paralog_suffix.fasta  2> /dev/null


	        # Can also remove supercontigs with paralogs, by getting a list of ids from the above file (includes paralogs),
	        # then use list to retrieve seqs from the supercontigs file (paralogs will not be retrieved):
	        cat ${sampleId}_all_single_copy_AND_paralogous_genes_PLUS_paralog_suffix.fasta \
	        | grep '>' | sed 's/^>//' \
	        > ${sampleId}_all_single_copy_AND_paralogous_genes_PLUS_paralog_suffix_ids_only.txt
	        
	        seqtk subseq ${sampleId}_all_genes_supercontig.fasta \
	        ${sampleId}_all_single_copy_AND_paralogous_genes_PLUS_paralog_suffix_ids_only.txt \
	        > ${sampleId}_all_genes_supercontig_no_paralogs.fasta

	        # Report on genes with paralogs removed:
	        echo "Number of genes: " `cat ${sampleId}_all_genes.fasta | grep '>' | wc -l ` > ${sampleId}_gene_and_paralog_summary.txt
	        echo "Number of gene sequences plus all paralog sequences:" `cat ${sampleId}_all_single_copy_AND_paralogous_genes_PLUS_paralog_suffix.fasta | grep '>' | wc -l ` >> ${sampleId}_gene_and_paralog_summary.txt
	        echo "Number of genes without paralogs:" `cat ${sampleId}_all_genes_supercontig_no_paralogs.fasta | grep '>' | wc -l ` >> ${sampleId}_gene_and_paralog_summary.txt

	        # Also report the length of the contigs before and after removing paralogs:
	 		#### NB - this doeesnt work because I'm comparing exon only and supercontigs - sp don't bother!!!! - remove code
	        #echo "Sum length of genes:" `fastalength ${sampleId}_all_genes.fasta | awk '{sum += $1} END {print sum}' ` >> ${sampleId}_gene_and_paralog_summary.txt
			#echo "Sum length of genes without paralogs:" `fastalength ${sampleId}_all_genes_supercontig_no_paralogs.fasta | awk '{sum += $1} END {print sum}' ` >> ${sampleId}_gene_and_paralog_summary.txt

	        # [Outside this script, could also report number of genes affected by paralogs for all samples in the set.
	        #  At same time could also get number of times each gene is found across all samples with 'sort | uniq -c:''
	        #  cat Sample_*/*_all_genes_supercontig_no_paralogs.fasta | grep '>' | awk -F '-' '{print $2}' | sort | uniq -c | wc -l  ]

	        # HybPiper cleanup - remvoves the spades dir (Sample_/$sampleId/$geneName/$geneName_spades)
	        cleanup.py $sampleId
		fi
	elif [[ $hybSeqProgram == 'hybpiper2' ]]; then

		echo Using HybPiper version 2 ...
		$exePrefix hybpiper assemble  --cpu $cpu $bwa\
		$targetFileFlag $targetsFile \
		-r ${sampleId}_R*_trimmomatic.fastq  \
		--cov_cutoff 4 \
		--prefix ${sampleId} \
		$unpairedFastqFile \
		> ${sampleId}_hybpiper_assemble.log 2>&1
		# NB - if a DNA targets file is supplied without specifying the --bwa option, the targets will be translated and the blastx  option will proceed.
		#      Also found that when using --targetfile_aa flag with a DNA targets file (no --bwa flag), the targets gets translated.  

		if [[ -s ${sampleId}/genes_with_seqs.txt ]]; then

			# Need to create a small file containing the sample name for the following steps:
			echo ${sampleId} > namelist.txt

			# Get stats:
			hybpiper stats $targetFileFlag $targetsFile  gene  namelist.txt > ${sampleId}_hybpiper_stats.log 2>&1

			# Retrieve all the recovered sequence types for this sample:
			# NB - it seems that --targetfile_aa and --targetfile_dna can both accept the opposite residue type - OK 
			hybpiper retrieve_sequences dna $targetFileFlag $targetsFile --single_sample_name $sampleId > ${sampleId}_hybpiper_retrieve_sequences_dna.log 2>&1 # Creates file ${sampleId}_FNA.fasta
			hybpiper retrieve_sequences aa $targetFileFlag $targetsFile --single_sample_name $sampleId > ${sampleId}_hybpiper_retrieve_sequences_aa.log 2>&1 # Creates file ${sampleId}_FAA.fasta
			hybpiper retrieve_sequences intron $targetFileFlag $targetsFile --single_sample_name $sampleId > ${sampleId}_hybpiper_retrieve_sequences_intron.log 2>&1 # Creates file ${sampleId}_introns.fasta
			hybpiper retrieve_sequences supercontig $targetFileFlag $targetsFile --single_sample_name $sampleId > ${sampleId}_hybpiper_retrieve_sequences_supercontig.log 2>&1 # Creates file ${sampleId}_supercontigs.fasta

			# Get paralogs:
			hybpiper paralog_retriever namelist.txt $targetFileFlag $targetsFile
			# Outputs go to:
			# paralogs_no_chimeras/gene001_paralogs_no_chimeras.fasta - seems like the file to use 
			# paralogs_all/gene001_paralogs_all.fasta
			# The fasta header format is: >[sampleId].[main|0|1|2|etc]	- the geneId is not present so need to add it

			#Create sample file for the paralogs_no_chimeras:
			for file in paralogs_no_chimeras/*_paralogs_no_chimeras.fasta; do
				geneName=`basename $file | sed 's/_paralogs_no_chimeras.fasta//' `
	 			cat $file \
	 			| awk -v gene=$geneName '{if($1 ~ /^>/) {print $1 "-" gene} else {print $0}}' \
				> ${sampleId}_paralogs_no_chimeras.fasta
			done
			# Create sample file for the paralogs_all:
			for file in paralogs_all/*_paralogs_all.fasta; do
				geneName=`basename $file | sed 's/_paralogs_all.fasta//' `
	 			cat $file \
	 			| awk -v gene=$geneName '{if($1 ~ /^>/) {print $1 "-" gene} else {print $0}}' \
 				> ${sampleId}_paralogs_all.fasta
			done
		fi
	else
		echo "WARNING: If option -y was used, the HybPiper program was not recognised. The options are hybpiper[-bwa] or hybpiper2[-bwa]."
	fi

	# Converting to a tarball archive the main HybPiper results directory beneath the main <samplePrefix>_<sampleId> folder.
	# Reason: it's difficult for the file storage system to cope with lots of these directories 
	# containing many other directories and small files.
	if [[ -d $sampleId ]]; then
		tar -cf ${sampleId}.tar ${sampleId}
		gzip -f ${sampleId}.tar
		rm -fR $sampleId
	fi


	### 14.7.2023 - write to the PAFTOL db here - e.g.: see above for paftools
	### and doing before removing the unzipped files - might need to check whether script copes with this


	if [[ $stats == 'no' ]]; then
		# Remove the large fastq files:
		if [[ -s ${sampleId}_R1_trimmomatic.fastq ]]; then rm ${sampleId}_R1_trimmomatic.fastq; fi
		if [[ -s ${sampleId}_R1_trimmomatic_unpaired.fastq.gz ]]; then rm ${sampleId}_R1_trimmomatic_unpaired.fastq.gz; fi
		if [[ -s ${sampleId}_R2_trimmomatic.fastq ]]; then rm ${sampleId}_R2_trimmomatic.fastq ${sampleId}_R2_trimmomatic_unpaired.fastq.gz; fi
		if [[ -s ${sampleId}_R1_R2_trimmomatic.log ]];then rm ${sampleId}_R1_R2_trimmomatic.log; fi
		if [[ -s ${sampleId}_R1_trimmomatic.log ]]; then rm ${sampleId}_R1_trimmomatic.log; fi
		# NB - not deleting this file here in case it is used in the future for the recovery stats: ${sampleId}_R1_R2_trimmomatic_unpaired.fastq
		#      Could also just get it remade in the stats clause
	fi
else
	echo "WARNING: If option -y was used, the Hyb-Seq program was not recognised. The options are paftools or hybpiper[2[-bwa]]'."
	#exit
fi


if [[ $stats != 'no' ]]; then

	echo "Collecting gene recovery stats..."

	if [[ $usePaftolDb != 'no' ]]; then

		echo "If using PaftolDB with Paftools (i.e. if hidden option -d is also being used), need to run Trimmomatic again as the previous run was only saved to /tmp/ - still to add Trimmomatic step here"
		# NB - in all other cases trimmomatic is being run again above before this clause is reached
		exit
	fi

	# Let's put the stats files into a separate directory and work in there:
	mkdir -p gene_recovery_stats
	cd gene_recovery_stats
	
	# Find the gene recovery file for indexing with BWA
	refFileName=''
	if [[ $refFilePathForStats == 'default' ]]; then
		# Try to find file in ../pwd, as if recoveries have just been done:
		if [[ -s ../${sampleId}.fasta ]]; then # For Paftools 
		 	# Need to copy the gene recovery file to pwd so that the BWA indices go to pwd(!):
			cp -p  ../${sampleId}.fasta  ${sampleId}.fasta
			refFileName=${sampleId}.fasta
			# Gene recovery file is now in pwd.
		elif [[ -s ../${sampleId}_all_genes.fasta ]]; then
			cp -p  ../${sampleId}_all_genes.fasta  ${sampleId}.fasta
			refFileName=${sampleId}_all_genes.fasta
			# Gene recovery file is  now in pwd.
		elif [[ -s ../${sampleId}_FNA.fasta ]]; then
			cp -p  ../${sampleId}_FNA.fasta  ${sampleId}_FNA.fasta
			refFileName=${sampleId}_FNA.fasta
			# Gene recovery file is  now in pwd.
		else
			echo "ERROR: Option -S selected but gene recovery fasta file not found or is empty. May need to use option -P. Stats cannot be calculated for sample: ${sampleId}."
			echo
			echo
			exit
		fi
	else
		# Try to use path given in option -P:
		# For Paftools output:
		if [[ -s $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}.fasta ]]; then
			# Need to copy the gene recovery file to pwd so that the BWA indices go to pwd(!):
			cp -p  $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}.fasta  ${sampleId}.fasta
			refFileName=${sampleId}.fasta
		elif [[ -s $refFilePathForStats/${sampleId}.fasta ]]; then
			# Need to copy the gene recovery file to pwd so that the BWA indices go to pwd:
			cp -p  $refFilePathForStats/${sampleId}.fasta  ${sampleId}.fasta
			refFileName=${sampleId}.fasta
		# For HybPiper output:
		elif [[ -s $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}_all_genes.fasta ]]; then
			cp -p  $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}_all_genes.fasta  ${sampleId}_all_genes.fasta
			refFileName=${sampleId}_all_genes.fasta
		elif [[ -s $refFilePathForStats/${sampleId}_all_genes.fasta ]]; then
			cp -p  $refFilePathForStats/${sampleId}_all_genes.fasta  ${sampleId}_all_genes.fasta
			refFileName=${sampleId}_all_genes.fasta
		# For HybPiper2 output:
		elif [[ -s $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}_FNA.fasta ]]; then
			cp -p  $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}_FNA.fasta  ${sampleId}_FNA.fasta
			refFileName=${sampleId}_FNA.fasta
		elif [[ -s $refFilePathForStats/${sampleId}_FNA.fasta ]]; then
			cp -p  $refFilePathForStats/${sampleId}_FNA.fasta  ${sampleId}_FNA.fasta
			refFileName=${sampleId}_FNA.fasta
		# For combining nuclear, plastid and barcode genes reference files (an internal PAFTOL project reference file)
		elif [[ -s $refFilePathForStats/${sampleId}_nc_pt_nr.fasta ]]; then
			# For the moment only have this file set up, all files in one directory:
			cp -p  $refFilePathForStats/${sampleId}_nc_pt_nr.fasta ${sampleId}_nc_pt_nr.fasta
			refFileName=${sampleId}_nc_pt_nr.fasta
		else 
			echo "ERROR: option -P - can't find the correct path to the gene recovery fasta file or it is empty. Stats cannot be calculated for sample: ${sampleId}."
			echo 
			echo
			exit
		fi
	fi

	echo "Found gene recovery fasta file: $refFileName"

	# Prepare to map reads for getting stats:
### 18.8.2023 - change to ../ for fastqs incl SE reads below then that's it I think except for removing fastq files
	bamFileWithDups=''
	bwa index $refFileName
	bwa mem -t $cpu $refFileName \
	../${sampleId}_R1_trimmomatic.fastq \
	../${sampleId}_R2_trimmomatic.fastq \
	> ${sampleId}_bwa_mem_with_dups.sam

	### Don’t forget to pipe in/out files as much as possible to save space
	### ALSO, uptodate samtools versions can output straight to bam - take a look - shoudl be able to go from bwa sam stright to samtools sort
	samtools view -bS ${sampleId}_bwa_mem_with_dups.sam > ${sampleId}_bwa_mem_with_dups.bam
	# Need to sort bam before indexing
	samtools sort ${sampleId}_bwa_mem_with_dups.bam > ${sampleId}_bwa_mem_with_dups_sort.bam
	samtools index  ${sampleId}_bwa_mem_with_dups_sort.bam
	bamFileWithDups=${sampleId}_bwa_mem_with_dups_sort.bam

	if [[ $hybSeqProgram == 'hybpiper'* ]]; then
		# Also need to map the single end reads file:
		bwa index $refFileName
		bwa mem -t $cpu $refFileName \
		../${sampleId}_R1_R2_trimmomatic_unpaired.fastq \
		> ${sampleId}_bwa_mem_with_dups_unpaired_reads.sam

		samtools view -bS ${sampleId}_bwa_mem_with_dups_unpaired_reads.sam > ${sampleId}_bwa_mem_with_dups_unpaired_reads.bam
		# Need to sort bam before indexing
		samtools sort ${sampleId}_bwa_mem_with_dups_unpaired_reads.bam > ${sampleId}_bwa_mem_with_dups_unpaired_reads_sort.bam
		# Merge sorted bam files:
		samtools merge -f  ${sampleId}_bwa_mem_with_dups_sort_merged.bam  ${sampleId}_bwa_mem_with_dups_sort.bam  ${sampleId}_bwa_mem_with_dups_unpaired_reads_sort.bam
		# Resort bam (just in case):
		samtools sort ${sampleId}_bwa_mem_with_dups_sort_merged.bam > ${sampleId}_bwa_mem_with_dups_unpaired_reads_sort_merged_resort.bam
		samtools index  ${sampleId}_bwa_mem_with_dups_unpaired_reads_sort_merged_resort.bam
		bamFileWithDups=${sampleId}_bwa_mem_with_dups_unpaired_reads_sort_merged_resort.bam

		# This file is only created in hybpiper mode:
		if [[ -s ../${sampleId}_R1_R2_trimmomatic_unpaired.fastq ]]; then rm ../${sampleId}_R1_R2_trimmomatic_unpaired.fastq; fi 
	fi


	# Before assessing read depth stats, remove duplicates from the mapped bam file:
	if [[ ! -d tmp ]]; then mkdir tmp; fi 	# Not sure if this is vital - maybe sample size dependant
	### Need to test whether the java -jar -Xmx${mem_gb}g flag is required - try to fix it independant of setting memory in Slurm
	### Example in docs: for 8.6GB/20GB file, run with 2GB (-Xmx2g) and 10 GB hard memory
	### So far not needed to set memory with PAFTOL data: ${sampleId}_largest _R1_trimmomatic.fastq file done to date = 4.9GB
	java -jar  -XX:ParallelGCThreads=$cpu -Djava.io.tmpdir=tmp $PICARD MarkDuplicates \
	INPUT=$bamFileWithDups \
	OUTPUT=${sampleId}_bwa_mem_sort.bam \
	METRICS_FILE=${sampleId}_bwa_mem_sort_markdup_metrics \
	REMOVE_DUPLICATES=true \
	ASSUME_SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
	TMP_DIR=tmp
	# Using the new command syntax e.g. from INPUT to -INPUT
	# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP - ulimit -n=1024 on Kew hatta cluster, ulimit -n=256 on macbook so setting to 100 seems safe enough
	# OPTICAL_DUPLICATE_PIXEL_DISTANCE - default=100 for unpatterned versions of the Illumina platform; for the patterned flow cell models, 
	# 2500 is more appropriate. The patterned flow cell models are: iSeq™ 100, NextSeq™ 1000/NextSeq™ 2000, HiSeq™ 3000/HiSeq™ 4000, HiSeq™ X and NovaSeq™ 6000
	# SORTING_COLLECTION_SIZE_RATIO - default=0.25 - not using for the moment: 
    #                          This number, plus the maximum RAM available to the JVM, determine the memory footprint
    #                          used by some of the sorting collections.  If you are running out of memory, try reducing
    #                          this number.  Default value: 0.25.
    # NB - removing duplicates after quality trimming. Duplicates are assessed at 5' end and trimming is more likely to be done towards the 3' end
    #      so 5' read2 (and read1) coordinate should be preserved for duplicate assessment

    ### NB - Picard commandline syntax is changing - new format would look like e.g.:	- just need add a dash char infront of each flag
    ### MarkDuplicates -INPUT 7210_bwa_mem_with_dups_sort.bam -OUTPUT 7210_bwa_mem_sort.bam -METRICS_FILE 7210_bwa_mem_sort_markdup_metrics -REMOVE_DUPLICATES true -ASSUME_SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 -TMP_DIR tmp

    if [[ ! -s ${sampleId}_bwa_mem_sort.bam ]]; then
    	echo "ERROR: Picard markduplicates output file doesn't exist or is empty: ${sampleId}_bwa_mem_sort.bam"
    	### 19.4.2023 - NB - if this file doesn't exist, then PICARD will already have failed and the script will exit then due to 'set -e' being used
    	exit
    fi

	
	# Getting contig recovery stats and adding each statistic to a file, one statistic per line.
	# Then in overall_gene_recovery_stats.sh script, the results can be converted to a table for 
	# all sample and genes.
	# Number of recovered genes:
	numbrRecoveredGenes=`cat $refFileName | grep '>' | wc -l `
	# Sum length of genes:
	if [[ $hybSeqProgram == 'hybpiper2' ]]; then
		# Remove strings of N's from the sequence line before counting the number of bases
		# NB - probably could use this code on all recovery methods; but will also remove a (very) 
		#      small number (?) of N's generated by the assembler - OK I think - Paftools overlapSerial produces some.
		cat $refFileName \
		| awk '{if($1 ~ /^>/) { print $0 } else { {gsub(/[Nn]/,"",$0)} {print $0} } }' \
		| grep -v ^$ \
		> $refFileName.Ns_removed_temp
		sumLengthOfGenes=`fastalength $refFileName.Ns_removed_temp | awk '{sum+=$1} END {print sum}' `
		#ls -l $refFileName.Ns_removed_temp
		###rm $refFileName.Ns_removed_temp  ### Delaying removing this file so I can compare it with N's included in the read stats
	else
		sumLengthOfGenes=`fastalength $refFileName | awk '{sum+=$1} END {print sum}' `
	fi
	echo "sampleId: $sampleId
numbrRecoveredGenes: $numbrRecoveredGenes
sumLengthOfGenes: $sumLengthOfGenes" > ${sampleId}_gene_recovery_stats.txt  # Also wipes out file contents from any previous run

	# Count number of all ambiguity codes:
	###numbrAmbiguityCodesInGenes=`cat $refFileName | grep -v '>' | grep -o '[RYMKSWHBDN]' | wc -l `
	###echo "numbrAmbiguityCodesInGenes: $numbrAmbiguityCodesInGenes" >> ${sampleId}_gene_recovery_stats.txt
	### For some reason doesn't work - may have to turn off set +u and +e - see below

	# Count the total number of raw fastq file reads (before any type of trimming this script does)
	numbrRawReads=`cat ../${sampleId}_trimmomatic.log | grep 'Input Read Pairs:' | awk '{print $4 *2}' ` 
	echo numbrRawReads: $numbrRawReads >> ${sampleId}_gene_recovery_stats.txt
	
	####################################
	# General stats on the BWA alignment
	####################################
	# Count the number of reads in the bam file just after mapping but before removing read duplicates:
	numbrTrimmedReadsInBamInclDups=`samtools view -c $bamFileWithDups `
	echo numbrTrimmedReadsInBamInclDups: $numbrTrimmedReadsInBamInclDups  >> ${sampleId}_gene_recovery_stats.txt

	percentTrimmedReads=`echo $numbrTrimmedReadsInBamInclDups $numbrRawReads | awk '{printf "%.1f", (($1 / $2) * 100)}' | awk '{print $1}' ` # Last awk produces the line return
	echo percentTrimmedReads: $percentTrimmedReads >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads in the bam file just after mapping and removing read duplicates.
	# NB - doesn't quite count all reads in the file.
	numbrTrimmedReadsInBam=`samtools view -c ${sampleId}_bwa_mem_sort.bam `
	echo numbrTrimmedReadsInBam: $numbrTrimmedReadsInBam  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads that are mapped to the reference with map quality >= 20:
	numbrMappedReads=`samtools view -c -F4 -q 20 ${sampleId}_bwa_mem_sort.bam `
	echo numbrMappedReads: $numbrMappedReads  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads mapped in proper pairs:
	numbrMappedReadsProperPairs=`samtools view -c -f2 -q 20 ${sampleId}_bwa_mem_sort.bam `
	echo numbrMappedReadsProperPairs: $numbrMappedReadsProperPairs  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads with an XA:Z flag (alternative hits).
	# NB - normally this flag is in column 16 or column 17 but not always.
	# 'XA:Z' might be in column 11 (quality) (?) so just check from 
	# column 12 onwards:
	numbrReadsWithAltHits=`samtools view -F4 -q 20 ${sampleId}_bwa_mem_sort.bam | cut -f12- | awk '$0 ~ /XA:Z/' | wc -l `
	# NB for some reason | grep 'XA:Z' didn't work!
	echo numbrReadsWithAltHits: $numbrReadsWithAltHits >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of unmapped reads:
	numbrUnmappedReads=`samtools view -c -f4 ${sampleId}_bwa_mem_sort.bam `
	echo numbrUnmappedReads: $numbrUnmappedReads  >> ${sampleId}_gene_recovery_stats.txt
	# Output the unmapped reads and convert to a fastq file.
	# (Need to re-sort bam by fastq record name for use with 'samtools fastq'.)
	# Not needed now - can select unmapped reads with samtools fastq:
	# samtools view -f4 ${sampleId}_bwa_mem_sort.bam > ${sampleId}_bwa_mem_sort_unmapped.bam

	# First, sort bam by fastq record id:
	samtools sort -n ${sampleId}_bwa_mem_sort.bam \
	| samtools fastq -f4 -1 ${sampleId}_bwa_mem_unmapped_R1.fastq.gz -2 ${sampleId}_bwa_mem_unmapped_R2.fastq.gz \
	-s ${sampleId}_bwa_mem_unmapped_single_ends.fastq.gz -N
	# samtools fastq -n - means that /1 and /2 are NOT added to output records - didn't work here
	# samtools fastq -N - means that /1 and /2 are ALWAYS added to fastq record ids
	# 	NB - not sure if this needs to be done - all singleton reads should be unique (would only be a problem if combining the R1 AND R2 read pairs)

	#######################
	# Reads on-target stats
	#######################
	# Reads per sample on-target with samtools view -L - need a bedfile:
	fastalength $refFileName | awk '{print $2 " 0 " $1}'  > ${sampleId}.bed
	# samtools view -F4 -L 10895.bed 10895_bwa_mem_sort.bam > 10895_bwa_mem_sort_st_-L.sam
	# 2,030,124 - total # reads = 2,784,802 = 72.9 % on-target
	numbrReadsOnTarget=`samtools view -c -F4 -q 20 -L ${sampleId}.bed ${sampleId}_bwa_mem_sort.bam `
	# -c counts the number of reads
	# -q mapping quality (NB -in other samtools program -Q is mapping quality!)
	# -F4 = INT of 4 = unmapped read but -F prints reads that are not INT=4 i.e. mapped - definitely NOT the same as samtools -c.
	# This value should be the same as $numbrMappedReads as there are no reference bases off target (?) - keep thinking the logic
	echo numbrReadsOnTarget: $numbrReadsOnTarget >> ${sampleId}_gene_recovery_stats.txt

	percentReadsOnTarget=`echo $numbrReadsOnTarget $numbrTrimmedReadsInBam | awk '{printf "%.1f", (($1 / $2) * 100)}' | awk '{print $1}' ` # Last awk produces the line return 
	echo percentReadsOnTarget: $percentReadsOnTarget >> ${sampleId}_gene_recovery_stats.txt

	###############################
	# Read coverage and depth stats
	###############################
	set +e	# NB - had a problem with script exiting half way through making samtools depth stats. Haven't got to the bottom of it but turning OFF these options here works for now.  
	set +u

	###################
	# samtools coverage - “meandepth” per gene (column 7 - includes positions with zero depth - see below)
	####################
	samtools coverage -q 20 -Q 20 ${sampleId}_bwa_mem_sort.bam > ${sampleId}_bwa_mem_sort_st_covrg.txt
	# Removed --reference $refFileName - I don't think it is required - not sure why you need to supply it - same for samtools depth
	# NB - It is necessary to implement a mapping quality threshold e.g. 20 would be OK - otherwise valeu a very high
	# -q 20 base quality - NBNB - on closer inspection there is an error in the samtools view command line docs - -q and -Q could be the other way round - Ok for now if I use 20 for each.
	# -Q 20 mapping quality

	# Can also visualise coverage across each contig as a histogram: 
	samtools coverage -q 20 -Q 20 -m ${sampleId}_bwa_mem_sort.bam > ${sampleId}_bwa_mem_sort_st_covrg_-m.txt 
	# To pick out a plot for a specific contig add: grep -A 11 '^contigId

	# Stats for all genes per sample:
	# coverage (column 6) - % bases covered by one or more reads per gene - might be an interesting value.
	# Mean coverage per sample:
	meanReadCovrg=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | awk '{sum+=$6} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
	echo meanReadCovrg: $meanReadCovrg >> ${sampleId}_gene_recovery_stats.txt
	# Median coverage per sample:
	medianPoint=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | awk 'END {printf "%.0f" , NR/2}' `	# Don't need to sort here!
	medianReadCovrg=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | sort -k6n | awk '{print $6}' | head -n $medianPoint | tail -n 1 `
	echo medianReadCovrg: $medianReadCovrg >> ${sampleId}_gene_recovery_stats.txt

	# 'meandepth' (column 7 -  includes positions with zero depth) per sample:
	meanReadDepth_min0x=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | awk '{sum+=$7} END {if(sum > 0) {print sum/NR} else {print "0"}}' `    # average
	echo "meanReadDepth_min0x_(samtools_coverage): $meanReadDepth_min0x" >> ${sampleId}_gene_recovery_stats.txt

	# Median “meandepth" per sample:
	medianPoint1=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | awk 'END {printf "%.0f", NR/2}' `
	medianReadDepth_min0x=`tail -n+2 ${sampleId}_bwa_mem_sort_st_covrg.txt | sort -k7n | awk '{print $7}' | head -n $medianPoint1 | tail -n 1 `
	echo "medianReadDepth_min0x_(samtools_coverage): $medianReadDepth_min0x" >> ${sampleId}_gene_recovery_stats.txt


	################
	# samtools depth - to obtain depth for read depth >= Xx
	################
	# The aim here is to count useful depth i.e. >1x or >4x is more informative - need to use samtools depth
	# The -a flag outputs all positions (including zero depth). We already have this above from samtools coverage - column 7) but not for ALL genes together.
	# So including the -a flag here then can assess from >=0x coverage and compare with samtools coverage - they are slightly different calculations!
	samtools depth -q 20 -q 20 -a ${sampleId}_bwa_mem_sort.bam > ${sampleId}_bwa_mem_sort_st_depth.txt

	# Mean read depth for bases with >= 0x depth across ALL genes:
	meanReadDepth_min0x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 0' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
	echo "meanReadDepth_min0x_(samtools_depth): $meanReadDepth_min0x" >> ${sampleId}_gene_recovery_stats.txt

	# Median read depth for bases with >= 0x depth across ALL genes:
	medianPoint2=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 0' | awk 'END {printf "%.0f" , NR/2}' `
	medianReadDepth_min0x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 0' | sort -k3n | awk '{print $3}' | head -n $medianPoint2 | tail -n 1 `
	echo "medianReadDepth_min0x_(samtools_depth): $medianReadDepth_min0x" >> ${sampleId}_gene_recovery_stats.txt

	# Mean read depth for bases with >= 1x depth across ALL genes:
	meanReadDepth_min1x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 1' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
	echo "meanReadDepth_min1x_(samtools_depth): $meanReadDepth_min1x" >> ${sampleId}_gene_recovery_stats.txt

	# Median read depth for bases with >= 1x depth across ALL genes:
	medianPoint3=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 1' | awk 'END {printf "%.0f" , NR/2}' `
	medianReadDepth_min1x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 1' | sort -k3n | awk '{print $3}' | head -n $medianPoint3| tail -n 1 `
	echo "medianReadDepth_min1x_(samtools_depth): $medianReadDepth_min1x" >> ${sampleId}_gene_recovery_stats.txt
	
	# Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepth_min4x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 4' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
	echo "meanReadDepth_min4x_(samtools_depth): $meanReadDepth_min4x" >> ${sampleId}_gene_recovery_stats.txt

	# Median read depth for bases with >= 4x depth  across ALL genes:
	medianPoint4=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 4' | awk 'END {printf "%.0f" , NR/2}' `
	if [[ $medianPoint4 -eq 0 ]]; then
		echo "medianReadDepth_min4x_(samtools_depth): 0" >> ${sampleId}_gene_recovery_stats.txt
		# head -n 0 gives error: head: illegal line count -- 0; only affects very poor quality samples, don't think it affects depth at >= 0 and >= 1.
	else
		medianReadDepth_min4x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 4' | sort -k3n | awk '{print $3}' | head -n $medianPoint4 | tail -n 1 `
		echo "medianReadDepth_min4x_(samtools_depth): $medianReadDepth_min4x" >> ${sampleId}_gene_recovery_stats.txt
	fi


	# Also comparing read depth at >=4x WITH duplicates NOT removed to read depth at >=4x without duplicates.
	samtools depth -q 20 -q 20 -a $bamFileWithDups > ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt

	# Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepthWithDups_min4x=`cat ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt | awk '$3 >= 4' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
	echo "meanReadDepthWithDups_min4x_(samtools_depth): $meanReadDepthWithDups_min4x" >> ${sampleId}_gene_recovery_stats.txt

	# Median read depth for bases with >= 4x depth  across ALL genes:
	medianPoint5=`cat ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt | awk '$3 >= 4' | awk 'END {printf "%.0f" , NR/2}' `
	if [[ $medianPoint4 -eq 0 ]]; then
		echo "medianReadDepthWithDups_min4x_(samtools_depth): 0" >> ${sampleId}_gene_recovery_stats.txt
	else
		medianReadDepthWithDups_min4x=`cat ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt | awk '$3 >= 4' | sort -k3n | awk '{print $3}' | head -n $medianPoint5 | tail -n 1 `
		echo "medianReadDepthWithDups_min4x_(samtools_depth): $medianReadDepthWithDups_min4x" >> ${sampleId}_gene_recovery_stats.txt
	fi


	# Number of recovered bases in ALL genes with read depth >= 4:
	# With duplicates removed:
	numbrBasesInAllGenes_ReadDepth_min4x=`cat ${sampleId}_bwa_mem_sort_st_depth.txt | awk '$3 >= 4' | wc -l `
	echo "NumbrBasesInAllGenes_ReadDepth_min4x_(samtools_depth): $numbrBasesInAllGenes_ReadDepth_min4x" >> ${sampleId}_gene_recovery_stats.txt
	# With duplicates:
	numbrBasesInAllGenes_ReadDepthWithDups_min4x=`cat ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt | awk '$3 >= 4' | wc -l `
	echo "NumbrBasesInAllGenes_ReadDepthWithDups_min4x_(samtools_depth): $numbrBasesInAllGenes_ReadDepthWithDups_min4x" >> ${sampleId}_gene_recovery_stats.txt


	#########################################
	# Further stats to do outside this script
	#########################################
	# 1. Read depth per gene across all samples --> results for 353 genes
	# 2. Read depth across all genes and samples - total mean and median values

	# Remove the large fastq files from any gene recovery method:
	if [[ -s ${sampleId}_R1_trimmomatic.fastq ]]; then rm ${sampleId}_R1_trimmomatic.fastq; fi
	if [[ -s ${sampleId}_R1_trimmomatic_unpaired.fastq.gz ]]; then rm ${sampleId}_R1_trimmomatic_unpaired.fastq.gz; fi
	if [[ -s ${sampleId}_R2_trimmomatic.fastq ]]; then rm ${sampleId}_R2_trimmomatic.fastq ${sampleId}_R2_trimmomatic_unpaired.fastq.gz; fi
	if [[ -s ${sampleId}_R1_R2_trimmomatic.log ]];then rm ${sampleId}_R1_R2_trimmomatic.log; fi
	if [[ -s ${sampleId}_R1_trimmomatic.log ]]; then rm ${sampleId}_R1_trimmomatic.log; fi
	# NB - ${sampleId}_R1_R2_trimmomatic_unpaired.fastq is only created in hybpiper mode and has already been removed above - OK

	# Finally, remove the largest files created for the stats (they take up TB of disk space if processing 1000's samples)
	if [[ -s ${sampleId}_bwa_mem_with_dups.sam ]]; then rm ${sampleId}_bwa_mem_with_dups.sam; fi
	if [[ -s ${sampleId}_bwa_mem_with_dups.bam ]]; then rm ${sampleId}_bwa_mem_with_dups.bam; fi
	if [[ -s ${sampleId}_bwa_mem_with_dups_sort.bam ]]; then rm ${sampleId}_bwa_mem_with_dups_sort.bam; fi
	if [[ -s ${sampleId}_bwa_mem_with_dups_sort.bam.bai ]]; then rm ${sampleId}_bwa_mem_with_dups_sort.bam.bai; fi
	if [[ -d tmp ]]; then rm -fR tmp; fi
	if [[ -s ${sampleId}_bwa_mem_sort.bam ]]; then rm ${sampleId}_bwa_mem_sort.bam; fi
fi
#####cd ../ # Back up to parent dir for next sample - 20.4.2020 - has no effect here now and not required anymore because looping through samples is done outside this script 
echo 
