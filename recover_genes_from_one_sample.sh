#!/bin/bash

##################################
# recover_genes_from_one_sample.sh

# Author: Paul Bailey
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
echo stats: $stats
echo refFilePathForStats: $refFilePathForStats
  

sampleId=`echo $line | cut -d ',' -f 1 `    		# I think this will be the PAFTOL_xxxxxx id - same as in the read file names
R1FastqFile=`echo $line | cut -d ',' -f 2 `
R2FastqFile=`echo $line | cut -d ',' -f 3 `
externalSequenceID=`echo $line | cut -d ',' -f 4 `	# For adding the external sequence Id to the paftol_da

  
if [[ ! -d ${samplePrefix}_$sampleId ]]; then mkdir ${samplePrefix}_$sampleId; fi
cd ${samplePrefix}_$sampleId
echo sampleId: $sampleId
ls $paftolDataSymlinksDir/$R1FastqFile
ls $paftolDataSymlinksDir/$R2FastqFile
pwd


#if [ $usePaftolDb != 'usePaftolDb' ]; then	- changed - now introducing data set type
if [[ $usePaftolDb == 'no' ]]; then
	
  			                                        #--nodelist=kppgenomics01.ad.kew.org  # mem normally set to 80000
	#  sbatch -J ${samplePrefix}_${sampleId}_fastqToGenes -p main -t 1-0:00 -c $cpu --mem=80000 -o ${samplePrefix}_${sampleId}_fastqToGenes.log   -e ${samplePrefix}_${sampleId}_fastqToGenes.log_err   --wrap "
	# RUNTIME: For 8 cpu, up to 2 mins; up to 18 GB mem (for 10 samples)
	#          For 4 cpu, up to 24 mins; up to 12 GB mem.    - Total # samples to date: ~2500 - would take 14h using 176 cpu (1 node) 
	#          For 2 cpu, up to 1h50', up to 10 GB mem (tested 40 samples)
	#          Conclusion: the above stats seems to suggest that running this step separately might be the most efficient
	pathToTrimmomatic=`which trimmomatic-0.39.jar `		# NB - 'which' requires the file to be executable!
	$exePrefix  java -jar $pathToTrimmomatic PE \
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
	# 8.3.2020 - changed to palidromic mode - confirm that the extra two parameters is Ok still for normal mode
  

	# Paftools requires unzipped fastq files (NB - if unzipped files already present, files will not be unzipped again (i.e gzip exits with an error), so I've included the -f flag to force uncompression!):
	#srun -J ${sampleId}_unzip_R1_R2 -n 1  gunzip -f ${sampleId}_R1_trimmomatic.fq.gz  ${sampleId}_R2_trimmomatic.fq.gz
	gunzip -f ${sampleId}_R1_trimmomatic.fastq.gz  ${sampleId}_R2_trimmomatic.fastq.gz
fi


### NB - I think that if I check file exists or is zero byte before running this is equivalent to make affect
### convert script to a makefile


usePaftolDbFlag=''
if [ $hybSeqProgram == 'paftools' ]; then

	#if [ $usePaftolDb != 'usePaftol' ]; then - To fit with other data types, changed to house the dataset acronym in $usePaftolDb 
	if [ $usePaftolDb != 'no' ]; then 
		usePaftolDbFlag='--usePaftolDb'

		# Remove the .gz ending from the file name for adding to the pafto_da database
		unzippedR1FastqFile=`basename -s .gz $paftolDataSymlinksDir/$R1FastqFile `
		unzippedR2FastqFile=`basename -s .gz $paftolDataSymlinksDir/$R2FastqFile `

		# Need to uncompress raw fastq files to get FastQCStats and upload into the paftol_da database:
		gunzip -f -c $paftolDataSymlinksDir/$R1FastqFile > $unzippedR1FastqFile
		gunzip -f -c $paftolDataSymlinksDir/$R2FastqFile > $unzippedR2FastqFile
### 6.11.2020 - still need to add a conditional for --sampleId - it can't be present but blank
### 			if $externalSequenceID is not null externalSequenceID='--sampleId $externalSequenceID'
###				then  can juast put the variable in the command and it will be empty or not --> almost; need an else clause to aslo add paftol info ()
###				need an extra flag in the wrapper for sampleId e.g. -e, also -d flag needs to have a dataOrigin value.
		export PYTHONPATH=$HOME/lib/python
		paftools addPaftolFastq  $unzippedR1FastqFile  $unzippedR2FastqFile \
		--fastqPath $paftolDataSymlinksDir \
		--dataOrigin $usePaftolDb \
		--sampleId $externalSequenceID
		# NB - the file name must be of this format e.g. PAFTOL_005853_R1.fastq
		# The fastqPath entered in the database consists of the path and filename e.g. $paftolDataSymlinksDir/$unzippedR1FastqFile
		# --sampleId=$sampleId - must be included for all data types, but in the case of paftol it is ignored (but there is no harm in always including it)

		### Add a new option to THIS program for uploading SRA e.g. -d PAFTOL, -d SRA

		echo "Exit status of paftools addpaftolFastq:" $?
		### NB - paftools recoverSeqs also runs FastQC stats against the raw reads but this is unnessary.
		###     Try to delete this step then I can just remove the raw fastqs here.

		# Paftools --usePaftolDb requires the targets filename WITHOUT the path
		# i.e. it insists that targets file is in the pwd (but I could change that):
		cp $targetsFile . 
		targetsFile=`basename $targetsFile `

		# Use of --usePaftolDb flag requires the use of the trimmomatic flag so have to have a separate paftoools recoverSeqs command here.
		export PYTHONPATH=$HOME/lib/python 			# I had to add this for the cluster ONLY - need to. check it is OK on Macbook, it should be.
		$exePrefix  paftools recoverSeqs \
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
		$usePaftolDbFlag \
		> ${sampleId}_overlapSerial.log 2>&1

		rm $targetsFile
		### 15.4.2020 - also remove PAFTOL_000948_R1.fastq and PAFTOL_000948_R2.fastq files to save space
	else 

		################
		# overlapRecover (OLC approach)
		################
		# RUNTIME: for 4 cpu, up to 16h21' + up to >>23GB mem (average mem = 14GB assessed over 142 samples; a few samples use up to 125GB mem in 8-11h); set to 80GB mem (or lower but will get more fails)
		# Was using:
		#srun -J ${sampleId}_overlapRecover -n 1  -o ${sampleId}_overlapRecover.log  -e ${sampleId}_overlapRecover.log_err \	- issue with srun
		export PYTHONPATH=$HOME/lib/python 			# I had to add this for the cluster ONLY - need to. check it is OK on Macbook, it should be.
		$exePrefix  paftools recoverSeqs \
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

		### 15.4.2020 - also remove ${sampleId}_R1_trimmomatic.fastq and ${sampleId}_R2_trimmomatic.fastq 7167_R1_trimmomatic_unpaired.fastq.gz 7167_R2_trimmomatic_unpaired.fastq.gz if they exist
	fi
elif [[ $hybSeqProgram == 'hybpiper'* ]]; then

	bwa=''
	if [[ $hybSeqProgram == 'hybpiper-bwa' ]];then 
		bwa='--bwa'
		echo Using HybPiper with the --bwa option...
	fi

	# First combine unpaired reads (both single end reads should have unique ids) - but won't I have the same problem as above?!
	gunzip -fc ${sampleId}_R1_trimmomatic_unpaired.fastq.gz ${sampleId}_R2_trimmomatic_unpaired.fastq.gz \
	> ${sampleId}_R1_R2_trimmomatic_unpaired.fastq

	$exePrefix reads_first.py --cpu $cpu $bwa \
	-b $targetsFile \
	-r ${sampleId}_R*_trimmomatic.fastq  \
	--cov_cutoff 4 \
	--prefix ${sampleId} \
	--unpaired ${sampleId}_R1_R2_trimmomatic_unpaired.fastq \
	> ${sampleId}_hybpiper.log 2>&1
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
else
	echo "WARNING: the Hyb-Seq program was not recognised. The options are 'paftools' or 'hybpiper' (without the quotes)."
	#exit
fi


if [[ $stats != 'no' ]]; then

	echo "Collecting gene recovery stats..."

	if [[ $usePaftolDb != 'no' ]]; then

		echo If using PaftolDB with Paftools, need to run Trimmomatic again as the previous run was only saved to /tmp/ - still to add
	fi	

	
	# Find the gene recovery file for indexing with BWA
	refFileName=''
	if [[ $refFilePathForStats == 'default' ]]; then
		# Try to find file in pwd, as if recoveries have just been done:
		if [[ -s ${sampleId}.fasta ]]; then 
			refFileName=${sampleId}.fasta
			# Gene recovery file is in pwd.
		elif [[ -s ${sampleId}_all_genes.fasta ]]; then
			refFileName=${sampleId}_all_genes.fasta
			# Gene recovery file is in pwd.
		else
			echo "ERROR: Option -S selected but gene recovery fasta file not found or is empty. May need to use option -P. Stats cannot be calculated for sample: ${sampleId}."
			echo
			echo
			exit
		fi
	else
		# Try to use path given in option -P:
		if [[ -s $refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}.fasta ]]; then 
			refFileName=$refFilePathForStats/${samplePrefix}_${sampleId}/${sampleId}.fasta
			# Copy the gene recovery file to pwd so that the BWA indices go to pwd:
			cp -p  $refFileName ${sampleId}.fasta
		elif [[ -s $refFilePathForStats/${sampleId}.fasta ]]; then
			refFileName=$refFilePathForStats/${sampleId}.fasta
			# Copy the gene recovery file to pwd so that the BWA indices go to pwd:
			cp -p  $refFileName ${sampleId}.fasta
		else 
			echo "ERROR: option -P - can't find the correct path to the gene recovery fasta file or is empty. Stats cannot be calculated for sample: ${sampleId}."
			echo 
			echo
			exit
		fi
	fi

	echo "Found gene recovery fasta file: $refFileName"


	# Prepare to map reads for getting stats:
	bwa index ${sampleId}.fasta
	bwa mem -t $cpu ${sampleId}.fasta \
	${sampleId}_R1_trimmomatic.fastq \
	${sampleId}_R2_trimmomatic.fastq \
	> ${sampleId}_bwa_mem_with_dups.sam
### NB - WHAT HAPPENS WHEN USING HYBPIPER? NEED TO ALSO MAP THE UNPAIRED READS
### Don’t forget to pipe in/out files as much as possible to save space
	samtools view -bS ${sampleId}_bwa_mem_with_dups.sam > ${sampleId}_bwa_mem_with_dups.bam
	# Need to sort bam before indexing
	samtools sort ${sampleId}_bwa_mem_with_dups.bam > ${sampleId}_bwa_mem_with_dups_sort.bam
	samtools index  ${sampleId}_bwa_mem_with_dups_sort.bam

	# Before assessing read depth stats, remove duplicates from the mapped bam file:
	if [[ ! -d tmp ]]; then mkdir tmp; fi 	# Not sure if this is vital - maybe sample size dependant
	### Need to test whether the java -jar -Xmx${mem_gb}g flag is required - try to fix it independaqnt of settign memory in Slurm
	java -jar  -XX:ParallelGCThreads=$cpu -Djava.io.tmpdir=tmp $PICARD MarkDuplicates \
	-INPUT=${sampleId}_bwa_mem_with_dups_sort.bam \
	-OUTPUT=${sampleId}_bwa_mem_sort.bam \
	-METRICS_FILE=${sampleId}_bwa_mem_sort_markdup_metrics \
	-REMOVE_DUPLICATES=true \
	-ASSUME_SORT_ORDER=coordinate \
	-VALIDATION_STRINGENCY=LENIENT \
	-MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
	-TMP_DIR=tmp
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


	> ${sampleId}_gene_recovery_stats.txt  # Wiping out file contents from any previous run
	
	####################################
	# General stats on the BWA alignment
	####################################
	# Count the number of reads in the bam file just after mapping but before removing read duplicates:
	numbrTrimmedReadsWithDups=`samtools view -c ${sampleId}_bwa_mem_with_dups_sort.bam `
	echo numbrTrimmedReadsWithDups: $numbrTrimmedReadsWithDups  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads in the bam file just after mapping and removinf read duplicates.
	# NB - doesn't quite count all reads in the file.
	numbrTrimmedReads=`samtools view -c ${sampleId}_bwa_mem_sort.bam `
	echo numbrTrimmedReads: $numbrTrimmedReads  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads that are mapped to the reference with map quality >= 20:
	numbrMappedReads=`samtools view -c -F4 -q 20 ${sampleId}_bwa_mem_sort.bam `
	echo numbrMappedReads: $numbrMappedReads  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads mapped in proper pairs:
	numbrMappedReadsProperPairs=`samtools view -c -f2  ${sampleId}_bwa_mem_sort.bam `
	echo numbrMappedReadsProperPairs: $numbrMappedReadsProperPairs  >> ${sampleId}_gene_recovery_stats.txt

	# Count the number of reads with an XA:Z flag (alternative hits).
	# NB - normally this flag is in column 16 or column 17 but not always.
	# 'XA:Z' might be in column 11 (quality) (?) so just check from 
	# column 12 onwards:
	numbrReadsWithAltHits=`samtools view -F4 -q 20 ${sampleId}_bwa_mem_sort.bam | cut -f12- | awk '$0 ~ /XA:Z/' | wc -l `
	# NB for some reason | grep 'XA:Z' didn't work!
	echo numbrReadsWithAltHits: $numbrReadsWithAltHits >> ${sampleId}_gene_recovery_stats.txt


	#######################
	# Reads on-target stats
	#######################
	# Reads per sample on-target with samtools view -L - need a bedfile:
	fastalength ${sampleId}.fasta | awk '{print $2 " 0 " $1}'  > ${sampleId}.bed
	# samtools view -F4 -L 10895.bed 10895_bwa_mem_sort.bam > 10895_bwa_mem_sort_st_-L.sam
	# 2,030,124 - total # reads = 2,784,802 = 72.9 % on-target
	numbrReadsOnTarget=`samtools view -c -F4 -q 20 -L ${sampleId}.bed ${sampleId}_bwa_mem_sort.bam `
	# -c counts the number of reads
	# -q mapping quality (NB -in other samtools program -Q is mapping quality!)
	# -F4 = INT of 4 = unmapped read but -F prints reads that are not INT=4 i.e. mapped - definitely NOT the same as samtools -c.
	# This value should be the same as $numbrMappedReads as there are no reference bases off target (?) - keep thinking the logic
	echo numbrReadsOnTarget: $numbrReadsOnTarget >> ${sampleId}_gene_recovery_stats.txt
	

	###############################
	# Read coverage and depth stats
	###############################
	set +e	# NB - had a problem with script exiting half way through making samtools depth stats. Haven't got to the bottom of it but turning OFF these options here works for now.  
	set +u

	###################
	# samtools coverage - “meandepth” per gene (column 7 - includes positions with zero depth - see below)
	####################
	samtools coverage -q 20 -Q 20 ${sampleId}_bwa_mem_sort.bam > ${sampleId}_bwa_mem_sort_st_covrg.txt
	# Removed --reference ${sampleId}.fasta - I don't think it is required - not sure why you need to supply it - same for samtools depth
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
	samtools depth -q 20 -q 20 -a ${sampleId}_bwa_mem_with_dups_sort.bam > ${sampleId}_bwa_mem_with_dups_sort_st_depth.txt
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


	#########################################
	# Further stats to do outside this script
	#########################################
	# 1. Read depth per gene across all samples --> results for 353 genes
	# 2. Read depth across all genes and samples - total mean and median values 
fi
#####cd ../ # Back up to parent dir for next sample - 20.4.2020 - has no effect here now and not required anymore because looping through samples is done outside this script 
echo
echo