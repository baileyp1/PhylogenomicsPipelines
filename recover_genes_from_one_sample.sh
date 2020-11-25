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
#echo exePrefix: "$exePrefix"	# check that variable is not split on space


  

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
###				then  can juast put the variable in the command and it will be empy or not
###				need an extra flag in the wrapper for sampleId e.g. -e, also -d flag needs to have a dataOrigin value.
		export PYTHONPATH=$HOME/lib/python
		paftools addPaftolFastq  $unzippedR1FastqFile  $unzippedR2FastqFile \
		--fastqPath $paftolDataSymlinksDir \
		--dataOrigin $usePaftolDb \
		--sampleId $externalSequenceID
		# NB - the file name must be of this format e.g. PAFTOL_005853_R1.fastq
		# The fastqPath entered in the database consists of the path and filename e.g. $paftolDataSymlinksDir/$unzippedR1FastqFile
		# --sampleId=$sampleId - must be included for all data types, except paftol (but there is no harm in always including it)

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
	print "ERROR: the Hyb-Seq program was not recognised. The options are 'paftools' or 'hybpiper' (without the quotes)."
	exit
fi

cd ../ # Back up to parent dir for next sample - 20.4.2020 - has no effect here now and not required anymore because looping through samples is done outside this script 
echo
echo