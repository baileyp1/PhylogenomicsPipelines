#!/bin/bash

############################
# various_tasks_in_python.sh
#
# Author: Paul Bailey
#
# Copyright (c) 2025 The Board of Trustees of the Royal Botanic Gardens, Kew
#
# Purpose: code to perform various (routine) tasks with Bash (with explanations for learning purposes).
#		   Tasks fall within discrete subroutines - see the documentation within each method for more info 
#
# Usage: various_tasks_in_bash.py <method_name> <option1 e.g. infile>  <option2 e.g. outfile_prefix>  <option3>  <option4> etc
#
# Procedure for adding a subroutine:
# 1. write a brief method decription here below
# 2. Add extra $1-9 variables to an option* variable, as required
# 3. Write the subroutine, placing it at the end of all the subroutines
# 4. Copy the brief method description to within the subroutine itself and extend the documentation there
# 4. Add a clause in the main code section for the method - see bottom of this file
#
# wget_sra_download():
#	Downloads SRA fastq files 
#
# retrieve_targets():
#	retrieves gene orthologs corresponding to a set of targets genes from gene coding sequences or a transcriptome assembly 
#
# *** Next method here ***
# 
############################

if [ "$#" -lt 1 ]; then
	echo 'ERROR: you need to specify a bash function to use!'
	echo 'List of main functions available:'
	echo '1. wget_sra_download'
	echo '2. retrieve_targets'
	exit 1
fi

method=$1	# method name

if [ "$#" -ge 2 ]; then
	option1=$2	# Often the main infile	
else
	echo "ERROR: having no infile or some input is set up to exit with an error at the moment"
	exit 1
fi

option2=''					# Defined these vars here so I can test whether option2 onwards is empty or not; actually not required to be predefined in bash
option3=''					
option4=''
option5=''
option6=''
option7=''
if [[ "$#" -ge 3 ]]; then	# Testing whether there are 3 or more paramters (otherwise script crashes - not tested this bash code here)
	option2=$3
fi
if [[ "$#" -ge 4 ]]; then	
	option3=$4
fi
if [[ "$#" -ge 5 ]]; then	
	option4=$5
fi
if [[ "$#" -ge 6 ]]; then	
	option5=$6
fi
if [[ "$#" -ge 7 ]]; then	
	option6=$7
fi


###################
wget_sra_download()	{
###################
# 	Purpose: Downloads SRA fastq files - pair end or single end - by manually creating the URL from the predictable 
#            FTP structure of the fastq file location in the SRA at ENA
#	         e.g. URL example: http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/006/SRR1653336/SRR1653336.fastq.gz 
# 			 Formats: 1. URL example: http://ftp.sra.ebi.ac.uk/vol1/fastq/<accession-prefix>/<00-last-digit-of-full-accession>/<full-accession>/
#			          2. fastq file name: <accession_number>_[12].fastq.gz
#
#			 Takes about 4 minutes to download each fastq file. If fastq file already exists, it will be renamed and a fresh download started.
#			 If fastq file(s) are downloaded unsucessfully (incomplete file), the wget log file hangs around showing the number of bytes downloaded. 
#			
#
# 	Usage: 1. To download a single accession:
#			  various_tasks_in_bash.sh wget_sra_download <accession_number> > wget_SRA_download.log 2>&1 &
#			  e.g. various_tasks_in_bash.sh wget_sra_download SRR14570809 > wget_SRA_download.log 2>&1 &
#
#		   2. To download multiple accessions together:
#			  cat <accn_list_file> | while read accn; do \
#				various_tasks_in_bash.sh wget_sra_download $accn > wget_SRA_download.log 2>&1
#			  done

### To do - UPTOHERE 10.2.2025
### READY to add method to recovery pipeline
### --> Review how SRA fastq files are removed afterwards - is this already set up in the main pipeline?
### --> Could have an option to download fastq files in main pipeline then exit!!
### Add documentation for this function in the PP readme and cmd line help; mention in the README that you need wget installed
###		--> THEN merge this various_tasks_in_bash.sh code with script of same name in PhyloPipeline repo and commit to PP repo
### Look at wget to see whether it can assess download - doesn't seem to but could look at the exit signal 
#			  

# 	Usage example (via Slurm): sbatch -J sra_download -p long -c 1 -n 1  --mem=5000  -o wget_SRA_download.log  -e wget_SRA_download.err  --wrap  "
#							   various_tasks_in_bash.sh wget_sra_download <accn_list_file> "
#							   NB - probably best not to parallelize these downloads!
#	Note: 1. it is strongly recommended to downoad files into a fresh/empty directory because the script tests whether files are absent or not
#		  2. the wget_SRA_download_[PE|SE]_fastqs.log files are appended to
#		  3. There is a change that what I think in the log e.g. [157921666/157921666] is bytes downloaded out of the total is not what I think it is.
#		  		BUT the info for the remote file size does seem available but only in the verbose wget output but I just haven't been able to test I can use it.
#				e.g.:
#			    wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14570809/SRR14570809_1.fastq.gz
#				--2025-02-08 13:45:05--  http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/009/SRR14570809/SRR14570809_1.fastq.gz
#				Output is:
#					Resolving ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)... 193.62.193.165
#					Connecting to ftp.sra.ebi.ac.uk (ftp.sra.ebi.ac.uk)|193.62.193.165|:80... connected.
#					HTTP request sent, awaiting response... 200 OK
#					Length: 157921666 (151M) [application/x-gzip]			<-- HERE is the file size in bytes
#					Saving to: 'SRR14570809_1.fastq.gz.2'
#			 I could try to use the wget --spider to get the info first without downloading if it turn out my way of checking doesn't work i.e. using e.g. [157921666/157921666]:
#			 Onlione e.g.'s:
#			 size_bytes=$(wget -S "${url}" --start-pos=500G 2>&1 | grep Content-Length | cut -d: -f2)
#			 Try: wget --spider --force-html -i bookmarks.html
#
# 	Notes:
#   1. If wget finds an existing file it will by default append .1, .2 etc to the fastq file and the subsequent file checks will not work.
#	   So this script will rename the file existing file so this doesn't happen
#	   Note: could look at wget options e.g. -r -p -nc - which can overwrite a file - but what I'm doing is probably OK as it leaves the
#	   old file around for comparison
#	2. The wget_SRA_download.log file shows in square brackets how many bytes were downloaded out of how many e.g. [157921666/157921666] 
#	   Matches the downloaded file size so using these values to check download success until I know how to download the md5sum for each accession
#	3. The URLs used here will download archive-generated fastq files from ENA but sometimes these files are not available for an accession.
#	   Haven't noticed this yet
#	   Alternative locations are:
#	   1. for submitted reads files:
#		  ftp://ftp.sra.ebi.ac.uk/vol1/run/<accession-prefix>/<full-accession>/
#		  e.g. ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR164/ERR164407/
#	   2. SRA read files - I think these require using the SRA Toolkit
#
#	   However, so far the problem has been because the inner folder doesn't exist or is different to described
#      i.e. <00-last-digit-of-full-accession> is not always true e.g. for SRR10237411:
#	   e.g.1 wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/011/SRR10237411/SRR10237411_1.fastq.gz - works
#		    wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR102/001/SRR10237411/SRR10237411_1.fastq.gz - doesnt work!
#	   So an extra search level needs to be added to find a PE read1 file, an SE file, then PE read 2 file
#      e.g.2 accession fastq file exists if you miss out the 00[0-9] folder level:
#	   		wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098689/SRR098689.fastq.gz - works
#	   		wget http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/009/SRR098689/SRR098689.fastq.gz - doesn t work
#	4. The downloaded fastq file has the original time stamp from ENA.
#	

#	ENA documentation:
#	https://www.ebi.ac.uk/ena/browser/downloading-data
#	https://ena-docs.readthedocs.io/en/latest/retrieval/file-download/sra-ftp-structure.html
# 	https://ena-docs.readthedocs.io/en/latest/faq/archive-generated-files.html#archive-generated-files
###################

	echo
	echo "Command: $0 $@"
	echo

	accn=$1
	outerFolder=`echo $accn | awk '{ print substr($0, 1, 6) }' `
	innerFolder00X=`echo $accn | awk '{ print "00" substr($0, length($0)) }' `	# The inner folder is defined in the ENA docs as being 00 followed by the last accession digit,
	innerFolder0XX=`echo $accn | awk '{ print "0" substr($0, length($0)-1) }' `	# BUT the inner folder is sometimes the last two digits of the accession number 
	echo INFO: $accn
	echo INFO: outerFolder: $outerFolder
	echo INFO: "Inner_folder_(starting_assumption): $innerFolder00X"
	# URL example: http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/006/SRR1653336/SRR1653336.fastq.gz
	urlPrefix='http://ftp.sra.ebi.ac.uk/vol1/fastq'
	
	if [[ -s ${accn}_1.fastq.gz ]]; then # if file already exists before attempting download, rename it before redownloading
		mv ${accn}_1.fastq.gz ${accn}_1.fastq.gz0
		echo "WARNING: ${accn}_1.fastq.gz file already exists, now renamed it to ${accn}_1.fastq.gz0 and will download afresh."
	fi

	echo URL_link: $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}_1.fastq.gz
	time wget --no-verbose -o ${accn}_1.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}_1.fastq.gz
	if [[ -s ${accn}_1.fastq.gz ]]; then 
		fastq_integrity_test  ${accn}_1.fastq.gz  ${accn}_1.fastq_wget_SRA_download.log
	else
		echo "WARNING: unable find an R1 fastq file name for accession $accn in the default location, will try this location:"
		# $innerFolder is sometimes the last two digits of the accession number:
		echo URL_link: $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}_1.fastq.gz
		time wget --no-verbose -o ${accn}_1.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}_1.fastq.gz
		if [[ -s ${accn}_1.fastq.gz ]]; then 
			fastq_integrity_test  ${accn}_1.fastq.gz  ${accn}_1.fastq_wget_SRA_download.log
		else
			echo "WARNING: still unable find an R1 fastq file name for accession $accn. Will remove the inner folder (absent in older (?) accessions):" 
			# URL example: http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR098/DRR098846/DRR098846_1.fastq.gz
			echo URL_link: $urlPrefix/$outerFolder/$accn/${accn}_1.fastq.gz
			time wget --no-verbose -o ${accn}_1.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$accn/${accn}_1.fastq.gz
			if [[ -s ${accn}_1.fastq.gz ]]; then
				fastq_integrity_test  ${accn}_1.fastq.gz  ${accn}_1.fastq_wget_SRA_download.log
			else
				echo "WARNING: unable find an R1 fastq file name for accession $accn. Will attempt to download a single end file with this filename format: ${accn}.fastq.gz"
				echo URL_link: $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}.fastq.gz
				time wget --no-verbose -o ${accn}.fastq_wget_SRA_download.log  $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}.fastq.gz
				if [[ -s ${accn}.fastq.gz ]]; then
					fastq_integrity_test  ${accn}.fastq.gz  ${accn}.fastq_wget_SRA_download.log
				else
					echo "WARNING: unable find an R1 fastq file name for accession $accn. Will attempt to download a single end file with this filename format: ${accn}.fastq.gz"
					echo URL_link: $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}.fastq.gz
					time wget --no-verbose -o ${accn}.fastq_wget_SRA_download.log  $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}.fastq.gz
					if [[ -s ${accn}.fastq.gz ]]; then
						fastq_integrity_test  ${accn}.fastq.gz  ${accn}.fastq_wget_SRA_download.log
					else
						echo "WARNING: Still unable to find an R1 fastq filename or a single end filename for accession $accn in the default location, will try this location:"
						# $innerFolder is sometimes absent in the accession URL (it seems to be for older accessions):
						echo URL_link: $urlPrefix/$outerFolder/$accn/${accn}.fastq.gz
						time wget --no-verbose -o ${accn}.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$accn/${accn}.fastq.gz
						if [[ -s ${accn}.fastq.gz ]]; then
							fastq_integrity_test  ${accn}.fastq.gz  ${accn}.fastq_wget_SRA_download.log
						else
							echo "ERROR: unable to find R1 fastq filename or a single end filename for accession $accn"
						fi
						exit # SE file so don t need to continue to get an R2 file
					fi
					exit
				fi
				exit
			fi
		fi
	fi

	if [[ -s ${accn}_2.fastq.gz ]]; then # if file already exists before attempting download, rename it before redownloading
		mv ${accn}_2.fastq.gz ${accn}_2.fastq.gz0
		echo "WARNING: ${accn}_2.fastq.gz file already exists, Have renamed it to ${accn}_2.fastq.gz0 and will download afresh."
	fi
	echo URL_link: $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}_2.fastq.gz
	time wget --no-verbose -o ${accn}_2.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$innerFolder00X/$accn/${accn}_2.fastq.gz
	if [[ -s ${accn}_2.fastq.gz ]]; then
		fastq_integrity_test  ${accn}_2.fastq.gz  ${accn}_2.fastq_wget_SRA_download.log
	else
		echo "WARNING: unable find an R2 fastq file name for accession $accn in the default location, will try this location:"
		# $innerFolder is sometimes the last two digits of the accession number:
		echo URL_link: $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}_2.fastq.gz
		time wget --no-verbose -o wget_SRA_download_PE_fastqs.log $urlPrefix/$outerFolder/$innerFolder0XX/$accn/${accn}_2.fastq.gz
		if [[ -s ${accn}_2.fastq.gz ]]; then
			fastq_integrity_test  ${accn}_2.fastq.gz  ${accn}_2.fastq_wget_SRA_download.log
		else
			echo "WARNING: still unable find an R2 fastq file name for accession $accn. Will remove the inner folder (absent in older (?) accessions):" 
			# URL example: http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR098/DRR098846/DRR098846_2.fastq.gz
			echo URL_link: $urlPrefix/$outerFolder/$accn/${accn}_2.fastq.gz
			time wget --no-verbose -o ${accn}_2.fastq_wget_SRA_download.log $urlPrefix/$outerFolder/$accn/${accn}_2.fastq.gz
			if [[ -s ${accn}_2.fastq.gz ]]; then
				fastq_integrity_test  ${accn}_2.fastq.gz  ${accn}_2.fastq_wget_SRA_download.log
			else
				echo "ERROR: unable to find an R2 fastq file name for accession. $accn"
			fi
		fi
	fi
}


######################
fastq_integrity_test()	{
######################
# 	Purpose: tests the integrity of a fastq file downloaded from ENA SRA.
#			
#	An internal function for wget_sra_download()
#
# 	Usage: fastq_integrity_test  <fastq_file_name>  <fastq_file_name>_wget_SRA_download.log>
#		   e.g. fastq_integrity_test  ${accn}_1.fastq.gz  ${accn}_1.fastq_wget_SRA_download.log 
#
# 	Input parameters: 
# 	$1 = <fastq_file_name>
#	$2 = <fastq_file_name>_wget_SRA_download.log
######################

	# Example wget output after download: 2024-11-13 12:33:24 URL:http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR098/DRR098846/DRR098846_1.fastq.gz [12964810382/12964810382] -> "DRR098846_1.fastq.gz" [1]
	# Identifying line with square brackets to get the number of bytes downloaded and bytes to download then storing the latter value.
	# NB - might not be using the perfect line identifier.
	bytesDownloadable=0
	bytesDownloadable=`tail -n 1 $2 |  grep  '[][] ->' | awk '{print $4}' | sed 's/[][]//g' | awk -F '/' '{print $2}' `
	fastqR1FileSize=0
	fastqR1FileSize=`ls -l $1 | awk '{print $5}' `
	if [[ $bytesDownloadable -gt 0 && $fastqR1FileSize -gt 0 && $bytesDownloadable -eq $fastqR1FileSize ]]; then
		echo "INFO: Download successful: $fastqR1FileSize of $bytesDownloadable bytes downloaded for $1"
		rm $2
	else
		echo "ERROR: Download unsuccessful: $fastqR1FileSize of $bytesDownloadable bytes downloaded for $1"
		# Log file stays around
	fi
}


#################
retrieve_targets()	{
#################
# 	Purpose: To retrieve a set of reference gene orthologs from gene DNA coding sequences (derived from an annotated genome)
#            or gene contigs from a transcriptome assembly
#
#	Description: Carries out the same task as Paftools retrieveTargets program, except it can use either blastn or tblastn to perform 
#				 the homologous searches, with tblastn being more sensitive.
#			
# 	Usage: export pathToScript=<path to required script, various_tasks_in_python.py # Work around for set up with PhylogenomicsPipelines respository 
#																					# so that this script can still be run as a standalone one 
#		   various_tasks_in_bash.sh retrieveTargets <fasta_file_of_reference_targets> <fasta_file_of_gene_sequences_to_search> \
#		   <sampleId> \
#		   <blast_program> \
#		   <residue_type_nucl_or_prot>

#          Example:
#		   export pathToScript=/Users/pba10kg/Documents/bin 
#		   various_tasks_in_bash.sh retrieve_targets \
# 		   /Users/pba10kg/Documents/workarea/1_data_references_etc/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected.fasta \
#		   GCA_024733475.1_NTU_Sgrande_1.0_cds_from_genomic_modified.fna \
#		   GCA_024733475.1 \
#		   tblastn \
# 		   nucl \
#		   4
# 		   > ${sampleId}_retrieve_targets.tblastn.log 2>&1 &
#
#	Subroutine input parameters: $1 = <fasta_file_of_reference_targets> 
#					  			 $2 = <fasta_file_of_gene_sequences_to_search>
#					  			 $3 = <sampleId>
#		   			  			 $4 = <blast_program>
#		   			  			 $5 = <residue_type_nucl_or_prot>
#		   			  			 $6 = <cpu> 
#					  			[also consider to add the post-BLAST filtering value being used here]
#
#	Additional software required:
#	BLAST
#	Exonerate (fastatranslate)
#	various_tasks_in_python.py retrieve_targets_magic 
# 		Python modules:
#		from Bio import SeqIO

#
#	Notes:
#	1. tblastn requires the query sequence input files to be translated
#   2. One issue with the approach is that only the 1st HSP of the top hit is evaluated which might be an issue if the gene 
#      translation goes out of frame. This should be OK but the top HSP might be shorter than another one but still have made 
#	   it to the top of the list e.g. if it just contains one domain of the gene. One solution is to use Exonerate or parse 
#	   the blast outputs to get the hit across the whole gene if present (Kevin's script does this for the unannotated gnomes; 
#	   also look at Captus extract)
#	3. Would it be worth to take the best of all 6 translations? Actually, no - tblastn will reveal hits from all frames
#	4. NB - Diamond sequence searcher only uses blastx and blastp - so maybe I should adapt for use with blastx
#	5. Could allow user to select the evalue and % id filters
#
#
#
### UPTOHERE 13.5.2025 - things to do
### Now check with Slurm on Gruffalo (Measure how much RAM is used)
	#####--exclude=node005,node010,node012,node002,node009,node007
### Also Test: ($hybSeqProgram != 'no' || $retrieveTargets != 'no') - line 240 in wrpper --> now check it's Ok for option -y
### Also there was an error on line ~592 w.r.t. --start-from-- - 
### Add stats sumHSPlength - done - now check

###	Also - Note to extract the busco genes and create an R plot of the pcIds - FIRST make a file with values ready for a histogram. (single column?)

### Added sort -4gr - make a note of this functionality in Linux notes

###	Prepare test data sets for PP repo
###	Once up to here report to Berta
### Also test with Angio353_v2 interim and mega353
### Repeat and complete stats +/- filtering, blast vs tblastn and with paftools and Captus

### Remove the files no longer required consider to put stats per line into the log file
### Also print out the by-row table header in csv format and put it in the log file
### ALSO, still need to check the fasta file header  - see lines ~438-441


### Next week Tue 13th May onwards / Future tasks for later:
### Add in my thoughts on how to get coord for all HSPs of a top hit - another dict is required for this I think 
###		{gene}{top query id}{genome contig}{lowest start} - start and end filled in in both parts of the condtional .. if .. and  else
###										   {highest end}
### 	Investigate Diamond to make it quicker - would enable testing on Macbook
### 	Use via HybPiper-2.2.0; can only use blastx so would have to transfer to use blastx - so delay for now
###		Probably is Ok but would need to check for hits in the other direction??? Should check for that anyway?

### Quickly investigate again to pull out genomes from ncbi genome website --> start method pseudocode to do this --> on hold 

### Try out Captus again - can use the extract command - could be method 2!!!

### Investigate Diamond to make it quicker - would enable testing on Macbook
### 	Use via HybPiper-2.2.0; can only use blastx so would have to transfer to use blastx - so delay for now
###		Probably is Ok but would need to check for hits in the other direction??? Should check for that anyway?

### Is it worth reporting multiple contigs for the same gene? Code is in place for that but need to finish it off:
### Need to create an array and append each row - actually, again just need to keep the best hit for each subj hit 
### This would be worthwhile doing as these would be the ||ogs!
###		See HybPiper option --paralog_min_length_percentage for length filter of 0.75

### There are 619 bad hits involving 73 genes so it shows that these ref targs do have significant homologs outside the 
### orthologous clade (they are < 55% identical so cannot be paralogs)
### more  ../testing_orig_angios353_my_subroutine/GCA_024733475.1_tblastn.sort-k13g_bad_hits.tab | awk '{print $1}' | awk -F '-' '{print $2}' | sort | wc -l
### THINK further - would I need to assess all contigs coming from each query seq - I think so - 
### makes it more difficult to code - at the moment it's by chance or not that I can see them with different query hits
### but it might not be accessing all the subject hits in the data if the paralog is never the top hit!!!

### Review whether any more input checks are required 

#################

	set -e
	set -u
	set -o pipefail
	shopt -s failglob

	echo "Command: $0 $@"
	echo
	sampleId=$3
	blastProgram=$4
	residueType=$5 	# nucl or prot (to agree with blast options)
	cpu=$6
	echo sampleId: $sampleId
	echo Blast program: $blastProgram
	echo residueType: $residueType
	# The sample fasta file needs to be copied to the sample folder, otherwise BLAST creates
	# the database files in the original location:
	cp $2 .
	geneSeqsToSearch=`basename $2`
	echo File of gene sequences to search: $geneSeqsToSearch
	echo

	# Check that the fasta file is zipped or not and decompress as required.
	# makeblastdb requires unzipped files.
	if [[ $geneSeqsToSearch == *'.gz' ]]; then
		geneSeqsToSearch=`basename -s .gz $2`
		gunzip -cf ${geneSeqsToSearch}.gz > $geneSeqsToSearch
	fi


	### NB - FIRST will need to examine the fasta header line to make sure that BLAST can use it:
	#1. pipe clean - just replace with an underscore
	#2. check fasta id has less than 50 chars for the main id!!! - see sygenium notes - exit if so with error
	#3. check all seqs have fasta records in them - what did I mean here


	echo "Making the BLAST db index of the gene sequences to search..."
	makeblastdb \
	-in $geneSeqsToSearch \
	-dbtype nucl \
	-parse_seqids
	# NB: -parse_seqids enables parsing of seqIds for retrieval by a blast utility
	# NBNB - maximum string length for the id = 50 (prior to first white space I think)
	# -dbtype nucl or prot \    - actually for tblastn you still need to create the index from DNA seqs
	# NB - when indexing aa seqs with -dbtype set to nucl, program doesn't crash, just ignores aa chars + blast search just gives zerto output. 

	# Translate the blast query sequences if using tblastn:
	queryFile=$1
	if [[ $blastProgram == 'tblastn' ]];then
		fastatranslate -F 1 $1 > ${sampleId}_queries.pep
		queryFile=${sampleId}_queries.pep
	fi

	echo "Blasting with $blastProgram..."
	$blastProgram -db $geneSeqsToSearch \
	-num_threads $cpu \
	-query $queryFile \
	-num_alignments 1 \
	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" \
	-out ${sampleId}_${blastProgram}.tab
	# Other parameters:
	# -num_alignments = Number of database sequences to show alignments for. Default = 250
	# -evalue = Expectation value (E) threshold for saving hits. Default = 10 - can also assess afterwards?
	#
	# BLAST query = from reference targets file
	# BLAST subject = from set of gene coding sequences or transcriptomes


	# Next sort by evalue across all query hits, then by % id, then by length
	# 'sort -k11g -k3gr [-k4gr - actually still to test]' does this which is pretty cool if you ask me - no need use use Python Pandas here
	# NB - a sequence can have evalue of 0.0 but pc id can still vary a lot so hits need to be sorted by % as well
	#	   so that the best subject sequence is chosen. It's less important for other evalues as they tend to be unique.  
	#
	# Once sorted, the top hit HSP for each gene will appear first in the list so if the same 
	# database sequence hits multiple genes, the best matching gene will be selected first.
	# Paftools retrieveTargets wasn't doing this for some reason.
	# Also need an evalue cut off: 0.05, 0.01 are options; hybpiper assemble uses 0.0001 for blastx hits, probably safer if we are demanding orthologs, not just homologs
	# Also using a cut off for pcid; hybpiper assemble exonerate cut off for exonerate hits = 55 - so use that for now.
	# Will also sort by length as well - would complete the logic of selecting best hit but it is very unlikely to result in any further improvement
	sort -k13g -k3gr -k4gr ${sampleId}_${blastProgram}.tab | awk '$13 <= 0.0001' | awk '$3 >= 55' > ${sampleId}_${blastProgram}.sort-k13g-k3gr-k4gr.tab
	
	# For testing with no filtering (also for testing +/- sorting):
	###sort -k13g -k3gr -k4gr ${sampleId}_${blastProgram}.tab > ${sampleId}_${blastProgram}.sort-k13g.tab 

	# Bad hits:
	# NB - don't need to use the '$11 > 0.0001' filter, only whether sequence is < 55% id, then we will also see hits with evalues > 0.0001.
	sort -k13g -k3gr -k4gr ${sampleId}_${blastProgram}.tab | awk '$3 < 55' > ${sampleId}_${blastProgram}.sort-k13g-k3gr-k4gr_bad_hits.tab
	

	# Now rearrange the fasta header of the gene coding or transcriptome sequence so the main id
	# corresponds to the reference target gene name from the query hit:
	# Usage: various_tasks_in_python retrieve_targets_magic  <blast_output.tab>  <fasta_file_for_blast_db>
	$pathToScript/various_tasks_in_python.py retrieve_targets_magic $sampleId ${sampleId}_${blastProgram}.sort-k13g-k3gr-k4gr.tab  $geneSeqsToSearch $blastProgram
	# Output file  of fasta records: <sampleId>.fasta

	# Check whether the same subject sequence is hitting the same reference target gene (incorrect orthology):
	numbrMultiGeneHits=`cat ${sampleId}.fasta | grep '>' | awk '{print $3}' | sort |uniq -c | awk '$1 > 1' | wc -l `
	#echo numbrMultiGeneHits: $numbrMultiGeneHits
	echo
	echo "INFO: BLAST hits filtered out with evalue of > 0.0001 and % id of < 55"
	echo
	if [[ $numbrMultiGeneHits -gt 0 ]]; then
		echo 'WARNING: Hits to the same contig (BLAST subject) found for more than one gene:'
		echo `cat ${sampleId}.fasta | grep '>' | awk '{print $3}' | sort | uniq -c | awk '$1 > 1' `
		# Now use printf to get a list of subj hits to use with grep:
		grepRegex=`cat ${sampleId}.fasta | grep '>' | awk '{print $3}' | sort |uniq -c | awk '$1 > 1' | awk '{printf $2 " "}' `
		echo $grepRegex | sed 's/ /\|/g'
		# Finally print out the fasta header for each gene affected
		grep "$grepRegex" ${sampleId}.fasta
	fi


	# Basic stats across all genes:
	numbrRecoveredGenes=`cat ${sampleId}.fasta | grep '>' | wc -l `
	sumLengthOfGenesWithNs=`fastalength ${sampleId}.fasta | awk '{sum+=$1} END {print sum}' `
	# Also removing strings of N's from the sequence line before counting the number of bases:
	cat ${sampleId}.fasta \
	| awk '{if($1 ~ /^>/) { print $0 } else { {gsub(/[Nn]/,"",$0)} {print $0} } }' \
	| grep -v ^$ \
	> ${sampleId}.fasta.Ns_removed_temp
	sumLengthOfGenes=`fastalength ${sampleId}.fasta.Ns_removed_temp | awk '{sum+=$1} END {print sum}' `
	rm ${sampleId}.fasta.Ns_removed_temp
	sumLengthOfHSPs=`cat ${sampleId}.fasta | grep '>' | awk '{print $5}' | sed 's/lenHSP=//' | awk '{sum+=$1} END {print sum}' `
	avPcIdAcrossTopHSP=`cat ${sampleId}.fasta | grep '>' | awk '{print $4}' | sed 's/pcid=//' | awk '{sum+=$1} END {if(sum > 0) {print sum/NR} else {print "0"}}' `
	minPcIdAcrossTopHSP=`cat ${sampleId}.fasta | grep '>' | awk '{print $4}' | sed 's/pcid=//' | sort -n | head -n 1 `
	maxPcIdAcrossTopHSP=`cat ${sampleId}.fasta | grep '>' | awk '{print $4}' | sed 's/pcid=//' | sort -n | tail -n 1 `
	# Average % coverage across top HSP against the query (reference target) genes:
	avPcHSPCovrgToQueryLen=`cat ${sampleId}.fasta | grep '>' | awk '{print $5 " " $6}' | sed 's/lenHSP=//' | sed 's/qlen=//' | awk '{print ($1 / $2) * 100}' | awk '{sum+=$1} END {if(sum > 0) {print sum/NR} else {print "0"}}' `
	medianPoint=`cat ${sampleId}.fasta | grep '>' | awk 'END {printf "%.0f", NR/2}' ` 
	medianPcHSPCovrgToQueryLen=`cat ${sampleId}.fasta | grep '>' | awk '{print $5 " " $6}' | sed 's/lenHSP=//' | sed 's/qlen=//' | awk '{print ($1 / $2) * 100}' | sort -n | head -n $medianPoint | tail -n 1 `
	minPcHSPCovrgToQueryLen=`cat ${sampleId}.fasta | grep '>' | awk '{print $5 " " $6}' | sed 's/lenHSP=//' | sed 's/qlen=//' | awk '{print ($1 / $2) * 100}' | sort -n | head -n 1 `
	# Will be some values > 100% if gaps are present in the HSP:
	maxPcHSPCovrgToQueryLen=`cat ${sampleId}.fasta | grep '>' | awk '{print $5 " " $6}' | sed 's/lenHSP=//' | sed 's/qlen=//' | awk '{print ($1 / $2) * 100}' | sort -n | tail -n 1 `
	# Average % coverage across the top query against all subject genes:
	avPcQueryCovrgToSubjectLen=`cat ${sampleId}.fasta | grep '>' | awk '{print $6 " " $7}' | sed 's/qlen=//' | sed 's/slen=//' | awk '{print ($1 / $2) * 100}' | awk '{sum+=$1} END {if(sum > 0) {print sum/NR} else {print "0"}}' `

echo "sampleId: $sampleId
numbrRecoveredGenes: $numbrRecoveredGenes
sumLengthOfGenesWithNs: $sumLengthOfGenesWithNs
sumLengthOfGenes: $sumLengthOfGenes
sumLengthOfHSPs: $sumLengthOfHSPs
avPcIdAcrossTopHSP: $avPcIdAcrossTopHSP
minPcIdAcrossTopHSP (min % allowed, 55%): $minPcIdAcrossTopHSP
maxPcIdAcrossTopHSP: $maxPcIdAcrossTopHSP
avPcHSPCovrgToQueryLen: $avPcHSPCovrgToQueryLen
medianPcHSPCovrgToQueryLen: $medianPcHSPCovrgToQueryLen
minPcHSPCovrgToQueryLen: $minPcHSPCovrgToQueryLen
maxPcHSPCovrgToQueryLen (might be > 100% if gaps present): $maxPcHSPCovrgToQueryLen
avPcQueryCovrgToSubjectLen: $avPcQueryCovrgToSubjectLen" > ${sampleId}_stats.txt

echo "$sampleId $numbrRecoveredGenes $sumLengthOfGenesWithNs $sumLengthOfGenes $avPcIdAcrossTopHSP $minPcIdAcrossTopHSP $maxPcIdAcrossTopHSP \
$avPcHSPCovrgToQueryLen $medianPcHSPCovrgToQueryLen $minPcHSPCovrgToQueryLen $maxPcHSPCovrgToQueryLen $avPcQueryCovrgToSubjectLen \
" > ${sampleId}_stats_by_row.txt


	# Remove the files no longer required:
	###GCA_024733475.1_NTU_Sgrande_1.0_cds_from_genomic_modified.fna.n* and GCA_024733475.1_NTU_Sgrande_1.0_cds_from_genomic_modified.fna
	###GCA_024733475.1_tblastn.tab
	###GCA_024733475.1_queries.pep
	### geneSeqsToSearch ???

}
#################


###########
# Main code
###########
if [[ $method == 'wget_sra_download' ]]; then
 	wget_sra_download $option1
elif [[ $method == 'retrieve_targets' ]]; then
 	retrieve_targets $option1 $option2 $option3 $option4 $option5 $option6

# elif method == '... - other methods here

else
	echo 'ERROR: you need to specify a correctly named bash function to use!'
	echo 'List of main functions available:'
	echo '1. wget_sra_download'
	echo '2. retrieve_targets'
	exit 1
fi
