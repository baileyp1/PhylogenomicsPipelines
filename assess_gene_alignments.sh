#!/bin/bash

###########################
# assess_gene_alignments.sh 

# Author:   Paul Bailey

# Purpose: to calcualte stats from the collective outputs of make_gene_trees.sh script

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
###########################
#set -e 			# NB - issue with creating ${fileNamePrefix}_summary_of_sample_quality.txt when using any of these set commands; trying to do my own checks anyway
#set -u
#set -o pipefail
shopt -s failglob


#echo Inside slurm_make_gene_trees.sh script:  
#echo SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID

fractnAlnCovrg=$1
fractnMaxColOcc=$2
fractnSpecies=$3
totalNumbrSamples=$4
fileNamePrefix=$5
geneFile=$6
treeTipInfoMapFile=$7		# NB - just testing if this filename was submitted, then will add tree tip info to species tree
option_u=$8
seqType=$9
alnFileForTreeSuffix=${10}

### Potential additional (hidden) parameters to set:
###minColOccToTolerate
###minSumMaxColOccPercentToTolerate=
#minSumOfContigLengthsToTolerate=		HIEDDNE vlue for the minSumOfContigLengthsToTolerate set it at '20'
#maxSumOfContigLengths=


# Convert $emptyMatchStateFractn and $fractnSpecies to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `
fractnSpecies_pc=`awk -v FRACTN=$fractnSpecies 'BEGIN{printf "%.0f", FRACTN * 100}' `
fractnMaxColOcc_pc=`awk -v FRACTN=$fractnMaxColOcc 'BEGIN{printf "%.0f", FRACTN * 100}' `
# Minimum number of samples to tolerate for including into Astral:
numbrSamplesThreshold=`awk -v FRACTN=$fractnSpecies -v numbrSamples=$totalNumbrSamples 'BEGIN{printf "%.0f", FRACTN * numbrSamples}' `
# Above awk code is zero proof - can have 0 * 100 - returns zero


echo fractnAlnCovrg to use: $fractnAlnCovrg
echo fractnMaxColOcc: $fractnMaxColOcc
echo fractnSpecies to use: $fractnSpecies
echo numbrSamples: $totalNumbrSamples
#echo numbrSamplesThreshold: $numbrSamplesThreshold
echo fileNamePrefix: $fileNamePrefix
echo geneFile: $geneFile
echo treeTipInfoMapFile: $treeTipInfoMapFile
echo option_u: $option_u
echo "seqType (in assess script): " $seqType
echo "alnFileForTreeSuffix (in assess script): " $alnFileForTreeSuffix


if [[ $geneFile != 'use_genewise_files' ]]; then 
	######################################
	# Sum length of all contigs per sample
	######################################
	echo "Sample NumbrGenes SumOfContigLengths" > ${fileNamePrefix}_summary_of_sample_quality.txt
	for file in *modified.fasta; do
		# 29.5.2020 - now taking the sample name from the fasta header of the *_modified.fasta sample-wise files.
		#sample=$(echo $file| sed 's/_modified.fasta//')
		sample=$(head -n 1 $file | awk '{print $2}')
		len_bp=$(fastalength $file  2>/dev/null | awk '{sum +=$1} END {print sum}')
		# Some files may be empty, if so print a zero:
		if [ -z $len_bp ]; then len_bp=0; fi 	# Need to test whether variable is empty when using set cmds; doesn't seem to work though, only works when all set cmd are removed!!
		numbrGenes=$(cat $file | grep '>' | wc -l)
		echo $sample " " $numbrGenes " " $len_bp
	done | sort -k3n \
	>> ${fileNamePrefix}_summary_of_sample_quality.txt
	#exit

	if [[ -s $treeTipInfoMapFile && $option_u == 'yes' ]]; then
		# Prepare the mapfile again to add in the sum length values onto the tree tips, if requested
		addTreeTipInfoFromTable_PlusL.py \
		tree_tip_info_mapfile.txt \
		${fileNamePrefix}_summary_of_sample_quality.txt \
		tree_tip_info_mapfile_plusL.txt
		# Now move new file to original mapfile created to avoid changing any other code:
		cp -p tree_tip_info_mapfile_plusL.txt  tree_tip_info_mapfile.txt
		### 27.5.2020 - copy not move for now to check outputs
	fi
	#exit
 
	maxSumOfContigLengths=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt | sort -k3n | tail -n 1 | awk '{print $3 " (sample " $1 ")" }' `
	minSumOfContigLengthsToTolerate=`echo $maxSumOfContigLengths | awk '{printf "%.0f", 0.20 * $1}' `

	# Now using median length to determine the 20% cut off values for poor quality samples to remove:
	medianPoint=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt | sort -k3n | awk 'END {printf "%.0f" , NR/2}' `
    medianGeneLength=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt | sort -k3n | awk '{print $3}' | head -n $medianPoint | tail -n 1 `
    minSumOfMedianGeneLengthsToTolerate=`echo $medianGeneLength | awk '{printf "%.0f", 0.20 * $1}' `
    # Now inserting above minSum into code just below.


	poorQualitySamples=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfMedianGeneLengthsToTolerate '$3 <= minSumOfContigLengthsToTolerate {print $0}' \
	| sort -k3n | wc -l | sed 's/ //g' `		# NB - Macbook Darwin wc -l print 7 space characters infront of the number. Just tidying up print with sed 's/ //g' in this script

	okQualitySamples=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfMedianGeneLengthsToTolerate '$3 > minSumOfContigLengthsToTolerate {print $0}' \
	| sort -k3n | wc -l | sed 's/ //g' `


	# If there are any poor quality samples, print their ids to a file (to be used for further investigation):
	if [ $poorQualitySamples -ge 1 ]; then
		tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
		| awk -v minSumOfContigLengthsToTolerate=$minSumOfMedianGeneLengthsToTolerate '$3 <= minSumOfContigLengthsToTolerate {print $1}' \
		> ${fileNamePrefix}_poor_quality_samples.txt
	fi

	# Now print a list of ok/good quality samples (to be used for repeating the species trees):
	tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfMedianGeneLengthsToTolerate '$3 > minSumOfContigLengthsToTolerate {print $1}' \
	> ${fileNamePrefix}_ok_quality_samples.txt
fi





###################################
# Amount of overlap between samples
###################################
# Sum Length of the longest gene seq in each alignment (AFTER filtering for maxColOcc > 120/150 and >3 seqs in file):
sumLenLongestGene=`cat *_aln_summary.txt | grep '^lenLongestGene:' | awk '{sum+=$2} END {print sum}' `

# Sum Length of the longest gene seq in each alignment after trimming to remove very rare insertions:
sumLenLongestGeneAfterTrim=`cat *_aln_summary.txt | grep '^lenLongestGeneAfterTrim:' | awk '{sum+=$2} END {print sum}' `

# Sum of MedianGeneLength in each alignment (before filtering/trimming):
sumMedianGeneLength=`cat *_aln_summary.txt | grep '^medianGeneLength:' | awk '{sum+=$2} END {print sum}' `
### 24.8.2020 - not prited this out to summary file yet


# Sum length of all columns with maxColOcc (same for all values of $fractnAlnCov) 
# NB - was filtering on min maxColOcc ($1) but now removed - filtering is inlcuded in the $numbrOverlapColsForSpeciesTree value.
sumMaxColOcc=`cat *_aln_summary.txt | grep '^maxColOcc:' | awk '{sum+=$2} END {print sum}' `


# Sum length of all parsimonious columns with maxColOcc (same for different values of $fractnAlnCov)
																	### 29.9.2020 - Not using this var now
sumMaxParsCols=`cat *_aln_summary.txt | grep '^maxParsCols:' | awk -v sumMaxColOcc=$sumMaxColOcc '{sum+=$2} END {print sum}' `
  

# Total number of residues in the maxColOcc area for all samples (prior to any filtering):
totalResiduesInMaxColOccRegion=`cat *_aln_summary.txt | grep '^totalResiduesInMaxColOccRegion:' | awk '{sum+=$2} END {print sum}' `


# Total number of ALL residues in the aln for tree (i.e. after filtering):
totalBasesInAlnForTree=`for file in *.$alnFileForTreeSuffix; do fastalength $file 2>/dev/null ; done | awk '{sum+=$1} END {print sum}' `


if [[ $seqType == 'protein' ]]; then
	totalNumbrSeqs=`cat  *.protein.fasta | grep '>' | wc -l `
	totalNumbrSeqsWithStops=`cat *.protein.log | grep '^TotalNumbrSeqsWithStops:' | awk '{sum+=$2} END {print sum}' `
	echo totalNumbrSeqs:  $totalNumbrSeqs
	echo totalNumbrSeqsWithStops: $totalNumbrSeqsWithStops
	if [[ $totalNumbrSeqs -eq 0 || $totalNumbrSeqsWithStops -eq 0 ]]; then
		percentSeqsWithStops=0
		echo percentSeqsWithStops: $percentSeqsWithStops
	else
		percentSeqsWithStops=`awk -v totalNumbrSeqsWithStops=$totalNumbrSeqsWithStops -v totalNumbrSeqs=$totalNumbrSeqs 'BEGIN{printf "%.2f", (totalNumbrSeqsWithStops/totalNumbrSeqs) * 100}' `
		echo percentSeqsWithStops1: $percentSeqsWithStops
	fi
fi


echo 'Table for gene alignments showing numbers per gene for:
1. longest recovered gene length (LenLongestGene)
2. longest recovered gene length after trimming to remove very rare insertions (LenLongestGeneTrim)
3. length of common overlap (MinColOcc, table ordered on this value, shortest to longest).
4. Number of parsimonious columns in length of common overlap (ParsimCols)
5. number of samples (NumberSamples, after filtering sequences by coverage)' > ${fileNamePrefix}_summary_gene_recovery.txt
echo 'GeneId LenLongestGene LenLongestGeneTrim MedianGeneLength MinColOcc ParsimCols Ratio:Col4/Col5 NumberSamples' | column -t >> ${fileNamePrefix}_summary_gene_recovery.txt
for file in *.$alnFileForTreeSuffix; do
 	gene=`echo $file | sed "s/.$alnFileForTreeSuffix//" `
 	numbrSamples=`cat $file | grep '>' | wc -l `
	lenLongestGene=`fastalength  $file 2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `				                 # The longest recovered gene in the aln
    maxColOcc=`fastalength  ${gene}_${seqType}_aln_AMAS_trim_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `     # minimum column occupancy (aln columns) - NB 24.8.2020 - can't i also get this from the summary.txt file??? And other values - woudl be easier
    maxColOccP=`fastalength  ${gene}_${seqType}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' ` # minimum column occupancy (parsimonious sites ONLY)
	lenLongestGeneAfterTrim=`fastalength ${gene}.$alnFileForTreeSuffix 2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `
	medianGeneLength=`cat ${gene}_aln_summary.txt | grep '^medianGeneLength:' | awk '{print $2}' `
	ratio=`echo  $medianGeneLength  $maxColOcc | awk '{printf "%.2f", $1/$2}' `		# Might be an indicator of overall effectiveness of gene recovery for samples submitted
	### 27.8.2020 - I think this is sometiems divisible by zero which is fatal!!!!
    echo "$gene $lenLongestGene $lenLongestGeneAfterTrim $medianGeneLength $maxColOcc $maxColOccP $ratio $numbrSamples"| column -t
done | sort -k5n| column -t >> ${fileNamePrefix}_summary_gene_recovery.txt
### NB - 12.10.2020 - the lenLongestGeneAfterTrim is now same as lenLongestGene!!!! Can get lenLongestGene from the summary.txt file which is before trimming!





# Number of gene alignments containing $fractnSpecies_pc % of samples (NB - also after trimming for rare insertions):
numbrGeneAlnsForSpeciesTree=`for file in *.$alnFileForTreeSuffix; do
 	numbrSamples=$(cat $file | grep '>' | wc -l)
    echo $numbrSamples
done \
| awk -v numbrSamplesThreshold=$numbrSamplesThreshold  '$1 >= numbrSamplesThreshold' \
| wc -l | sed 's/ //g'`


# Total number of alignment columns in the area of common overlap for these genes (NB - AFTER all filtering and trimming, also after trimming for rare insertions):
numbrOverlapColsForSpeciesTree=`for file in *.$alnFileForTreeSuffix; do
	gene=$(echo $file | sed "s/.$alnFileForTreeSuffix//")
 	numbrSamples=$(cat $file | grep '>' | wc -l)
    #lenLongestGene=$(fastalength  $file | sort -n | tail -n 1 | awk '{print $1}')				                 		   # The longest recovered gene in the aln
    maxColOcc=$(fastalength  ${gene}_${seqType}_aln_AMAS_trim_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}')   # maximum column occupancy (aln columns)
    echo $numbrSamples  " " $maxColOcc
    ####echo 2000  "test" $maxColOcc
done \
| awk -v numbrSamplesThreshold=$numbrSamplesThreshold  '$1 >= numbrSamplesThreshold {sum+=$2} END {print sum}' `
### NBNB - 11.10.2020 - I think I DON'T also need to filter by maxColOcc because these alignment files 
###        only exist if they pass maxColOcc filtering! Keeping checking logic though




### 20.8.2020 - Need to re-work this section 
### Still get the same stats but better now to redo the protein lists in the species tree script i think - the file might not be used
### 22.5.2021 - Also, if fasta file has < 4 seqs it should not enter this list beacuse trees are not beign made !!!!

################################################################################################################################
# Create list of gene alns with > ${fractnAlnCovrg_pc} gene coverge AND containing more than the $fractnSpecies % of the samples
# Used to concatenate gene alns for supermatrix methods AND to gather the stats below
# NB - also checking whether ther are > 3 species in each gene tree (this check was already made before building each gene tree
#      but the alignment file still exists so is removed here below.
################################################################################################################################
for file in *.$alnFileForTreeSuffix; do
 	gene=`echo $file | sed "s/.$alnFileForTreeSuffix//" `
 	numbrSamples=`cat $file | grep '>' | wc -l `;
    echo $gene " " $numbrSamples
done \
| awk -v alnFileForTreeSuffix=$alnFileForTreeSuffix -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc}  '$2 >= numbrSamplesThreshold && $2 > 3 {print $1 "." alnFileForTreeSuffix}' \
> mafft_dna_alns_fasta_file_list.txt

# For protein file list:
cat mafft_dna_alns_fasta_file_list.txt | sed 's/dna/protein/' > mafft_protein_alns_fasta_file_list.txt


# Count the number of samples in the gene trees being used and place a sorted list into a file:
### 6.9.2020 - if i now do this on the aln_for_tree.fasta file, then I don't need to generate the file list above here in this script!!
numbrSamplesInTrees=`cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sort -u | wc -l | sed 's/ //g'`
cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sed 's/>//' | sort -u > mafft_dna_alns_fasta_samples_in_tree.txt
### NB - 11.10.2020 - what happens if this is protein trees or if dna is not beign used?!!


# Print out a list of all samples submitted to the whole workflow:
#cat $geneFile | awk '{print $2}' | grep -v  '^$'| sort -u > samples_submitted.txt
# Altered above line to now use a set of files that is common to both species-wise and dna-wise datasets, i.e. *_dna.fasta:
cat *_dna.fasta | grep '>' | sed 's/^>//' | sort -u > samples_submitted.txt


# Number of samples not present in ASTRAL or RAxML trees:
samplesNotInSpeciesTrees=$(comm -23 samples_submitted.txt  mafft_dna_alns_fasta_samples_in_tree.txt | wc -l | sed 's/ //g')





echo "##########################################
Summary statistics for the gene alignments
##########################################

Filtering parameters set:
fractnAlnCovrg: $fractnAlnCovrg (user set)
fractnMaxColOcc: $fractnMaxColOcc (fixed)
fractnSpecies: $fractnSpecies (user set)

Total number of samples: $totalNumbrSamples" > ${fileNamePrefix}_summary_stats.txt

if [[ $geneFile != 'use_genewise_files' ]]; then 
	echo "
Largest sum length of all genes in bases for a sample found (that sample in brackets): $maxSumOfContigLengths
Median sum length of all genes per sample in bases: $medianGeneLength
Number of poor quality samples (<= 20 % of the median sum length of all genes per sample): $poorQualitySamples
Number of ok/good quality samples (> 20 % of the median sum length of all genes per sample): $okQualitySamples" >> ${fileNamePrefix}_summary_stats.txt
fi

echo "
Sequence type assessed: $seqType
Sum Length of the longest gene sequence in each alignment: $sumLenLongestGene
Sum Length of the longest gene sequence in each alignment after trimming rare inserts: $sumLenLongestGeneAfterTrim
Sum Length of the median gene sequence in each alignment: $sumMedianGeneLength
Total number of homologous columns occupied by >= $fractnMaxColOcc_pc % of samples: $sumMaxColOcc (area of common overlap)
Total number of homologous + parsimonious columns occupied by >= $fractnMaxColOcc_pc % of samples: $sumMaxParsCols
Total number of residues in area of common overlap: $totalResiduesInMaxColOccRegion " >> ${fileNamePrefix}_summary_stats.txt

if [[ $seqType == 'protein' ]]; then
	echo " 
Number of sequences with STOP codons (total number of seqs): $totalNumbrSeqsWithStops ($totalNumbrSeqs) " >> ${fileNamePrefix}_summary_stats.txt
fi

echo "
Number of gene alignments containing >= $fractnSpecies_pc % of the samples after filtering by coverage (will be used in the species tree): $numbrGeneAlnsForSpeciesTree
Total number of alignment columns in the area of common overlap for these genes (after any gene filtering): $numbrOverlapColsForSpeciesTree
Total number of residues in all these genes (after any gene filtering): $totalBasesInAlnForTree

After filtering genes by coverage AND/OR by number of samples in each gene tree:
Number of samples that will be present in final ASTRAL and RaxML trees (total samples submitted): $numbrSamplesInTrees ($totalNumbrSamples)
Number of samples that will NOT be in the species trees: $samplesNotInSpeciesTrees
" >> ${fileNamePrefix}_summary_stats.txt


if [ $numbrSamplesInTrees -ne $totalNumbrSamples ]; then
	echo "These samples will not be present in final species tree(s): "  >> ${fileNamePrefix}_summary_stats.txt
	comm -23 samples_submitted.txt  mafft_dna_alns_fasta_samples_in_tree.txt >> ${fileNamePrefix}_summary_stats.txt
	echo >> ${fileNamePrefix}_summary_stats.txt
fi

echo 'Number times each sample will appear in the gene trees (Format: number_of_gene_trees sample).' >> ${fileNamePrefix}_summary_stats.txt
echo '(Low numbers suggest low gene recovery for that sample)' >> ${fileNamePrefix}_summary_stats.txt
cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sed 's/>//' | sort | uniq -c | sort -k1n >> ${fileNamePrefix}_summary_stats.txt


rm samples_submitted.txt mafft_dna_alns_fasta_samples_in_tree.txt





