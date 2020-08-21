#!/bin/bash

###########################
# assess_gene_alignments.sh 

# Author:   Paul Bailey

# Purpose: to calcualte stats from the collective outputs of make_gene_trees.sh script

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
treeTipInfoMapFile=$7
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


#echo fractnAlnCovrg to use: $fractnAlnCovrg
#echo fractnSpecies to use: $fractnSpecies
#echo numbrSamples: $totalNumbrSamples
#echo numbrSamplesThreshold: $numbrSamplesThreshold
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


	poorQualitySamples=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfContigLengthsToTolerate '$3 <= minSumOfContigLengthsToTolerate {print $0}' \
	| sort -k3n | wc -l | sed 's/ //g' `		# NB - Macbook Darwin wc -l print 7 space characters infront of the number. Just tidying up print with sed 's/ //g' in this script

	okQualitySamples=`tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfContigLengthsToTolerate '$3 > minSumOfContigLengthsToTolerate {print $0}' \
	| sort -k3n | wc -l | sed 's/ //g' `


	# If there are any poor quality samples, print their ids to a file (to be used for further investigation):
	if [ $poorQualitySamples -ge 1 ]; then
		tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
		| awk -v minSumOfContigLengthsToTolerate=$minSumOfContigLengthsToTolerate '$3 <= minSumOfContigLengthsToTolerate {print $1}' \
		> ${fileNamePrefix}_poor_quality_samples.txt
	fi

	# Now print a list of ok/good quality samples (to be used for repeating the species trees):
	tail -n+2 ${fileNamePrefix}_summary_of_sample_quality.txt \
	| awk -v minSumOfContigLengthsToTolerate=$minSumOfContigLengthsToTolerate '$3 > minSumOfContigLengthsToTolerate {print $1}' \
	> ${fileNamePrefix}_ok_quality_samples.txt
fi





###################################
# Amount of overlap between samples
###################################
# Sum Length of the longest gene seq in each alignment (AFTER filtering for maxColOcc > 120/150 and >3 seqs in file):
sumLenLongestGene=`cat *_aln_summary.log | grep '^lenLongestGene:' | awk '{sum+=$2} END {print sum}' `

# Sum Length of the longest gene seq in each alignment after trimming to remove very rare insertions:
sumLenLongestGeneAfterTrim=`cat *_aln_summary.log | grep '^lenLongestGeneAfterTrim:' | awk '{sum+=$2} END {print sum}' `


# Sum length of all columns with maxColOcc (same for all values of $fractnAlnCov) 
# NB - still need to filter on min maxColOcc to tolerate.
sumMaxColOcc=`cat *_aln_summary.log | grep '^maxColOcc:' | awk '$2 >= 90 {sum+=$2} END {print sum}' `


# Sum length of all parsimonious columns with maxColOcc (same for different values of $fractnAlnCov) 
sumMaxParsCols=`cat *_aln_summary.log | grep '^maxParsCols:' | awk '$2 >= 90 {sum+=$2} END {print sum}' `
### NB - should not need to filter on min maxParsCols to tolerate BUT
### there is one value of 144 which I don't understand - they should arleady be filtered   


# Total number of bases in the maxColOcc area for all samples:
totalBasesInMaxColOccRegion=`for file in *.$alnFileForTreeSuffix; do fastalength $file 2>/dev/null ; done | awk '{sum+=$1} END {print sum}' `
										 ### 20.8.2020 - was file before trim 0.003 


echo 'Table for gene alignments showing numbers per gene for:
1. longest recovered gene length (LenLongestGene)
2. longest recovered gene length after trimming to remove very rare insertions (LenLongestGeneTrim)
3. length of common overlap (MinColOcc, table ordered on this value, shortest to longest).
4. Number of parsimonious columns in length of common overlap (ParsimCols)
5. number of samples (NumberSamples, after filtering sequences by coverage)' > ${fileNamePrefix}_summary_gene_recovery.txt
echo 'GeneId LenLongestGene LenLongestGeneTrim MinColOcc ParsimCols Ratio:Col3/Col4 NumberSamples' | column -t >> ${fileNamePrefix}_summary_gene_recovery.txt
for file in *.$alnFileForTreeSuffix; do
 	gene=`echo $file | sed "s/.$alnFileForTreeSuffix//" `
 	numbrSamples=`cat $file | grep '>' | wc -l `
	lenLongestGene=`fastalength  $file 2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `				                 # The longest recovered gene in the aln
    maxColOcc=`fastalength  ${gene}_${seqType}_aln_AMAS_trim_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `     # minimum column occupancy (aln columns)
    maxColOccP=`fastalength  ${gene}_${seqType}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' ` # minimum column occupancy (parsimonious sites ONLY)
	lenLongestGeneAfterTrim=`fastalength ${gene}.$alnFileForTreeSuffix 2>/dev/null | sort -n | tail -n 1 | awk '{print $1}' `
	ratio=`echo  $lenLongestGeneAfterTrim  $maxColOcc | awk '{printf "%.2f", $1/$2}' `								         # Might be an indicator of overall effectiveness of gene recovery for samples submitted
    echo "$gene $lenLongestGene $lenLongestGeneAfterTrim $maxColOcc $maxColOccP $ratio $numbrSamples"| column -t
done | sort -k4n| column -t >> ${fileNamePrefix}_summary_gene_recovery.txt





# Number of gene alignments containing $fractnSpecies_pc % of samples (NB - also after trimming for rare insertions):
numbrGeneAlnsForSpeciesTree=`for file in *.$alnFileForTreeSuffix; do
 	numbrSamples=$(cat $file | grep '>' | wc -l)
    echo $numbrSamples
done \
| awk -v numbrSamplesThreshold=$numbrSamplesThreshold  '$1 >= numbrSamplesThreshold' \
| wc -l | sed 's/ //g'`


# Total number of alignment columns in the area of common overlap for these genes (NB - also after trimming for rare insertions):
numbrOverlapColsForSpeciesTree=`for file in *.$alnFileForTreeSuffix; do
	gene=$(echo $file | sed "s/.$alnFileForTreeSuffix//")
 	numbrSamples=$(cat $file | grep '>' | wc -l)
    #lenLongestGene=$(fastalength  $file | sort -n | tail -n 1 | awk '{print $1}')				                 		   # The longest recovered gene in the aln
    maxColOcc=$(fastalength  ${gene}_${seqType}_aln_AMAS_trim_${fractnMaxColOcc}.fasta  2>/dev/null | sort -n | tail -n 1 | awk '{print $1}')   # maximum column occupancy (aln columns)
    echo $numbrSamples  " " $maxColOcc
    ####echo 2000  "test" $maxColOcc
done \
| awk -v numbrSamplesThreshold=$numbrSamplesThreshold  '$1 >= numbrSamplesThreshold {sum+=$2} END {print sum}' `





### 20.8.2020 - Need to re-work this section 
### Still get the same stats but better now to redo the protein lists in the species tree script i think 

################################################################################################################################
# Create list of gene alns with > ${fractnAlnCovrg_pc} gene coverge AND containing more than the $fractnSpecies % of the samples
# Used to concatenate gene alns for supermatrix methods AND to gather the stats below
################################################################################################################################
for file in *.$alnFileForTreeSuffix; do
 	gene=`echo $file | sed "s/.$alnFileForTreeSuffix//" `
 	numbrSamples=`cat $file | grep '>' | wc -l `;
    echo $gene " " $numbrSamples
done \
| awk -v alnFileForTreeSuffix=$alnFileForTreeSuffix -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc}  '$2 >= numbrSamplesThreshold  {print $1 "." alnFileForTreeSuffix}' \
> mafft_dna_alns_fasta_file_list.txt

# For protein file list:
cat mafft_dna_alns_fasta_file_list.txt | sed 's/dna/protein/' > mafft_protein_alns_fasta_file_list.txt


# Count the number of samples in the gene trees being used and place a sorted list into a file:
numbrSamplesInTrees=`cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sort -u | wc -l | sed 's/ //g'`
cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sed 's/>//' | sort -u > mafft_dna_alns_fasta_samples_in_tree.txt


# Print out a list of all samples submitted to this tool:
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
Number of poor quality samples (<= 20% of the largest sum length of all genes per sample): $poorQualitySamples
Number of ok/good quality samples (> 20% of the largest sum length of all genes per sample): $okQualitySamples" >> ${fileNamePrefix}_summary_stats.txt
fi


echo "
Sequence type assessed: $seqType
Sum Length of the longest gene sequence in each alignment: $sumLenLongestGene
Sum Length of the longest gene sequence in each alignment after trimming rare inserts: $sumLenLongestGeneAfterTrim
Total number of homologous columns occupied by >= $fractnMaxColOcc_pc % of samples: $sumMaxColOcc (area of common overlap)
Total number of homologous + parsimonious columns occupied by >= $fractnMaxColOcc_pc % of samples: $sumMaxParsCols
Total number of residues in area of common overlap: $totalBasesInMaxColOccRegion

Number of gene alignments containing >= $fractnSpecies_pc % of the samples after filtering by coverage (will be used in the species tree): $numbrGeneAlnsForSpeciesTree
Total number of alignment columns in the area of common overlap for these genes: $numbrOverlapColsForSpeciesTree

After filtering genes by coverage AND/OR by number of samples in each gene tree:
Number of samples that will be present in final ASTRAL and RaxML trees (total samples submitted): $numbrSamplesInTrees ($totalNumbrSamples)
Number of samples that will NOT be in the species trees: $samplesNotInSpeciesTrees
" >> ${fileNamePrefix}_summary_stats.txt


if [ $numbrSamplesInTrees -ne $totalNumbrSamples ]; then
	echo "These samples will not be present in final ASTRAL or RAxML trees: "  >> ${fileNamePrefix}_summary_stats.txt
	comm -23 samples_submitted.txt  mafft_dna_alns_fasta_samples_in_tree.txt >> ${fileNamePrefix}_summary_stats.txt
	echo >> ${fileNamePrefix}_summary_stats.txt
fi

echo 'Number times each sample will appear in the gene trees (Format: number_of_gene_trees sample).' >> ${fileNamePrefix}_summary_stats.txt
echo '(Low numbers suggest low gene recovery for that sample)' >> ${fileNamePrefix}_summary_stats.txt
cat mafft_dna_alns_fasta_file_list.txt | xargs cat | grep '>' | sed 's/>//' | sort | uniq -c | sort -k1n >> ${fileNamePrefix}_summary_stats.txt


rm samples_submitted.txt mafft_dna_alns_fasta_samples_in_tree.txt
