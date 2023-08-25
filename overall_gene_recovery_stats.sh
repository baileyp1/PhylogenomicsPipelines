#!/bin/bash

##################################
# recover_genes_from_one_sample.sh

# Author: Paul Bailey

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
##################################
#set -e
#set -u
#set -o pipefail
shopt -s failglob


### Need to bring in the Sample prefix; might vary from 'Sample*'
### Yes, at the moment script is assuming sample prefix is 'Sample' (option -p default)


# Required for the per gene stats (readDepthPerGene_st_depth.txt):
geneList=/data/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt
### NBNB - maybe use the fasta file so that i get the length of each gene



# Print ALL stats per sample on a single line for producing scatterplots and a table for importing into Excel
# (printf prints out all values from each sample on a single line.)
for file in Sample_*/*_gene_recovery_stats.txt; do line=`cat $file | awk '{printf $2 " "}' `; echo "$file $line"; done > gene_recovery_stats_per_sample.txt
### NBNB - don't need $file variable in line now as I have the sampleID in file itself now 



# Read depth per gene across all samples(from Samtools depth)
# Collecting average of average and median of median values:
# Output file header:
echo gene medianReadDepth_min0x medianReadDepth_min1x meanReadDepth_min4x medianReadDepth_min4x meanReadDepthWithDups_min4x \
medianReadDepthWithDups_min4x numbrBasesInAllSamples_ReadDepth_min0x numbrBasesInAllSamples_ReadDepth_min4x numbrBasesInAllSamples_ReadDepthWithDups_min4x > readDepthPerGene_st_depth.txt

### NB - could also collect # samples per gene from the recovered fasta file - need to think about this one
### Also shoudl report the number of samples being used to geenrate the stats accross all samples

cat $geneList | \
while read gene; do
	# Read depth after removing duplicate reads:
	depth=0
	# Number of recovered bases in ALL samples with read depth >= 0 (i.e. total of all bases):
	numbrBasesInAllSamples_ReadDepth_min0x=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | wc -l `

	# Mean read depth for bases with >= 0x depth  across ALL genes:
	###meanReadDepth1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
    
    medianPoint1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint1 | tail -n 1 `
	
    depth=1
    # Mean read depth for bases with >= 4x depth  across ALL samples:
	###meanReadDepth2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
    
    medianPoint2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint2 | tail -n 1 `
    
    depth=4
     # Mean read depth concatfor bases with >= 4x depth  across ALL samples:
	meanReadDepth3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average

    ### May need to add for very poor quality samples where some genes might have hardly any coverage:
    ###if [[ $medianPoint4 -eq 0 ]]; then
	###	echo "medianReadDepth_min4x_(samtools_depth): 0" >> ${sampleId}_gene_recovery_stats.txt
    medianPoint3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint3 | tail -n 1 `
    # Read depth without removing duplicate reads:
    depth=4
    # Mean read depth for bases with >= 4x depth across ALL samples:
	meanReadDepth4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene |  awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average

    medianPoint4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint4 | tail -n 1 `

    # Number of recovered bases in ALL samples with read depth >= 4:
    # Minus read duplicates:
	numbrBasesInAllSamples_ReadDepth_min4x=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | wc -l `
	# Plus read duplicates:
	numbrBasesInAllSamples_ReadDepthWithDups_min4x=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | wc -l `
	

	# Count number of all ambiguity codes across all samples per gene:
	### NB - this doesn't work - it is NOT per gene - need to use a gene-wise file so better done as part of the phylo step
	###numbrAmbiguityCodesInGene=`cat Sample_*/*.fasta | grep -v '>' | grep -o '[RYMKSWHBDN]' | wc -l `


    # These commands take a long time so restricting mean value calculations to >=4x depth only for now:
    echo $gene $medianValue1 $medianValue2 $meanReadDepth3 $medianValue3 $meanReadDepth4 \
$medianValue4 $numbrBasesInAllSamples_ReadDepth_min0x $numbrBasesInAllSamples_ReadDepth_min4x $numbrBasesInAllSamples_ReadDepthWithDups_min4x    
done >> readDepthPerGene_st_depth.txt
# NB - for the mediumPoint value accross all samples, there will be very likely all genes present (in the case of Angiosperms353 set, MediumPoint=353/2), but per sample unlikely to have all genes - OK



### Print out histograms per gene:
#cat $geneList | \
###cat /science/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt | \
###while read gene; do 
###	cat Sample_*/*_bwa_mem_sort_st_covrg_-m.txt | grep -A 12 "^$gene "
###done > bwa_mem_sort_st_covrg_-m_per_gene_all_genes.txt  						
