#!/bin/bash

##################################
# recover_genes_from_one_sample.sh

# Author: Paul Bailey
##################################
#set -e
#set -u
#set -o pipefail
shopt -s failglob


### Need to bring in the Sample prefix


# Required for the per gene stats:
geneList=/science/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt
### NBNB - maybe use the fasta file so that i get the length of each gene



# Print ALL stats per sample on a single line for producing scatterplots and a table for importing into Excel
# (printf prints out all values from each sample on a single line.)
for file in Sample_*/*_gene_recovery_stats.txt; do line=`cat $file | awk '{printf $2 " "}' `; echo "$file $line"; done > gene_recovery_stats_per_sample.txt
### NBNB - don't need $file variable in line now as I have the sampleID in file itself now 



# Read depth per gene across all samples(from Samtools coverage)
# Collect average of average and median of median values:
### Create a header line here for medianReadDepthPerGene_st_depth.txt

### NB - could also collect # samples per gene from the recovered fasta file - need to think about this one

cat $geneList | \
while read gene; do
	# Read depth after removing duplicate reads:
	depth=0
	# Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepth1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
    
    medianPoint1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint1 | tail -n 1 `
	
    depth=1
    # Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepth2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average
    
    medianPoint2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint2 | tail -n 1 `
    
    depth=4
     # Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepth3=`cat Sample_*/*__bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average

    ### May need to add for very poor quality samples where some genes might have hardly any coverage:
    ###if [[ $medianPoint4 -eq 0 ]]; then
	###	echo "medianReadDepth_min4x_(samtools_depth): 0" >> ${sampleId}_gene_recovery_stats.txt
    medianPoint3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint3 | tail -n 1 `
    # Read depth without removing duplicate reads:
    depth=4
    # Mean read depth for bases with >= 4x depth  across ALL genes:
	meanReadDepth4=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | awk '{sum+=$3} END {if(sum > 0) {print sum/NR} else {print "0"}}' `	# average

    medianPoint4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint4 | tail -n 1 `

    # Number of recovered bases in ALL genes with read depth >= 4:
    # Minus read duplicates:
	numbrBasesInAllGenes_ReadDepth_min4x=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | wc -l
	# Plus read duplicates:
	numbrBasesInAllGenes_ReadDepthWithDups_min4x=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | awk -v depth=$depth '$3 >= depth' | wc -l
	

	# Count number of all ambiguity codes across all samples per gene:
	numbrAmbiguityCodesInGene=`cat Sample_*/*.fasta | grep -v '>' | grep -o '[RYMKSWHBDN]' | wc -l `


    echo $gene $numbrAmbiguityCodesInGene meanReadDepth1 $medianValue1 meanReadDepth2 $medianValue2 meanReadDepth3 $medianValue3 meanReadDepth4 $medianValue4 $numbrBasesInAllGenes_ReadDepth_min4x $numbrBasesInAllGenes_ReadDepthWithDups_min4x
done > readDepthPerGene_st_depth.txt
# NB - for the mediumPoint value accross all samples, there will be very likely all genes present (in the case of Angiosperms353 set, MediumPoint=353/2), but per sample unlikely to have all genes - OK



### Print out histograms per gene:
#cat $geneList | \
###cat /science/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt | \
###while read gene; do 
###	cat Sample_*/*_bwa_mem_sort_st_covrg_-m.txt | grep -A 12 "^$gene "
###done > bwa_mem_sort_st_covrg_-m_per_gene_all_genes.txt  						
