#!/bin/bash

##################################
# recover_genes_from_one_sample.sh

# Author: Paul Bailey
##################################
#set -e
#set -u
#set -o pipefail
shopt -s failglob


geneList=/science/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt
### NBNB - maybe use the fasta file so that i get the length of each gene


# Read depth per gene across all samples(from Samtools coverage)
# Collect median of median values:
### Create a header line here for medianReadDepthPerGene_st_depth.txt

#### NB - could also collect # samples per gene from the recovered fasta file -- ooo thing about this one

cat $geneList | \
while read gene; do
	# Read depth after removing duplicate reads:
	depth=0
    medianPoint1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue1=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint1 | tail -n 1 `
    depth=1
    medianPoint2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue2=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint2 | tail -n 1 `
    depth=4
    ### May need to add for very poor quality samples where some genes might have hardly any coverage:
    ###if [[ $medianPoint4 -eq 0 ]]; then
	###	echo "medianReadDepth_min4x_(samtools_depth): 0" >> ${sampleId}_gene_recovery_stats.txt
    medianPoint3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint3 | tail -n 1 `
    # Read depth without removing duplicate reads:
    depth=4
    medianPoint4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint4 | tail -n 1 `

    echo $gene $medianValue1 $medianValue2 $medianValue3 $medianValue4
done > medianReadDepthPerGene_st_depth.txt
# NB - for the mediumPoint value accross all samples, there will be very likely all genes present (in the case of Angiosperms353 set, MediumPoint=353/2), but per sample unlikely to have all genes - OK


### NBNB - now i can just do this - gets each sample stats on a single line for producing scatterplots and a table for importing into Excel
### for file in Sample_*/*_gene_recovery_stats.txt; do line=`cat $file | awk '{printf $2 " "}' `; echo "$file $line"; done > gene_recovery_stats_per_sample.txt
### NBNB - don't need $file variable in line now as I have the sampkleID in file itself now 


### Print out histograms per gene:
#cat $geneList | \
###cat /science/projects/paftol/PAFTOL_additional_files/Angiosperms353_targetSequences_organism-gene_format_corrected_geneIds_ONLY.txt | \
###while read gene; do 
###	cat Sample_*/*_bwa_mem_sort_st_covrg_-m.txt | grep -A 12 "^$gene "
###done > bwa_mem_sort_st_covrg_-m_per_gene_all_genes.txt  						


# Per sample read depth across all loci
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min0x_(samtools_depth)' | sort -k2n | head -n 1083 | tail -n 1
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min0x_(samtools_depth):' | sort -k2n > gene_recovery_stats_medianReadDepth_min0x_samtools_depth.txt
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min1x_(samtools_depth)' | sort -k2n | head -n 1083 | tail -n 1
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min1x_(samtools_depth):' | sort -k2n > gene_recovery_stats_medianReadDepth_min1x_samtools_depth.txt
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min4x_(samtools_depth)' | sort -k2n | head -n 1083 | tail -n 1
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepth_min4x_(samtools_depth):' | sort -k2n > gene_recovery_stats_medianReadDepth_min4x_samtools_depth.txt
# cat 	Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepthWithDups_min4x_(samtools_depth)' | sort -k2n | head -n 1083 | tail -n 1
# cat   Sample_*/*_gene_recovery_stats.txt | grep 'medianReadDepthWithDups_min4x_(samtools_depth):' | sort -k2n > gene_recovery_stats_medianReadDepthWithDups_min4x_samtools_depth.txt


### For on-target assessemnt:

###cat   Sample_*/*_gene_recovery_stats.txt | grep 'numbrTrimmedReads:\|numbrTrimmedReadsInBam:' | sort -k2n > gene_recovery_stats_numbrTrimmedReadsInBam.txt

#cat   Sample_*/*_gene_recovery_stats.txt | grep 'numbrReadsOnTarget' | sort -k2n | head -n 1083 | tail -n 1
#cat   Sample_*/*_gene_recovery_stats.txt | grep 'numbrReadsOnTarget' | sort -k2n > gene_recovery_stats_numbrReadsOnTarget:.txt

#cat   Sample_*/*_gene_recovery_stats.txt | grep 'numbrMappedReads' | sort -k2n | head -n 1083 | tail -n 1
