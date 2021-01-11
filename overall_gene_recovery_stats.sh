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


# Read depth per gene across all samples
# Collect median of median values:
### Create a header line here for medianReadDepthPerGene_st_depth.txt
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
    medianPoint3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue3=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint3 | tail -n 1 `
    # Read depth without removing duplicate reads:
    #depth=4	# 
    #medianPoint4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    #medianValue4=`cat Sample_*/*_bwa_mem_with_dups_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint3 | tail -n 1 `

    echo $gene $medianValue1 $medianValue2 $medianValue3 # $medianValue4
done > medianReadDepthPerGene_st_depth.txt
# NB - for the mediumPoint value accross all samples, there will be very likely all genes present (in the case of Angiosperms353 set, MediumPoint=353/2), but per sample unlikely to have all genes - OK
