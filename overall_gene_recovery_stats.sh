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
depth=4
cat $geneList | \
while read gene; do
    medianPoint=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | awk 'END {printf "%.0f " , NR/2}' `
    medianValue=`cat Sample_*/*_bwa_mem_sort_st_depth.txt | grep ^$gene | awk -v depth=$depth '$3 >= depth' | sort -k3n | awk '{print $3}' | head -n $medianPoint | tail -n 1 `
    echo $gene $medianValue
done > medianReadDepthPerGene_min${depth}x_st_depth.txt
"