#!/bin/bash

##########################################
# make_species_trees_pipeline_functions.sh

# Author:   Paul Bailey

# Copyright © 2023 The Board of Trustees of the Royal Botanic Gardens, Kew
###########################################
### Consider to drop these settings if they caue issues - yes they do! They also affect main script - makes sense!!!
#set -e				# if a command finishes with a non-zero status, exit immediately
#set -u				# treat unset variables as an error and exit immediately
#set -o pipefail		# ?
shopt -s failglob	# if e.g. *.fasta doesn’t expand to anything, exit immediataly


gtm()	{
	###########
    # Function: Guilde Tree Merger (GTM)
    #			Builds the gene tree from alignment subsets derived from the main gene sequence alignment
    #			Procedure is summarized in a three step process as follows:
    #			1. Build a guide tree from the main gene sequence alignment file
    #			2.  

    # Input parameters:
    # $1 = residue type: dna, aa or codon 
    # $2 = filtered alignment filename in fasta format: <geneID>.<residue_type>.aln.for_tree.fasta
	# $3 = output file prefix including the path to any directory
	# $4 = phylogeny program to use ($phyloProgramDNA or $phyloProgramPROT)
    # $5 = raxmlng_model e.g. raxmlngModel='GTR+G' - for DNA; for protein raxmlngModel='JTT+G'
    # $6 = iqtree2Model - currently set to JTT+F+G
    # $7 = iqtree2_seq_type  iqTree2SeqType='DNA'; for protein iqTree2SeqType='AA'
    # $8 = fasttreeFlags='-nt -gtr' - for DNA; for protein fasttreeFlags='' - NB - this last flag needs to be last in case it is blank - it needs to be blank for protein analysis)

    # NB - 29.11.2023 - need to work towards simplifying the specification of the models of evolution/residue name for the makeGeneTree and gtm methods 																													

    # NB - only importing some variables into this function. All other variables are globally available and do not change during running of this script
    # NBNB - 28.11.2023 - this might not be true now I'm working in a different script to these variables - are they inherited? Check! Yes they are! e.g. 
    # 		     			$exePrefix is visible to this script - so it seems I couls also call a method from the other file as well!! - try it out 
    #						NB - function variables $1, $2, $3 etc are local to the script parameters or function - so can have nested functions with their own $1, $2 etc      
    ###########
    
    echo 'Creating subalignments for use in the guide tree merger (GTM) method'

    > ${gene}.dna.gtm_stats.txt	# Compiling stats in this file rather than the log file
    echo "#######################" >> ${3}/${gene}.dna.gtm_stats.txt
    echo "Stats on the GTM Method" >> ${3}/${gene}.dna.gtm_stats.txt
    echo "#######################" >> ${3}/${gene}.dna.gtm_stats.txt
   

    ### Testing import of the gene alignment file name
	echo $2
	### As this is an internal standard file to pipeline, can easy get the gene name again:
### Yes I don't need these prints anymore but should organise the models/residue type inputs better 
	geneId=`echo $2 | awk -F '.' '{print $1}' `  		#### NBNB - I think $gene is inherited from the make_gene_trees script so no need to redefine here - not sure now, STILL CHECK thouh
	echo geneId: $geneId
	### Could also obtain the residue type from the gene name as well, then create all the info I need concerning models/residue type internally to the makeGeneTree method
	residueType=`echo $2 | awk -F '.' '{print $2}' ` 
	echo residueType: $residueType
	

	echo "###########################################################" >> ${3}/${gene}.dna.gtm_stats.txt 
	echo "Step 1: create a full guide tree against the full alignment" >> ${3}/${gene}.dna.gtm_stats.txt
	echo "###########################################################" >> ${3}/${gene}.dna.gtm_stats.txt
	# Provides the backbone for breaking up the alignments as required and/or for the for GTM step
	# Runtime: took 11mins on the smallest tree:
	echo  exePrefix is: $exePrefix
	$exePrefix fasttree $fasttreeFlags  $2  > ${3}/${gene}.${1}.guide_gene_tree.nwk

	numbrSeqsInAlnForTree=`cat $2 | grep '>' | wc -l ` 
	numbrTipsInGuideTree=`nw_labels -I ${3}/${gene}.${1}.guide_gene_tree.nwk | wc -l`
	
	echo "Number of seqs in alignment used for the guide tree: $numbrSeqsInAlnForTree" >> ${3}/${gene}.dna.gtm_stats.txt
	echo "Number of tips in guide tree: $numbrTipsInGuideTree" >> ${3}/${gene}.dna.gtm_stats.txt
	echo "Note: the two above numbers are very likely to be smaller than the number of seqs in the original alignment from EMMA if any filtering has been applied" >> ${3}/${gene}.dna.gtm_stats.txt
	nw_labels -I ${3}/${gene}.${1}.guide_gene_tree.nwk | sort > ${3}/${gene}.${1}.guide_gene_tree.labels_sort.txt

	nw_rename  ${3}/${gene}.${1}.guide_gene_tree.nwk tree_tip_info_mapfile.txt | nw_order -c n - | nw_topology - \
	| nw_display -s -w 1500 -v 15  -b 'visibility:hidden' -I r -l 'font-size:small;font-family:Arial' \
	-i 'font-size:x-small;font-family:Arial;stroke:blue' - \
	> ${3}/${geneId}.dna.guide_gene_tree.html
	### Will need to think about access to tree_tip_info_mapfile.txt in the protein/codon folders

	### To do here:
	### Assuming all labels are unique - they might not be
	### Next step would be to call the makeGeneTree function to run fastatree --> then finally move the function to this script
	### Still keep the gtm method as a development area
	### Need to be able to also call NJ? Can I call the make_genes_trees.sh script in here as well OR are the functions already accessable in here - yes I believe they are
	###		first should set up the NJ program - see notes and lit     


    echo "######" >> ${3}/${gene}.dna.gtm_stats.txt
    echo "Step 2" >> ${3}/${gene}.dna.gtm_stats.txt
    echo "######" >> ${3}/${gene}.dna.gtm_stats.txt 
    # Trying both methods in ||el
    ########

    ### 28.12.2024 - not using step 2a - duplicate samples in subtrees
	# echo "###################################################################" >> ${gene}.dna.emma_and-or_gtm_stats.txt 
	# echo "Step 2a: Using the subset decomposition output alignments from EMMA" >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# echo "###################################################################" >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# mkdir -p ${geneId}_emma/sub-alignments-trees
	# # ### NBNB - sometimes there are *_query_[123].est.aln.fasta files which appear similr to the zero files - these will cause funny numbers - see:
	# ### alignert.py getSubproblemArguments() method - need to work out which one to use
	# for aln in ${geneId}_emma/sub-alignments/subset_*_query_0.est.aln.fasta; do

	# 	### Temporary code to remove duplicate sequences
	# 	alnWithDups=`cat $aln | grep '>' | sort | uniq -c | awk '$1 > 1' | wc -l`
	# 	echo "subalignment with duplicate seqs: $aln;  $alnWithDups duplicates"    
	# 	### Write a hash-based python method to remove the dup seqs; look if dup seqs are identical
	# 	### Might be quicker just to send Chengze an email about this!!!!!!!

	# 	echo $aln
	# 	subAlnNumbr=`echo $aln | awk -F '_' '{print $3}' `;
	# 	subAlnSubNumbr=`echo $aln | awk -F '_' '{print $5}' `;
	# 	fasttree -nt -gtr $aln > ${geneId}_emma/sub-alignments-trees/${geneId}.subset_${subAlnNumbr}_query_${subAlnSubNumbr}.dna.gene_tree.nwk
	# 	### Next step is to call the existing method with fasttree/iqtree2
	# done
	# ### Actually the above trees also contain the backbone seqs twice so they need removing.
	# ### Could run TreeShrink on individual trees!!
	# ### Need to Count the number of trees made an flag a warning if not!!!!!!!!

	# ### Finally, count number of unique sample names in all the subtrees, sort and compare with full guide tree --> list them
	# GT_vs_EMMAsubAlns=`comm -23 ${3}/${geneId}.${1}.guide_gene_tree.labels_sort.txt  ${geneId}_emma/sub-alignments/${geneId}_all_subsets_labels_sort_uniq-c.txt | wc -l `
	# echo "Number labels specific to the full guide tree (if any) compared to EMMA subalignments: $GT_vs_EMMAsubAlns"  >> ${gene}.dna.emma_and-or_gtm_stats.txt

	# if [[ $GT_vs_EMMAsubAlns -gt 0 ]]
	# 	echo "List of labels specific to the full guide tree (if any) compared to EMMA subalignments:"  >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# 	comm -23 ${3}/${geneId}.${1}.guide_gene_tree.labels_sort.txt  ${geneId}_emma/sub-alignments/${geneId}_all_subsets_labels_sort_uniq-c.txt >> ${gene}.dna.emma_and-or_gtm_stats.txt

	# 	# Remove samples, if any from the guide tree.
	# 	nw_prune ${3}/${geneId}.${1}.guide_gene_tree.nwk \
	# 	`comm -23 ${3}/${geneId}.${1}.guide_gene_tree.labels_sort.txt  ${geneId}_emma/sub-alignments/${geneId}_all_subsets_labels_sort_uniq-c.txt` \
	# 	> ${3}/${geneId}.${1}.guide_gene_tree.pruned.nwk
	# 	### Assign pruned.nwk to a variable
	# else 
	# 	echo "No need to prune guide tree w.r.t. EMMA subalns"
	# fi

	# ### Step 2a next step:
	# ### Count the number of subtrees made and flag a warning if not the correct number.


	echo "###################################################################### " >> ${3}/${gene}.dna.gtm_stats.txt
	echo "Step 2b: Alternative method to break up tree (Smirnov and Warnow, 2020)" >> ${3}/${gene}.dna.gtm_stats.txt
	echo "###################################################################### " >> ${3}/${gene}.dna.gtm_stats.txt
	################################################
	### Use the build_subsets_from_tree.py (NJMerge data repository) and mentioned in Smirnov and Warnow (2020) GTM paper to break up a full backbone tree
	mkdir -p build_subsets_from_tree
### 12.12.2023 - Sort out the PATH - same solution as EMMA
	python $NJMERGETOOLS/build_subsets_from_tree.py \
	-n 120 \
	-t ${3}/${geneId}.dna.guide_gene_tree.nwk \
	-o ${3}/build_subsets_from_tree/${geneId}_subset
	### Had to bring in bisect_tree() method into build_subsets_from_tree.py to avoid errors coming from pasta.pastaalignerjob which I don't need
	### Then script ran but the outputs appear to be without 74 samples (gene 6886)!! and 64 samples (gene 6557) but this is the same number as each subset - is that a clue?
	### It means that the gtm.py script will not work but who knows - NOW Investigate here which ones they are - 28.12.2023 - I think I solved this but not sure what the issue was.
	# Script outputs a list of tip labels.
	numbrSubsets=`ls ${3}/build_subsets_from_tree/${geneId}_subset-*-outof-*.txt | wc -l `
	echo "Number of subsets created: $numbrSubsets" >> ${3}/${gene}.dna.gtm_stats.txt
	numbrTipsInAllSubsets=`cat ${3}/build_subsets_from_tree/${geneId}_subset*.txt | sort | wc -l `
	echo "Number of tips in all subsets: $numbrTipsInAllSubsets" >> ${3}/${gene}.dna.gtm_stats.txt 
	cat ${3}/build_subsets_from_tree/${geneId}_subset*.txt | sort > ${3}/build_subsets_from_tree/${geneId}_all_subsets_labels_sort.txt


	labelsSpecficToGuideTree=`comm -23 ${3}/${geneId}.${1}.guide_gene_tree.labels_sort.txt  ${3}/build_subsets_from_tree/${geneId}_all_subsets_labels_sort.txt | wc -l `
	echo "Number labels specific to the full guide tree(if any) compared to NJMerge subsets - for gene $gene: $labelsSpecficToGuideTree" >> ${3}/${gene}.dna.gtm_stats.txt
	if [[ $labelsSpecficToGuideTree -gt 0 ]]; then 
		echo "List of labels specific to the full guide tree (if any) compared to NJMerge subsets: " >> ${3}/${gene}.dna.emma_and-or_gtm_stats.txt
		comm -23 ${3}/${geneId}.${1}.guide_gene_tree.labels_sort.txt  ${3}/build_subsets_from_tree/${geneId}_all_subsets_labels_sort.txt >> ${3}/${gene}.dna.gtm_stats.txt
	fi


	# Extract sequences from the EMMA alignment and build a tree
	for subset in ${3}/build_subsets_from_tree/${geneId}_subset-*-outof-*.txt; do

		subsetFileName=`basename -s .txt  $subset `

		seqtk subseq $2 $subset > ${3}/build_subsets_from_tree/${subsetFileName}.fasta

		#subAlnNumbr=`echo $aln | awk -F '_' '{print $3}' `
		#subAlnSubNumbr=`echo $aln | awk -F '_' '{print $5}' `
		# Used for testing subtree building:
		#$exePrefix fasttree -nt -gtr ${3}/build_subsets_from_tree/${subsetFileName}.fasta > ${3}/build_subsets_from_tree/${subsetFileName}.nwk
		if [[ $4 == 'gtm:fasttree'* ]]; then 
			makeGeneTree $1 "${3}/build_subsets_from_tree/${subsetFileName}.fasta" ${3}/build_subsets_from_tree 'fasttree' $5 $6 $7 $8
		elif [[ $4 == 'gtm:iqtree2-B1000-nm1000'* ]]; then 
			makeGeneTree $1 "${3}/build_subsets_from_tree/${subsetFileName}.fasta" ${3}/build_subsets_from_tree 'iqtree2-B1000-nm1000' $5 $6 $7 $8
		fi

		if [[ ! -s ${3}/build_subsets_from_tree/${gene}_${1}_gene_tree_USE_THIS.nwk ]]; then
			echo "WARNING: During the GTM method, a subset tree failed to be built for gene $gene: ${3}/build_subsets_from_tree/${geneId}_subset-*-outof-*.txt"
			echo "WARNING: Therefore the final merged tree should fail to be created."
		else
			# The makeGeneTree function renames trees from all phylo methods to: ${gene}_${1}_gene_tree_USE_THIS.nwk
			# So need to rename Newick file to a subtree name: 
			mv ${3}/build_subsets_from_tree/${gene}_${1}_gene_tree_USE_THIS.nwk  ${3}/build_subsets_from_tree/${subsetFileName}.nwk
			### Q: what happens if there is a file from a previous run? This shouldn't happen unless a run is interupted then phylo method fails on the second run
			if [[ -s ${3}/build_subsets_from_tree/${gene}.dna.aln_iqtree.ckp.gz ]]; then	# For iqtree analysis. To activate checkpointing for subtrees, would need to alter how files are named
				rm ${3}/build_subsets_from_tree/${gene}.dna.aln_iqtree.ckp.gz				# My own version of checkpointing will still work - ok
																							# Removing the iqtree.ckp.gz file allows the tree to be repeated and old results overwritten
			fi
		fi
	done

	numbrSubsetTrees=`ls ${3}/build_subsets_from_tree/${gene}_subset-*-outof-*.nwk | wc -l `
	echo "Number of subset trees created: $numbrSubsetTrees" >> ${3}/${gene}.dna.gtm_stats.txt
	# Measuring the subset sizes for the gene:
	avgSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_subset-*-outof-*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | awk '{sum+=$1} END {print sum/NR}'  `
	minSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_subset-*-outof-*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | sort -n | head -n 1 `
	maxSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_subset-*-outof-*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | sort -n | tail -n 1 `

	echo "Min/avg/max number of seqs in subset(s): ${minSubsetSize}/${avgSubsetSize}/${maxSubsetSize}" >> ${3}/${gene}.dna.gtm_stats.txt 

	if [[ $numbrSubsets	-ne $numbrSubsetTrees ]]; then
			echo "ERROR: During the GTM method, one or more subset trees failed to be created for gene $gene"
	fi

	

	### 28.12.2024 - not using step 3a - duplicate samples in subtrees
	# echo "###################################################" >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# echo "Step 3a: run GTM on the EMMA subtrees to merge them" >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# echo "###################################################" >> ${gene}.dna.emma_and-or_gtm_stats.txt
	# ### NB - there's always a temporary option to prune the full guide tree to remove ones not in the sub sets and continue!

	# #python /Users/pba10kg/Documents/ProgramFiles/GTM/gtm.py \
	# #--start ${3}/${gene}.${1}.guide_gene_tree.nwk \
	# #--trees ${geneId}_emma/sub-alignments-trees/${geneId}.subset_*_query_*.dna.gene_tree.nwk \
	# #--output ${geneId}.dna.GTM_gene_tree_USE_THIS.nwk

	# ### Switch in the pruned tree  - I think the numbers still don't match - yes - it's due to the alignments having been filtered
	# #python /Users/pba10kg/Documents/ProgramFiles/GTM/gtm.py \
	# #--start 6557.dna.guide_gene_tree.pruned.nwk \
	# #--trees 6557_emma/sub-alignments-trees/6557.subset_*_query_*.dna.gene_tree.nwk \
	# #--output 6557.dna.GTM_gene_tree_USE_THIS.nwk

	# #python /Users/pba10kg/Documents/ProgramFiles/GTM/gtm.py \
	# #--start 6886.dna.guide_gene_tree.nwk \
	# #--trees 6886_emma/sub-alignments-trees/6886.subset_*_query_*.dna.gene_tree.nwk \
	# #--output 6886.dna_GTM_gene_tree_USE_THIS.nwk
	# ### 1.12.2023 - This worked, despite there being duplicate samples in the guide tree - yes there are and they come out in different places
	# ### Results:
	# ### 1. The tree looks very different from the original!!!
	# ### 2. Should get dups removed now
	# ### 3. NBNB - for gene 6886 this sample is in the orig guide tree but NOT in any subset aln BUT it does appear in hte GTM tree
	# ### 	so it must be coming again from the orig guide tree?????:
	# #>9299
	# #ATGGCCTGGAGATCGCAGTTGTCTAAGAACCTGAAAGAGGTTCGGATTCTCTTCTGCCAAACATCCCCTGCAAGTGCTGAAGC---T-A--GG------------------------------A-GA---ACCTTTGTCGAGAAGAATTACAAGGAGCTCAAGACATTGAACCCGAAGTTGCCAATCTTGATTCGCGAATGTCGTGGGGTGGAACCTCAGCTTTGGGCTAGATAT



	######################################################################################
	# Step 3b: run GTM on the other subtree method to merge them (build_subsets_from_tree) 
	######################################################################################
	python $GTM \
	--start ${3}/${gene}.${1}.guide_gene_tree.nwk \
	--trees ${3}/build_subsets_from_tree/${gene}_subset-*-outof-*.nwk \
	--output ${3}/${gene}_dna_gene_tree_USE_THIS.nwk

	numbrTipsInGTM_Tree=`nw_labels -I ${3}/${gene}_dna_gene_tree_USE_THIS.nwk | wc -l`
	echo "Number of tips in final GTM tree: $numbrTipsInGTM_Tree" >> ${3}/${gene}.dna.gtm_stats.txt

	# NB - this tree image is not required to be made here - it gets done in the createGeneAlignmentAndTreeImages subfunction:
	nw_rename  ${3}/${gene}_dna_gene_tree_USE_THIS.nwk  tree_tip_info_mapfile.txt | nw_order -c n - | nw_topology - \
	| nw_display -s -w 1500 -v 15  -b 'visibility:hidden' -I r -l 'font-size:small;font-family:Arial' \
	-i 'font-size:x-small;font-family:Arial;stroke:blue' - \
	> ${3}/${gene}_dna_gene_tree_USE_THIS.taxonomy.html
	

	if [[ -s ${3}/${gene}_dna_gene_tree_USE_THIS.nwk ]]; then

		# Calculate Robinson-Foulds between the guide and GTM trees:
		RF_guide_vs_GTM=`python $NJMERGETOOLS/compare_trees.py  ${3}/${gene}.${1}.guide_gene_tree.nwk  ${3}/${gene}_dna_gene_tree_USE_THIS.nwk | awk '{print $6}' `
		echo "Robinson-Foulds distance (guide vs GTM tree): $RF_guide_vs_GTM" >> ${3}/${gene}.dna.gtm_stats.txt

		if [[ $numbrTipsInGuideTree -ne $numbrTipsInGTM_Tree ]]; then
			echo "WARNING: GTM method was not able to create a gene tree or expected number of tree tips is wrong - skipping tree building of $2"
			exit 0	# zero allows Slurm to continue with the dependancies; NB - also exits from make_gene_trees.sh
		fi
	fi
}





