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
	echo "Note: the two above numbers are very likely to be smaller than the number of seqs in the original alignment from EMMA if any filtering has been applied"
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
	echo "Number of subsets: $numbrSubsets" >> ${3}/${gene}.dna.gtm_stats.txt
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
		makeGeneTree $1 $2 $3 'fasttree' $5 $6 $7 $8
		###makeGeneTree $1 $2 $3 'iqtree2-B1000-nm1000' $5 $6 $7 $8

		if [[ ! -s ${3}/build_subsets_from_tree/${subsetFileName}.nwk ]]; then
			echo "WARNING: During the GTM method, a subset tree failed to be built for gene $gene: ${3}/build_subsets_from_tree/${geneId}_subset-*-outof-*.txt"
			echo "WARNING: Therefore the final merged tree should fail to be created."
		fi
	done

	# Measuring the subset sizes for the gene:
	avgSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | awk '{sum+=$1} END {print sum/NR}'  `
	minSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | sort | head -n 1 `
	maxSubsetSize=`for subsetTree in ${3}/build_subsets_from_tree/${gene}_*.nwk; do
				   nw_labels -I $subsetTree | wc -l
				   done | sort | tail -n 1 `

	echo "Min/avg/max number of seqs in subset size: ${minSubsetSize}/${avgSubsetSize}/${maxSubsetSize}" >> ${3}/${gene}.dna.gtm_stats.txt 

	

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
	

	#######################
	### UPTOHERE 30.12.2023
	#######################
	### ASAP get either iqtree, raxml or raxml ng working - might be good to test the recent paper on the ML plateaux (raxml)
	###		TEST function call
	###			Need to sort out options -q and -r if value is gtm:
	###				gtm:fasttree:100  --> need to break this up so that can use extra  

	###		Commit to Git:
	###		1.In gene_tree script . $pathToScripts/make_species_trees_pipeline_functions.sh # Sourcing a file of bash functions
	###		2.Set up basic working of this script
	###		2.Add call to the GTM module

	###		FIRST set up on KewHPC FastTree (as a control) --> iqtree --> raxml(-ng)
	###			Once I know how long each subset takes for the longest gene and can estimate a time/RAM for a subset size, then this work will be going well 
	###		Think about how to specify size of subset
	###		THEN consolidate and plan next dev steps incl sample list + aa alns		

	### 	Also look at other phylo papers for ideas I've forgotten about and my dev notes
	###		THEN set up on Gruffalo for future repeats e.g. with gene tree backbones


	### Further work on the GTM function:
	### Go through method and tidy up again
	### Step 2a and 3a - using EMMA subsets - need to remove the dups from the EMMA subsets --> see if trees look any better - DELAY this now
	### Make sure I record how I've modified the scripts and make a note - e.g. line 154 - 158 - add notes to install.sh script
	###		Add to bash_profile like so:
	###		export NJMERGETOOLS=/Users/pba10kg/Documents/ProgramFiles/NJMerge_repository_Erin_and_Tandy_2018/tools
	###		export GTM=/Users/pba10kg/Documents/ProgramFiles/GTM/gtm.py
	###		ALSO: add in code to check that the two scripts are installed at the start of the method, else report an error!!!! 
	
	
	### Compile message to Chengze and plan to update the Google doc - looks at the other papers in there
	### 1. Confirm what the extra subset_*_query_[12345].est.aln.fasta - are those the subaln problem?
	### 2. Have another go at removing the 'queries' value - just try removing the whole clause --> ask
	### Duplicate samples

	### Now I'm using outgroup labels, the tree images are NOT being made - need to check first what
	### what happens when no labels are present + relogic

	###-->Note to try and makeGeneTree() method more flexible w.r.t. model of evolution/residue type
	### NB - note that the GTM will not work for paralogs instended for ASTRAL-PRO the way I intended to use them
	### 	Mihgt afterall need to indicate paraologs first but the pipeline would need changing - think

	### NB - I know why the EMMA subalns don't add up to the guide treee: I'm uisng the subalns before any filtering,
	###		 so you would have to run the filtering on all the EMMA sub alns - possible - but easiest to use my second approach.
	###		 Actually I think I could just compare a subaln ID list with the guide tree ID list and remove seqs in the sub aln that no longer exist in the guide tree.
	### 	 Hey, I think this solves the issue of duplicate seqs in the EMMA subalns

	### NB - note to check to test RAxML+CAT model with no bootstrapping for the guide tree - using 10 starting trees - it might not be too slow
	### Note to think about doing the concatenated aln in the same way, use FastTree or Astral as starting tree, then divide a concatenated aln up into  
	### subalns and trees using a proper ML method - might want to also consider NJst as the starting tree if Fasttree becomes too slow https://doi.org/10.1093/sysbio/syr027

	###Next action 7:—>Implement filtering and trimming for concatenated aln by using existing methods moved to this script - required by two scripts 
	
	###	I could use that R primer from EI to shade the subsets a different colour on the guide tree so I can track where they are in the full tree!!!!!
	###		https://www.earlham.ac.uk/articles/plotting-phylogenetic-trees-r-alternating-clade-highlights?utm_source=Earlham+Institute+Digital+Updates&utm_campaign=eb7d13f07d-Mailchimp_Nov23&utm_medium=email&utm_term=0_3c2e0987cf-eb7d13f07d-457393009
	### 	Plotting phylogenetic trees in R: alternating clade highlights - A guide to highlighting clades in phylogenetic trees in R 23 November 2023 Rowena Hill

	### Note to calculate Robinson-Foulds between each gene tree and the species tree - has to be done after the species tree is made:
	###		Use compare_trees.py
	###		Try and output to the existing tree stats file

	### Remmber that I can force nw_reroot to use the ingroup to get the tree rooted on the outgroup - see plastid work - alter createGeneAlignmentAndTreeImages function to do this.
	###		NB - need to catch the error or whether file is empty depending on what happens
	###		Sometimes no labels given might be in the tree - what happens if some are missing?

	### Look into FastTree MPI - Chengze has it working - compare between 1cpu fasttree and ||el version
}





