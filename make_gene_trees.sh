#!/bin/bash

####################
# make_gene_trees.sh

# Author: Paul Bailey
####################
geneId=$1
geneFile=$2
fractnAlnCovrg=$3
phyloProgramDNA=$4
phyloProgramPROT=$5
fractnMaxColOcc=$6
cpuGeneTree=$7
mafftAlgorithm="$8"
exePrefix="$9"
alnProgram="${10}"
dnaSelected="${11}"
proteinSelected="${12}"
codonSelected="${13}"
filterSeqs1="${14}"
pathToScripts="${15}"
maxColOccThreshold="${16}"
filterSeqs2="${17}"
trimAln1="${18}"
trimAln2="${19}"
treeshrink="${20}"

# Convert $emptyMatchStateFractn to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `


echo Inside make_gene_trees.sh script:
echo SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID
echo SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID


geneId=`echo $geneId | cut -d ',' -f 1 `		#### 20.8.2019 - does this do anything now? Required only for ${geneId}_aln_summary.log but i coudl change that.
gene=`echo $geneId | tr -d '\n' `				# Need to remove line return, if present after the first field! Was for the sample- to gene-wise conversion 
### Might also be good to have the family
#genus=`echo $line | cut -d ',' -f 2 `
#species=`echo $line | cut -d ',' -f 3 `
#R1FastqFile=`echo $line | cut -d ',' -f 4 `
#R2FastqFile=`echo $line | cut -d ',' -f 5 `


echo GeneId: $geneId
echo fractnAlnCovrg: $fractnAlnCovrg
echo fractnAlnCovrg_pc: $fractnAlnCovrg_pc
echo geneFile: $geneFile
echo mafftAlgorithm: $mafftAlgorithm
echo alnProgram: $alnProgram
echo dnaSelected: $dnaSelected
echo proteinSelected: $proteinSelected
echo codonSelected: $codonSelected
echo phyloProgramDNA: $phyloProgramDNA
echo phyloProgramPROT: $phyloProgramPROT
echo maxColOccThreshold: $maxColOccThreshold
echo filterSeqs1: $filterSeqs1
echo filterSeqs2: $filterSeqs2
echo trimAln1: $trimAln1
echo trimAln2: $trimAln2
echo treeshrink: $treeshrink


filterSeqs1()	{
	###########
    # Function: filters sequences from an alignment of any residue type
    #           depending on their coverage across well occupied columns,
    #			as determined in this function. Also filter sequences from DNA
    #			and codon alignments using the results of filtering the protein 
    #			alignment, if a protein alignment has been selected.
    #			Also creates filtering stats and prepares files for stats.

    # Input parameters:
    # $1 = alignment filename in fasta format
    # $2 = residue type for AMAS.py option -d: dna, aa  (NB - AMAS trim -d dna also seems to work with aa!) 
    # $3 = maxColOcc threshold - 90 for DNA, 30 for protein
    # $4 = seqLength threshold - 84 for DNA, 28 for protein - deprecated - will rely on a minimum overall seq length to tolerate after all filtering and trimming
    #
    # NB - only importing variables into this function if they vary depending on the residue type of $1,
    #      all other variables a globally available and do not change during running of this script
    ###########
    if [[ $2 == 'protein' ]]; then 
    	residueType=aa 
    else
    	residueType=$2
    fi

    # Empty this file so contents is removed from the last run and new content can be appended:
    > ${gene}_filtered_out_alignments.txt


    # Create gene alignment summary file or empty an existing one:
	>${geneId}_aln_summary.log
	echo "gene Id: ${geneId}
residueType: $2" >> ${geneId}_aln_summary.log

	# Find the length of the alignment to calculate the coverage of each sequence:
	# NB - this code works because seqtk also counts dashes - c.f. fastalength which does not.
	# NBNB - putting seq on a single line then can count length of line i.e seq (incl. dashes)
	# NBNBNB - /dev/fd/0 - have to explicitly redirect into seqtk on Macbook Darwin.
	# NBNBNBNB - Not using this full aln value now - using length of the longest gene in the aln (see just below).
	alnLength=`cat $1 | seqtk seq -l 0 /dev/fd/0 | tail -n1 | awk '{print length($1)}' `
	#echo Alignment length: $alnLength


# Using fastalength to determine the longest gene in the alignment, rather than the full length of each alignment (as above)
# This is probably a better messure for the total aln length esp for a phylogeny subset e.g. single family; total aln length is probably too strict.
#lenLongestGene=`fastalength  ${gene}_mafft_dna_aln.fasta | sort -n | tail -n 1 | awk '{print $1}' `
#echo lenLongestGene: $lenLongestGene  >> ${geneId}_aln_summary.log
### NB - 28.1.2020 - now calculating this AFTER filtering for $maxColOcc and for > 3 seqs so it is directly comparable to lenLongestGeneAfterTrim
### 10.10.2019 - NBNB - the issue with this approach to calculating the longest gene is that some genes are much bigger than the majority of the alignment.
### BUT - I think the issue occurs with only a small number of seqs so I could calculate the median value for the top say 10 of seqs and use that value.
### Could print out the following stats here: longest gene, avrg and median of top 10%.
### NB - 13.1.2019 - this stat also counts longest gene for genes < ~100 bp maxColOcc - won;t make musch differecne but could move the calculation down to within the conditional!!


# Can now easily identify the number of columns occupied (ideally) by ALL seqs or at least a high proportion e.g. 70%.
# These columns will be the most informative ones for phylogeny, i.e. ones that are mostly occupied (being generous with 70 %):
###AMAS.py trim -f fasta -d $2 -t ${fractnMaxColOcc} -i $1 -o ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta
AMAS.py trim -f fasta -d $residueType -t ${fractnMaxColOcc} -i $1 -o ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta
###maxColOcc=`fastalength  ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
maxColOcc=`fastalength  ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
# NB - some seqs will have zero length seqs but the fastalength warning can go to /dev/null; using this command to check the occupancy but nothing more 
echo maxColOcc: $maxColOcc  >> ${geneId}_aln_summary.log

### Old method and ideas:
# $maxColOcc needs to be at least 150 bp so that even the shortest seqs are ~100bp over the region.
# 150bp is length of the shortest target seqs (except for a few - looking at Angiosperm353 and plastid target seqs)
# Was also going to fail a gene if the $lenLongestGene : $maxColOcc ratio > 3x - would suggest that there is a likely issue with obtaining overlapping seqs:
# ratioLenOcc=`echo $length $maxOcc | awk '{printf "%.0f", $1/$2}' `
# However, then I'm relying on $lenLongestGene again and the are some issues with that (see above)   
# Just testing $maxColOcc should be OK.

### New change - 25.9.2020 - allow $maxColOcc to work down to zero so all genes can be included if desired.
# This completes the logic of this filtering approach: filtering genes at $fractnAlnCovrg across region of $maxColOcc, 
# the latter being filtered by $maxColOccThreshold set by user.
# NB - if there are no overlapping columns found for a gene, then this filtering is not possibe AND 
# the awk command will crash at $1/LENG if $maxColOcc and $3 are both zero
# So now allowing all seqs for these genes through where $maxColOcc==0 (they tend to be small genes),
# just need to set $maxColOcc to 1 for awk to work.
if [[ $maxColOcc -eq 0 ]]; then
	maxColOcc = 1
fi
if [[ $maxColOcc -ge $3 ]]; then

	# Calculate the number of parsimonious columns (-p option) for the columns found from $maxColOcc:
	AMAS.py trim -f fasta -d dna -t $fractnMaxColOcc -p -i ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta -o ${gene}_${2}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta
	maxParsCols=`fastalength ${gene}_${2}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
	echo maxParsCols: $maxParsCols  >> ${geneId}_aln_summary.log

	# Count total number of residues in the MCOT region here:
	totalResiduesInMaxColOccRegion=`fastalength ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta 2>/dev/null | awk '{sum+=$1} END {print sum}' `
	echo totalResiduesInMaxColOccRegion: $totalResiduesInMaxColOccRegion >> ${geneId}_aln_summary.log

	# Now find the lengths of each sequence in bases only (not dashes!)
	# and create a list of sequences that are >= $fractnAlnCovrg across $maxColOcc region
	# NB - this code works because, luckily for this case, fastalength ignores dash characters (unlike seqtk used above)!	
	fastalength ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta 2>/dev/null \
	| awk -v LENG=$maxColOcc -v FRACTN=$fractnAlnCovrg '{if( ($1/LENG) >= FRACTN) {print $2} }' \
	> ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
	# NB - 25.9.2020 - removed seqLenThreshold from awk line: | awk -v LENG=$maxColOcc -v FRACTN=$fractnAlnCovrg -v seqLenThreshold=$4 '{if( ($1/LENG) >= FRACTN && $1 >= seqLenThreshold) {print $2} }' \
	# Note if .fasta file is empty, fastalength still exits normally.

### 19.8.2020 - consider to do this: else if < $maxColOcc just remove seqs below seqLenThreshold abd continue with same file as above but no filtering with $fractnAlnCovrg

	# Get seqs with ideal coverage, except, if this seq list file made just above 
	# only has 3 or less seqs then there is also no point in making a gene tree.
	numbrSeqs=`cat ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt | wc -l `
	echo numbrSeqs for tree = $numbrSeqs
	if [ "$numbrSeqs" -gt 3 ]; then 		# Conditonal here not essential - but don't know how AMAS.py behaves with empty input files so leaving in here
											# Actually now using to print out filename if it fails
											### 27.8.2020 - can improve this check; if possible leave it to the main code to not build tree if < 3 seqs after all filtering checks 
											### For now have just repeated the identical warning message at this stage as well
											### 8.9.2020 - thinking further, this conditional should be removed but still check if numbrSeqs < 3 so I can print out alns
											###				Also still print out alns with maxColOcc < 30 separately
### 28.9.2020 - now I'm thinking that I shoud lstill keep it:
### 1. saves producing files that don't get usec for tree 
### 2. don't then have to check the files when trimming as they will be empty... - check this!
		seqtk subseq -l $alnLength \
		$1 \
		${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
		> ${gene}.${2}.aln.after_filter1.fasta
		# Old output file name: ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
		# NB - if .txt and .fasta are empty seqtk subseq still exits normally.

		# Record stats (NB - before trimming)
		lenLongestGene=`fastalength  $1 | sort -n | tail -n 1 | awk '{print $1}' `
		echo lenLongestGene: $lenLongestGene  >> ${geneId}_aln_summary.log

		# Median value of sequence length for current gene:
		medianPoint=`fastalength $1 | awk 'END {printf "%.0f" , NR/2}' `
		medianGeneLength=`fastalength $1 | awk '{print $1}' | sort -n | head -n $medianPoint | tail -n 1 `
		echo medianGeneLength: $medianGeneLength  >> ${geneId}_aln_summary.log


		# Extra steps here to do if DNA and/or codon is also selected.
		# Prepare protein and codon alns using result of the filtering on the protein aln using this list: 
		# ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
		echo "proteinSelected == 'yes' && dnaSelected == 'yes'": $proteinSelected $dnaSelected
		if [[ $proteinSelected == 'yes' && $dnaSelected == 'yes' ]]; then
			echo Also preparing the filtered dna aln...
			seqtk subseq -l $alnLength \
			${gene}.dna.aln.fasta \
			${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
			> ${gene}.dna.aln.after_filter1.fasta
			# Old output file name: ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
### DON'T NEED TO DO THIS STEP HERE NOW - DONE IN TRIMALN2 FUNCTION - REMOVE
			# AMAS.py trim -t 0.003 \
			# -f fasta \
			# -d dna \
			# -i ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
			# -o ${gene}.dna.aln.for_tree.fasta
			### 12.8.2020 - changed the filename:
			###-o ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta

			### if aligning the high occupancy columns do above for DNA and protein here - add them to a new directory
			### if protein and DNA only, need to add above with a condtional for high occ columns
		fi
		if [[ $codonSelected == 'yes' ]]; then
			echo Also preparing the filtered codon aln...
			seqtk subseq -l $alnLength \
			codonAln/${gene}.codon.aln.fasta \
			${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
			> codonAln/${gene}.codon.aln.after_filter1.fasta
			# Old output file name: codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
### DON'T NEED TO DO THIS STEP HERE NOW - DONE IN TRIMALN2 FUNCTION - REMOVE
			# AMAS.py trim -t 0.003 \
			# -f fasta \
			# -d dna \
			# -i codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
			# -o codonAln/${gene}.codon.aln.for_tree.fasta
			### 12.8.2020 - changed the filename:
			###-o codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta

			### if aligning the high occupancy columns do above for DNA here - add them to a new directory
		fi
	else
		# Printing out the alignments with < 4 sequences for assessing further (can concatenate them together)
		# NB - this file is produced only for the protein seqs, even if DNA is selected; but is produced for DNA if DNA only has been selected.
		echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
		echo $1 >> ${gene}_filtered_out_alignments.txt
	fi
else
	# Printing out the alignments with little sequence overlap for assessing further (can concatenate them together)
	echo $1 >> ${gene}_filtered_out_alignments.txt 
fi 
} # End of filterSeqs1 function


filterSeqs2()	{
	###########
   	# Function: filters sequences based on assessing sequence coverage length 
   	#           based on reference target length - specifically median length
   	#           of all sequences in set. Also filter sequences from DNA
    #			and codon alignments using the results of filtering the protein 
    #			alignment, if a protein alignment has been selected.
   	#
   	# Input parameters:
   	# $1 = input alignment fasta file
   	# $2 = residue type
   	# $3 = percent minimum sequence length to tolerate
   	# $4 = output filtered alignment fasta file
   	###########

   	# Create gene alignment summary file or empty an existing one:
	>${geneId}_aln_summary.log
	echo "gene Id: ${geneId}
	residueType: $2" >> ${geneId}_aln_summary.log

   	# Find the median value of the sequence lengths for current gene set (NB - fastalength ignores dash chars):
	medianPoint=`fastalength $1 | awk 'END {printf "%.0f" , NR/2}' `
	medianGeneLength=`fastalength $1 | awk '{print $1}' | sort -n | head -n $medianPoint | tail -n 1 `
	echo medianGeneLength: $medianGeneLength  >> ${geneId}_aln_summary.log

	# Use the median length to find sequences that are >= $fractnAlnCovrg of the median length
	# and create a list of them.
 	# NB - this code works because, luckily for this case, fastalength ignores dash characters!
 	# NB - can pipe identifier list into seqtk subseq! Saves making an output file, EXCEPT do need a list file to mirror output for DNA and codon below!
	fastalength $1 2>/dev/null \
	| awk -v LENG=$medianGeneLength -v minPercent=$filterSeqs2 '{if( (($1/LENG) * 100) >= minPercent) {print $2} }' \
	| seqtk subseq -l 0 $1 /dev/fd/0 > "$4"


	# Still need to get list file to mirror output for DNA and codon below, if selected:
	fastalength $1 2>/dev/null \
	| awk -v LENG=$medianGeneLength -v minPercent=$filterSeqs2 '{if( (($1/LENG) * 100) >= minPercent) {print $2} }' \
	> ${gene}.${2}.aln.after_filter2.txt
	

	### UPTOHERE 29.9.2020
	### Now think about the stats:
	### 1. might want to put the numbrSeqs in the logs - before and after
	###numbrSeqs=`cat  <fasta file>  | wc -l `
	###echo numbrSeqs for tree = $numbrSeqs


	# Extra steps here to do if DNA and/or codon is also selected.
	# Prepare protein and codon alns using result of the filtering on the protein aln using this list: 
	# ${gene}.${2}.aln.after_filter2.txt
	echo "proteinSelected == 'yes' && dnaSelected == 'yes'": $proteinSelected $dnaSelected
	if [[ $proteinSelected == 'yes' && $dnaSelected == 'yes' ]]; then
		echo Also preparing the filtered dna aln...
		seqtk subseq -l 0 \
		${gene}.dna.aln.fasta \
		${gene}.${2}.aln.after_filter2.txt \
		> ${gene}.dna.aln.after_filter2.fasta
	fi
	if [[ $codonSelected == 'yes' ]]; then
		echo Also preparing the filtered codon aln...
		seqtk subseq -l $alnLength \
		codonAln/${gene}.codon.aln.fasta \
		${gene}.${2}.aln.after_filter2.txt \
		> codonAln/${gene}.codon.aln.after_filter2.fasta
		fi
}


trimAln1()	{
   	###########
   	# Function: trims alignment columns to remove low occupancy sites with optrimAl
   	#			(Shee et al ... Polkorny 2020 https://doi.org/10.3389/fpls.2020.00258).
   	#			optrimAl optimizes the gap threshold value in trimAl to obtain the 
   	#			highest proportion of parsimony-informative characters.
   	#			optrimAl scripts, PASTA_taster.sh and optrimAl.R are here:
   	#			https://github.com/keblat/bioinfo-utils/blob/master/docs/advice/scripts/optrimAl.txt
   	#			Altered the code to fit this pipeline - see scripts.
   	# 			Could have transferred the code into this subroutine but didn't.			 
   	#
   	# Input parameters:
   	# $1 = gene name
    # $2 = input fasta file
   	# $3 = residue type for AMAS.py option -d: dna, aa  (NB - AMAS trim -d dna also seems to work with aa!) 
   	# $4 = output file name
   	# $5 = $pathToScripts
   	###########
   	$pathToScripts/PASTA_taster.sh $1 $2 $3 $4 $5
}


trimAln2()	{
   	###########
   	# Function: trims alignment columns to remove low occupancy sites, 
   	#           % occupancy ($3) defined by user
   	#
   	# Input parameters:
   	# $1 = input fasta file
   	# $2 = residue type for AMAS.py option -d: dna, aa  (NB - AMAS trim -d dna also seems to work with aa!)  
   	# $3 = maximum limit of fraction to trim at (e.g. 0.003 = 0.3 % of residues - would remove very rare inserts)
   	# $4 = output file name AFTER trimming.

   	# Note: There are many rare inserts when looking at lots of samples across the Angiosperms!
	# 		Removing these columns is meant to speed up the phylogeny programs, as stated in Fasttree documentation.
	# 		0.001 is equivalent to 1 residue in every 1000 seqs (0.1%) - could use 0.003 for this purpose. 
   	###########
   	if [[ $2 == 'protein' ]]; then 
    	residueType=aa 
    else
    	residueType=$2
    fi
	AMAS.py trim -t $3 \
	-f fasta \
	-d $residueType \
	-i $1 \
	-o $4
	# Stats on the number of bases removed (it will be considerable for large trees with diverse taxa)
	lenLongestGeneAfterTrim=`fastalength $4 | sort -n | tail -n 1 | awk '{print $1}' `
	if [[ -s ${geneId}_aln_summary.log ]]; then
		echo lenLongestGeneAfterTrim: $lenLongestGeneAfterTrim  >> ${geneId}_aln_summary.log
	else
		echo gene Id: $geneId
		echo lenLongestGeneAfterTrim: $lenLongestGeneAfterTrim  > ${geneId}_aln_summary.log
	fi
}


makeGeneTree()	{
	###########
    # Function: makes the gene tree from the available methods in this function 
    # 			from a DNA alignment
 	#
 	# Input parameters:
 	# $1 = residue type: dna, aa or codon 
    # $2 = filtered alignment filename in fasta format
    # $3 = output file prefix including the path to any directory
    # $4 = raxmlng_model e.g. raxmlngModel='GTR+G' - for DNA; for prtoein raxmlngModel='JTT+G'
    # $5 = iqtree2_seq_type  iqTree2SeqType='DNA'; for protein iqTree2SeqType='AA'
    # $6 = fasttreeFlags='-nt -gtr' - for DNA; for protein fasttreeFlags='' - NB - this last flag needs to be last in case it is blank - it needs to be blank for protein analysis)
    ###########

    raxmlngModel=$4
    iqTree2SeqType=$5
    fasttreeFlags=$6

### 5.9.2020 - incorrect logic when different programs have been seelected for DNA and protein:
    echo phyloProgramDNA: $phyloProgramDNA
    echo phyloProgramPROT: $phyloProgramPROT
### I think you need to do the following:
### programToUse= <bring in to method as a varaible e.g. $phyloProgramDNA
### Then it's  if [[ $programToUse == 'fasttree' ]]; then


    if [[ "$phyloProgramDNA" == 'fasttree' || "$phyloProgramPROT" == 'fasttree' ]]; then
		###srun -J ${gene}_make_tree -n 1 \
		echo echo phyloProgramDNA: $phyloProgramDNA
		echo Running fasttree on the alignment...
		$exePrefix fasttree $fasttreeFlags \
    	$2 \
		> ${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
	elif [[ "$phyloProgramDNA" == 'raxml-ng' || "$phyloProgramPROT" == 'raxml-ng' ]]; then
		# RAxML-NG with DNA (uses conventional bs support - how fast is it? it is slow!):
		### RAXML-NG will crash if there is not enough data to ||elize so need to use 1 cpu for small # seqs,
		### and also aln length - keep an eye on this - test larger dataet with 2, 4 cpu
		echo
		echo Running raxml-ng on the DNA alignment...
		cpuGeneTree=1		# NB - at the moment keeping cpu to 1 because very small trees can crash if # cpu is higher - can do ||elisation other ways with RAxML though 
		### 22.10.2020 - was --threads auto{$cpuGeneTree} \; now testing --threads auto{MAX} \ - doesn't work! Says ERROR: Invalid number of threads: %s auto{1}, please provide a positive integer number! but I was!
		### Also the --workers auto option doesn't exist!!!! \
		$exePrefix raxml-ng --threads $cpuGeneTree \
		--redo \
		--all \
		--msa $2 \
		--model $raxmlngModel \
		--prefix ${3}/${gene}.${1}.aln \
		--seed 2 \
		--bs-metric fbp,tbe \
		--bs-trees 100
### NB - changed to bs of 10 - also try: --bs-trees     autoMRE{N} default N=1000
		# NB - it seems that the input file fasta header lines mustn't contain any space-separated fields, just an identfier next to the '>' char!
		# NB --prefix is the output filename prefix; you can also add a path and the output files go into the specified directory
		# --outgroup - can supply a list of outgroups
		# --redo - overwrites files with the output prefix

		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln.raxml.supportFBP \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm ${3}/${gene}.${1}l.aln.raxml.supportFBP
	elif [[ "$phyloProgramDNA" == 'iqtree2' || "$phyloProgramPROT" == 'iqtree2' ]]; then
		echo
		echo Running IQ-Tree on the DNA alignment with these options: -B 1000 ... 
		$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
		-redo \
		--seqtype $iqTree2SeqType \
		-s $2 \
		--prefix ${3}/${gene}.${1}.aln_iqtree \
		-B 1000
		# Command notes:
		# --prefix		you can also add a path and the output files go into the specified directory
		# -B 			ultrafast bootstrap option (NB - the minimum number you can set is 1000!)
		# -alrt 		Removed -alrt 1000 option but it is very fast to compute so no real need - BUT it does add an extra value at the node I think
		# -b			standard nonparametric bootstrap
		# -fast			builds just 2 starting trees and there is no iteration with NNI to try and find higher likelihoods 

		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln_iqtree.contree \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm  ${3}/${gene}.${1}.aln_iqtree.contree
	elif [[ "$phyloProgramDNA" == 'iqtree2-fast-b100' || "$phyloProgramPROT" == 'iqtree2-fast-b100' ]]; then
		echo
		echo Running IQ-Tree on the DNA alignment with these options: -fast, -b 100, -m GTR+F+G7 ...
		$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
		-redo \
		--seqtype $iqTree2SeqType \
		-s $2 \
		--prefix ${3}/${gene}.${1}.aln_iqtree \
		-b 100 \
		-fast \
		-m GTR+F+G7
		
		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln_iqtree.contree \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm  ${3}/${gene}.${1}.aln_iqtree.contree
	elif [[ "$phyloProgramDNA" == 'iqtree2-alrt' || "$phyloProgramPROT" == 'iqtree2-alrt' ]]; then
		echo
		echo Running IQ-Tree on the DNA alignment with these options: -alrt, -m GTR+F+G7 ...
		$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
		-redo \
		--seqtype $iqTree2SeqType \
		-s $2 \
		--prefix ${3}/${gene}.${1}.aln_iqtree \
		-alrt \
		-m GTR+F+G7
		
		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln_iqtree.contree \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm  ${3}/${gene}.${1}.aln_iqtree.contree
	elif [[ "$phyloProgramDNA" == 'iqtree2-B1000-nstep40-nm50' || "$phyloProgramPROT" == 'iqtree2-B1000-nstep40-nm50' ]]; then
		echo
		echo Running IQ-Tree on the DNA alignment with these options: -B 1000, -nstep 40, -nm 50, -m GTR+F+G7 ...
		$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
		-redo \
		--seqtype $iqTree2SeqType \
		-s $2 \
		--prefix ${3}/${gene}.${1}.aln_iqtree \
		-B 1000 \
		-nstep 40 \
		-nm 50 \
		-m GTR+F+G7
		
		# Rename final tree file to a clearer name:			
		cp -p ${3}/${gene}.${1}.aln_iqtree.contree \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm  ${3}/${gene}.${1}.aln_iqtree.contree

#### 17.10.2020 -testing fast bootstraps with raxml!!!!

	# if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng'  || "$phyloProgramPROT" == 'iqtree2' ]]; then
 #     echo Running RAxML on the protein supermatrix...
 #    # RAxML won't run if files already exists from a previous run so remove them here: 
 #    ### Could look for a --force option
 #    if [ -a RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
 #    if [ -a RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
 #    if [ -a RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
 #    if [ -a RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
 #    if [ -a RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
 #    if [ -a ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced ]; then rm ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced; fi

 #    $exePrefix raxmlHPC-PTHREADS-SSE3 -T $cpu \
 #    -f a \
 #    -x 12345 \
 #    -p 12345 \
 #    -# 100 \
 #    -m PROTCATJTT \
 #    -s ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
 #    -n ${fileNamePrefix}__raxmlHPC-PTHREADS-SSE



	else 
		echo "Phylogeny program for gene trees from $1 sequences not detected - will not make gene trees from $1 sequence."
	fi

	### Other gene tree programs could go here
	### NB - will need to confirm whether program outputs the gene tree node support values as a fraction or a percent
	### So far, FASTTREE=fraction,0-1; raxml-ng=integer(0-100); iqtree=integer(0-100)

} # End of makeDnaGeneTree function


filterShortSeqs()	{
	###########
    # Function: filters short sequences from the gene set as defined by fasta file in $2.
 	#
 	# Input parameters:
 	# $1 = input fasta file
    # $2 = minimum length to tolerate
    # $3 = output fasta file
    #
    # Note: has to be done after all trimming steps (otherwise length might reduce further below the minimum value)
    # Note: these short sequences are removed after reporting the alignment stats - so the aln stats may be based on a few more genes - small point
    ###########

    # Can pipe identifier list into seqtk subseq! Saves an output file. 
	fastalength $1 | awk -v minLength=$2 '$1 >= minLength {print $2}' \
	| seqtk subseq -l 0 $1 /dev/fd/0 > "$3"
	numbrSeqs=`cat ${3} | grep '>' | wc -l`
}










############
# Main code:
############
dnaFastaFileForAln="${gene}_dna.fasta"	# Name of the gene-wise DNA fasta input file, ready for aligning; two other input files are possible and dealt with below
if [[ $geneFile != 'use_genewise_files' ]]; then

	echo
	echo
	echo #########################################
	echo 'Organising sequences on a per gene basis (from data input on a per sample basis)...'
	echo #########################################

	# Description of code just below that does the organising of sequences into a gene file:
	# For each gene:
	# 1. find specific gene from all species in turn (NB - sequence is on a single line 
	#    which is why this code can work - bash line length limit much much longer than required here)
	#    NB - using ">" ensures that I don't also grep out a species id: grep -A1 ">$gene" \
	#	 Improvement: current grep now matches whole words so can have gene identifiers
	#                 e.g. 50, 506 and 5064 and they won't pick each other out (BUT these chars count as word endings: '.'  '/' )   
	# 2. print seqs for current gene to a separate file
	# 3. reorganise the fasta header so that species info is next to '>' and the gene id is removed (required to run in coalescence)
	###gene=`echo $geneId | tr -d '\n' `	# Need to remove line return, if present after the first field! Now moved furtehr up.
	#echo Gene: $gene
	cat $geneFile \
	| grep -Ew -A1 ">$gene" \
	| awk '{if($1 ~ /^>/) {print ">" $2} else {print $0} }' \
	| grep -v '^--$' \
	> $dnaFastaFileForAln
	# NB - grep -v '^--$' removes output between each set of rows that grep produces.

	# It is possible that this file ${gene}_dna.fasta will be empty if option -a is set but gene-wise files have been input
	# OR if there are no sequences for this gene
	### Consider to copy the below elif to here as well to capture this - it just means that when no gene seqs exist for
	### a gene the pipeline will exit - this is probablsbly a good thing to happen, no?
	### he decision is to print an error ONLY or also to exit - I'm just ignoring the gene if # seqs < 3 so should do that here as well 
	### -OK this is easy!!!   Also for lien 195 if [ "$numbrSeqs" -gt 3 ]; then --> need a WARNING in an else clause!!!!!!!!!



# If an input fasta file name doesn't exist then the 'modified.fasta' created above will not exist.
# Testing whether the gene-wise file exists here: ${gene}_dna.fasta, 4802_after_treeshrink_dna.fasta or ${gene}_after_filterSeqs_dna.fasta
#	None of them existing triggers the error.
elif [[ ! -s $dnaFastaFileForAln ]]; then

	if [[ -s ${gene}_after_treeshrink_dna.fasta ]]; then 
		dnaFastaFileForAln=${gene}_after_treeshrink_dna.fasta
		echo File for aligning: ${gene}_after_treeshrink_dna.fasta
	elif [[ -s ${gene}_after_filterSeqs_dna.fasta ]]; then
		dnaFastaFileForAln=${gene}_after_filterSeqs_dna.fasta
		echo File for aligning: ${gene}_after_filterSeqs_dna.fasta
	else
		echo "WARNING: the input gene-wise fasta file for this gene does not exist or is empty: $gene - skipping alignment of this gene.
(It indicates that there are no samples for this gene (possibly after a filtering step),
or in gene-wise mode (option -G), the gene list is incompatible with the input gene-wise fasta files.)"
		# NB - acknowledged error above so OK to exit with zero.
		# Anyway it has to be zero to satisfy Slurm --dependancy afterok:$jobId parameter;
		# Only works if $jobId exit code is 0, otherwise would need to use --dependancy afterany:$jobId     
		exit 0
	fi
fi



echo ###############################
echo 'Making the gene alignments...'
echo ###############################
dnaAlnToUse=''		# Stores the alignment file from the chosen aligner - is then used to update/set the filename to use subsequently after EACH filtering or trimming step  
proteinAlnToUse=''
codonAlnToUse=''
if [[ $dnaSelected == 'yes' ]]; then
    if [[ "$alnProgram" == 'mafft' ]]; then 	# If aligners can be set to auto residue detect, can use a generic subR - or brign in a variable.
    	echo
        echo Creating a DNA alignment with MAFFT...
		$exePrefix mafft --thread $cpuGeneTree \
		$mafftAlgorithm \
		--reorder \
		--preservecase \
		$dnaFastaFileForAln \
		> ${gene}.dna.aln.fasta
		# Possible run modes:
		# 1. Progressive method (fast - up to 5k seqs) --retree 1                           			FFT-NS-1/NW-NS-1
		# 2. Progressive method (fast - up to 5k seqs) --retree 2 (better than retree 1)    			FFT-NS-2/NW-NS-2
		# 3. Iterative refinement method (slower, more accurate) --maxiterate 1000 - use this if poss   Builds on FFT-NS-2 --> FFT-NS-i - refinement repeated until no more improvement in the WSP score is made or the number of cycles reaches 1,000.         
		# 4. L-INS-i, E-INS-i, G-INS-i — Iterative refinement methods using WSP and consistency scores - not sure if I need this level (slowest, most accurate)
		# NB - if using Slurm srun with mafft AND fasttree with -c set to > 1, you need to pin the task down with -n flag to 1 task, 
		# otherwise it spawns > 1 runs of the same task.
		### 25.7.2020 - look at docs to see best/max number of thread worthwhile to use - hard code so it only uses max # thread (== to cpu?)
        dnaAlnToUse=${gene}.dna.aln.fasta
    elif [[ "$alnProgram" == 'upp' ]]; then
    	echo
   		echo Creating a DNA alignment with UPP...
   		# Before running check whether pasta has already been run and delete previous files (they can't be overwritten!):
   		### NB - had an issue with removing all files when only testing one of them, in the end my own error I think.
   		### If it happens again (UPP will complain about overwriting files) can use this conditional instead to check all the files in the set at once:
   		### if ls *.dna.upp* 2>&1 >/dev/null; then echo exists; ls -l *.dna.upp* ; fi
   		if [[ -f ${gene}.dna.upp_pasta.fasta ]]; then
   			rm ${gene}.dna.upp_pasta.fasta  ${gene}.dna.upp_pasta.fasttree  ${gene}.dna.upp_alignment.fasta
   			# Don't think this file always exists (for very small datasets):
   			if [[ -f ${gene}.dna.upp_insertion_columns.txt ]]; then rm ${gene}.dna.upp_insertion_columns.txt; fi
   		fi
		run_upp.py -x $cpuGeneTree -M -1 -s $dnaFastaFileForAln -M -1 -o ${gene}.dna.upp
		# Other options to consider:
		# UPP(Fast): run_upp.py -s input.fas -B 100. - what's the -B option??!!
		# -m [dna|rna|amino]
		# -x <cpus> - runs a ||elised version of UPP
		# -M <median_full_length> - specify the median full length of the gene for selecting a set of backbone seqs that are within 25% of that value
		#						    would have to bring in an external file of size - also useful for filterSeqs option 2 - actually '-M -1' will use the median full length of the seqs  
		# Ouput files:
		# _pasta.fasta 	- backbone aln (?)
		# _pasta.Fasttree - backbone tree (?)
		# _alignment.fasta 	- main aln but also contains regions outside the HMM columns which remain unaligned (so do not use!!)
		# _alignment_masked.fasta - masked aln where non-homologous sites in the query set are removed (This is the file to use!!)

		# NB - UPP might not be able to align small gene sets (e.g. from the test data set) - will report this and skip this gene.
		#      UPP doesn’t seem to align two seqs but it does align three seqs!
    	#	   With the -M -1 flag it couldn’t align 5 genes but sucess depends on how complete the sequences are.
		if [[ ! -s ${gene}.dna.upp_alignment.fasta ]]; then 
			echo "ERROR: UPP was not able to align this gene set - skipping alignment of $dnaFastaFileForAln"
			exit 0	# zero allows Slurm to continue with the dependancies
		fi
		mv ${gene}.dna.upp_alignment_masked.fasta ${gene}.dna.aln.fasta 
       	dnaAlnToUse=${gene}.dna.aln.fasta
    fi
fi

if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then 
 	# NB - The protein fasta headers will contain  \[translate(1)\] - fastatranslate (v2.4.x) adds this string to the header.
	# It makes raxml-ng crash so remove it here:
	fastatranslate -F 1  $dnaFastaFileForAln \
	| sed 's/ \[translate(1)\]//' \
	> ${gene}.protein.fasta

	# Detect STOP codons and create STOPS stats, then switch to use file containing 0 or 1 STOP codons:
	$pathToScripts/various_tasks_in_python.py detect_stops ${gene}.protein.fasta  ${gene}.protein
	cp ${gene}.protein.0or1_STOP.fasta ${gene}.protein.fasta
	if [[ ! -s ${gene}.protein.0or1_STOP.fasta ]]; then 
	 	echo "WARNING: after checking for sequences with many STOP codons, this gene set is now empty - skipping alignment of $dnaFastaFileForAln"
	 	exit 0	# zero allows Slurm to continue with the dependancies
	fi

 	if [[ "$alnProgram" == 'mafft' ]]; then 	# If aligners can be set to auto residue detect, then could use a generic subR - or bring in a variable.
		echo
		echo Creating a protein alignment with MAFFT...
		$exePrefix mafft --thread $cpuGeneTree \
		$mafftAlgorithm \
		--reorder \
		--preservecase \
		${gene}.protein.fasta \
		> ${gene}.protein.aln.fasta
		proteinAlnToUse=${gene}.protein.aln.fasta
	elif [[ "$alnProgram" == 'upp' ]]; then
		echo
   		echo Creating a protein alignment with UPP...
   		if [[ -f ${gene}.protein.upp_pasta.fasta ]]; then
   			rm ${gene}.protein.upp_pasta.fasta ${gene}.protein.upp_pasta.fasttree ${gene}.protein.upp_alignment.fasta
   			if [[ -f ${gene}.protein.upp_insertion_columns.txt ]]; then rm ${gene}.protein.upp_insertion_columns.txt; fi
   		fi
		run_upp.py -x $cpuGeneTree -M -1 -m amino -s ${gene}.protein.fasta -o ${gene}.protein.upp
		if [[ ! -s ${gene}.protein.upp_alignment.fasta ]]; then 
			echo "ERROR: UPP was not able to align this gene set - skipping alignment of ${gene}.protein.fasta"
			exit 0
		fi
		mv ${gene}.protein.upp_alignment_masked.fasta ${gene}.protein.aln.fasta
		proteinAlnToUse=${gene}.protein.aln.fasta
	fi

	if [[ $codonSelected == 'yes' ]]; then
		echo
		echo Creating a DNA alignment guided by the protein alignment...
		if [[ ! -d codonAln ]]; then mkdir codonAln; fi
		### Not tested yet - need to check protein fasta header is identical to dna header.
    	pal2nal.pl \
    	-output fasta \
    	${gene}.protein.aln.fasta \
    	$dnaFastaFileForAln \
    	> codonAln/${gene}.codon.aln.fasta
    	### From previous notes:
    	### 11.8.2018 - noticed that where there are small repeats for which one is an insertion, 
    	### mafft or pal2nal can misplace repeat seqs that are not actually themselves repeated in the seq in question - see sg312 as an example.
    	# This file is almost ready for phylogeny and PAML dN/dS analysis.
    	codonAlnToUse=codonAln/${gene}.codon.aln.fasta
	fi
fi





echo ########################################################
echo 'Options for filtering and/or trimming the alignments...'
echo ########################################################
dnaAlnForTree=$dnaAlnToUse				# Variables for applying the correct alignment filename ready 
										# for building the gene trees after the filter and trimming steps.
										# Just need to assign the output filename at the end after each filtering step.
										# These variables ARE required because in the case of filterSeqs1 some
										# gene alignments can be filtered out of the tree building!
proteinAlnForTree=$proteinAlnToUse
codonAlnForTree=$codonAlnToUse
if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then
	if [[ $filterSeqs1 != 'no' ]]; then 
		echo "Filter sequences option 1 - assessing the protein aln for filtering (and the DNA aln if DNA has also been selected)."
		# maxColOccThreshold is required in filterSeqs1 method and in tree making conditionals below.
		# For protein need to divide by 3:
		maxColOccThreshold=`echo $maxColOccThreshold | awk -v maxColOccThreshold=$maxColOccThreshold 'BEGIN{printf "%.0f", maxColOccThreshold/3}' `
		echo maxColOccThreshold for protein: $maxColOccThreshold
							# was 90, now testing 30; seq covrg across region was 84, now testing 28
		# Function parameters: input_fasta_file, residue_type, maxColOccThreshold, minimum seq to tolerate (deprecated)
### 7.10.2020 - I think now you coudl also use *AlnForTree here as well - more logical!!!! Ditto for filterSeqs2
    	filterSeqs1 $proteinAlnToUse protein $maxColOccThreshold 9	# NB - unlike the other functions below not specifying outfile here because we need 1 to 3 output files, protein and/or DNA and/or codon
    	# Alternative for returning filename from function - to stdout (can only echo the filename though):
    	# filteredSeqListToUse=$(filterSeqs1 ${gene}.dna.aln.fasta  dna 90 84)
    	# OR submit the output filename to the function, then I have more control in this code 
    	# here if I want to change filenames and variables - now dong this in the other functions below

		proteinAlnForTree=${gene}.protein.aln.after_filter1.fasta	# Filtered output file - NB - this file will not exist if it doesn't pass the $naxColOcc filtering - OK!
    	echo maxColOcc: $maxColOcc
    	echo numbrSeqs: $numbrSeqs

    	if [[ $dnaSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs1 function when protein has been selected.
			dnaAlnForTree=${gene}.dna.aln.after_filter1.fasta    		
    		echo  dnaAlnForTree: $dnaAlnForTree
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs1 function when protein has been selected.
			codonAlnForTree=codonAln/${gene}.codon.aln.after_filter1.fasta
    		echo  codonAlnForTree: $codonAlnForTree
   		fi
   	elif [[ $filterSeqs2 != 'no' ]]; then
   		echo "Filter sequences option 2 - assessing the protein aln for filtering (and the DNA and codon alns if they have also been selected)."
   		# Function parameters: input_fasta_file, residue_type (REQUIRED HERE???????), percent minimum sequence length to tolerate, outfile_name
   		filterSeqs2 $proteinAlnToUse protein  $filterSeqs2  ${gene}.protein.aln.after_filter2.fasta
   		proteinAlnForTree=${gene}.protein.aln.after_filter2.fasta # Filtered output file

		if [[ $dnaSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs2 function when protein has been selected.
			dnaAlnForTree=${gene}.dna.aln.after_filter2.fasta    		
    		echo  dnaAlnForTree: $dnaAlnForTree
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs2 function when protein has been selected.
			codonAlnForTree=codonAln/${gene}.codon.aln.after_filter2.fasta
    		echo  codonAlnForTree: $codonAlnForTree
   		fi


   	####################################
	# Other filtering steps can go here.
	# NB - need to register the output file in the make_species_trees_pipeline.sh script
	#      to check whether the files need to be deleted from a previous run - see line ~588
	####################################


    else
		echo "No filtering selected."
		# NB - *AlnForTree variables contain the raw aln, ready for trimming 
		#      if those options are selected, otherwise raw aln is ready for 
		#      small seq filtering before tree building.
    fi # End of filtering options


    ##################
    # Trimming options
    ##################
    # NB - conditionals below ensures that trimming is only done after a realignment step following Treeshrink and/or filtering.      
	if [[ "$trimAln1" == 'yes' && "$filterSeqs1" == 'no' && "$filterSeqs2" == 'no' && "$treeshrink" == 'no' ]]; then
   		echo "Trim sequences option 1 (opTrimAL) - assessing the protein aln for trimming (and the DNA and codon alns if they have also been selected)."
### NB - I think this trim option needs to be done before trimming low occupancy columns - but need to run opTrimAL first before I know this for certain
   		# Function parameters: gene_name input_aln_fasta_file, residue_type, outfile_name, $pathToScripts
   		trimAln1 $gene  $proteinAlnForTree  protein  ${gene}.protein.aln.after_trim1.fasta  $pathToScripts
   		proteinAlnForTree=${gene}.protein.aln.after_trim1.fasta

		if [[ $dnaSelected == 'yes' ]]; then
			# Seqs already prepared in filterSeqs* functions when protein has been selected.
			# Still need to trim them here:
			trimAln1 $gene  $dnaAlnForTree  dna  ${gene}.dna.aln.after_trim1.fasta  $pathToScripts
			dnaAlnForTree=${gene}.dna.aln.after_trim1.fasta		#### NB - coudl define variable and use it to call in function first!
    		echo  dnaAlnForTree: $dnaAlnForTree
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs* functions when protein has been selected.
    		# Still need to trim them here:
#### NBNB - need to check the format of codonAlnForTree is consistent i.e. has path or not!!! Have added the path above
			codonAlnForTree=codonAln/${gene}.codon.aln.after_trim1.fasta	# define output file
			trimAln1 $gene  $codonAlnForTree  dna  $codonAlnForTree  $pathToScripts
    		echo  codonAlnForTree: $codonAlnForTree
   		fi
    fi
    #if [[ "$trimAln2" != 'no' && "$filterSeqs1" == 'no' && "$filterSeqs2" == 'no' && "$treeshrink" == 'no' ]]; then
    # 13.10.2020 - the above conditional would ensure that trimAln2 would only operate after the second re-alignment step,
    # but need now to use it before building the gene tree (should be much faster to build tree without them):
    if [[ "$trimAln2" != 'no' ]]; then
   		echo "Trim sequences option 2 - assessing the protein aln for trimming (and the DNA and codon alns if they have also been selected)."
   		# Function parameters: input_fasta_file, residue_type_for_AMAS, maximum_fraction_limit_to_trim  fraction_to_trim_at, outfile_name
### NB - if $proteinAlnForTree is empty after filtering, AMAS will exit with an error - Ok I think for now - BUT
### COULD ADD A CONDITIONAL ROUND THE FUNCTION To ONLY TRIM IF FILE IS > zero byte BUT STILL ASSIGN A NAME TO VARIABLE
### Not that important as most datasets should still have at least one fasta record
   		trimAln2 $proteinAlnForTree protein $trimAln2 ${gene}.protein.aln.after_trim2.fasta
   		proteinAlnForTree=${gene}.protein.aln.after_trim2.fasta
   		echo  proteinAlnForTree: $proteinAlnForTree

		if [[ $dnaSelected == 'yes' ]]; then
			# Seqs already prepared in filterSeqs* functions when protein has been selected.
			# Still need to trim them here:
			trimAln2 $dnaAlnForTree dna $trimAln2 ${gene}.dna.aln.after_trim2.fasta
			dnaAlnForTree=${gene}.dna.aln.after_trim2.fasta		#### NB - coudl define variable and use it to call in function first!
    		echo  dnaAlnForTree: $dnaAlnForTree
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSeqs* functions when protein has been selected.
    		# Still need to trim them here:
#### NBNB - need to check the format of codonAlnForTree is consistent i.e. has path or not!!! Have added the path above
			codonAlnForTree=codonAln/${gene}.codon.aln.after_trim2.fasta
			trimAln2 $codonAlnForTree dna $trimAln2 $codonAlnForTree
    		echo  codonAlnForTree: $codonAlnForTree
   		fi


   	###################################
	# Other trimming steps can go here.
	# NB - need to register the output file in the make_species_trees_pipeline.sh script
	#      to check whether the files need to be deleted from a previous run - see line ~588
	###################################


	#else 
		#echo "No trimming selected."  ### NB - misleading - may have done trim1 but not trimAln2 step!
		# NB - *AlnForTree variables contain either the raw aln or filtered/trimmed aln file,
		#	   ready for small seq filtering before tree building.
	fi


elif [[ $dnaSelected == 'yes' ]]; then
	# NB - most explanations of the workflow logic are presented in the above protein conditional (it's the same)
    if [[ $filterSeqs1 != 'no' ]]; then 
    	echo "Filter sequences option 1 - assessing the DNA aln for filtering"
    	###############maxColOccThreshold=15	# was 90, now testing 30, then 15; seq covrg across region was 84, now testing 27, then 14
    	filterSeqs1 $dnaAlnToUse dna $maxColOccThreshold 14
### ACTUALLY NOW USING UNTRIMMED ALN ANY WAY SO NO NEED TO ASSIGN    
    	dnaAlnForTree=${gene}.dna.aln.after_filter1.fasta
    	echo maxColOcc: $maxColOcc
    	echo numbrSeqs: $numbrSeqs

### COPYING THE TRIMMING CLAUSES FROM THE PROTEIN CLAUSE ABOVE
### UPTOHERE 8.10.2020 transferring code to this clause
	elif [[ $filterSeqs2 != 'no' ]]; then
   		echo "Filter sequences option 2 - assessing the DNA aln for filtering"
   		# Function parameters: input_fasta_file, residue_type (REQUIRED HERE???????), percent minimum sequence length to tolerate, outfile_name
   		filterSeqs2 $dnaAlnToUse  dna  $filterSeqs2  ${gene}.dna.aln.after_filter2.fasta
   		dnaAlnForTree=${gene}.dna.aln.after_filter2.fasta # Filtered output file


   	####################################
	# Other filtering steps can go here.
	# NB - need to register the output file in the make_species_trees_pipeline.sh script
	#      to check whether the files need to be deleted from a previous run - see line ~588
	####################################


    else
		echo "No filtering selected."
		# NB - *AlnForTree variables contain the raw aln, ready for trimming 
		#      if those options are selected, otherwise raw aln is ready for 
		#      small seq filtering before tree building.
    fi # End of filtering options


    ##################
    # Trimming options
    ##################
    # NB - conditionals below ensures that trimming is only done after a realignment step following Treeshrink and/or filtering.      
	if [[ "$trimAln1" == 'yes' && "$filterSeqs1" == 'no' && "$filterSeqs2" == 'no' && "$treeshrink" == 'no' ]]; then
   		echo "Trim sequences option 1 (opTrimAL) - assessing the DNA aln for trimming"
### NB - I think this trim option needs to be done before trimming low occupancy columns - but need to run opTrimAL first before I know this for certain
   		# Function parameters: gene_name input_aln_fasta_file, residue_type, outfile_name, $pathToScripts
   		trimAln1 $gene  $dnaAlnForTree  dna  ${gene}.dna.aln.after_trim1.fasta  $pathToScripts
   		dnaAlnForTree=${gene}.dna.aln.after_trim1.fasta
    fi
    #if [[ "$trimAln2" != 'no' && "$filterSeqs1" == 'no' && "$filterSeqs2" == 'no' && "$treeshrink" == 'no' ]]; then
    if [[ "$trimAln2" != 'no' ]]; then
   		echo "Trim sequences option 2 - assessing the DNA aln for trimming"
   		# Function parameters: input_fasta_file, residue_type_for_AMAS, maximum_fraction_limit_to_trim  fraction_to_trim_at, outfile_name
   		trimAln2 $dnaAlnForTree dna $trimAln2 ${gene}.dna.aln.after_trim2.fasta
   		dnaAlnForTree=${gene}.dna.aln.after_trim2.fasta
   		echo  dnaAlnForTree: $dnaAlnForTree
   	fi


   	###################################
	# Other trimming steps can go here.
	# NB - need to register the output file in the make_species_trees_pipeline.sh script
	#      to check whether the files need to be deleted from a previous run - see line ~588
	###################################
fi


#########################################################
# Note on output alignment filenames before tree building
#########################################################
# At the end of all filtering and trimming steps, the file entering the 
# makeGeneTree subroutine, assess_gene_alignments.sh and run_treeshrink_and_realign.sh
# has to have the same name. Currently this is: ${gene}.protein.aln.for_tree.fasta i.e.
# it needs to be known upfront.
# So just need to move/copy the final output filename and variable to the above name after
# adding any other filtering or trimming step - now doing this below just before tree building.
#########################################################


echo
echo
echo ############################
echo 'Making the gene tree(s)...'
echo ############################
#if [[ $maxColOcc -ge 90 ]]; then 	 	# $maxColOcc created in filterSeqs1 function - not making trees if $maxColOcc is below specified threshold
										# NB - this conditional is only appropriate with the filterSeqs1 function, so insteaad assessing whether a 
										#      filtered fasta file has been created from any filtering function and making tree if it exists. Note that
										#      if $maxColOcc is increased then any previously run files will also be picked up, so best to run in a new
										#	   directory.
### UPDATE THERSE NOTES TO EXPLAIN FINAL LOGIC
										#		Also now thinking of removing *aln.for_tree.fasta and the *.nwk file before running wrapper script - 25.9.2020 - I think I did this - check
										### 10.8.2020 - alternatively, could have a condisitonal to set $maxColOcc to 0 if maxColOcc is not used if filterSeqs is OFF
										###             need to think further - I think you set it to 91 not 0 then all seqs will go through!
###echo dnaAlnForTree: $dnaAlnForTree
###echo proteinAlnForTree: $proteinAlnForTree
if [[ -s $dnaAlnForTree || -s $proteinAlnForTree ]]; then
	###if [ "$numbrSeqs" -gt 3 ]; then 	# $numbrSeqs created in filterSeqs1 function - now calculating just before building tree - better!
		if [[ $dnaSelected == 'yes' ]]; then
			echo dnaAlnForTree: $dnaAlnForTree
			# Function parameters: input_fasta_file, minimum_seq_length_to_tolerate, output_fasta_file_for_tree_building					
			filterShortSeqs $dnaAlnForTree 80 ${gene}.dna.aln.for_tree.fasta	
			echo numbrSeqs: $numbrSeqs
			# Function parameters: residue_type, input_fasta_file, out_dir, raxmlng_model, iqtree2_seq_type, fasttree_flags (NB - this last flag needs to be last - it needs to be blank for protein analysis)
### Still need to confirm file/variable input
			if [ "$numbrSeqs" -gt 3 ]; then 	# Check goes here after ALL filtering steps
				makeGeneTree dna ${gene}.dna.aln.for_tree.fasta '.' 'GTR+G' 'DNA' '-nt -gtr'
			else
				echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
			fi
		fi
		if [[ $codonSelected == 'yes' ]]; then
			filterShortSeqs $codonAlnForTree 80 ${gene}.codon.aln.for_tree.fasta
			echo numbrSeqs: $numbrSeqs
### Still need to confirm file/variable input
			if [ "$numbrSeqs" -gt 3 ]; then
				makeGeneTree codon ${gene}.codon.aln.for_tree.fasta 'codonAln' 'GTR+G' 'DNA' '-nt -gtr'
			else
				echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
			fi
		fi
		if [[ $proteinSelected == 'yes' ]]; then
### Still need to confirm file/variable input
			echo proteinAlnForTree: $proteinAlnForTree
			filterShortSeqs $proteinAlnForTree 28 ${gene}.protein.aln.for_tree.fasta
			echo numbrSeqs: $numbrSeqs
			if [ "$numbrSeqs" -gt 3 ]; then
				makeGeneTree protein ${gene}.protein.aln.for_tree.fasta '.' 'JTT+G' 'AA' ''
			else
				echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
			fi
		fi
	###else
	###	echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
	#####fi # end of block testing $numbrSeqs > 3

fi # end of block testing $maxColOcc threshold
sleep 3		# Helps slurm not to 'fall over' on the Cluster
