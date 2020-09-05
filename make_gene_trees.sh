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



filterSequences()	{
	###########
    # Function: filters sequences from an alignment of any residue type
    #           depending on their coverage across well occupied columns,
    #			as determined in this function.
    #			Also creates filtering stats and prepares files for stats.

    # Input parameters:
    # $1 = alignment filename in fasta format
    # $2 = residue type for AMAS.py option -d: dna, aa  (NB - AMAS trim -d dna also seems to work with aa!) 
    # $3 = maxColOcc threshold - 90 for DNA, 30 for protein
    # $4 = seqLength threshold - 84 for DNA, 28 for protein
    # NB - only importing variables into this function if they vary depending on the residue type of $1,
    #      all other variables a globally available and do not change during running of this script
    ###########
    if [[ $2 -eq 'protein' ]]; then 
    	residueType=aa 
    else
    	residueType=$2
    fi

    # Empty this file so contents is removed from the last run and new content can be appended:
    > filtered_out_alignments.txt


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


# $maxColOcc needs to be at least 150 bp so that even the shortest seqs are ~100bp over the region.
# 150bp is length of the shortest target seqs (except for a few - looking at Angiosperm353 and plastid target seqs)
# Was also going to fail a gene if the $lenLongestGene : $maxColOcc ratio > 3x - would suggest that there is a likely issue with obtaining overlapping seqs:
# ratioLenOcc=`echo $length $maxOcc | awk '{printf "%.0f", $1/$2}' `
# However, then I'm relying on $lenLongestGene again and the are some issues with that (see above)   
# Just testing $maxColOcc should be OK:
if [[ $maxColOcc -ge $3 ]]; then

	# Calculate the number of parsimonious columns (-p option) for the columns found from $maxColOcc:
	AMAS.py trim -f fasta -d dna -t $fractnMaxColOcc -p -i ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta -o ${gene}_${2}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta
	maxParsCols=`fastalength ${gene}_${2}_aln_AMAS_trim_-p_${fractnMaxColOcc}.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
	echo maxParsCols: $maxParsCols  >> ${geneId}_aln_summary.log


	# Now find the lengths of each sequence in bases only (not dashes!)
	# and create a list of sequences that are > $lenLongestGene:
	# NB - this code works because, luckily for this case, fastalength ignores dash characters (unlike seqtk used above)!
	#### NB 9.12.2019 - but in the case of $maxColOcc do/should I need to ignore dash chars????? I think it needs to be theobsolute length plus
	###			31.1.2020 - no I think it's OK as it is.		
	fastalength ${gene}_${2}_aln_AMAS_trim_${fractnMaxColOcc}.fasta 2>/dev/null \
	| awk -v LENG=$maxColOcc -v FRACTN=$fractnAlnCovrg -v seqLenThreshold=$4 '{if( ($1/LENG) >= FRACTN && $1 >= seqLenThreshold) {print $2} }' \
	> ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
	# Note if .fasta file is empty, fastalength still exits normally.

### 19.8.2020 - consider to do this: else if < $maxColOcc just remove seqs below seqLenThreshold abd continue with same file as above but no filtering with $fractnAlnCovrg

	# Get seqs with ideal coverage, except, if this seq list file made just above 
	# only has 3 or less seqs then there is also no point in making a gene tree.
	numbrSeqs=`cat ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt | wc -l `
	echo numbrSeqs for tree = $numbrSeqs
	if [ "$numbrSeqs" -gt 3 ]; then 		# Conditonal here not essential - but don't know how AMAS.py behaves with empty input files so leaving in here
											# Actually now using to print out filename if it fails
											### 27.8.2020 - can improve this check; if possible leave it to the main code to not build tree if < 3 seqs after all filtering checks 
											### For now have just repeated the identical error message at this stage as well
### 3.9.2020 - now consider to drop seqs < x residues long
### BUT wait - go throught logic of this next section  and work out what residue tyep I'm workign on
### It's Ok - just need to use the residueType var in AMAS -d flag but it seems to work with it beign DNA
		seqtk subseq -l $alnLength \
		$1 \
		${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
		> ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
		# NB - if .txt and .fasta are empty seqtk subseq still exits normally.

		# Record stats (NB - before trimming)
		lenLongestGene=`fastalength  $1 | sort -n | tail -n 1 | awk '{print $1}' `
		echo lenLongestGene: $lenLongestGene  >> ${geneId}_aln_summary.log

		# Median value of sequence length for current gene:
		medianPoint=`fastalength $1 | awk 'END {printf "%.0f" , NR/2}' `
		medianGeneLength=`fastalength $1 | awk '{print $1}' | sort -n | head -n $medianPoint | tail -n 1 `
		echo medianGeneLength: $medianGeneLength  >> ${geneId}_aln_summary.log


		# Trim columns to remove rare inserts, say ones filled only 0.1% with residues.
		# There are many rare inserts when looking at lots of samples across the Angiosperms!
		# Removing these columns is meant to speed up the phylogeny programs, as stated in Fasttree documentation.
		# 0.001 is equivalent to 1 residue in every 1000 seqs (0.1%). - 9.5.2020 - change to 0.003 
		AMAS.py trim -t 0.003 \
		-f fasta \
		-d dna \
		-i ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
		-o ${gene}.${2}.aln.for_tree.fasta
		### 12.8.2020 - changed the filename:
		###-o ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta
		### NB - 23.5.2020 - switched in this file now for steps below and also for preparing the Newick files for Astral and supermatrix method.
		

		# Stats on the number of bases removed (it will be considerable for large trees with diverse taxa)
		lenLongestGeneAfterTrim=`fastalength ${gene}.${2}.aln.for_tree.fasta | sort -n | tail -n 1 | awk '{print $1}' `
		echo lenLongestGeneAfterTrim: $lenLongestGeneAfterTrim  >> ${geneId}_aln_summary.log


		# Extra steps here to do if DNA and/or codon is also selected.
		# Prepare protein and codon alns using result of the filtering on the protein aln using this list: 
		# ${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
		echo "proteinSelected == 'yes' && dnaSelected == 'yes'": $proteinSelected $dnaSelected
		if [[ $proteinSelected == 'yes' && $dnaSelected == 'yes' ]]; then
			echo entered here for making trimmed dna aln
			seqtk subseq -l $alnLength \
			${gene}.dna.aln.fasta \
			${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
			> ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
			AMAS.py trim -t 0.003 \
			-f fasta \
			-d dna \
			-i ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
			-o ${gene}.dna.aln.for_tree.fasta
			### 12.8.2020 - changed the filename:
			###-o ${gene}_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta

			### if aligning the high occupancy columns do above for DNA and protein here - add them to a new directory
			### if protein and DNA only, need to add above with a condtional for high occ columns
		fi
		if [[ $codonSelected == 'yes' ]]; then
			seqtk subseq -l $alnLength \
			codonAln/${gene}.codon.aln.fasta \
			${gene}_${2}_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
			> codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
			AMAS.py trim -t 0.003 \
			-f fasta \
			-d dna \
			-i codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
			-o codonAln/${gene}.codon.aln.for_tree.fasta
			### 12.8.2020 - changed the filename:
			###-o codonAln/${gene}_codon_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta

			### if aligning the high occupancy columns do above for DNA here - add them to a new directory
		fi
	else
		# Printing out the alignments with < 4 sequences for assessing further (can concatenate them together)
		# NB - this file is produced only for the protein seqs, even if DNA is selected; but is produced for DNA if DNA only has been selected.
		echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
		echo $1 >> filtered_out_alignments.txt
	fi
else
	# Printing out the alignments with little sequence overlap for assessing further (can concatenate them together)
	echo $1 >> filtered_out_alignments.txt 
fi 
} # End of filterSequences function


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


    if [[ "$phyloProgramDNA" == 'fastttree' || "$phyloProgramPROT" == 'fasttree' ]]; then
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
		$exePrefix raxml-ng --threads $cpuGeneTree \
		--redo \
		--all \
		--msa $2 \
		--model $raxmlngModel \
		--prefix ${3}/${gene}.${1}.aln \
		--seed 2 \
		--bs-metric fbp,tbe \
		--bs-trees 100
		# NB --prefix is the output filename prefix; you can also add a path and the output files go into the specified directory
		# --outgroup - can supply a list of outgroups
		# --redo - overwrites files with the output prefix

		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln.raxml.supportFBP \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm ${3}/${gene}.${1}l.aln.raxml.supportFBP
	elif [[ "$phyloProgramDNA" == 'iqtree2' || "$phyloProgramPROT" == 'iqtree2' ]]; then
		echo
		echo Running IQ-Tree on the DNA alignment...
		$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
		-redo \
		--seqtype $iqTree2SeqType \
		-s $2 \
		--prefix ${3}/${gene}.${1}.aln_iqtree \
		-B 1000
		### 9.6.2020 - removign this option for the moment incase it taeks up time: -alrt 100
		# --prefix - you can also add a path and the output files go into the specified directory

		# Rename final tree file to a clearer name:
		cp -p ${3}/${gene}.${1}.aln_iqtree.contree \
		${3}/${gene}_${1}_gene_tree_USE_THIS.nwk
		rm  ${3}/${gene}.${1}.aln_iqtree.contree
	else 
		echo "Phylogeny program for gene trees from $1 sequences not detected - will not make gene trees from $1 sequence."
	fi

	### Other gene tree programs could go here

} # End of makeDnaGeneTree function










# Main code:
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
		echo "ERROR: the input gene-wise fasta file for this gene does not exist or is empty: $gene - skipping alignment of this gene.
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
dnaAlnToUse=''		# Stores the alignment file from the chosen aligner - useful for specifying the filename later in script 
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
   			rm ${gene}.dna.upp_pasta.fasta  ${gene}.dna.upp_pasta.fasttree  ${gene}.dna.upp_alignment_masked.fasta
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
		# _alignment.fasta 	- main aln
		# _alignment_masked.fasta - masked aln where non-homologous sites in the query set are removed

		# NB - UPP might not be able to align small gene sets (e.g. from the test data set) - will report this and skip this gene.
		#      UPP doesn’t seem to align two seqs but it does align three seqs!
    	#	   With the -M -1 flag it couldn’t align 5 genes but sucess depends on how complete the sequences are.
		if [[ ! -s ${gene}.dna.upp_alignment.fasta ]]; then 
			echo "ERROR: UPP was not able to align this gene set - skipping alignment of $dnaFastaFileForAln"
			exit 0	# zero allows Slurm to continue with the dependancies
		fi
		mv ${gene}.dna.upp_alignment.fasta ${gene}.dna.aln.fasta 
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


 	if [[ "$alnProgram" == 'mafft' ]]; then 	# If aligners can be set to auto residue detect, can use a generic subR - or brign in a variable.
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
   			rm ${gene}.protein.upp_pasta.fasta ${gene}.protein.upp_pasta.fasttree ${gene}.protein.upp_alignment_masked.fasta
   			if [[ -f ${gene}.protein.upp_insertion_columns.txt ]]; then rm ${gene}.protein.upp_insertion_columns.txt; fi
   		fi
		run_upp.py -x $cpuGeneTree -M -1 -m amino -s ${gene}.protein.fasta -o ${gene}.protein.upp
		if [[ ! -s ${gene}.protein.upp_alignment.fasta ]]; then 
			echo "ERROR: UPP was not able to align this gene set - skipping alignment of ${gene}.protein.fasta"
			exit 0
		fi
		mv ${gene}.protein.upp_alignment.fasta ${gene}.protein.aln.fasta
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
filteredDnaAlnToUse=''
filteredProteinAlnToUse='' 
filteredCodonAlnToUse=''
trimmedDnaAlnToUse=''
trimmedProteinAlnToUse='' 
trimmedCodonAlnToUse=''
dnaAlnForTree=$filteredDnaAlnToUse				# Variables for applying the correct alignment filename ready for building the gene trees.
												# Just need to assign the output filename at the end after each filtering step.
												# These variables are not strictly required now but will keep using them.
												### 24.8.2020 - not sure now why I need to assign $filtered*AlnToUse variables to therse vars here!!
proteinAlnFOrTree=$filteredProteinAlnToUse
codonAlnForTree=$filteredCodonAlnToUse
if [[ $proteinSelected == 'yes' || $codonSelected == 'yes' ]]; then
	if [[ $filterSeqs1 != 'no' ]]; then 
		echo "Filter sequences option 1 - assessing the protein aln for filtering (even if DNA has also been selected)."
		maxColOccThreshold=10	# Required in filterSequences method and in tree making conditionals below
								# was 90, now testing 30; seq covrg across region was 84, now testing 28
		# Function parameters: input_fasta_file, residue_type, maxColOccThreshold, minimum seq to tolerate
    	filterSequences $proteinAlnToUse protein $maxColOccThreshold 9
    	### Alternative for returning filename from function - to stdout (can only echo the filename though):
    	### filteredSeqListToUse=$(filterSequences ${gene}.dna.aln.fasta  dna 90 84)
    	### OR submit the output filename to the function, then I have more control in this code here if I want to change filenames and variables
    	proteinAlnForTree=${gene}.protein.aln.for_tree.fasta	# Filtered output file
    	echo maxColOcc: $maxColOcc
    	echo numbrSeqs: $numbrSeqs

    	if [[ $dnaSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSequences function when protein has been selected.
    		dnaAlnForTree=${gene}.dna.aln.for_tree.fasta
    		echo  dnaAlnForTree: $dnaAlnForTree
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		# Seqs already prepared in filterSequences function when protein has been selected.
    		codonAlnForTree=codonAln/${gene}.codon.aln.for_tree.fasta
    		echo  codonAlnForTree: $codonAlnForTree
   		fi

    ### filterSeqs2 option - assessing seq covrg length based on reference target length.

    ### trimAln() option 1 function here if selected - for trimming at 0.003 or higher where some datasets are longer than paftol
    ### Could use AMAS.py for this

    # alnTrim option 2 - this could be opTrimAI
	else
		echo "No filtering selected."
		# if no filtering or trimming has been selected assign filename to variable ready for tree building
		proteinAlnForTree=$proteinAlnToUse
		# Count $numbrSeqs in this file (already done for filtered file in filterSeqs1 function but not for the original aln):
		numbrSeqs=`cat $proteinAlnToUse | grep '>' | wc -l`
		echo numbrSeqs: $numbrSeqs

		# Copy the alignment file to a name used for subsequent scripts (has to be a fixed name unfortunately I think).
		# Also need to prepare the same for dna and codon alns if selected (same logic as above).
		cp -p $proteinAlnToUse ${gene}.protein.aln.for_tree.fasta
		if [[ $dnaSelected == 'yes' ]]; then
			dnaAlnForTree=$dnaAlnToUse	
    		cp -p $dnaAlnToUse ${gene}.dna.aln.for_tree.fasta
    	fi
    	if [[ $codonSelected == 'yes' ]]; then
    		codonAlnForTree=$codonAlnToUse	
    		cp -p $codonAlnToUse codonAln/${gene}.codon.aln.for_tree.fasta
   		fi	
    fi
elif [[ $dnaSelected == 'yes' ]]; then
    if [[ $filterSeqs1 != 'no' ]]; then 
    	echo "Filter sequences option 1 - assessing the DNA aln for filtering."
    	maxColOccThreshold=30	# was 90, now testing 30; seq covrg across region was 84, now testing 28
    	filterSequences $dnaAlnToUse dna $maxColOccThreshold 27
    	dnaAlnForTree=${gene}.dna.aln.for_tree.fasta
    	echo maxColOcc: $maxColOcc
    	echo numbrSeqs: $numbrSeqs

   	### filterSeqs2 option - assessing seq covrg length based on reference target length.

    ### trimAln() option 1 function here if selected - for trimming at 0.003 or higher where some datasets are longer than paftol
    ### Could use AMAS.py for this

    # alnTrim option 2 - this could be opTrimAI
    else
    	echo "No filtering selected."
		dnaAlnForTree=$dnaAlnToUse
		numbrSeqs=`cat $dnaAlnToUse | grep '>' | wc -l`
		echo numbrSeqs: $numbrSeqs
		cp -p $dnaAlnToUse ${gene}.dna.aln.for_tree.fasta
	fi
fi


#########################################################
# Note on output alignment filenames before tree building
#########################################################
# At the end of all filtering and trimming steps, the file entering the 
# makeGeneTree subroutine, assess_gene_alignments.sh and run_treeshrink_and_realign.sh
# has to have the same name. Currently this is: ${gene}.protein.aln.for_tree.fasta i.e.
# it needs to be known upfront.
# So just need to move/copy the final output filename and variable to the above name after
# adding any other filtering or trimming step
#########################################################


echo
echo
echo ############################
echo 'Making the gene tree(s)...'
echo ############################
#if [[ $maxColOcc -ge 90 ]]; then 	 	# $maxColOcc created in filterSequences function - not making trees if $maxColOcc is below specified threshold
										# NB - this conditional is only appropriate with the filterSequences function, so insteaad assessing whether a 
										#      filtered fasta file has been created from any filtering function and making tree if it exists. Note that
										#      if $maxColOcc is increased then any previously run files will also be picked up, so best to run in a new
										#	   directory.
										#		Also now thinking of removing *aln.for_tree.fasta and the *.nwk file before running wrapper script.
										### 10.8.2020 - alternatively, could have a condisitonal to set $maxColOcc to 0 if maxColOcc is not used if filterSeqs is OFF
										###             need to think further - I think you set it to 91 not 0 then all seqs will go through!
if [[ -s $dnaAlnForTree || -s $proteinAlnForTree ]]; then
	if [ "$numbrSeqs" -gt 3 ]; then 	# $numbrSeqs created in filterSequences function - what happens if filtering is OFF?!
		if [[ $dnaSelected == 'yes' ]]; then
			# Function parameters: residue_type, input_fasta_file, out_dir, raxmlng_model, iqtree2_seq_type, fasttree_flags (NB - this last flag needs to be last - it needs to be blank for protein analysis)
			echo dnaAlnForTree: $dnaAlnForTree
			makeGeneTree dna $dnaAlnForTree '.' 'GTR+G' 'DNA' '-nt -gtr'
		fi
		if [[ $codonSelected == 'yes' ]]; then
			makeGeneTree codon $codonAlnForTree 'codonAln' 'GTR+G' 'DNA' '-nt -gtr'
		fi
		if [[ $proteinSelected == 'yes' ]]; then 
			makeGeneTree protein $proteinAlnForTree '.' 'JTT+G' 'AA' ''
		fi
	else
		echo "WARNING: Not able to build a tree for this gene: $gene (less than four sequences)"
	fi # end of block testing $numbrSeqs > 3

fi # end of block testing $maxColOcc threshold

# If there are no .nwk trees built, need to report that and exit with message - reporting that in the next step (make_species_trees.sh)
sleep 3		# Helps slurm not to 'fall over' on the Cluster
