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

# Convert $emptyMatchStateFractn to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `


echo Inside slurm_make_gene_trees.sh script:
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


# Create gene alignment summary file or empty an existing one:
>${geneId}_aln_summary.log
echo gene Id: ${geneId} >> ${geneId}_aln_summary.log


if [[ $geneFile != 'use_genewise_files' ]]; then

	echo
	echo
	echo ###############################################
	echo 'Step 1 - Organising sequences on a per gene basis (from data input on a per sample basis)...'
	echo ###############################################

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
	> ${gene}_dna.fasta
	# NB - grep -v '^--$' removes output between each set of rows that grep produces.

	# It is possible that this file ${gene}_dna.fasta will be empty if option -a is set but gene-wise files have been input
	# OR if there are no sequences for this gene
	### Consider to copy the below elif to here as well to capture this - it just means that when no gene seqs exit for
	### a gene the pipeline will exit - this is probablsbly a good thing to happen, no?
	### he decision is to print an error ONLY or also to exit - I'm just ignoring the gene if # seqs < 3 so should do that here as well 
	### -OK this is easy!!!   Also for lien 195 if [ "$numbrSeqs" -gt 3 ]; then --> need a WARNING in an else clause!!!!!!!!!ß



# If an input fasta file name doesn't exist then the 'modified.fasta' created above will not exist.
# Testing whether file exists here:
elif [[ ! -s ${gene}_dna.fasta ]]; then
	echo "ERROR: this input gene-wise fasta file does not exist or is empty: ${gene}_dna.fasta 
It indicates that there are no samples with this gene, or in gene-wise mode (option -G), the gene list is incompatible with the input gene-wise fasta files"
	exit 1
fi


echo ########################################
echo 'Step 2 - Making the gene alignments...'
echo ########################################

#cpu=4	# NBNB - also needs changing in wrapper script - slurm command; 4cpu's work if srun is pinned down to -n 1.
		# Now using $cpuGeneTree													
# RUNTIME: up to 10 mins and 2MB mem (so far)	  -         --nodelist=kppgenomics01.ad.kew.org
# sbatch -J gene${gene}_make_aln_and_tree -p main  -t 0-2:00 -c $cpu --mem=1000 -o ${gene}_align_make_trees.log   -e ${gene}_align_make_trees.log_err   --wrap "

# NB - when using srun with mafft AND fasttree with -c set to > 1, you need to pin the task down with -n flag to 1 task, 
# otherwise it spawns > 1 runs of the same task.
###srun -J ${gene}_make_aln -n 1 \
echo Running mafft on the DNA alignment...
$exePrefix mafft --thread $cpuGeneTree \
$mafftAlgorithm \
--reorder \
--preservecase \
${gene}_dna.fasta \
> ${gene}_mafft_dna_aln.fasta
# Possible run modes:
# 1. Progressive method (fast - up to 5k seqs) --retree 1                           			FFT-NS-1/NW-NS-1
# 2. Progressive method (fast - up to 5k seqs) --retree 2 (better than retree 1)    			FFT-NS-2/NW-NS-2
# 3. Iterative refinement method (slower, more accurate) --maxiterate 1000 - use this if poss   Builds on FFT-NS-2 --> FFT-NS-i - refinement repeated until no more improvement in the WSP score is made or the number of cycles reaches 1,000.         
# 4. L-INS-i, E-INS-i, G-INS-i — Iterative refinement methods using WSP and consistency scores - not sure if I need this level (slowest, most accurate)


### Alignment with UPP
###run_upp.py -s ${gene}_dna.fasta -o ${gene}
# Other options to consider:
# UPP(Fast): run_upp.py -s input.fas -B 100
# -m [dna|rna|amino]
# -x <cpus> - runs a ||elised version of UPP
# Ouput files:
# _pasta.fasta 	- backbone aln (?)
# _pasta.Fasttree - backbone tree (?)
# _alignment.fasta 	- main aln
# _alignment_masked.fasta - masked aln where non-homologous sites in the query set are removed






# Find the length of the alignment to calculate the coverage of each sequence:
# NB - this code works because seqtk also counts dashes - c.f. fastalength which does not.
# NBNB - putting seq on a single line then can count length of line i.e seq (incl. dashes)
# NBNBNB - /dev/fd/0 - have to explicitly redirect into seqtk on Macbook Darwin.
# NBNBNBNB - Not using this full aln value now - using length of the longest gene in the aln (see just below).
alnLength=`cat ${gene}_mafft_dna_aln.fasta | seqtk seq -l 0 /dev/fd/0 | tail -n1 | awk '{print length($1)}' `
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
AMAS.py trim -f fasta -d dna -t ${fractnMaxColOcc} -i ${gene}_mafft_dna_aln.fasta -o ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta
maxColOcc=`fastalength  ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
# NB - some seqs will have zero length seqs but the fastalength warning can go to /dev/null; using this command to check the occupancy but nothing more 
echo maxColOcc: $maxColOcc  >> ${geneId}_aln_summary.log

# $maxColOcc needs to be at least 150 bp so that even the shortest seqs are ~100bp over the region.
# 150bp is length of the shortest target seqs (except for a few - looking at Angiosperm353 and plastid target seqs)
# Was also going to fail a gene if the $lenLongestGene : $maxColOcc ratio > 3x - would suggest that there is a likely issue with obtaining overlapping seqs:
# ratioLenOcc=`echo $length $maxOcc | awk '{printf "%.0f", $1/$2}' `
# However, then I'm relying on $lenLongestGene again and the are some issues with that (see above)   
# Just testing $maxColOcc should be OK:
if [[ $maxColOcc -ge 90 ]]; then

	# Calculate the number of parsimonious columns (-p option) for the columns found from $maxColOcc:
	AMAS.py trim -f fasta -d dna -t $fractnMaxColOcc -p -i ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta -o ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}_-p.fasta
	maxParsCols=`fastalength ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}_-p.fasta 2>/dev/null | sort -k1n | tail -n 1 | awk '{print $1}' `
	echo maxParsCols: $maxParsCols  >> ${geneId}_aln_summary.log


	# Now find the lengths of each sequence in bases only (not dashes!)
	# and create a list of sequences that are > $lenLongestGene:
	# NB - this code works because, luckily for this case, fastalength ignores dash characters (unlike seqtk used above)!
	#### NB 9.12.2019 - but in the case of $maxColOcc do/should I need to ignore dash chars????? I think it needs to be theobsolute length plus...   31.1.2020 - no I think it's OK as it is
	# 20.7.2020 - now also filtering out below a minimum length - 85 bases seems OK to me 		
	#fastalength ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta 2>/dev/null \
	#| awk -v LENG=$maxColOcc -v FRACTN=$fractnAlnCovrg '{if( ($1/LENG) >= FRACTN ) {print $2} }' \
	#> ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
	fastalength ${gene}_mafft_dna_aln_AMAS_${fractnMaxColOcc}.fasta 2>/dev/null \
	| awk -v LENG=$maxColOcc -v FRACTN=$fractnAlnCovrg '{if( ($1/LENG) >= FRACTN && $1 >= 85) {print $2} }' \
	> ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt



	echo
	echo
	echo ###################################
	echo 'Step 3 - Making the gene trees...'
	echo ###################################
	# Get seqs with ideal coverage, except, if this seq list file made just above 
	# only has 3 or less seqs then there is also no point in making a gene tree.
	numbrSeqs=`cat ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt | wc -l `
	echo numbrSeqs for tree = $numbrSeqs
	if [ "$numbrSeqs" -gt 3 ]; then

		seqtk subseq -l $alnLength \
		${gene}_mafft_dna_aln.fasta \
		${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
		> ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta

		# Record stats:
		lenLongestGene=`fastalength  ${gene}_mafft_dna_aln.fasta | sort -n | tail -n 1 | awk '{print $1}' `
		echo lenLongestGene: $lenLongestGene  >> ${geneId}_aln_summary.log


		# Trim columns to remove rare inserts, say ones filled only 0.1% with residues.
		# There are many rare inserts when looking at lots of samples across the Angiosperms!
		# Removing these columns is meant to speed up the phylogeny programs, as stated in Fasttree documentation.
		# 0.001 is equivalent to 1 residue in every 1000 seqs (0.1%). - 9.5.2020 - change to 0.003 
		AMAS.py trim -t 0.003 \
		-f fasta \
		-d dna \
		-i ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
		-o ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta
		### NB - 23.5.2020 - switched in this file now for steps below and also for preparing the Newick files for Astral and supermatrix method.
		

		# Stats on the number of bases removed (it will be considerable for large trees with diverse taxa)
		lenLongestGeneAfterTrim=`fastalength ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta | sort -n | tail -n 1 | awk '{print $1}' `
		echo lenLongestGeneAfterTrim: $lenLongestGeneAfterTrim  >> ${geneId}_aln_summary.log


		


### After any trimming need to end up with a common filename e.g.:
### ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trim.fasta
### Don't I just assign it to a varialbe and use that?!!!!!!
### 11.7.2020 - this is the point common to any trimming step before.
		if [ "$phyloProgramDNA" == 'fasttree' ]; then
			###srun -J ${gene}_make_tree -n 1 \
			echo
			echo Running fasttree on the DNA alignment...
			$exePrefix fasttree -nt \
    		-gtr \
    		${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
			> ${gene}_dna_gene_tree_USE_THIS.nwk
		elif [ "$phyloProgramDNA" == 'raxml-ng' ]; then

			# RAxML-NG with DNA (uses conventional bs support - how fast is it? it is slow!):
			### RAXML-NG will crash if there is not enough data to ||elize so need to use 1 cpu for small # seqs,
			### and also aln length - keep an eye on this - test larger dataet with 2, 4 cpu
			echo
			echo Running raxml-ng on the DNA alignment...
			$exePrefix raxml-ng --threads $cpuGeneTree \
			--redo \
			--all \
			--msa ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
			--model GTR+G \
			--prefix ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg \
			--seed 2 \
			--bs-metric fbp,tbe \
			--bs-trees 100
			# NB --prefix is the output filename prefix
			# --outgroup - can supply a list of outgroups
			# --redo - overwrites files with the output prefix

			# Rename final tree file to a clearer name:
			cp -p ${gene}_mafft_dna_aln_ovr60pc_aln_covrg.raxml.supportFBP \
			${gene}_dna_gene_tree_USE_THIS.nwk
			rm ${gene}_mafft_dna_aln_ovr60pc_aln_covrg.raxml.supportFBP
		elif [ "$phyloProgramDNA" == 'iqtree2' ]; then

			echo
			echo Running IQ-Tree on the DNA alignment...
			$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
			-redo \
			--seqtype DNA \
			-s ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
			--prefix ${gene}_mafft_dna_aln_iqtree \
			-B 1000
			### 9.6.2020 - removign this option for the moment incase it taeks up time: -alrt 100

			# Rename final tree file to a clearer name:
			cp -p ${gene}_mafft_dna_aln_iqtree.contree \
			${gene}_dna_gene_tree_USE_THIS.nwk
			rm  ${gene}_mafft_dna_aln_iqtree.contree
		else 
			echo "Phylogeny program for gene trees from DNA sequences not detected - will not make gene trees from DNA sequence."
		fi


		if [[ "$phyloProgramPROT" == 'fasttree' || "$phyloProgramPROT" == 'raxml-ng' || "$phyloProgramPROT" == 'iqtree2'  ]]; then

			echo ########################################
			echo 'Creating proteins and aligning them...'  
			echo ########################################
			# Mirroring the above commands for DNA - all conditionals already sorted on the DNA aln.
			### NB - 9.5.2020 - another possibility is to:
			### 1.seqtk subseq with the ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt file on ${gene}_protein.fasta
			### 2.THEN translate output and MAFFT aln on the filtered data set.
			### 3.BUT in the long run it would be better to align DNA against the protein aln so this step would then come first
			###   and filtering would be best done on the protein.

			# NB - The protein fasta headers will contain  \[translate(1)\] - fastatranslate (v2.4.x) adds this string to the header.
			# It makes raxml-ng crash so remove it here:
			fastatranslate -F 1  ${gene}_dna.fasta \
			| sed 's/ \[translate(1)\]//' \
			> ${gene}_protein.fasta
			### 27.1.2020 - should be removing seqs with > 1-3 STOP codons - there are some!!!!
			echo Running mafft on the protein alignment...
			$exePrefix mafft --thread $cpuGeneTree \
			$mafftAlgorithm \
			--reorder \
			--preservecase \
			${gene}_protein.fasta \
			> ${gene}_mafft_protein_aln.fasta
			
			# Now get the relevant protein seqs:
			### 9.5.2020 - might be safer to to use -l 0 - but check -l 0 means all seq in on one line
			seqtk subseq -l $alnLength \
			${gene}_mafft_protein_aln.fasta \
			${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt \
			>${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta
			### This file also needs to be trimmed which will result in slightly different versions w.r.t. DNA so can't use pal2nal on this.


			# Trim columns to remove rare inserts, say ones filled only 0.1% with residues.
			# There are many rare inserts when looking at lots of samples across the Angiosperms!
			# Removing these columns is meant to speed up the phylogeny programs, as stated in Fasttree documentation.
			# 0.001 is equivalent to 1 residue in every 1000 seqs (0.1%). - 9.5.2020 - change to 0.003 
			AMAS.py trim -t 0.003 \
			-f fasta \
			-d dna \
			-i ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.fasta \
			-o ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta
			### NB - 23.5.2020 - switched in this file now for steps below and also for preparing the Newick files for Astral and supermatrix method.


			# Make a tree with protein seqs if requested:
			if [ "$phyloProgramPROT" == 'fasttree' ]; then
				###  Gives a segmentation fault with -wag flag but the protein is default (JTT?)
				echo
				echo Running fasttree on the protein alignment...
				$exePrefix fasttree ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
				> ${gene}_protein_gene_tree_USE_THIS.nwk
			elif [ "$phyloProgramPROT" == 'raxml-ng' ]; then

				# RAxML-NG with protein (uses conventional bs support - how fast is it? it is slow!):
				echo
				echo Running raxml-ng on the protein alignment...
				$exePrefix raxml-ng --threads $cpuGeneTree \
				--redo \
				--all \
				--msa ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
				--model JTT+G \
				--prefix ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg \
				--seed 2 \
				--bs-metric fbp,tbe \
				--bs-trees 100

				# Rename final tree file to a clearer name:
				cp -p ${gene}_mafft_protein_aln_ovr60pc_aln_covrg.raxml.supportFBP \
				${gene}_protein_gene_tree_USE_THIS.nwk
				rm ${gene}_mafft_protein_aln_ovr60pc_aln_covrg.raxml.supportFBP
			elif [ "$phyloProgramDNA" == 'iqtree2' ]; then

				echo
				echo Running IQ-Tree on the protein alignment...
				$exePrefix iqtree2 -T AUTO -ntmax $cpuGeneTree \
				-redo \
				--seqtype AA \
				-s ${gene}_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta \
				--prefix ${gene}_mafft_protein_aln_iqtree \
				-B 100 \
				-alrt 100

				# Rename final tree file to a clearer name:
				cp -p ${gene}_mafft_protein_aln_iqtree.contree \
				${gene}_protein_gene_tree_USE_THIS.nwk
				rm  ${gene}_mafft_protein_aln_iqtree.contree
			else 
				echo "Phylogeny program for gene trees from protein sequences not detected - will not make gene trees from protein sequence."
			fi


			### Other gene tree programs/tests to go here

			### Treeshrink - can shrink the tree here but then need to repeat the aln and tree again - i.e.:
			### seqtk again on *_dna.fasta with   ${gene}_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg.txt
			### MAFFT aln + AMAS trim
			### redo tree but need to know which program to use (maybe this code can go into a new script)




		fi
	fi

fi # end of block testing $maxColOcc threshold

# If there are no .nwk trees built, need to report that and exit with message - reporting that in the next step (make_species_trees.sh)
sleep 3		# Helps slurm not to 'fall over' on the Cluster
