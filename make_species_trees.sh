#!/bin/bash
#SBATCH -J make_species_trees               # -c option: need to set cpu before sbatch calls so doing that with the sbatch call instead.
#SBATCH -n 1                                # Also setting these paramters in the sbatch call:  -p long, --mem 100000 
#SBATCH -o make_species_trees.log			# If this line is not set, default output name is slurm-%A_%a.out: %A = $SLURM_ARRAY_JOB_ID; %a = $SLURM_ARRAY_TASK_ID		 
#SBATCH -e make_species_trees.log			# NB - if I just specify the %a, then I don't get an accumulation of ouput files for each run of the script
												### Test whether you can specify redirection to .log file via 2>&1

shopt -s failglob 


echo Inside slurm_make_gene_trees.sh script:
echo SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID


fractnAlnCovrg=$1
fractnSpecies=$2
totalNumbrSamples=$3
fileNamePrefix=$4
geneFile=$5                     ### 27.1.2020 - don't need this var once stats code has been removed
cpu=$6               ### 28.1.2020 - added this paramter in but not usign it (at the moment - trying to make filenames more generic)
phyloProgramDNA=$7
phyloProgramPROT=$8
exePrefix="$9"
treeTipInfoMapFile=${10}

# Convert $emptyMatchStateFractn and $fractnSpecies to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `
fractnSpecies_pc=`awk -v FRACTN=$fractnSpecies 'BEGIN{printf "%.0f", FRACTN * 100}' `
# Minimum number of samples to tolerate for including into Astral:
numbrSamplesThreshold=`awk -v FRACTN=$fractnSpecies -v numbrSamples=$totalNumbrSamples 'BEGIN{printf "%.0f", FRACTN * numbrSamples}' `


echo Preparing to run species trees...
echo fractnAlnCovrg to use: $fractnAlnCovrg
echo fractnSpecies to use: $fractnSpecies
echo numbrSamples: $totalNumbrSamples
echo numbrSamplesThreshold: $numbrSamplesThreshold
echo exePrefix: $exePrefix
echo treeTipInfoMapFile: $treeTipInfoMapFile


#### 14.1.2019 - this step belwo can stay in this script
#### Need to use the trim0.01.fasta file for this step so I concatenate the correct files.


# Concatenate trees containing more than $fractnSpecies of samples for use with ASTRAL:
# NB - have made filenames generic so they can be derived from different phylo programs.
for file in *_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta; do
 	gene=`echo $file | sed "s/_mafft_dna_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta//" `
 	numbrSamples=`cat $file | grep '>' | wc -l `;
    echo $gene " " $numbrSamples
done \
| awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
'$2 >= numbrSamplesThreshold  {print $1 "_dna_gene_tree_USE_THIS.nwk"}' \
| xargs cat > ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo.nwk


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
    # Concatenate the protein gene trees as well (almost repeat of the above code):
    for file in *_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta; do
 	  gene=`echo $file | sed "s/_mafft_protein_aln_ovr${fractnAlnCovrg_pc}pc_aln_covrg_trimCols0.003.fasta//" `
 	  numbrSamples=`cat $file | grep '>' | wc -l `;
        echo $gene " " $numbrSamples
    done \
    | awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
    '$2 >= numbrSamplesThreshold  {print $1 "_protein_gene_tree_USE_THIS.nwk"}' \
    | xargs cat > ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo.nwk
fi


echo
echo
echo ##########################
echo Step 2 - running ASTRAL...
echo ##########################
echo
echo

# Before running Astral, remove clades with low bootstrap values (less than 15 %) from all trees at once:
nw_ed  ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo.nwk 'i & (b <= 0.15)' o > ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo_bs_less15pc_rmed.nwk
### I think I will have to assess the number of samples lost at this stage as well - yes - in which case this needs to be done in the stats script!!!! or could append results here
### Could be interesting - they should indicate difficult seqs - NBNB - does it just remove internal nodes or some tree tips as well?!

if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
    # Remove clades with low bootstrap values (15%) from the protein trees:
    nw_ed  ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo.nwk 'i & (b <= 0.15)' o > ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo_bs_less15pc_rmed.nwk
fi


### 8.5.2020 - LIKE FOR PROTEIN CREATE A CONDITIONAL HERE FOR DNA - NB - BUT SORT ASTRAL PROGAM PATH OUTSIDE CONDITIOANL
# Need to find the location of ASTRAL and supply the absolute path to the jar file.
# This assumes that the ASTRAL directory is already in $PATH. 
### NB - 11.12.2019 - could always set up an alias! e.g. astral='java -jar $PATH/astral.5.6.3.jar' 
pathToAstral=`which astral.5.7.3.jar `
echo Running Astral on the DNA gene trees...
$exePrefix java -jar $pathToAstral \
-i ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo_bs_less15pc_rmed.nwk \
-o ${fileNamePrefix}_astral_dna_species_tree.nwk

# Add tree tip info.
# NB - for this option, the $treeTipInfoMapFile must have been submitted BUT I'm re-formatting it --> 'tree_tip_info_mapfile.txt'!):
if [ -s $treeTipInfoMapFile ]; then 
    nw_rename -l  ${fileNamePrefix}_astral_dna_species_tree.nwk \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}_astral_dna_species_tree_USE_THIS.nwk
fi


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
    echo Running Astral on the protein gene trees...
    $exePrefix java -jar $pathToAstral \
    -i ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo_bs_less15pc_rmed.nwk \
    -o ${fileNamePrefix}_astral_protein_species_tree.nwk

    # Add tree tip info:
    if [[ -s $treeTipInfoMapFile ]]; then
        nw_rename -l ${fileNamePrefix}_astral_protein_species_tree.nwk \
        tree_tip_info_mapfile.txt \
        > ${fileNamePrefix}_astral_protein_species_tree_USE_THIS.nwk
    fi
fi


# Step2a - tree rooting
# Investigated outgroups - used NewickTools (there’s also GoTree and another tools I made a note of)
# Re-rooting with the whole outgroup (just need to specify 2 leaf labels for a clade):
# /Users/baileyp/Documents/ProgramFiles/newick-utils-1.6/src/nw_reroot \
# Apiaceae_astral.nwk   839_Aralia_cordata  5427_Pennantia_cunninghamii > Apiaceae_astral_outgroup_rerooted_by_8_species.nwk
# Re-rooting by one species
# /Users/baileyp/Documents/ProgramFiles/newick-utils-1.6/src/nw_reroot Apiaceae_astral.nwk   839_Aralia_cordata  > Apiaceae_astral_outgroup_rerooted_by_839_Aralia_cordata.nwk
# NB - can also root the tree in iToL on a specific branch.
# raxml-ng has a rooting option.


echo
echo
echo ###########################################
echo 'Step 3 - concatenating gene alignments...'
echo ###########################################
echo
echo
# Concatenate the alignment for 2 reasons:
# 1.easy browsing in Jalview - yes this works OK.
# 2.for running RAxML with a supermatrix tree.
AMAS.py concat  -c 1 \
-i `cat mafft_dna_alns_fasta_file_list.txt` \
--in-format  fasta \
-d dna \
--out-format fasta \
-t ${fileNamePrefix}__mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta
# NB- also creates a partitions file with coords for each gene which could be used to set different models.


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
    # Also concatenate protein alns:
    AMAS.py concat  -c 1 \
    -i `cat mafft_protein_alns_fasta_file_list.txt` \
    --in-format  fasta \
    -d aa \
    --out-format fasta \
    -t ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated_temp.fasta
    ### On cluster the protein fasta headers now contain _\[translate(1)\] - fastatranslate (v2.4.x) adds this string to the header,
    ### then AMAS adds it onto the required fasta record name - so remove them here:
    cat ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated_temp.fasta \
    | sed 's/_\[translate(1)\]//' \
    > ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta
    rm ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated_temp.fasta
fi


echo
echo
echo ####################################################
echo 'Step 4 - Running fasttree on the DNA supermatrix...'
echo ####################################################
echo
echo


$exePrefix fasttree -nt \
-gtr \
${fileNamePrefix}__mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
> ${fileNamePrefix}_fasttree_dna_species_tree.nwk

# Add tree tip info:
if [[ -s $treeTipInfoMapFile ]]; then 
    nw_rename -l  ${fileNamePrefix}_fasttree_dna_species_tree.nwk \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}_fasttree_dna_species_tree_USE_THIS.nwk
fi

 
if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
    echo Running fasttree on the protein supermatrix...
    $exePrefix fasttree \
    ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
    > ${fileNamePrefix}_fasttree_protein_species_tree.nwk

    # Add tree tip info:
    if [[ -s $treeTipInfoMapFile ]]; then 
        nw_rename -l  ${fileNamePrefix}_fasttree_protein_species_tree.nwk \
        tree_tip_info_mapfile.txt \
        > ${fileNamePrefix}_fasttree_protein_species_tree_USE_THIS.nwk
    fi
fi


echo
echo
echo #######################################################################
echo 'Step 5 - running RAxML on the concatenated alignment (supermatrix)...'
echo #######################################################################
echo
echo
# RAxML won't run if files already exists from a previous run so remove them here: 
if [ -a RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
if [ -a RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
if [ -a RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
if [ -a RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
if [ -a RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
if [ -a ${fileNamePrefix}__mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced ]; then rm ${fileNamePrefix}__mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced; fi


### For Slurm, add: srun -J raxml
# Runtime: took 1h8’, 283MB mem for Achariaceae aln
echo Running RAxML on the DNA supermatrix...
$exePrefix raxmlHPC-PTHREADS-SSE3 -T $cpu \
-f a \
-x 12345 \
-p 12345 \
-# 100 \
-m GTRGAMMA \
-s ${fileNamePrefix}__mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
-n ${fileNamePrefix}__raxmlHPC-PTHREADS-SSE
# This is the Newick  output file compatible with NewickTools):
# RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE

# Add tree tip info:
if [[ -s $treeTipInfoMapFile ]]; then 
    nw_rename -l  RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}_RAxML_100_bs_dna_species_tree_USE_THIS.nwk
fi


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' ]]; then
     echo Running RAxML on the protein supermatrix...
    # RAxML won't run if files already exists from a previous run so remove them here: 
    if [ -a RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bootstrap.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
    if [ -a RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bestTree.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
    if [ -a RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
    if [ -a RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitionsBranchLabels.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
    if [ -a RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE ]; then rm RAxML_info.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE; fi
    if [ -a ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced ]; then rm ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta.reduced; fi

    $exePrefix raxmlHPC-PTHREADS-SSE3 -T $cpu \
    -f a \
    -x 12345 \
    -p 12345 \
    -# 100 \
    -m PROTCATJTT \
    -s ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
    -n ${fileNamePrefix}__raxmlHPC-PTHREADS-SSE
    # This is the Newick  output file compatible with NewickTools):
    # RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE

    # Add tree tip info:
    if [[ -s $treeTipInfoMapFile ]]; then 
        nw_rename -l   RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE \
        tree_tip_info_mapfile.txt \
        > ${fileNamePrefix}_RAxML_100_bs_protein_species_tree_USE_THIS.nwk
    fi
fi










# To mid-point root the tree - e.g.:    I THINK THE TREE NEEDS TO REMAIN UNROOTED!!!!
# /usr/users/ga002/baileyp/ProgramFiles/gotree/gotree_amd64_linux reroot midpoint \
# -i ${fileNamePrefix}__RAxML-PTHREADS-SSE_100_bs_USE_THIS.nwk \
# > ${fileNamePrefix}__RAxML-PTHREADS-SSE_100_bs_USE_THIS_midpoint_rooted.nwk


### Also assessing the gene coverage per species across the concatenated aln.
### This is because so far, seqs have been removed on a per gene aln basis - but 
### there is still no idea whether any samples are still failing across many genes.
### This coverage may also be low for for species with low gene recovery.
### These species need to be listed - if their positions in the tree are unclear
### i.e. low bootstrap support, then this suggests:
### 1. samples might need to be regenerated and genes re-recovered
### 2. supermatrix tree could be generated without these samples
### Already have similar stats - but still these seqs should be removed I guess


### NB - this code is copy of make_gene_trees.sh script - needs converting to Python function (explore a Python one liner -c option:
###cat ${fileNamePrefix}_mafft_dna_alns__ovr70pc_acpg_ovr80pc_spgt__concatenated.out | seqtk seq -l 0 /dev/fd/0 | tail -n1 | awk '{print length($1)}'


# Calculate the full aln length:
###alnLength=`seqtk seq -l 0 ${fileNamePrefix}_mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr80pc_spgt__concatenated.fasta | tail -n1 | awk '{print length($1)}' `

### Calculte length of longets gene:
###lenLongestGene=`fastalength  ${fileNamePrefix}_mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr80pc_spgt__concatenated.fasta | sort -n | tail -n 1 | awk '{print $1}' `
###echo lenLongestGene: $lenLongestGene

# Now find the lengths of each sequence in bases only (not dashes!)
# and create a list of sequences that are > $lenLongestGene (used to use $fractnAlnCovrg):
# NB - this code works because, luckily for this case, fastalength ignores dash characters (unlike seqtk used above)!		
###fastalength ${fileNamePrefix}_mafft_dna_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr80pc_spgt__concatenated.fasta  \
###| awk -v LENG=$lenLongestGene -v FRACTN=$fractnAlnCovrg '{if( ($1/LENG) >= FRACTN ) {print $2} }' \
###> ${fileNamePrefix}_concat_dna_aln__ovr${fractnAlnCovrg_pc}pc_acpg_ovr80pc_spgt__concatenated.fasta


