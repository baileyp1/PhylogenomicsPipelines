#!/bin/bash
#SBATCH -J make_species_trees                               # -c option: need to set cpu before sbatch calls so doing that with the sbatch call instead.
#SBATCH -n 1                                                # Also setting these paramters in the sbatch call:  -p long, --mem 100000 
####SBATCH -o ${fileNamePrefix}_make_species_trees.log			# If this line is not set, default output name is slurm-%A_%a.out: %A = $SLURM_ARRAY_JOB_ID; %a = $SLURM_ARRAY_TASK_ID		 
####SBATCH -e ${fileNamePrefix}_make_species_trees.log			# NB - if I just specify the %a, then I don't get an accumulation of ouput files for each run of the script
												            ### Test whether you can specify redirection to .log file via 2>&
                                                            ### NBNB - moved the -o and -flags to the main script because you can't specify fileNamePrefix variable at this point in the script!

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
treeTipInfoMapFile=${10}        # NB - just testing if this filename was submitted, then will add tree tip info to species tree
dnaSelected=${11}
proteinSelected=${12}
codonSelected=${13}
collapseNodes=${14}
alnFileForTreeSuffix=${15}      # Should be generic for dna, protein or codon filenames!  aln.for_tree.fasta

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
echo collapseNodes: $collapseNodes


numbrLowSupportNodesThreshold=70    # For use with getTreeStats() - if using a program like fasttree that
                                    # outputs support values as fractions, need to convert percent value
                                    # beforehand then import into function.

> ${fileNamePrefix}_summary_tree_stats.txt  # Wiping out file contents from any previous run 


getTreeStats () {
    ###########
    # Function: for tree stats and tree comparisons

    # Input parameters:
    # $1 = NewickFile
    # $2 = numbrLowSupportNodesThreshold    - current set to bootstrap value of 0.7 or 70 depending on program used

    # NB - pretty sure that the Newick file must only hae a single value at each node i.e. not the Astral -t 2 outputs
    ###########

    newickTree=$1

    # Count number of low support nodes:
    numbrLowSupportNodes=`nw_ed -r $newickTree  "i & b < $2" s | wc -l | sed 's/ //g' `
    # Count the total number of support nodes with support values: 
    totalNumbrSupportNodes=`nw_ed -r $newickTree  "i & b > 0" s | wc -l | sed 's/ //g' `
    # Write summary stats:
    echo >> ${fileNamePrefix}_summary_tree_stats.txt
    echo "$newickTree - number of internal nodes with low support values (total number of internal nodes): ${numbrLowSupportNodes} (${totalNumbrSupportNodes})" >> ${fileNamePrefix}_summary_tree_stats.txt 

    outFilePrefix=`basename $newickTree .nwk`

    # Also prepare a list of sorted clades for comparing with other trees: 
    nw_ed $newickTree  "i & b >0" s | nw_topology -I - | \
    while read clade; do 
    echo $clade | sed 's/[)(;]//g' | tr ',' '\n' | sort |tr '\n' ' '; echo
    done | sort > ${outFilePrefix}_sorted_clade_list.txt
    ### 21.7.2020 - also bring in the USE_THIS .nwk file and repeat - required for annotation file
}


#### 14.1.2019 - this step belwo can stay in this script
#### Need to use the trim0.01.fasta file for this step so I concatenate the correct files.

### 20.8.2020 - need to put conditonals for each residue type:
### if [[ $dnaSelected == 'yes' ]]; then
###     seqType=dna   
###     make a dna species tree - can I make a subroutine for dna protein and codon and for supermatrix trees? 


# Concatenate trees containing more than $fractnSpecies of samples for use with ASTRAL:
# NB - have made filenames generic so they can be derived from different phylo programs.
echo dnaSelected: $dnaSelected
if [[ $dnaSelected == 'yes' ]]; then    ### TEMP SET UP - this doesn't woerk!!!!
    seqType=dna     ### Still to finish implementation
    echo Listing required files:
    ls *.${seqType}.$alnFileForTreeSuffix       # This gives the whole set
    for file in *.${seqType}.$alnFileForTreeSuffix; do
 	  gene=`echo $file | sed "s/.${seqType}.$alnFileForTreeSuffix//" `
 	  numbrSamples=`cat $file | grep '>' | wc -l `;
        echo $gene " " $numbrSamples
    done \
    | awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
    '$2 >= numbrSamplesThreshold  {print $1 "_dna_gene_tree_USE_THIS.nwk"}' \
    | xargs cat > ${fileNamePrefix}_${seqType}_gene_trees_for_coelescence_phylo.nwk
fi

### 24.8.2020 - improving the above selection of tree files.
### make_species_trees.sh line ~84ff - have had an error with gene dna_gene_tree_USE_THIS.nwk files not existing
### - Yes it is because when I’m tryign to use the gene aln file but it has been removed from the analysis due to that 
### it only has 2 seqs so the .nwk file doesn’t exist but the aln does (and that’s OK); at least occurs when there is no
### filtering (?) and now because of filtering small seqs just before the aln. The code is OK but to remvoe the error shoudl do:
### 1.just use the .nwk files to begin with. To coutn the # seqs in fasta - could use Newick Utils on the trees instead - 
### I still think i need to filter here unless I remove the fratnSpecies variable
### 2.implement removal of the for_tree aln andnwk files before start of analysis, then any previous runs can’t affect 
### the outcome of next run if run in same directory - better solution i think 


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' || "$phyloProgramPROT" == 'iqtree2' ]]; then
    seqType=protein
    # Concatenate the protein gene trees as well (almost repeat of the above code):
    for file in *.${seqType}.$alnFileForTreeSuffix; do
 	  gene=`echo $file | sed "s/.${seqType}.$alnFileForTreeSuffix//" `
 	  numbrSamples=`cat $file | grep '>' | wc -l `;
        echo $gene " " $numbrSamples
    done \
    | awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
    '$2 >= numbrSamplesThreshold  {print $1 "_protein_gene_tree_USE_THIS.nwk"}' \
    | xargs cat > ${fileNamePrefix}_${seqType}_gene_trees_for_coelescence_phylo.nwk
fi


echo
echo
echo ##########################
echo Step 2 - running ASTRAL...
echo ##########################
echo
echo

# Before running Astral, collapse clades with low bootstrap values (less than $collapseNodes) from all trees at once, if requested.
# NB - this step collapses nodes whose certainty is not clear by creating multifurcations with a parent node that is well supported.
# No taxa are removed! 
if [[ $collapseNodes != 'no' ]]; then
    
### 15.10.2020 - need a conditional for the DNA programs as for protein

    if [[ $phyloProgramDNA == 'fasttree' ]]; then
        # For fasttree (only, so far), need to convert $collapseNodes percent value to a fraction and use that in nw_ed. 
        collapseNodes=`echo $collapseNodes | awk 'fractn=$1/100 {print fractn}' `
        echo "\$CollapseNodes should now be a fraction for FASTTREE (option -L): $collapseNodes"

        # Also convert $numbrLowSupportNodesThreshold to a fraction for FASTTREE for use in getTreeStats function:
        numbrLowSupportNodesThreshold=`echo $numbrLowSupportNodesThreshold | awk 'fractn=$1/100 {print fractn}' `
        echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE (internal value): $numbrLowSupportNodesThreshold"
    fi
    collapseNodesD=$collapseNodes    # Required later for remembering the actual value for DNA (It might get changed for prtoein!!!)
    numbrLowSupportNodesThresholdD=numbrLowSupportNodesThreshold    # Same for this var
   

    nw_ed  ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo.nwk "i & (b < $collapseNodesD)" o > ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo_bs_less_${collapseNodesD}_rmed.nwk
    
   
    if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' || "$phyloProgramPROT" == 'iqtree2' ]]; then

        ### Not tested logic for this conditional yet - I think I need a separate varialbe for DNA and protein - done
        if [[ $phyloProgramDNA == 'fasttree' ]]; then
            # For fasttree (only, so far), need to convert $collapseNodes percent value to a fraction and use that in nw_ed. 
            collapseNodes=`echo $collapseNodes | awk 'fractn=$1/100 {print fractn}' `
            echo "\$CollapseNodes should now be a fraction for FASTTREE (option -L): $collapseNodes"

            # Also convert $numbrLowSupportNodesThreshold to a fraction for FASTTREE for use in getTreeStats function:
            numbrLowSupportNodesThreshold=`echo $numbrLowSupportNodesThreshold | awk 'fractn=$1/100 {print fractn}' `
         echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE (internal value): $numbrLowSupportNodesThreshold"
        fi 

        # Remove clades with low bootstrap values from the protein trees:
        nw_ed  ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo.nwk "i & (b < $collapseNodes)" o > ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo_bs_less_${collapseNodes}_rmed.nwk
    fi
fi


### 8.5.2020 - LIKE FOR PROTEIN CREATE A CONDITIONAL HERE FOR DNA - NB - BUT SORT ASTRAL PROGAM PATH OUTSIDE CONDITIOANL
# Need to find the location of ASTRAL and supply the absolute path to the jar file.
# This assumes that the ASTRAL directory is already in $PATH. 
### NB - iqtree11.12.2019 - could always set up an alias! e.g. astral='java -jar $PATH/astral.5.6.3.jar' 
pathToAstral=`which astral.5.7.3.jar `
echo Running Astral on the DNA gene trees...
$exePrefix java -jar $pathToAstral -t 2 \
-i ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo_bs_less_${collapseNodesD}_rmed.nwk \
-o ${fileNamePrefix}.dna.species_tree.astral_-t2.nwk
### Old name: -o ${fileNamePrefix}_astral_-t2_dna_species_tree.nwk
# Astral commands:
# -t 2 - outputs all branch annotations
# -t 16 - output a table file called freqQuad.csv - useful for finding high pp scores for the alternative probabilities

# Extract the pp1 scores from the main Astral tree:
cat ${fileNamePrefix}.dna.species_tree.astral_-t2.nwk \
| sed "s/'\[[fq=.\;0-9']\{1,\}pp1=//g" \
| sed "s/;[QCENp=.\;0-9-]\{1,\}\]':/:/g" \
> ${fileNamePrefix}.dna.species_tree.astral_pp1_value.nwk
### NB - 9.10.2020 - might want to remove all *.nwk files at the start of this script so they can't be used by a previous run.

getTreeStats ${fileNamePrefix}.dna.species_tree.astral_pp1_value.nwk $numbrLowSupportNodesThresholdD

# Add tree tip info.
# NB - for this option, the $treeTipInfoMapFile must have been submitted BUT I'm re-formatting it --> 'tree_tip_info_mapfile.txt'!):
if [ -s $treeTipInfoMapFile ]; then 
    nw_rename -l  ${fileNamePrefix}.dna.species_tree.astral_pp1_value.nwk \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}.dna.species_tree.astral_USE_THIS.nwk
fi

# Also running with -t option set to 16:
$exePrefix java -jar $pathToAstral -t 16 \
-i ${fileNamePrefix}_dna_gene_trees_for_coelescence_phylo_bs_less_${collapseNodesD}_rmed.nwk \
-o ${fileNamePrefix}.dna.species_tree.astral_-t16.nwk
mv freqQuad.csv ${fileNamePrefix}.dna.species_tree.astral_-t16_freqQuad.txt



if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng'  || "$phyloProgramPROT" == 'iqtree2' ]]; then
    echo Running Astral on the protein gene trees...
    $exePrefix java -jar $pathToAstral -t 2 \
    -i ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo_bs_less_${collapseNodes}_rmed.nwk \
    -o ${fileNamePrefix}.protein.species_tree.astral_-t2.nwk

    # Extract the pp1 scores from the main Astral tree:
    cat ${fileNamePrefix}.protein.species_tree.astral_-t2.nwk \
    | sed "s/'\[[fq=.\;0-9']\{1,\}pp1=//g" \
    | sed "s/;[QCENp=.\;0-9-]\{1,\}\]':/:/g" \
    > ${fileNamePrefix}.protein.species_tree.astral_pp1_value.nwk

    getTreeStats ${fileNamePrefix}.protein.species_tree.astral_pp1_value.nwk $numbrLowSupportNodesThreshold

    # Add tree tip info:
    if [[ -s $treeTipInfoMapFile ]]; then
        nw_rename -l ${fileNamePrefix}.protein.species_tree.astral_pp1_value.nwk \
        tree_tip_info_mapfile.txt \
        > ${fileNamePrefix}.protein.species_tree.astral_USE_THIS.nwk
    fi

    # Also running with -t option set to 16:
    $exePrefix java -jar $pathToAstral -t 16 \
    -i ${fileNamePrefix}_protein_gene_trees_for_coelescence_phylo_bs_less_${collapseNodes}_rmed.nwk \
    -o ${fileNamePrefix}.protein.species_tree.astral_-t16.nwk
    mv freqQuad.csv ${fileNamePrefix}.protein.species_tree.astral_-t16_freqQuad.txt
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


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng' || "$phyloProgramPROT" == 'iqtree2' ]]; then
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

getTreeStats ${fileNamePrefix}_fasttree_dna_species_tree.nwk $numbrLowSupportNodesThresholdD

# Add tree tip info:
if [[ -s $treeTipInfoMapFile ]]; then 
    nw_rename -l  ${fileNamePrefix}_fasttree_dna_species_tree.nwk \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}_fasttree_dna_species_tree_USE_THIS.nwk
fi

 
if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng'  || "$phyloProgramPROT" == 'iqtree2' ]]; then
    echo Running fasttree on the protein supermatrix...
    $exePrefix fasttree \
    ${fileNamePrefix}__mafft_protein_alns__ovr${fractnAlnCovrg_pc}pc_acpg_ovr${fractnSpecies_pc}pc_spgt__concatenated.fasta \
    > ${fileNamePrefix}_fasttree_protein_species_tree.nwk

    getTreeStats ${fileNamePrefix}_fasttree_protein_species_tree.nwk $numbrLowSupportNodesThreshold

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

getTreeStats RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE $numbrLowSupportNodesThresholdD

# Add tree tip info:
if [[ -s $treeTipInfoMapFile ]]; then 
    nw_rename -l  RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE \
    tree_tip_info_mapfile.txt \
    > ${fileNamePrefix}_RAxML_100_bs_dna_species_tree_USE_THIS.nwk
fi


if [[ "$phyloProgramPROT" == 'fasttree' ||  "$phyloProgramPROT" == 'raxml-ng'  || "$phyloProgramPROT" == 'iqtree2' ]]; then
     echo Running RAxML on the protein supermatrix...
    # RAxML won't run if files already exists from a previous run so remove them here: 
    ### Could look for a --force option
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

    getTreeStats RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE $numbrLowSupportNodesThreshold

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



