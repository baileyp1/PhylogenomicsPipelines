#!/bin/bash
#SBATCH -J make_species_trees                               # -c option: need to set cpu before sbatch calls so doing that with the sbatch call instead.
#SBATCH -n 1                                                # Also setting these paramters in the sbatch call:  -p long, --mem 100000 
####SBATCH -o ${fileNamePrefix}_make_species_trees.log			# If this line is not set, default output name is slurm-%A_%a.out: %A = $SLURM_ARRAY_JOB_ID; %a = $SLURM_ARRAY_TASK_ID		 
####SBATCH -e ${fileNamePrefix}_make_species_trees.log			# NB - if I just specify the %a, then I don't get an accumulation of ouput files for each run of the script
												            ### Test whether you can specify redirection to .log file via 2>&
                                                            ### NBNB - moved the -o and -flags to the main script because you can't specify fileNamePrefix variable at this point in the script!
#######################
# make_species_trees.sh

# Author: Paul Bailey

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
#######################
shopt -s failglob 


echo Inside slurm_make_species_trees.sh script:
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
dnaSelected=${11}               # NB - for these <seqType>Selected variables, this script should be able to process all three types together OR separately!
proteinSelected=${12}
codonSelected=${13}
collapseNodes=${14}
alnFileForTreeSuffix=${15}      # Should be generic for dna, protein or codon filenames!  aln.for_tree.fasta
option_u=${16}
speciesTreeProgram=${17}
outgroupRoot=${18}
speciesTreeSlurmMem=${19}
numbrBootstraps=${20}


echo Preparing to run species trees...
echo fractnAlnCovrg to use: $fractnAlnCovrg
echo fractnSpecies to use: $fractnSpecies
echo numbrSamples: $totalNumbrSamples
echo numbrSamplesThreshold: $numbrSamplesThreshold
echo exePrefix: $exePrefix
echo treeTipInfoMapFile: $treeTipInfoMapFile
echo dnaSelected: $dnaSelected
echo proteinSelected $proteinSelected
echo codonSelected: $codonSelected
echo collapseNodes: $collapseNodes
echo speciesTreeProgram: $speciesTreeProgram
echo "Option -U (species tree memory to use in mb): $speciesTreeSlurmMem"
echo "Number of bootstrap data sets to use for RAxML: $numbrBootstraps"

# Convert $emptyMatchStateFractn and $fractnSpecies to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `
fractnSpecies_pc=`awk -v FRACTN=$fractnSpecies 'BEGIN{printf "%.0f", FRACTN * 100}' `
# Minimum number of samples to tolerate for including into Astral (9.2.2022 - looks like the number is rounded down!)
numbrSamplesThreshold=`awk -v FRACTN=$fractnSpecies -v numbrSamples=$totalNumbrSamples 'BEGIN{printf "%.0f", FRACTN * numbrSamples}' `
# Above awk code is zero proof - can have 0 * 100 - returns zero


# Input checks for the sequence type option:
astralSelected=no   # Need to match whole word with grep -Ew to avoid getting astral when astralmp input!
astralmpSelected=no
fasttreeSelected=no
raxmlSelected=no                    
if [[ `echo $speciesTreeProgram | grep -Ew -o 'astral' ` == 'astral' ]];then
     astralSelected=yes
fi
if [[ `echo $speciesTreeProgram | grep -o 'astralmp' ` == 'astralmp' ]];then
     astralmpSelected=yes
fi
if [[ `echo $speciesTreeProgram | grep -o 'fasttree' ` == 'fasttree' ]];then
     fasttreeSelected=yes
fi
if [[ `echo $speciesTreeProgram | grep -o 'raxml' ` == 'raxml' ]];then
     raxmlSelected=yes
fi
if [[ $astralSelected == 'no' &&  $astralmpSelected == 'no' && $fasttreeSelected == 'no' && $raxmlSelected == 'no' ]]; then
    echo "ERROR: No phylogeny program (option -s) was entered or recognised."; exit
fi
### Might also want to check option -s in the wrapper script so error is reported eariler.


if [[ -s $treeTipInfoMapFile && $option_u == 'yes' ]]; then
    cp -p tree_tip_info_mapfile_plusL.txt  tree_tip_info_mapfile.txt
fi


numbrLowSupportNodesThreshold=95    # For use with getTreeStats() - if using a program like fasttree that
                                    # outputs support values as fractions, need to convert percent value
                                    # beforehand then import into function - 23.10.2020 - changed logic slightly - too complicated to use this value

# 
echo "Tree statistics
---------------
Support nodes threshold set to: $numbrLowSupportNodesThreshold
ASTRAL: local posterior probability (considered to be well supported if >= 95% (Sayyari et al., 2016) Molecular Biology and Evolution https://doi.org/10.1093/molbev/msw079
IQTREE2: ultrafast bootstrap support (considered to be well supported if >= 95% (Minh et al., 2013) Molecular Biology and Evolution https://doi.org/10.1093/molbev/mst024
---------------"> ${fileNamePrefix}_summary_tree_stats.txt  # Wiping out file contents from any previous run 



getTreeStats () {
    ###########
    # Function: for tree stats and tree comparisons

    # Input parameters:
    # $1 = NewickFile
    # $2 = numbrLowSupportNodesThreshold    - value is global so doesn't need importing - clearer though
    # $3 = tree type e.g. astral, fasttree

    # The Newick file contains the identifiers from the gene trees, not info/taxonomy from the option -t in file 

    # NB - pretty sure that the Newick file must only hae a single value at each node i.e. not the Astral -t 2 outputs
    ###########

    newickTree=$1
    bootstrapThreshold=$2
    treeType=$3
    echo bootstrapThreshold: $bootstrapThreshold


### NB - 23.10.2020 - the logic of $bootstrapThreshold (and possiby collapseNodes) is wrong and doesn't work for anythign other than fasttree gene trees i think 
### Better to just ask here whether tree is astral or fastatree then preapre bootstrapThreshold as a fraction
    if [[ $treeType == 'astral'* || $treeType == 'fasttree' ]]; then 
       # Convert $numbrLowSupportNodesThreshold to a fraction:
        bootstrapThreshold=`echo $bootstrapThreshold | awk 'fractn=$1/100 {print fractn}' `
        echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE and ASTRAL (internal value): $bootstrapThreshold"  
    fi

    # Count number of low support nodes:
### not sure why I'm using -r!!
    numbrLowSupportNodes=`nw_ed -r $newickTree  "i & b >= $bootstrapThreshold" s | wc -l | sed 's/ //g' `
### NBNB - 15.10.2020 - this gives '1' even when the file is empty!!!#
    # Count the total number of support nodes with support values:
    ### 23.10.2020 - Astral tree seems to have two root nodes which get printed by "i" (i.e. get 2 more nodes than expected, only can notice with v small trees)
    ### and can't be ignored by "!r".
    ### Could remove them by: "i & (b != 0)" s | wc -l
    ### But it might also remove internal nodes with a zero bootstrap - small point. 
    ### Just changing to use "i" only
    #totalNumbrSupportNodes=`nw_ed -r $newickTree  "i & b > 0" s | wc -l | sed 's/ //g' `
    totalNumbrSupportNodes=`nw_ed -r $newickTree  "i" s | wc -l | sed 's/ //g' `
    # Write summary stats:
    echo >> ${fileNamePrefix}_summary_tree_stats.txt
    if [[ $totalNumbrSupportNodes > 0 ]]; then
        echo "$newickTree - number of internal nodes with high support values >= $bootstrapThreshold (total number of internal nodes): ${numbrLowSupportNodes} (${totalNumbrSupportNodes})" >> ${fileNamePrefix}_summary_tree_stats.txt 
    fi

    outFilePrefix=`basename $newickTree .nwk`

    # Also prepare a list of sorted clades for comparing with other trees:
    nw_ed $newickTree  "i & b >0" s | nw_topology -I - | \
    while read clade; do 
    echo $clade | sed 's/[)(;]//g' | tr ',' '\n' | sort |tr '\n' ' '; echo
    done | sort > ${outFilePrefix}_sorted_clade_list.txt
    ### 21.7.2020 - also bring in the USE_THIS .nwk file and repeat - required for annotation file
}



makeSpeciesTree () {
    ###########
    # Function: makes a SINGLE species tree from the available methods in this function 
    #           from a sequence alignment
    #
    # Input parameters:
    # $1 = residue type: dna, aa or codon
    # $2 = input file (list of Newick trees or a concatenated alignment in fasta format)
    # $3 = output file prefix including the path to any directory
    # $4 = program to use: astral, astralmp, fasttree, raxml
    # $5 = raxml_model e.g. raxmlModel='GTR+G' - for DNA; for protein raxmlModel='JTT+G'
    ### UPTOHERE
    # $6 = iqtree2_seq_type  iqTree2SeqType='DNA'; for protein iqTree2SeqType='AA'
    # $7 = fasttreeFlags='-nt -gtr' - for DNA; for protein fasttreeFlags='' - NB - this last flag needs to be last in case it is blank - it needs to be blank for protein analysis) 
    
    # NB - only importing variables into this function if they have different input values.
    #      All other variables are globally available and do not change during running of this script.
    ###########

    residueType=$1
    infile=$2
    outFilePrefix=$3
    programToUse=$4
    raxmlModel=$5
    iqTree2SeqType=$6
    fasttreeFlags=$7

    if [[ "$programToUse" == 'astral' ]]; then
        echo programToUse: $programToUse
        echo AstralInFile: $infile
        echo Running Astral on the DNA gene trees...
        # NB - location of ASTRAL must be presented to the java -jar command as a environment variable
        # containing the absolute path to the jar file.
        # NB - 22.10.2020 - added the -Xmx12000m (12G memory) to try and increase speed of ASTRAL. OK for 353 gene trees and 3500 samples, may need to increase for larger data sets
        $exePrefix java -Xmx12000m -jar $ASTRAL -t 2 \
        -i $infile \
        -o ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astral_-t2.nwk
        # Astral commands:
        # -t 2 - outputs all branch annotations
        # -t 16 - output a table file called freqQuad.csv - useful for finding high pp scores for the alternative probabilities

        # Extract the pp1 scores from the main Astral tree:
        cat ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astral_-t2.nwk \
        | sed "s/'\[[fqEaN=.\;0-9'-]\{1,\}pp1=//g" \
        | sed "s/;[QCENp=.\;0-9-]\{1,\}\]':/:/g" \
        > ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astral_pp1_value.nwk
        # Also running with -t option set to 16:
        #$exePrefix java -Xmx12000m -jar $ASTRAL -t 16 \
        #-i $inFile --output ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astral_-t16.nwk
        ### 21.7.2021 - NB - in this script the -o or --output flag doesn't work for ASTRAL or ASTRAL-MP!
        ### It works from the command line though so no idea what's going on - error is:
        ### Error: -o does not exist.
        ### Error: Unexpected argument: ./test.dna.species_tree.astralmp_-t16.nwk
        #mv ${outFilePrefix}/freqQuad.csv ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astral_-t16_freqQuad.txt
    elif [[ "$programToUse" == 'astralmp' ]]; then
        echo programToUse: $programToUse
        echo AstralInFile: $infile
        echo Running ASTRAL-MP on the DNA gene trees...
        ### Could also add a fixed value to Slurm memory OPTION -U so astral-mp will still work if slurm mem is set to zero:
        ### 4.8.2021 - if $speciesTreeSlurmMem == 0 , then give it some memory
        ### let "speciesTreeSlurmMem= $speciesTreeSlurmMem + 12000" - might not be optimal though, but then user should specify the amount of memory to use
        ###     Check this was enough mem for the release 1.0 tree ~ 3100 tree tips
        ###     echo statement about memory using e.g. "Defualt memory for ASTRAL-MP is 12GB. Use option -U to adjust memory up or down."

        ### 4.8.2021 - Also if -$cpu == 1 increase to 2 otherwise prgoram crashes
       
        $exePrefix java -Xmx${speciesTreeSlurmMem}m $ASTRALMPLIB -jar $ASTRALMP -C -T $cpu -t 2 \
        -i $infile \
        -o ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astralmp_-t2.nwk

        cat ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astralmp_-t2.nwk \
        | sed "s/'\[[fqEaN=.\;0-9'-]\{1,\}pp1=//g" \
        | sed "s/;[QCENp=.\;0-9-]\{1,\}\]':/:/g" \
        > ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astralmp_pp1_value.nwk
        # Also running with -t option set to 16:
        #$exePrefix java -Xmx${speciesTreeSlurmMem}m $ASTRALMPLIB -jar $ASTRALMP -C -T $cpu -t 16 \
        #-i $inFile \
        #-o ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astralmp_-t16.nwk
        #mv ${outFilePrefix}/freqQuad.csv ${outFilePrefix}/${fileNamePrefix}.${residueType}.species_tree.astralmp_-t16_freqQuad.txt
    elif [[ "$programToUse" == 'fasttree' ]]; then
        echo programToUse: $programToUse
        echo Running fasttree on the concatenated alignment...
        $exePrefix fasttree $fasttreeFlags \
        $infile \
        > $outFilePrefix/${fileNamePrefix}.${residueType}.species_tree.fasttree.nwk
    elif [[ "$programToUse" == 'raxmlq' ]]; then
        echo programToUse: $programToUse
        echo "Running RAxML on the concatenated alignment using a partitioned analysis of each gene.
        Number of samples: ${totalNumbrSamples}; model of evolution: $RAxML_ModelOfEvolution"
        # RAxML won't run if files already exists from a previous run so remove them here: 
        if [ -a RAxML_bootstrap.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bootstrap.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bestTree.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bestTree.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bipartitions.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitions.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bipartitionsBranchLabels.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitionsBranchLabels.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_info.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_info.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a ${infile}.reduced ]; then rm ${infile}.reduced; fi

        $exePrefix raxmlHPC-PTHREADS-SSE3 -T $cpu \
        -f a \
        -x 12345 \
        -p 12345 \
        -# $numbrBootstraps \
        -m $raxmlModel \
        -s $infile \
        -n ${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE \
        -q ${fileNamePrefix}.${residueType}.alns.partitions.txt
        # This is the Newick  output file compatible with NewickTools):
        # RAxML_bipartitions.${fileNamePrefix}__raxmlHPC-PTHREADS-SSE
        # -m - was using GTRGAMMA but GTRCAT is ~4x quicker - should not use it if have < 50 seqs in dataset.
        #      Also, shouldn't need to use alpha (GAMMA) parameter for really small data sets but can't use
        #      GTR model without GAMMA with RAxML; can use CAT with -V option but then shouldn't use CAT model
        #      for small data sets! RAxML gives an error: 'WARNING the alpha parameter...' - OK for now.
        # -q Added -q option for partitioning genes:
        #    If, e.g.,-m GTRGAMMA is used, individual alpha-shape parameters, GTR rates, and empirical 
        #    base frequencies will be estimated and optimized for each partition. Note - RAxML doesn't have 
        #    other DNA models so can't specify them in the partitions file
        #    Also note - you can not assign different models of rate heterogeneity to different partitions,
        #    i.e., it will be either CAT , GAMMA , GAMMAI etc. for all partitions, as specified with -m.
    elif [[ "$programToUse" == 'raxml' ]]; then
        echo programToUse: $programToUse
        echo "Running RAxML on the concatenated alignment WITHOUT using a partitioned analysis of each gene.
        Number of samples: ${totalNumbrSamples}; model of evolution: $RAxML_ModelOfEvolution"
        # RAxML won't run if files already exists from a previous run so remove them here: 
        if [ -a RAxML_bootstrap.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bootstrap.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bestTree.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bestTree.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bipartitions.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitions.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_bipartitionsBranchLabels.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_bipartitionsBranchLabels.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a RAxML_info.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE ]; then rm RAxML_info.${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE; fi
        if [ -a ${infile}.reduced ]; then rm ${infile}.reduced; fi

        $exePrefix raxmlHPC-PTHREADS-SSE3 -T $cpu \
        -f a \
        -x 12345 \
        -p 12345 \
        -# $numbrBootstraps \
        -m $raxmlModel \
        -s $infile \
        -n ${fileNamePrefix}.${residueType}.raxmlHPC-PTHREADS-SSE
    fi
}





############
# Main code:
############
if [[ "$astralSelected" == 'yes' || "$astralmpSelected" == 'yes' || "$astralproSelected" == 'yes' ]]; then
    echo "Preparing to run ASTRAL..."
    if [[ $dnaSelected == 'yes' ]]; then
        seqType=dna
        echo dnaSelected: $dnaSelected   

        # Prepare to concatenate trees containing more than $fractnSpecies of samples (AND with species > 3 ) for use with ASTRAL:
        # NB - have made filenames generic so they can be derived from different phylo programs.
        echo Listing required gene set:
        ls *.${seqType}.$alnFileForTreeSuffix       # This gives the whole set
        for file in *.${seqType}.$alnFileForTreeSuffix; do
            gene=`echo $file | sed "s/.${seqType}.$alnFileForTreeSuffix//" `
            numbrSamples=`cat $file | grep '>' | wc -l `;
            echo $gene " " $numbrSamples
        done \
        | awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
        '$2 >= numbrSamplesThreshold && $2 > 3 {print $1 "_dna_gene_tree_USE_THIS.nwk"}' \
        > ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt
        # Now concatenate the gene trees into a single file:
        cat ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt | xargs cat > ${fileNamePrefix}.${seqType}.gene_trees_set.nwk
    fi
    ### NBNB - 2nd -v fractnAlnCovrg_pc not used here anymore! Can delete... - same for protein below
    ###        Also, the above can be made into a function as it can be the same for protein

    ### 24.8.2020 - improving the above selection of tree files.
    ### make_species_trees.sh line ~84ff - have had an error with gene dna_gene_tree_USE_THIS.nwk files not existing
    ### - Yes it is because when I’m tryign to use the gene aln file but it has been removed from the analysis due to that 
    ### it only has 2 seqs so the .nwk file doesn’t exist but the aln does (and that’s OK); at least occurs when there is no
    ### filtering (?) and now because of filtering small seqs just before the aln. The code is OK but to remvoe the error shoudl do:
    ### 1.just use the .nwk files to begin with. To coutn the # seqs in fasta - could use Newick Utils on the trees instead - 
    ### I still think i need to filter here unless I remove the fratnSpecies variable
    ### 2.implement removal of the for_tree aln and nwk files before start of analysis, then any previous runs can’t affect 
    ### the outcome of next run if run in same directory - better solution i think - 12.7.2021 - I think this is what I'm doign now.
    ###     except still need to remove genes with < 3 samples (also now doing above) --> summarize thes points and leave a short note.

    if [[ $proteinSelected == 'yes' ]]; then
        seqType=protein
        # Prepare to concatenate the protein gene trees as well (almost repeat of the above code):
        for file in *.${seqType}.$alnFileForTreeSuffix; do
 	      gene=`echo $file | sed "s/.${seqType}.$alnFileForTreeSuffix//" `
 	      numbrSamples=`cat $file | grep '>' | wc -l `;
            echo $gene " " $numbrSamples
        done \
        | awk -v numbrSamplesThreshold=$numbrSamplesThreshold -v fractnAlnCovrg_pc=${fractnAlnCovrg_pc} \
        '$2 >= numbrSamplesThreshold && $2 > 3 {print $1 "_protein_gene_tree_USE_THIS.nwk"}' \
        > ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt
        # Now concatenate the gene trees into a single file:
        cat ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt | xargs cat > ${fileNamePrefix}.${seqType}.gene_trees_set.nwk
    fi
    ### 7.5.2022 - codon clause as well here? Shoudl make into a function

    # Before running Astral, collapse clades with low bootstrap values (less than $collapseNodes) from all trees at once, if requested.
    # NB - this step collapses nodes whose certainty is not clear by creating multifurcations with a parent node that is well supported.
    # No taxa are removed!
    dnaAstralInFile=${fileNamePrefix}.dna.gene_trees_set.nwk
    proteinAstralInfile=${fileNamePrefix}.protein.gene_trees_set.nwk
    dnaGeneTreesSetFilenames=${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt
    proteinGeneTreesSetFilenames=${fileNamePrefix}.${seqType}.gene_trees_set_filenames.txt
    if [[ $collapseNodes != 'no' ]]; then
        if [[ $dnaSelected == 'yes' ]]; then
            seqType=dna
            if [[ $phyloProgramDNA == 'fasttree' ]]; then
                # For fasttree (only, so far), need to convert $collapseNodes percent value to a fraction and use that in nw_ed. 
                collapseNodes=`echo $collapseNodes | awk 'fractn=$1/100 {print fractn}' `
                echo "\$CollapseNodes should now be a fraction for FASTTREE (option -L): $collapseNodes"

                # # Also convert $numbrLowSupportNodesThreshold to a fraction for FASTTREE for use in getTreeStats function:
                # numbrLowSupportNodesThreshold=`echo $numbrLowSupportNodesThreshold | awk 'fractn=$1/100 {print fractn}' `
                # echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE (internal value): $numbrLowSupportNodesThreshold"
            fi
            collapseNodesD=$collapseNodes    # NB - $collapseNodes is also required for the protein trees - will not need converting to a fraction if RAxML is used for protein aln so need temporary variable here.
            #numbrLowSupportNodesThresholdD=$numbrLowSupportNodesThreshold    # Same for this var
            nw_ed  $dnaAstralInFile "i & (b < $collapseNodesD)" o > ${fileNamePrefix}.${seqType}.gene_trees_set.collapse${collapseNodesD}.nwk
            # Also create list of filenames for the collapsed trees for use in creating an archive file:
            cat $dnaGeneTreesSetFilenames \
            | xargs ls | sed "s/.nwk$/.collapse${collapseNodesD}.nwk/" > ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.collapse${collapseNodesD}.txt
            # Also require separate gene tree files for the Kew Tree of Life Explorer for generating the gene tree tarball (see below):
            cat $dnaGeneTreesSetFilenames | \
            while read file; do
                filePrefix=`basename -s .nwk $file `
                cat $file | nw_ed  /dev/fd/0  "i & (b < $collapseNodesD)" o > ${filePrefix}.collapse${collapseNodesD}.nwk 
            done
            # The set of 'collapsed' files now required for ASTRAL:
            dnaAstralInFile=${fileNamePrefix}.${seqType}.gene_trees_set.collapse${collapseNodesD}.nwk
            dnaGeneTreesSetFilenames=${fileNamePrefix}.${seqType}.gene_trees_set_filenames.collapse${collapseNodesD}.txt
        fi

        if [[ $proteinSelected == 'yes' ]]; then
            seqType=protein
            ### Not tested logic for this conditional yet - I think I need a separate varialbe for DNA and protein - done
            if [[ $phyloProgramPROT == 'fasttree' ]]; then
                # For fasttree (only, so far), need to convert $collapseNodes percent value to a fraction and use that in nw_ed. 
                collapseNodes=`echo $collapseNodes | awk 'fractn=$1/100 {print fractn}' `
                echo "\$CollapseNodes should now be a fraction for FASTTREE (option -L): $collapseNodes"

            #    # Also convert $numbrLowSupportNodesThreshold to a fraction for FASTTREE for use in getTreeStats function:
            #    numbrLowSupportNodesThreshold=`echo $numbrLowSupportNodesThreshold | awk 'fractn=$1/100 {print fractn}' `
            # echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE (internal value): $numbrLowSupportNodesThreshold"
            fi 

            # Remove clades with low bootstrap values from the protein trees:
            nw_ed  $proteinAstralInfile "i & (b < $collapseNodes)" o > ${fileNamePrefix}.${seqType}.gene_trees_set.collapse${collapseNodes}.nwk
            cat $proteinGeneTreesSetFilenames \
            | xargs ls | sed "s/.nwk$/.collapse${collapseNodes}.nwk/" > ${fileNamePrefix}.${seqType}.gene_trees_set_filenames.collapse${collapseNodes}.txt
             cat $proteinGeneTreesSetFilenames | \
            while read file; do
                filePrefix=`basename -s .nwk $file `
                cat $file | nw_ed  /dev/fd/0  "i & (b < $collapseNodes)" o > ${filePrefix}.collapse${collapseNodes}.nwk 
            done
            proteinAstralInfile=${fileNamePrefix}.${seqType}.gene_trees_set.collapse${collapseNodes}.nwk
            proteinGeneTreesSetFilenames=${fileNamePrefix}.${seqType}.gene_trees_set_filenames.collapse${collapseNodes}.txt
        fi
        # else
        #     ### NB - 15.10.2020 - need to set these new variables for the rest of the code if collapse nodes is not set!!!
        #     ### The outfile name is not ideal if collspe node is not set - best to set a variable to deal with this
        #     collapseNodesD=$collapseNodes    # Var not used/needed !!
        #     # Need to mirror what's happening above; required later for remembering the actual value for DNA (It might get changed for protein!!!)
        #     if [[ $phyloProgramDNA == 'fasttree' ]]; then
        #         numbrLowSupportNodesThresholdD=`echo $numbrLowSupportNodesThreshold | awk 'fractn=$1/100 {print fractn}' `
        #         echo "\$numbrLowSupportNodesThreshold should now be a fraction for FASTTREE (internal value): $numbrLowSupportNodesThreshold"
        #     fi
        ### Does this need a codon clause here as well???? Yes and better still, convert to a function for use with DNA, protein and codon trees         
    fi
    echo
    echo
    echo #################
    echo Running ASTRAL...
    echo #################
    echo
    echo
    if [[ $dnaSelected == 'yes' ]]; then
        if [[ "$astralSelected" == 'yes' ]]; then
            # Function parameters: residue_type, input_file, outfile_prefix, program, program-specifIC paramters
### SHOULD REALLY TRY TO CONVERT PROGRAM PARAMTERS TO A GENERIC STRIGN FOR USE WITH ANY PRGORAM  
            makeSpeciesTree dna $dnaAstralInFile '.' astral 'GTR+G' 'DNA' '-nt -gtr'    # 30.4.2022 - Can these last params be removed? - misleading
            # Add tree tip info.
            # NB - for this option, the $treeTipInfoMapFile must have been submitted BUT I'm re-formatting it --> 'tree_tip_info_mapfile.txt'!):
            if [ -s $treeTipInfoMapFile ]; then
                nw_rename -l  ${fileNamePrefix}.dna.species_tree.astral_pp1_value.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.dna.species_tree.astral_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.dna.species_tree.astral_pp1_value.nwk $numbrLowSupportNodesThreshold astral
        elif [[ "$astralmpSelected" == 'yes' ]]; then 
            makeSpeciesTree dna $dnaAstralInFile '.' astralmp 'GTR+G' 'DNA' '-nt -gtr'
            # Add tree tip info.
            if [ -s $treeTipInfoMapFile ]; then
                nw_rename -l  ${fileNamePrefix}.dna.species_tree.astralmp_pp1_value.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.dna.species_tree.astralmp_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.dna.species_tree.astralmp_pp1_value.nwk $numbrLowSupportNodesThreshold astral
        fi
        ### Astral-Pro can go here
    fi
    if [[ $proteinSelected == 'yes' ]]; then
        if [[ "$astralSelected" == 'yes' ]]; then
            makeSpeciesTree protein "$proteinAstralInfile" '.' astral 'GTR+G' 'DNA' '-nt -gtr'
            if [ -s $treeTipInfoMapFile ]; then
                nw_rename -l  ${fileNamePrefix}.protein.species_tree.astral_pp1_value.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.protein.species_tree.astral_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.protein.species_tree.astral_pp1_value.nwk $numbrLowSupportNodesThreshold astral
        elif [[ "$astralmpSelected" == 'yes' ]]; then 
            makeSpeciesTree protein $proteinAstralInFile '.' astralmp 'GTR+G' 'DNA' '-nt -gtr'
            if [ -s $treeTipInfoMapFile ]; then
                nw_rename -l  ${fileNamePrefix}.protein.species_tree.astralmp_pp1_value.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.protein.species_tree.astralmp_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.protein.species_tree.astralmp_pp1_value.nwk $numbrLowSupportNodesThreshold astral
        fi
        ### Astral-Pro can go here
    fi
    ### 7.5.2022 - shoudln't there be a codon step here?????
fi # End of Astral step


# Creating a zipped tarball for the gene alignments and trees files. 
# These archives can then be moved around easily and checked with a single checksum.
# Doing this before the concatenated alignment trees which will take ages!  
if [[ $dnaSelected == 'yes' ]]; then
    tar -czf ${fileNamePrefix}.dna.aln.for_tree.fasta.tar.gz  *dna.aln.for_tree.fasta
    # Newick files may not exist if only building a concatenated alignment tree:
    if [[ -s $dnaAstralInFile ]]; then
        # Trees will be collapsed or not according to the use of option -L; $dnaAstralInFile variable 
        # already contains file with either tree type, with nodes collapsed or not.
        # Using the file of filnames here to create the archive, rather than the file containing all 
        # the trees (and printing each tree per line) because very big trees will start to print to > 1 line 
        # once ARG_MAX is reached)
        tar -czf ${dnaGeneTreesSetFilenames}.tar.gz `cat $dnaGeneTreesSetFilenames`
    fi
    # Also reating a zipped tarball for the original unaligned gene-wise fasta files.
    if ls ../*dna.fasta >/dev/null 2>&1; then
        tar -czf ${fileNamePrefix}.dna.fasta.tar.gz ../*dna.fasta
    fi
fi
if [[ $proteinSelected == 'yes' ]]; then
    tar -czf  ${fileNamePrefix}.protein.aln.for_tree.fasta.tar.gz  *protein.aln.for_tree.fasta
    if [[ -s $proteinAstralInFile ]]; then
        tar -czf ${proteinGeneTreesSetFilenames}.tar.gz `cat $proteinGeneTreesSetFilenames`
    fi
    if ls ../*protein.fasta >/dev/null 2>&1; then
        tar -czf ${fileNamePrefix}.protein.fasta.tar.gz ../*protein.fasta
    fi
fi
#### NEED TO ADD CODON CONDITIONAL AS WELL HERE


if [[ "$fasttreeSelected" == 'yes' || "$raxmlSelected" == 'yes' ]]; then
    echo 'Concatenating gene alignments...'
    # Concatenate the alignment for 2 reasons:
    # 1.easy browsing in Jalview - yes this works OK, slow for very big trees though
    # 2.for running RAxML with a supermatrix tree.

    if [[ $dnaSelected == 'yes' ]]; then
        AMAS.py concat  -c 1 \
        -i `cat mafft_dna_alns_fasta_file_list.txt` \
        --in-format  fasta \
        -d dna \
        --out-format fasta \
        -t dna.alns.concatenated.fasta
        # NB - also creates a partitions.txt file with coords for each gene - preparing for use with RAxML:
        cat partitions.txt | awk -F '_' '{print "DNA, " $2}' > ${fileNamePrefix}.dna.alns.partitions.txt

        # Compressing file for easier transfer of large alignments:
        gzip -c dna.alns.concatenated.fasta \
        > ${fileNamePrefix}.dna.alns.concatenated.fasta.gz
    fi
    if [[ $proteinSelected == 'yes' ]]; then
        # Also concatenate protein alns:
        AMAS.py concat  -c 1 \
        -i `cat mafft_protein_alns_fasta_file_list.txt` \
        --in-format  fasta \
        -d aa \
        --out-format fasta \
        -t protein.alns.concatenated_temp.fasta
        # NB - also creates a partitions.txt file with coords for each gene - preparing for use with RAxML:
        cat partitions.txt | awk -F '_' '{print "JTT, " $2}' > ${fileNamePrefix}.protein.alns.partitions.txt
        ### NB - model is hard coded at the moment

        ### On cluster the protein fasta headers now contain _\[translate(1)\] - fastatranslate (v2.4.x) adds this string to the header,
        ### then AMAS adds it onto the required fasta record name - so remove them here:
        cat protein.alns.concatenated_temp.fasta \
        | sed 's/_\[translate(1)\]//' \
        > protein.alns.concatenated.fasta
        rm protein.alns.concatenated_temp.fasta

         # Compressing file for easier transfer of large alignments:
        gzip -c protein.alns.concatenated.fasta \
        > ${fileNamePrefix}.protein.alns.concatenated.fasta.gz
    fi
    #### NEED TO ADD CODON CONDITIONAL AS WELL HERE


    ### IF filtewrSeqs1 IS ON
    ###     use filterSeqs() function here
    ### NB TRIMMING SHOUDL ALREADY HAVE BEEN DONE



    if [[ "$fasttreeSelected" == 'yes' ]]; then
        if [[ $dnaSelected == 'yes' ]]; then
            makeSpeciesTree dna dna.alns.concatenated.fasta '.' fasttree 'not_required_here' 'not_required_here' '-nt -gtr'
            # Add tree tip info:
            if [[ -s $treeTipInfoMapFile ]]; then 
                nw_rename -l  ${fileNamePrefix}.dna.species_tree.fasttree.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.dna.species_tree.fasttree_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.dna.species_tree.fasttree.nwk $numbrLowSupportNodesThreshold fasttree
        fi
        if [[ $proteinSelected == 'yes' ]]; then
            makeSpeciesTree protein protein.alns.concatenated.fasta '.' fasttree 'not_required_here' 'not_required_here' ''
            # Add tree tip info:
            if [[ -s $treeTipInfoMapFile ]]; then 
                nw_rename -l  ${fileNamePrefix}.protein.species_tree.fasttree.nwk \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.protein.species_tree.fasttree_USE_THIS.nwk
            fi
            getTreeStats ${fileNamePrefix}.protein.species_tree.fasttree.nwk $numbrLowSupportNodesThreshold fasttree
        fi
        #### NEED TO ADD CODON CONDITIONAL AS WELL HERE
    fi
    if [[ "$raxmlSelected" == 'yes' ]]; then
        if [[ $dnaSelected == 'yes' ]]; then
            RAxML_ModelOfEvolution=GTRCAT
            phyloProgramSwitch=raxml
            if [[ $totalNumbrSamples -lt 200 ]]; then   # 21.4.2020 - set to 200, only really need speed up with larger trees 
                RAxML_ModelOfEvolution=GTRGAMMA
            fi
            # if [[ $totalNumbrSamples -lt 20 ]]; then    # Maybe have to adjust this if number samples over 20 give this error: 'WARNING the alpha parameter with a value of' 
            #     RAxML_ModelOfEvolution=   [ use GTR only ]  - will need to use IQTREE2 or raxml-ng; can use RAxML -V option with the CAT model but number of sample here is low!
            #     THEN: makeSpeciesTree dna dna.alns.concatenated.fasta '.' [IQTREE2 or raxml-ng - I think IQTREE is better] $RAxML_ModelOfEvolution 'not_required_here' 'not_required_here'
            #     ALSO: need to include the option for partitioning here as well
            # else USE existing command below - or use the phyloProgramSwitch variable above and just have ONE call to makeSpeciesTrees()
            ### NB - need to check that the chosen program is installed here
            ### For turning on partitioning or not:     [ also add to protein cmd ]     
            ### if [[ $speciesTreeProgram == *'raxmlq'* ]];then
            ###    makeSpeciesTree dna dna.alns.concatenated.fasta '.' raxmlq $RAxML_ModelOfEvolution 'not_required_here' 'not_required_here'
            ### else
            makeSpeciesTree dna dna.alns.concatenated.fasta '.' raxml $RAxML_ModelOfEvolution 'not_required_here' 'not_required_here'
            # Add tree tip info:
            if [[ -s $treeTipInfoMapFile ]]; then 
                nw_rename -l  RAxML_bipartitions.${fileNamePrefix}.dna.raxmlHPC-PTHREADS-SSE \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.dna.species_tree.raxml_${numbrBootstraps}bs_USE_THIS.nwk
            fi
            getTreeStats RAxML_bipartitions.${fileNamePrefix}.dna.raxmlHPC-PTHREADS-SSE $numbrLowSupportNodesThreshold astral
        fi
        if [[ $proteinSelected == 'yes' ]]; then
            RAxML_ModelOfEvolution=PROTCATJTT
            if [[ $totalNumbrSamples -lt 200 ]]; then    # 21.4.2020 - set to 200, only really need speed up with larger trees 
                RAxML_ModelOfEvolution=PROTGAMMAJTT
            fi
            makeSpeciesTree protein protein.alns.concatenated.fasta '.' raxml $RAxML_ModelOfEvolution 'not_required_here' 'not_required_here'
            # Add tree tip info:
            if [[ -s $treeTipInfoMapFile ]]; then 
                nw_rename -l   RAxML_bipartitions.${fileNamePrefix}.protein.raxmlHPC-PTHREADS-SSE \
                tree_tip_info_mapfile.txt \
                > ${fileNamePrefix}.protein.species_tree.raxml_${numbrBootstraps}bs_USE_THIS.nwk
            fi
            getTreeStats RAxML_bipartitions.${fileNamePrefix}.protein.raxmlHPC-PTHREADS-SSE $numbrLowSupportNodesThreshold astral   
        fi
        #### NEED TO ADD CODON CONDITIONAL AS WELL HERE
    fi
fi