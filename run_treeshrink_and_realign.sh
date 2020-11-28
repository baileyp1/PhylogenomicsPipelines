#!/bin/bash

###################
# run_treeshrink.sh

# Purpose: runs treeshrink on existing gene tree Newick files, then realigns surviving sequences 
#          again in gene-wise mode via the master script, make_species_trees_pipeline.sh.
#
# Usage: an internal pipeline script

# Author: Paul Bailey
###################

shopt -s failglob

numbrSamples="$1"
phyloProgramDNA="$2"
phyloProgramPROT="$3"
sampleTableFile="$4"
geneListFile="$5"
fractnAlnCovrg="$6"
fractnMaxColOcc="$7"
fractnSamples="$8"
mafftAlgorithm="$9"
cpuGeneTree="${10}"
partitionName="${11}"
pathToScripts="${12}"
geneTreesOnly="${13}"
dnaSelected="${14}"
proteinSelected="${15}"
codonSelected="${16}"
treeshrink="${17}"
filterSeqs1="${18}"
alnProgram="${19}"
maxColOccThreshold="${20}"
filterSeqs2="${21}"
trimAln1="${22}"
trimAln2="${23}"
collapseNodes="${24}"
fileNamePrefix="${25}"
cpu="${26}"
geneTreeSlurmMem="${27}"
speciesTreeSlurmMem="${28}"
geneTreeSlurmTime="${29}"
speciesTreeSlurmTime="${30}"
option_u="${31}"


echo "$numbrSamples"
echo "$phyloProgramDNA"
echo "$phyloProgramPROT"
echo "$sampleTableFile"
echo "$geneListFile"
echo "$fractnAlnCovrg"
echo "$fractnMaxColOcc"
echo "$fractnSamples"
echo "$mafftAlgorithm"
echo "$cpuGeneTree"
echo "$partitionName"
echo "$pathToScripts"
echo maxColOccThreshold: $maxColOccThreshold
echo "geneTreesOnly: $geneTreesOnly"
echo "filterSeqs1: $filterSeqs1"
echo "filterSeqs2: $filterSeqs2"
echo "trimAln1: $trimAln1"
echo "trimAln2: $trimAln2"
echo "collapseNodes: $collapseNodes"
echo "fileNamePrefix: $fileNamePrefix"
echo "cpu: $cpu"
echo "geneTreeSlurmMem: $geneTreeSlurmMem"
echo "speciesTreeSlurmMem: $speciesTreeSlurmMem"
echo "geneTreeSlurmTime: $geneTreeSlurmTime"
echo "speciesTreeSlurmTime: $speciesTreeSlurmTime"
echo "option_u: $speciesTreeSlurmTime"

# Convert $emptyMatchStateFractn to a percent for use in the output files:
fractnAlnCovrg_pc=`awk -v FRACTN=$fractnAlnCovrg 'BEGIN{printf "%.0f", FRACTN * 100}' `


runTreeShrink() {
    ###########
    # Function: prepares files (DNA or protein), runs TreeShrink and 
    #           set up gene-wise files for running make_species_trees_pipeline.sh 
    # $1 = $seqType
    # $2 = $phyloProgramMethod  - not using now
    # All global variables are also available so no need to bring them in as such - keep confirming...
    ###########

    echo "Entered Treeshink function...."


    # To use Treeshrink, need to prepare a directory to hold each tree and alignment per gene.
    # Current aln name: *.${1}.aln.for_tree.fasta
    # Current tree name: *_${1}_gene_tree_USE_THIS.nwk
    if [[ ! -d treeshrink_${1}_gene_trees ]]; then mkdir treeshrink_${1}_gene_trees; fi
    for file in  ../*_${1}_gene_tree_USE_THIS.nwk; do
        # Also need to remove the relative path at the front!
        gene=`echo $file | sed "s/_${1}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
        # echo "gene: $gene"

        if [[ ! -d treeshrink_${1}_gene_trees/$gene  ]]; then mkdir treeshrink_${1}_gene_trees/$gene ; fi

        # Copy files into the gene directory - need to give them a name common to all the gene tree and aln files:
        cp -p $file  treeshrink_${1}_gene_trees/$gene/${1}_gene_tree_USE_THIS.nwk
        cp ../${gene}.${1}.aln.for_tree.fasta  treeshrink_${1}_gene_trees/$gene/${1}_gene_tree_aln.fasta 
    done

    bParameter=20
    run_treeshrink.py -b $bParameter -i treeshrink_${1}_gene_trees  -t ${1}_gene_tree_USE_THIS.nwk -a ${1}_gene_tree_aln.fasta -O ${1}_gene_tree
    # Summary of usage:
    # Three modes: 'per-gene', 'all-genes', 'per-species' - auto will select per-species 
    #              unless there are rare species (i.e. a species that occurs in less than 20 gene trees),
    #              in which case 'all-genes' mode is used.
    #              Will apply to some species - need to understand the difference between the two modes  
    # -q quantiles threshold - default = 0.05 - can use multiple values e.g. -q "0.05 0.10" - might be useful to do
    # -k and -s control the maximum number of species that can be removed
    # -b controls if a species is removed if it has a significant impact on the tree diameter (maximum distance between any two leaves of the tree) - set by default to 5%
    #    NB - documentation says: a higher value, like -b 20 may make more sense for your dataset if you want to be more conservative. We suggest exploring this option.
    #         i.e. at -b 20 the species would have to have a higher distance before being considered for removal

    # Output filenames:
    # 1. dna_gene_tree_RS_shrunk_0.05.txt   - contains a tab-separated list of sample names removed
    #                                       - I think the 0.05 refers to option -q = 0.05 by default - can't confirm until I change this option
    ### NB - there is also another file dna_gene_tree_USE_THIS_shrunk_RS_0.05.txt - don't know what this is - empty so far.
    # 2. dna_gene_tree_tree_shrunk_0.05.nwk  -  NB - these trees don't have bootstrap values so file size appears smaller
    # 3. dna_gene_tree_aln_shrunk_0.05.fasta

    # TreeShrink results
    # Count number of samples removed from any tree:
    numbrSeqsRemoved=`cat treeshrink_${1}_gene_trees/*/${1}_gene_tree_RS_shrunk_0.05.txt | sed 's/\t/\n/g' | sort -u | wc -l`
    echo "TreeShrink results
==================
-b parameter = $bParameter
Number of samples in data set before TreeShrink: $numbrSamples
Number of samples removed from any gene tree: $numbrSeqsRemoved " > ${fileNamePrefix}_treeshrink_results.txt
    # Only print if > 0 samples have been removed:
    if [[ $numbrSeqsRemoved -ge 1 ]]; then
        echo "NumberOfGeneTrees RemovedSampleFrom
`cat treeshrink_${1}_gene_trees/*/${1}_gene_tree_RS_shrunk_0.05.txt | sed 's/\t/\n/g' | grep -v '^$' | sort | uniq -c | sort -k1nr` " >> ${fileNamePrefix}_treeshrink_results.txt
    fi

 
    # Move the TreeShrunk alignment files to *_dna.fasta for each gene after unaligning them by removing all '-' gaps:
    ### NB - will treeshrink always produce a file even if it has no contents? - look at code
    # Using .nwk file as a guide i.e. you only create alignments again for existing gene trees (i.e. ignoring any gene trees that have been filtered for whatever reason)
    # for file in  ../*_${1}_gene_tree_USE_THIS.nwk; do
    #     gene=`echo $file | sed "s/_${1}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
    #     cat treeshrink_${1}_gene_trees/$gene/${1}_gene_tree_aln_shrunk_0.05.fasta \
    #     | awk '{ if($1 ~ /^>/) {print $1} else { gsub(/-/,"",$1);print $1 }  }' \
    #     > ${gene}_after_treeshrink.fasta
    # done
    # 21.8.2020 - Now having to use an alternative approach to obtain the gene-wise DNA sequence file.
    # Still using the Newick file to obtain a list of seq names but then extract sequences from
    # the original genewise file.
    # (This approach is now necessary after TreeShrink with protein because we need to start again from DNA but the above code
    # has removed gaps from a protein aln - no good! Could go looking for the DNA aln but 'dna' might not have been seelected.)
    ### NBNB - this means that there is no point in using Treeshrink -a flag, but could just leave as is.
    for file in  ../*_${1}_gene_tree_USE_THIS.nwk; do
        gene=`echo $file | sed "s/_${1}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
        # Get the leaf labels from the TreeShrunk Newick file (could also use the equivalent aln file) then extract from the unaligned starting DNA file:
        ### NB - 21.10.2020 - just realised I have this wrong - need to use the TS output to get the retained labels:
        ### nw_labels -I $file > ${gene}_${1}_tree_leaf_labels.txt
        nw_labels -I treeshrink_${1}_gene_trees/$gene/${1}_gene_tree_tree_shrunk_0.05.nwk > ${gene}_${1}_tree_leaf_labels.txt
        seqtk subseq -l 0 \
        ../${gene}_dna.fasta \
        ${gene}_${1}_tree_leaf_labels.txt \
        > ${gene}_after_treeshrink.fasta
    done
}


reAlignSeqs()   {
    ###########
    # Function: realigns sequences from gene-wise files that have passed through a filtering, trimming or
    #           a TreeShrink step. 
    # $1 = $seqType
    # $2 = $phyloProgramMethod
    # $3 = log file suffix and used for the input fasta files
    # All global variables are also available so no need to bring them in as such - keep confirming...
    ###########

    # Now ready to redo gene trees using gene-wise alignments.
    # First need to alter option -t, -g and any other future option 
    # containing an input file to adjust a relative path to the new pwd,
    # or not if it is a full path.
    # NB - need to use a fresh variable otherwise relative path is extended after every call of reAlignSeqs() !

    # Ammend option -t path if required and selected:
    uOption=''
    if [[ $sampleTableFile != 'no' ]]; then
        if [[ $sampleTableFile != /* ]]; then 
            sampleTableFileForReAln="../$sampleTableFile"
        else
            sampleTableFileForReAln=$sampleTableFile
        fi

        if [ $option_u == 'yes' ]; then
            uOption='-u'
            # Prepare the mapfile again to add in the sum length values onto the tree tips, if requested
            # For the realignment need to get this file into the new directory!
            # addTreeTipInfoFromTable_PlusL.py \
            # ../tree_tip_info_mapfile.txt \
            # ../${fileNamePrefix}_summary_of_sample_quality.txt \
    #### Actually just need to copy it across here
            cp -p tree_tip_info_mapfile_plusL.txt .
            # NB - still need to copy this file to 'tree_tip_info_mapfile.txt' in the species
            # script because it gets prepared again at start of the realn step.
        fi
    fi
    ls -l $sampleTableFileForReAln

    # Amend option -g path if necessary (mandatory option so always have to check):
    if [[ $geneListFile != /* ]]; then 
        geneListFileForReAln="../$geneListFile"
    else
        geneListFileForReAln=$geneListFile
    fi
    ls -l $geneListFileForReAln
    # exit


    # Now run parent script again in gene-wise mode on a uniquely named a set of files.
    # If user has selected to build gene trees ONLY, then add that option flag.
    iOption=''
    if [ $geneTreesOnly == 'yes' ]; then iOption='-i'; fi
    
    trimAlnOption1=''
    trimAlnOption2=''
    if [ $trimAln1 != 'no' ]; then trimAlnOption1='-J'; fi
    if [ $trimAln2 != 'no' ]; then trimAlnOption2="-K $trimAln2"; fi
    #echo trimAlnOption1: $trimAlnOption1  
    #echo trimAlnOption2: $trimAlnOption2
    echo collapseNodes: $collapseNodes
    collapseNodesOption=''
    if [ $collapseNodes != 'no' ]; then collapseNodesOption="-L $collapseNodes"; fi


    ### At the moment these files don't exist in /after_treeshrink_USE_THIS_dna/ dir
    ### which exists by this stage - but they do exist in the above dir so will copy 
    ### them over here so can now make a concatenated aln and tree:
    cp ../mafft_*_alns_fasta_file_list.txt .


    # Note the quotes around variables with spaces! BUT $2 must not be quoted!
    #echo " # For checking option values that need to be quoted (contain spaces)
    $pathToScripts/make_species_trees_pipeline.sh $iOption $trimAlnOption1 $trimAlnOption2 $collapseNodesOption $uOption \
    -p $fileNamePrefix \
    -G \
    -D "$1" \
    -A $alnProgram \
    -M "$mafftAlgorithm" \
    -t $sampleTableFileForReAln \
    -g $geneListFileForReAln \
    -m $fractnMaxColOcc \
    -O $maxColOccThreshold \
    $2 \
    -C $cpuGeneTree \
    -c $cpu \
    -Q $partitionName \
    -R $geneTreeSlurmMem \
    -U $speciesTreeSlurmMem \
    -V $geneTreeSlurmTime \
    -W $speciesTreeSlurmTime \
    *_${3}.fasta \
    > make_species_trees_pipeline_${3}.log 2>&1 #"
}





# Main code:
if [[ $treeshrink == 'yes' ]]; then
    
    echo
    echo
    echo ###################################
    echo Running TreeShrink on gene trees...
    echo ###################################
    echo
    echo

    if [[ $dnaSelected == 'yes' ]]; then
        #echo "Entered TreeShrink for DNA..."
        seqType=dna         # required for file and directory names     ### new var required now?
        if [[ ! -d after_treeshrink_USE_THIS_$seqType ]]; then mkdir after_treeshrink_USE_THIS_$seqType; fi
        cd after_treeshrink_USE_THIS_$seqType
        pwd
        # Function only knows about one alignment filename (the one after all filtering/trimming, if any)
### There is an issue with codon alns - need to move into the codonAln dir - just need to inlcude 
### filoe path varible like I've done in makeGeneTree function
        runTreeShrink $seqType
        phyloProgramsToUse="-q $phyloProgramDNA"
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_treeshrink"
        cd ../  # Back to main dir to apply TreeShrink to protein seqs if required
    fi
    if [[ $proteinSelected == 'yes' ]]; then
        #echo "Entered TreeShrink for protein..."
        seqType=protein
        if [[ ! -d after_treeshrink_USE_THIS_$seqType ]]; then mkdir after_treeshrink_USE_THIS_$seqType; fi
        cd after_treeshrink_USE_THIS_$seqType
        pwd
        runTreeShrink $seqType
        phyloProgramsToUse="-r $phyloProgramPROT"
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_treeshrink"
        cd ../  # Back to main dir to apply TreeShrink to codon seqs if required
    fi
    if [[ $codonSelected == 'yes' ]]; then
        #echo "Entered TreeShrink for codon..."
        seqType=codon
        if [[ ! -d after_treeshrink_USE_THIS_$seqType ]]; then mkdir after_treeshrink_USE_THIS_$seqType; fi
        cd after_treeshrink_USE_THIS_$seqType
        pwd
### NB - see above for change to make for path to codn dir 
        runTreeShrink $seqType
        phyloProgramsToUse="-q $phyloProgramDNA"
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_treeshrink"
        cd ../
    fi  
elif [[ $filterSeqs1 != 'no' || $filterSeqs2 != 'no' ]]; then
    echo
    echo
    echo ##########################################################
    echo Re-aligning gene alignments before making species trees...
    echo ##########################################################
    echo
    echo
    # Realign sequences if Treeshrink not set but seqs have been filtered - working directory: after_reAlnFilterSeqs_USE_THIS
    # First re-create the -D flag value for residue type(s):
    # For re-running the gene trees, need to prepare which residue type(s) have been selected:
    seqType=''
    phyloProgramsToUse=''
    if [[ $dnaSelected == 'yes' ]]; then
        seqType="$seqType dna"
        phyloProgramsToUse="$phyloProgramsToUse -q $phyloProgramDNA"
    fi
    if [[ $proteinSelected == 'yes' ]]; then
        seqType="$seqType protein"
        phyloProgramsToUse="$phyloProgramsToUse -r $phyloProgramPROT"
    fi
    if [[ $codonSelected == 'yes' ]]; then
        seqType="$seqType codon"
        phyloProgramsToUse="$phyloProgramsToUse -q $phyloProgramDNA"
    fi
    echo seqType for filterSeqs option: $seqType
    echo phyloProgramsToUse for filterSeqs option: $phyloProgramsToUse
    if [[ ! -d after_reAlnFilterSeqs_USE_THIS ]]; then mkdir after_reAlnFilterSeqs_USE_THIS; fi
    cd after_reAlnFilterSeqs_USE_THIS
    # Copy the (filtered) alignment DNA files to *_after_filterSeqs.fasta for each gene after unaligning them by removing all '-' gaps:
    # Using .nwk file as a guide i.e. you only create alignments again for existing gene trees 
    # (i.e. ignoring any gene trees that have been filtered for whatever reason)

### Need to determine which residue type to loop through for two reasons:
### 1. just need to know which dataset(s) exist: use protein if it exists, then codon, then DNA
###    Prioritizing protein if protein has been selected (should give better alns), then codon, then DNA.
### 2. Some trees may filter slightly differently between the residue types with filterSeqs1 but really need to choose one option
###    here, otherwise it gets too complicated because I need to start with single set of fasta files for next iteration of the script
###    So prioritizing protein if protein has been selected.
###    (should give better alns), then codon, then DNA.
###    If you want to use DNA no matter what, then can still set the -D option to 'dna' only.

### NBNB - only can start realigning seq with script using DNA seq, therefore need to get a list of fasta records from the right residue file,
###        otherwise at some point there will be a zero byte file as input which won't work! Done - now checking...
### Will also need to bring in the path so I can grab the previuous codon file.
    if [[ `echo $seqType | grep -o 'protein' ` == 'protein' ]]; then
        treeType=protein
        for file in  ../*_${treeType}_gene_tree_USE_THIS.nwk; do
            gene=`echo $file | sed "s/_${treeType}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
            # Get the leaf labels and extract from the unliagned starting DNA file:
            nw_labels -I $file > ${gene}_${treeType}_tree_leaf_labels.txt
            seqtk subseq -l 0 \
            ../${gene}_dna.fasta \
            ${gene}_${treeType}_tree_leaf_labels.txt \
            > ${gene}_after_filterSeqs.fasta
        done
        echo Before reAlnSeqs:  seqType: $seqType, mafftAlgorithm: $mafftAlgorithm 
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_filterSeqs"
    elif [[ `echo $seqType | grep -o 'codon' ` == 'codon' ]]; then
        treeType=codon
        for file in  ../*_${treeType}_gene_tree_USE_THIS.nwk; do
            gene=`echo $file | sed "s/_${treeType}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
            # Get the leaf labels and extract from the unliagned starting DNA file:
            nw_labels -I $file > ${gene}_${treeType}_tree_leaf_labels.txt
            seqtk subseq -l 0 \
            ../${gene}_dna.fasta \
            ${gene}_${treeType}_tree_leaf_labels.txt \
            > ${gene}_after_filterSeqs.fasta
        done
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_filterSeqs"
     elif [[ `echo $seqType | grep -o 'dna' ` == 'dna' ]]; then 
        treeType=dna
        for file in  ../*_${treeType}_gene_tree_USE_THIS.nwk; do
            gene=`echo $file | sed "s/_${treeType}_gene_tree_USE_THIS.nwk//" | sed "s/^\.\.\///" `
            # Get the leaf labels and extract from the unliagned starting DNA file:
            nw_labels -I $file > ${gene}_${treeType}_tree_leaf_labels.txt
            seqtk subseq -l 0 \
            ../${gene}_dna.fasta \
            ${gene}_${treeType}_tree_leaf_labels.txt \
            > ${gene}_after_filterSeqs.fasta
        done
        reAlignSeqs "$seqType" "$phyloProgramsToUse" "after_filterSeqs"
    fi
fi
