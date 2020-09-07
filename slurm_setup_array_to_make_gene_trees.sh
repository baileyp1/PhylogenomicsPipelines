#!/bin/bash
#SBATCH -J make_gene_trees
#SBATCH -n 1
#SBATCH -o make_gene_trees-%a.log		# If this line is not set, default output name is slurm-%A_%a.out: %A = $SLURM_ARRAY_JOB_ID; %a = $SLURM_ARRAY_TASK_ID		 
#SBATCH -e make_gene_trees-%a.log 		# NB - if I just specify the %a, then I don't get an accumulation of ouput files for each run of the script


geneFile=$1
listFile=$2
fractnAlnCovrg=$3
pathToScripts=$4
phyloProgramDNA=$5
phyloProgramPROT=$6
fractnMaxColOcc=$7
cpuGeneTree=$8
mafftAlgorithm="$9"
exePrefix="${10}"
alnProgram="${11}"
dnaSelected="${12}"
proteinSelected="${13}"
codonSelected="${14}"
filterSeqs1="${15}"
pathToScripts="${16}"

echo Inside Slurm array script, listFile: $listFile
echo Inside Slurm array script, fractnAlnCovrge: $fractnAlnCovrg


# When --array is specified, any commands in this script are executed and repeated for each element of the array
# e.g. these lines will appear in each slurm-%A_%a.out file for each element of the array. 
echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
echo "SLURM_ARRAY_TASK_COUNT: " $SLURM_ARRAY_TASK_COUNT
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID		# Should be a different value reported in each eio_pipeline-%a.out file
echo


# Need to associate different input filenames or info for different samples to each sbatch array element.
# The bash 'mapfile' command places each line of the input file to a Bash array called SAMPLELIST:
### NB - the only thing I don't understand is that the array is created again and again for each sbatch array element which is OK but redundant.
mapfile -t SAMPLELIST < $listFile
echo GeneId: ${SAMPLELIST[$SLURM_ARRAY_TASK_ID]}
echo


# Finally input the line of current sample into the worker script.
### Confirm I don't need srun here - if I used it maybe I could get the memory used.
exePrefix="/usr/bin/time -v "		# this time command gets the RSS memory
$exePrefix $pathToScripts/make_gene_trees.sh \
${SAMPLELIST[$SLURM_ARRAY_TASK_ID]} \
$geneFile \
$fractnAlnCovrg \
$phyloProgramDNA \
$phyloProgramPROT \
$fractnMaxColOcc \
$cpuGeneTree \
"$mafftAlgorithm" \
"$exePrefix" \
"$alnProgram" \
$dnaSelected \
$proteinSelected \
$codonSelected \
"$filterSeqs1" \
$pathToScripts \
$maxColOccThreshold
### 16.3.2020 - could remove this $$exePrefix - but it might be useful to see whether sub-processes get added to the mem used for this script.
### 3.6.2020 - yes they do - the time and mem from this command is the total of tiem within script AND max mem used within the child processes - good