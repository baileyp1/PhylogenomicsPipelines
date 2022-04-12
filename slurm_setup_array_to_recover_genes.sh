#!/bin/bash
#SBATCH -J recover_genes
#SBATCH -n 1
#SBATCH -o recover_genes-%A-%a.log	# If this line is not set, default output name is slurm-%A_%a.out: %A = $SLURM_ARRAY_JOB_ID; %a = $SLURM_ARRAY_TASK_ID		 
#SBATCH -e recover_genes-%A-%a.log

csvFile=$1					### This should be samples names - I think it also can contain the full line
targetsFile=$2
paftolDataSymlinksDir=$3
adapterFasta=$4
samplePrefix=$5
cpu=$6
pathToScript=$7
exePrefix="$8"
hybSeqProgram=$9
usePaftolDb=${10}
stats=${11}
refFilePathForStats=${12} 
echo Inside Slurm array script, csvFile: $csvFile
echo Inside Slurm array script, exePrefix: $exePrefix


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
mapfile -t "SAMPLELIST" < $csvFile
echo Sample info: "${SAMPLELIST[$SLURM_ARRAY_TASK_ID]}"
echo


# Finally input the line of current sample into the worker script.
### Confirm I don't need srun here - if I used it maybe I could get the memory used.
$pathToScript/recover_genes_from_one_sample.sh "${SAMPLELIST[$SLURM_ARRAY_TASK_ID]}"  $targetsFile  $paftolDataSymlinksDir  $adapterFasta  $samplePrefix  $cpu  "$exePrefix"  $hybSeqProgram $usePaftolDb $stats $refFilePathForStats 
sleep 1		# Not yet sure whether sleep positioned here makes the script sleep between suubmitting samples