#!/bin/bash
# PASTA_taster.sh - run this first

# This script tests if the pasta is gluten-free.
# JK. This is a wrapper around optrimAl that generates trimmed alignments using trimAl 1.2 (Capella-Gutierrez et al, 2009) at successively stricter trimming thresholds then summarises statistics for these trimmed alignments using AMAS (Borowiec, 2016).
# It calls the optrimAl script, which then returns all alignment files trimmed to an optimum threshold defined as yielding the maximum proportion of parsimony informative characters but losing no more data than one median absolute deviation above the median data loss across the entire range of trimming thresholds being tested.
# Alignments that lose more than a set cap of data (default 30% in the script) after optimal trimming are not returned.
# Theoretically, I prefer this approach because it considers the amount of missing data in each data set and avoids excessive trimming, instead of setting an arbitrary fixed gap threshold, which DOES result in loss of informativeness in some data sets.
# Empirically, this approach has NOT been tested.
# If you find any bugs, please modify accordingly and let me know (totally optional). I will probably not have any time to troubleshoot in the near future (maybe when I'm finally retired). Good luck! ZQ

# This script produces ALOT of output.
# Alignment files (e.g. *.aln) returned to the working directory are the optimally trimmed alignments.
# overlost.txt lists the alignments where data loss exceeded the cap.
# dldp*.png are graphs showing the proportion of parsimony informative characters and data loss at each trimming threshold, as well as the selected trimming threshold, for each alignment.
# dldp*.csv are the raw data from which the graphs are produced.
# summary*.txt are the summary statistics produced by AMAS.
# Directories named with the specified trimming threshold values (e.g. 0.1) should be deleted immediately once done with analysis as they take up ALOT of space.

# Provide a text file with desired trimming threshold values, one per line (cutoff_trim.txt).
# Make sure to set the working directory and trimAl path correctly, that optrimal.R and cutoff_trim.txt are in the same directory, update the file name pattern for the alignment where required and provide a set of trimming thresholds (any number of thresholds from 0 to 1, one threshold per line, must include 0 and 1) in the cutoff_trim.txt input file.
# My working directory in this case was ‘~/zq/working/optrimal’, my trimAl path was ‘~/zq/bin/trimAl/source/trimal’ and my file name pattern was 'g*' so just change those accordingly.
# This script WILL generate non-fatal errors where alignments are missing - check if these alignments were intentionally omitted or went missing for some other reason.

# Paul B. - setup:
gene=$1
alignmentFile=$2
residueType=$3      
outAlnFile=$4
pathToScripts=$5
#echo $gene  $alignmentFile  $residueType  

if [[ ! -d ${gene}.${residueType}.aln.optrimal ]]; then mkdir ${gene}.${residueType}.aln.optrimal; fi
cd ${gene}.${residueType}.aln.optrimal

# Create the cutoff_trim.txt file:
echo "0
0.05
0.10
0.15
0.20
0.25
0.30
0.35
0.40
0.45
0.50
0.55
0.60
0.65
0.70
0.75
0.80
0.85
0.90
0.95
1.00" > cutoff_trim.txt


while read cutoff_trim
do
    #cd ~/zq/working/optrimal
    #mkdir $cutoff_trim
    if [[ ! -d $cutoff_trim ]]; then mkdir $cutoff_trim; fi

    ###for alignment in g*  # Paul B - I would like to process each alignment file separately!
    for alignment in ../$alignmentFile # Could remove this loop, OK for now
    do
        ls -l $alignment
        ###~/zq/bin/trimAl/source/trimal -in ${alignment}/output_alignment.fasta -out ${cutoff_trim}/${alignment}.aln -htmlout ${cutoff_trim}/${alignment}.htm -gt $cutoff_trim
        trimal -in $alignment -out $cutoff_trim/${gene}.${residueType}.aln.optrimal.fasta -htmlout $cutoff_trim/${gene}.${residueType}.aln.optrimal.htm -gt $cutoff_trim

        # check if alignment was trimmed to extinction by trimAl
        #if grep ' 0 bp' ${cutoff_trim}/${alignment}.aln
        if grep ' 0 bp' $cutoff_trim/${gene}.${residueType}.aln.optrimal.fasta
        then
            #rm -f ${cutoff_trim}/${alignment}.aln
            rm -f $cutoff_trim/${gene}.${residueType}.aln.optrimal.fasta
        fi
    done
    #cd ~/zq/working/optrimal/${cutoff_trim}
    #python3 ~/zq/bin/AMAS-master/amas/AMAS.py summary -f fasta -d dna -i *.aln
    AMAS.py summary -f fasta -d $residueType -i ${cutoff_trim}/${gene}.${residueType}.aln.optrimal.fasta -o summary_${cutoff_trim}.txt
    ####mv summary.txt ../summary_${cutoff_trim}.txt

done < cutoff_trim.txt

###xvfb-run Rscript –vanilla optrimal.R
Rscript $pathToScripts/optrimal.R

# Paul B. - Moving chosen alignment file to pwd of main pipeline so it can be picked up there for next step:
mv ${gene}.${residueType}.aln.optrimal.fasta ../${gene}.${residueType}.aln.after_trim1.fasta

# Stats for chosen file from dldp_<gene_name>.dna.aln.optrimal.fasta.csv file:
### Loop through csv, chop at , chars
### But then not sure how to work out which cut off file has been chosen!!!
### Once have info, print to log file


