#!/usr/bin/env python

############################
# various_tasks_in_python.py
#
# Author: Paul Bailey
#
# Copyright (c) 2025 The Board of Trustees of the Royal Botanic Gardens, Kew
#
# Purpose: code to perform various (routine) tasks with Python (with explanations for learning purposes).
#		   Tasks fall within discrete methods - see the docstring of each method for more info 
#
# Usage: an internal pipeline script with a simple interface (no argparse)
# 		 various_tasks_in_python.py <method_name> <option1 e.g. infile>  <option2 e.g. outfile_prefix>  <option3>  <option4> etc
#
# Procedure for adding a method:
# 1. write a method decription here below
# 2. Add extra sys.argv[] variables to an 'option*' variable', as required
# 3. Write the method at the bottom of the methods section.
#	 Import any modules within each method rather than outside
# 4. Copy the method description within the docstring of method itself and elaborate
#	 Make a note of the Python version required e.g. python 2, 3 etc
# 4. Add a clause in the main code section for the method - see bottom of this file
#
#
# detect_stops()
#	Detects STOP codons in a protein aligment, removes sequences with > 1 STOP codon and create stats
#
# orderAlnByTreeOrder()
#	Orders a sequence alignment by the order in a Newick tree file
#
#retrieve_targets_magic():
#	An internal function for retrieve_targets() in various_tasks_in_bash.sh
#
# *** Next method here ***
#
############################
from __future__ import print_function
import sys
import os
import csv
import re
from Bio import SeqIO					# https://biopython.org/DIST/docs/api/;  http://biopython.org/DIST/docs/tutorial/Tutorial.html;  https://biopython.org/wiki/SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


if len(sys.argv) == 1:
	print("ERROR: specify one of more method to use")
	print('List of popular functions available:')
	print('1. detect_stops')
	print('2. orderAlnByTreeOrder')
	print('3. retrieve_targets_magic')
	exit()
if len(sys.argv) >= 2:
	method = sys.argv[1]	# method name
if len(sys.argv) > 2:
	option1 = sys.argv[2]	# Often the main infile	
else:
	print("ERROR: having no infile is set up to exit with an error at the moment")
	exit()	

option2 = ''					# Defined these vars here so I can test whether option2 onwards is empty or not
option3 = ''					
option4 = ''
if len(sys.argv) > 3:			# Testing whether sys.argv has 3 or more elements (otherwise script crashes)
	option2 = sys.argv[3]
if len(sys.argv) > 4:
	option3 = sys.argv[4]
if len(sys.argv) > 5:
	option4 = sys.argv[5]



# Methods section:
def detect_stops(infile, outfilePrefix):
	'''
	Detects STOP codons in a protein aligment, removes sequences with > 1 STOP codon and create stats
	Assumes that a STOP codon is denoted by a '*' char.

	Usage: various_tasks_in_python.py detect_stops <protein_aln_infile>  <outfile_prefix>
	Usage example: /Users/pba10kg/Documents/ProgramFiles/PhylogenomicsPipelines/various_tasks_in_python.py  detect_stops  4848.protein.fasta  4848.protein

	'''

	# Output files:	
	outfile = outfilePrefix + '.0or1_STOP.fasta'
	fh = open(outfile, "w")
	outfile1 = outfilePrefix + '.ovr1_STOP.fasta'
	fh1 = open(outfile1, "w")
	outfile2 = outfilePrefix + '.stats.log'
	fh2 = open(outfile2, "w")


	stopCountr = {}		# Dict to count seqs with x number of stops, each hash key = x number of stops per seq
	for record in SeqIO.parse(infile, "fasta"):		# returns a SeqRecord object, includes a Seq object called seq
		#print(record.id, "\n", record.seq, len(record))
		#print('Number of stops in seq: ', record.seq.count('*'))
		numbrStops = record.seq.count('*')

		if numbrStops in stopCountr:
			stopCountr[numbrStops] += 1
		else:
			stopCountr[numbrStops] = 1
		#print 'stopCountr[', numbrStops, ']:', stopCountr[numbrStops]

		# If sequence has zero or one stops, print to file (for use in the phylogeny).
		# NB - one stop might be the real stop codon so it should be allowed through.
		#      In current work I don't think the sequences have the end STOP codon and 
		#      in any case one STOP and maybe one or two more in a sequence seem to be 
		#      in frame aligned and well with many other sequences, so should be allowed through.
		### Might want to increase or decrease the number of stops in seqs to print.
		### 30.4.2022 - a better idea  might be to let everything through then assess % id of seqs instead

		# UPP refuses to align a sequences containing STOP chars so removing them 
		# before printing seqs to file - actually converting them to 'X' chars (for now).
		# (c.f. MAFFT strips any '?' and '*' chars out of the alignment, retains 'X' chars though)
		seqStopsRemoved = re.sub('\*', 'X', str(record.seq))	# NB - record.seq is a Seq object but the full sequence can be returned as a python string, using e.g. str(my_seq).	
		recordNoStops = '>' + record.id + '\n' + seqStopsRemoved + '\n'

		if numbrStops == 0 or numbrStops == 1:
			
			#print('Print record 0/1 ', record)
			#SeqIO.write(record, fh, "fasta")	# Need to use a file handle here if printing to file multiple times, otherwise previous write gets overwritten, not appended
			fh.write(str(recordNoStops))
		else:
			#print('Print record ovr1 ', record)
			#SeqIO.write(record, fh1, "fasta")
			fh1.write(str(recordNoStops))
		

	# Print STOP stats to summary file:
	#fh2.write('NumbrStops NumbrSeqs' + "\n")
	totalNumbrSeqsWithStopsCountr = 0
	for key in sorted(stopCountr):
		if key == 0:
			fh2.write('numbrSeqsWith_' + str(key) + '_Stops: ' + str(stopCountr[key]) + "\n")
		else:
			totalNumbrSeqsWithStopsCountr += stopCountr[key]
			fh2.write('numbrSeqsWith_' + str(key) + '_Stops: ' + str(stopCountr[key]) + "\n")
	fh2.write('TotalNumbrSeqsWithStops: ' + str(totalNumbrSeqsWithStopsCountr) + "\n")
	fh.close()
	fh1.close()
	fh2.close()


def bioSeqIOLoop(infile):
	'''
	Template use of Bio.SeqIO to manipulate sequence records from a fasta file
	'''

	for record in SeqIO.parse(infile, "fasta"):     # returns a SeqRecord object, includes a Seq object called seq
		#print(record.id, "\n", record.seq, len(record) )
		print(record.description)	# This variable appears to be the full fasta file header
		# Parse the organism name from the fasta header (2nd field in record.description):
		headrFields = record.description.split()
		print('headrFields:', headrFields[1])


def orderAlnByTreeOrder(infile, infile1):
	'''

	***** 7.5.2025 - REMEMBER TO UPDATE A STABLE VERSION TO BIN FOLDER *****

	Orders a sequence alignment by the order in a Newick tree file

	Usage: various_tasks_in_python.py orderAlnByTreeOrder <alnfile> <Newick_file_ordered_tip_list>
	Usage example: various_tasks_in_python.py orderAlnByTreeOrder  6636.protein.aln.for_tree.fasta  6636.protein.guide_gene_tree_labels_test_temp.nwk 

	Ordered tip list can be found for the input e.g. nw_labels -I <Newick_file> > <Newick_file_ordered_tips.txt>

	Output: fasta records in STDOUT ordered by the tree input

	Tested with Python3

	Possible bug: if the tree labels contain single quote characters (which they do coming out of the GTM program, then
	              this method will crash. Consider to remove them like so:
	              row = row.strip(\'\'\').rstrip('\n') - but without the backslashes - NOT TESTED
	              Removing them outside this method for now. 

	'''

	fastaFileIndex = SeqIO.index(infile, "fasta")
	#print(fastaFileIndex["10081"])
	#print(fastaFileIndex["10081"].id)
	#print(fastaFileIndex["10081"].seq)
	#exit()
	with open(infile1, "r") as fh:
		for row in fh:
			row = row.rstrip('\n')
			print(">" + fastaFileIndex[row].id)
			print(fastaFileIndex[row].seq)


def retrieve_targets_magic(sampleId, blast_output_file, fasta_file_for_blast_db, blastProgram):
	'''

	***** 7.5.2025 - REMEMBER TO UPDATE STABLE VERSION TO $HOME/BIN FOLDER *****

	Python3+
	Purpose: Takes BLAST results of hits to the reference targets and organises the fasta header line 
			 main id to refer to the reference gene and the sample id in HybPiper fasta format

	An internal function for retrieve_targets() in various_tasks_in_bash.sh

	Usage: various_tasks_in_python retrieve_targets_magic  <blast_output.tab>  <fasta_file_for_blast_db> 

	Format of blast_output_file: qaccver saccver pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore"
	Example: BEFC-5406	JAMXDD010000010.1_cds_KAI6674735.1_2641	95.286	297	14	0	203	499	499	1	891	1005	0.0	592

	Output file  of fasta records: <sampleId>.fasta
	'''

	# Modules required (for when I put them within each method):
	#from Bio import SeqIO


	# Index the  <fasta_file_for_blast_db> file:
	fastaFileBlastDBDict = SeqIO.index(fasta_file_for_blast_db, "fasta")
	#print(vars(fastaFileBlastDBDict))
	#print(dir(fastaFileBlastDBDict))
	#for key in fastaFileBlastDBDict:
	#	print(fastaFileBlastDBDict[key])
	#	print(fastaFileBlastDBDict[key].seq)  
		#exit()
	#print(fasta_file_for_blast_db['JAMXDD010000010.1_cds_KAI6672684.1_590'])  # use any record ID
	# #print record_dict["JAMXDD010000005.1_cds_KAI6692379.1_20089"].id 
	# #print record_dict["Q39056"].annotations['gene_name']
	# #exit() 


	# Main output file:
	outfile1=f'{sampleId}_{blastProgram}.fasta'
	outfile2=f'{sampleId}_{blastProgram}_trimmed.fasta'
	outfile3=f'{sampleId}_{blastProgram}_paralogs.fasta'
	fh1 = open(outfile1, "w")
	# Trimmed output file:
	fh2 = open(outfile2, "w")
	# Paralogs output file
	fh3 = open(outfile3, "w")


	geneHitDict = {}	# Stores the top matching gene coding or transcritome hit from the BLAST subject
	geneAltHitDict = {} # Stores any other hits found after the top hit
	with open(blast_output_file, "r") as fh:
		for row in fh:
			row = row.rstrip('\n')
			# Split the tsv row, add geneName to key and row to a nested dict with gene coding seq/transcriptome id as key
			rowArray = row.split('\t')
			(refSpecies, refGeneName) = rowArray[0].split('-')
			subjectHspStart=int(rowArray[9]) # int() required for numerical comparison later
			subjectHspEnd=int(rowArray[10])
			qlen = rowArray[8]
			slen = rowArray[11]
			#print(refGeneName)
			if refGeneName not in geneHitDict:
				#print(refGeneName)
				geneHitDict[refGeneName] = {}
				# Indicate that gene coding/transcript sequence has been selected:
				row = row + '\tselected'
				geneHitDict[refGeneName][rowArray[1]] = {}
				geneHitDict[refGeneName][rowArray[1]]['row'] = row
				# Also store the start and end coordinates for trimming the sequences at 5' and 3' ends:
				geneHitDict[refGeneName][rowArray[1]]['lowest_start'] = subjectHspStart
				geneHitDict[refGeneName][rowArray[1]]['highest_end'] = subjectHspEnd
				#print(geneHitDict)
				#print seq to the main outfile in new fasta header format: >geneName <reftarg> <original CDS/transcriptomeId pcid lenHSP qlen slen evalue
				if blastProgram == 'tblastn' or blastProgram == 'blastx':
					# slen is in bases! Makes sense as the input is DNA
					slen = round(int(rowArray[11]) / 3) # Dividing by 3 to fit the rest of the numbers (amino acids)

				#print(refGeneName + " " + rowArray[1])	
				#print(fastaFileBlastDBDict[rowArray[1]].id)
				#print(fastaFileBlastDBDict[rowArray[1]].seq)
				#print(fastaFileBlastDBDict[rowArray[1]].description)
				#lineToWrite = '>' + sampleId + "-" + refGeneName + " " + rowArray[0] + " " + rowArray[1] + " " \
#+ 'pcid=' + rowArray[2] + ' lenHSP=' + rowArray[3] + ' qlen=' + rowArray[8] + ' slen=' + str(slen) + ' evalue=' + rowArray[12] + "\n" + fastaFileBlastDBDict[rowArray[1]].seq + "\n"
				# Now writing as an f-string (much easier syntax! Can also use it directly in a print statement.)
				lineToWrite = f'>{sampleId}-{refGeneName} query={rowArray[0]} subject={rowArray[1]} \
pcid={rowArray[2]} lenTopHSP={rowArray[3]} qlen={rowArray[8]} slen={str(slen)} evalue={rowArray[12]}\n{fastaFileBlastDBDict[rowArray[1]].seq}\n'
				fh1.write(str(lineToWrite))

				
				#print(f"INFO: inital values {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['lowest_start'])}")
				#print(f"INFO: initial values {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['highest_end'])}")
			else:
				# Already found a top hit for this gene.
				# Gathering coordinates for all HSPs in the top hit:
				if rowArray[1] in geneHitDict[refGeneName]:
###27.5.2025 - consider an extra dict level for this step
###if  rowArray[0]/target ref in geneHitDict[refGeneName]
###		if rowArray[1] in geneHitDict[refGeneName][targetRef]:	- still thinking it is Ok as it is the same gene and same blast subject, just not same query gene ortholog
					# Test for lower/higher start/end subject HSP coordinates:
### 28.5.2025 - could also check here whether seq hit is in the rev ori!! - actually do in loop below 
					if geneHitDict[refGeneName][rowArray[1]]['lowest_start'] > subjectHspStart:
						print(f"INFO - before: lower start found for {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['lowest_start'])}; new: {str(subjectHspStart)}")
						geneHitDict[refGeneName][rowArray[1]]['lowest_start'] = subjectHspStart
						print(f"INFO - after: lower start found for {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['lowest_start'])}; new: {str(subjectHspStart)}")
					if geneHitDict[refGeneName][rowArray[1]]['highest_end'] < subjectHspEnd:
						print(f"INFO - before: lower end found for {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['highest_end'])}; new: {str(subjectHspEnd)}")
						geneHitDict[refGeneName][rowArray[1]]['highest_end'] = subjectHspEnd
						print(f"INFO - after: lower end found for {refGeneName} - {rowArray[1]}: existing: {str(geneHitDict[refGeneName][rowArray[1]]['highest_end'])}; new: {str(subjectHspEnd)}")
				# else:
				# 	# Not already seen an alt hit for this gene.
				# 	if rowArray[1] not in geneAltHitDict[refGeneName]:
				# 		geneAltHitDict[refGeneName] = {}
				# 		geneAltHitAltDict[refGeneName][rowArray[1]] = {}
				# 		geneAltHitDict[refGeneName][rowArray[1]]['row'] = row
				# 		# Also store the start and end coordinates for trimming the sequences at 5' and 3' ends:
				# 		geneAltHitDict[refGeneName][rowArray[1]]['lowest_start'] = subjectHspStart
				# 		geneAltHitDict[refGeneName][rowArray[1]]['highest_end'] = subjectHspEnd
				# 	else:
				# 		# Subject already come across so need to collect the lowest/highest start/end coordinates for
				#		# for the remaining HSPs for this alternative hit.

						### Consider to use this method here for use also above: find low_start_high_end()
						### But first, find out whether you can create a method inside another method, if not it needs to go outside!
						# if geneAltHitDict[refGeneName][rowArray[1]]['lowest_start'] > subjectHspStart:
						# 	print(f"INFO - before: lower start found for {refGeneName} - {rowArray[1]}: existing: {str(geneAltHitDict[refGeneName][rowArray[1]]['lowest_start'])}; new: {str(subjectHspStart)}")
						# 	geneAltHitDict[refGeneName][rowArray[1]]['lowest_start'] = subjectHspStart
						# 	(f"INFO - after: lower start found for {refGeneName} - {rowArray[1]}: existing: {str(geneAltHitDict[refGeneName][rowArray[1]]['lowest_start'])}; new: {str(subjectHspStart)}")
						# if geneAltHitDict[refGeneName][rowArray[1]]['highest_end'] < subjectHspEnd:
						# print(f"INFO - before: lower end found for {refGeneName} - {rowArray[1]}: existing: {str(geneAltHitDict[refGeneName][rowArray[1]]['highest_end'])}; new: {str(subjectHspEnd)}")
						# geneAltHitDict[refGeneName][rowArray[1]]['highest_end'] = subjectHspEnd
						# print(f"INFO - after: lower end found for {refGeneName} - {rowArray[1]}: existing: {str(geneAltHitDict[refGeneName][rowArray[1]]['highest_end'])}; new: {str(subjectHspEnd)}")




				### else if already seen gene - suject gene contig before,
				###		Work out if coords need updating for start and end


    # Loop through the geneHitDict which is now filled with top only hits and print out the trimmed version:
	for geneName in geneHitDict:
		print(geneName, "   ", len(geneHitDict[geneName])) 	# just testing whether geneHitDict is an ordered dict - yes it is! I think all dicts are now ordered above a certain python version (3.6+?)
		print(geneHitDict[geneName])

		for subject in geneHitDict[geneName]:
			#print(dir(fastaFileBlastDBDict[subject]))
			#print(vars(fastaFileBlastDBDict[subject]))
			testString = str(fastaFileBlastDBDict[subject].seq)
			print(f'Untrimmed length: {len(testString)}')
			print(f'Untrimmed seq: {testString}')
			lowestStart = geneHitDict[geneName][subject]['lowest_start'] - 1 # might want to check value is not -ve
			highestEnd = geneHitDict[geneName][subject]['highest_end']		# # might want to check value is greater than the end the subject - coudl also reverse the seq!!!!

			### 28.5.2025 - if coords in geneHitDict[geneName][subject]['row'] show query/subj top hit is on opposite strand
			### do not print to fasta  file but print a warning instead

			print(f'lowestStart: {lowestStart}; highestEnd: {highestEnd}')
			print(f'Trimmed length: {len(testString[lowestStart:highestEnd])}')
			print(f'Trimmed seq: {testString[lowestStart:highestEnd]}')
			trimmedSeq = testString[lowestStart:highestEnd]


			rowArray2 = geneHitDict[geneName][subject]['row'].split('\t')
			qlen2 = rowArray2[8]
			slen2 = rowArray2[11]
			lineToWrite2 = f'>{sampleId}-{geneName} query={rowArray2[0]} subject={rowArray2[1]} \
pcid={rowArray2[2]} lenTopHSP={rowArray2[3]} qlen={rowArray2[8]} slen={str(slen2)} evalue={rowArray2[12]}\n{trimmedSeq}\n'
			fh2.write(str(lineToWrite2))
	print(f'End of {geneName}')

	### UPTOHERE 22.5.2025
	### Now tidy up the above code
	### Also have other notes in K67.P3 to move to here
	### Will need to take an average for pcid= and lenHSP=; avpcid= and lenHSPs - or just say it is the result from the top HSP!!!
	###		Take an average of each two values before storing???? biased calc? only a few rounds so probably not important
	###		This might not be good if the start and end coords are from different target refs!!!!
	###		Consider to have an extra level to store speciesId-geneId of target ref!! - 27.5.2025 - still not sure
	###		Taking an average is not appropraite if I'm not using the same query target.
	###			Instead just rename the fastaheader from lenAllHSPs/lenHSP to lenTopHSP - also check stats.txt file
	###			Then in logs have info on filtering method - lowest start highest end subject coords using all availabel hits to targt refs of the same gene 
	### Finally get blast to assess the top 10 hits but I happily I think the code would stay the same.
	###		Actually now uncertain aobut this 
	### 27.5.2025 - test some alignments to see whether trimming looks good
	### NBNB - HORVU3Hr1G038140.15 is in the antisense direction so seq is zero length - need to trap this error
	### Check that gene 5562 that was repeat - both repeats should be present now.
	### Also check the numbers on the fasta header when using tblastn - check lenHSP is now lenAllHSPs!
	### -->Note in the docs that the trimmed seqs are all in frame 1 - very handy; NB - TBLASTN might have had to insert residues into the 
	###		sequence but I think they are just to fit with the query seq, not that they are added to keep the seq in frame - need to confirm this 
	### 27.5.2025 - Would be good to check that no fasta records have no seq afterwards
	### 	For now just print a  warniogn aobut seqs in opposite ori and do a cehck in bash for empty line | grep ^$ - then name file
	### 	if empty line detected,
	###		cat  
	###		mv *.trimmed.fasta *.trimmed.fasta_before_rming_empty_lines
	###		cat *.trimmed.fasta_before_rming_empty_lines | grep -B 1 -v ^$ - actuually this won't work or is dangerous!!!!
	### 28.5.2025 - in program notes say that the ref targets reallt need to be in sense ori - program can't do anything about antisense ori 
	### but it can detect and eliminate sequences whose hits are in the opposite orientation!
	### 1.6.2025 - think about this further - might be useful to be able to detect and revcom seqs if necessary e.g. for rRNA - it shoudl be easy to do
	### 			ands report a warning!!!
	### 6.6.2025 - w.r.t. extracting genes from unannotated genomes, it might be possible to use all subject coords then remove regions falling outide the hit coords
	###				I think this could work - is this how HybPiper does it?
	

	### Then loop through the geneAltHitDict (still to be made) to print out the paralogous genes:
	for geneName in geneHitDict:
		pass # still to devel
	
	fh1.close()
	fh2.close()
	fh3.close()


# Main code:
if method == 'detect_stops':
	detect_stops(option1, option2)

elif method == 'orderAlnByTreeOrder':
	orderAlnByTreeOrder(option1, option2)

elif method == 'retrieve_targets_magic':
	retrieve_targets_magic(option1, option2, option3, option4)

else:
	print('ERROR: you need to specify an existing Python method to use!')
	print('List of popular functions available:')
	print('1. detect_stops')
	print('2. orderAlnByTreeOrder')
	print('3. retrieve_targets_magic')


