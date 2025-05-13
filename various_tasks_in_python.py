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

	***** 7.5.2025 - REMEMBER TO UPDATE STABLE VERSION TO BIN FOLDER *****

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
	outfile=sampleId + '.fasta'
	fh1 = open(outfile, "w")


	geneHitDict = {}	# Stores the top matching gene coding or transcritome hit from the BLAST subject
	with open(blast_output_file, "r") as fh:
		for row in fh:
			row = row.rstrip('\n')
			# Split the tsv row, add geneName to key and row to a nested dict with gene coding seq/transcriptome id as key
			rowArray = row.split('\t')
			(refSpecies, refGeneName) = rowArray[0].split('-')
			#print(refGeneName)
			if refGeneName not in geneHitDict:
				#print(refGeneName)
				geneHitDict[refGeneName] = {}
				# Indicate that gene coding/transcript sequence has been selected:
				row = row + '\tselected'
				geneHitDict[refGeneName][rowArray[1]] = row
				#print(geneHitDict)
				#print seq to the main outfile in new fasta header format: >geneName <reftarg> <original CDS/transcriptomeId pcid lenHSP qlen slen evalue
				qlen = rowArray[8]
				slen = rowArray[11]
				if blastProgram == 'tblastn' or blastProgram == 'blastx':
					# slen is in bases! Makes sense as the input is DNA
					slen = round(int(rowArray[11]) / 3) # Dividing by 3 to fit the rest of the numbers (amino acids)


				#print(refGeneName + " " + rowArray[1])	
				#print(fastaFileBlastDBDict[rowArray[1]].id)
				#print(fastaFileBlastDBDict[rowArray[1]].seq)
				#print(fastaFileBlastDBDict[rowArray[1]].description)
				#print('>' + refGeneName + " " + rowArray[0] + " " + rowArray[1] + "\n" + fastaFileBlastDBDict[rowArray[1]].seq + "\n")
				#lineToWrite = '>' + sampleId + "-" + refGeneName + " " + rowArray[0] + " " + rowArray[1] + " " \
#+ 'pcid=' + rowArray[2] + ' lenHSP=' + rowArray[3] + ' qlen=' + rowArray[8] + ' slen=' + str(slen) + ' evalue=' + rowArray[12] + "\n" + fastaFileBlastDBDict[rowArray[1]].seq + "\n"
				# Now writing as an f-string (much easier syntax!)
				lineToWrite = f'>{sampleId}-{refGeneName} query={rowArray[0]} subject={rowArray[1]} \
pcid={rowArray[2]} lenHSP={rowArray[3]} qlen={rowArray[8]} slen={str(slen)} evalue={rowArray[12]}\n{fastaFileBlastDBDict[rowArray[1]].seq}\n'
				fh1.write(str(lineToWrite))

				# Now store the gene coding seq/transcriptome hit in a separate hash for testing
				# whether the same id appears 
			  
			else:

				### 8.5.2025 - if collecting start end coord of top hit, can do that here I think

				# Now add all the remaining subject hits for the same gene (but different ref target))
				### This is not an essential step but would store the different hits for comparison.
				### Instead have improved the sort step in the bash script so probably don't need to know about these.
				if rowArray[1] not in geneHitDict[refGeneName]:
					geneHitDict[refGeneName][rowArray[1]] = {}
					geneHitDict[refGeneName][rowArray[1]] = row
				else:
					geneHitDict[refGeneName][rowArray[1]] = row
					### This is still incorrect - I need to append rows to an array that have the same subj hit.
					### As it is, the last subj overwriting any previous hit, including the selected one!!!
					### Actually, also onlt need to store the top row hit for each subj hit come across
					### These could be paralogs which we might want.
					###		A more thorough extension of this is to assess the top 10 hits but I happily I think
					###		the code would stay the same.
	####close(fh1)

    # Loop through the geneHitDict which is now filled with selected and failed hits:
	#for geneName in geneHitDict:
		#print(geneName) 	# just testing whether geneHitDict is an ordered dict - yes it is! I think all dicts are now ordered above a certain python version (3.6+?)

		
		#if len(geneHitDict[geneName]) > 0:
			#print(geneHitDict[geneName])
    		###	print assosicated stats to a file for each selected hit
    		### THEN: If geneHitDict[geneName] contains more than 2 entries, print to a warnings file with the output stats 



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


