#!/usr/bin/env python

############################
# various_tasks_in_python.py

# Purpose: methods to perform various easy tasks in Python rather than in bash
#
# Usage: an internal pipeline script with a simple interface (no argparse)
#        various_tasks_in_python.py <method_name> <option1 e.g. infile>  <option2 e.g. outfile_prefix>  <option3>  <option4> etc
#
# Author: Paul Bailey
#
# detect_stops
#	Detects STOP codons in a protein aligment, removes sequences with > 1 STOP codon and create stats
#
# orderAlnByTreeOrder():
#	Orders a sequence alignment by the order in a Newick tree file
#
# *** Next method here ***
#
# Copyright (c) 2024 The Board of Trustees of the Royal Botanic Gardens, Kew
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
	outfile2 = outfilePrefix + '.log'
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


# Main code:
if method == 'detect_stops':
	detect_stops(option1, option2)

elif method == 'orderTableByTreeTips':
	orderTableByTreeTips(option1, option2)

elif method == 'orderAlnByTreeOrder':
	orderAlnByTreeOrder(option1, option2)

else:
	print('ERROR: you need to specify an existing Python method to use!')


