#!/usr/bin/env python

############################
# various_tasks_in_python.py

# Purpose: methods to perform various easy tasks in Python rather than in bash
#
# Usage: an internal pipeline script
#        various_tasks_in_python.py <method_name> <infile> <outfile_prefix>
#
# Method 1: detect_stops
# Method 2: etc

# Author: Paul Bailey

# Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew
############################
import sys
import os
import csv
import re
from Bio import SeqIO					# https://biopython.org/DIST/docs/api/;  http://biopython.org/DIST/docs/tutorial/Tutorial.html;  https://biopython.org/wiki/SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#import __future__ [as name]


method = sys.argv[1]
infile = sys.argv[2]
outfilePrefix = sys.argv[3]


def detect_stops(infile, outfilePrefix):
	''' Detects STOP codons in a protein aligment, remove sequences with > 1 STOP codon and create stats.
	    Assumes that a STOP codon is denoted by a '*' char.

		Input parameters: filename of sequence records in fasta format and an outfile prefix to use
		Assumes that a STOP codon is denoted by a '*' char.

		Test command and input data: /Users/pba10kg/Documents/ProgramFiles/PhylogenomicsPipelines/various_tasks_in_python.py detect_stops 4848.protein.fasta  4848.protein
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
		print record.id, "\n", record.seq, len(record)
		print 'Number of stops in seq: ', record.seq.count('*')
		numbrStops = record.seq.count('*')

		if numbrStops in stopCountr:
			stopCountr[numbrStops] += 1
		else:
			stopCountr[numbrStops] = 1
		#print 'stopCountr[', numbrStops, ']:', stopCountr[numbrStops]

		# If sequence has zero or one stops, print to file (for use in the phylogeny).
		# NB - one stop might be the real stop codon so it should be allowed through.
		#      In current work I don;t think the sequences have the end STOP codon and 
		#      in any case one STOP and maybe one or two more in a sequence seem to be 
		#      in frame aligned well with many other sequences, so should be allowed through.
		### Might want to increase or decrease the number of stops in seqs to print:

		# UPP refuses to align a sequences containing STOP chars so removing them 
		# before printing seqs to file - actually converting them to 'X' chars (for now).
		# (c.f. MAFFT strips any '?' and '*' chars out of the alignment, retains 'X' chars though)
		seqStopsRemoved = re.sub('\*', 'X', str(record.seq))	# NB - record.seq is a Seq object but the full sequence can be returned as a python string, using e.g. str(my_seq).	
		recordNoStops = '>' + record.id + '\n' + seqStopsRemoved + '\n'

		if numbrStops == 0 or numbrStops == 1:
			
			print 'Print record 0/1 ', record 
			#SeqIO.write(record, fh, "fasta")	# Need to use a file handle here if printing to file multiple times, otherwise previous write gets overwritten, not appended
			fh.write(str(recordNoStops))
		else:
			print 'Print record ovr1 ', record 
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


# Main code:
if method == 'detect_stops':
	detect_stops(infile, outfilePrefix)


