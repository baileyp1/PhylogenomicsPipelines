#!/usr/bin/env python

import sys
import os
import csv
from Bio import SeqIO
import re
#import __future__ [as name]


sampleTableFile = sys.argv[1]
treeTipInfoMapfile = sys.argv[2]


csvDict = {}
with open(sampleTableFile, 'r') as csvfile:
	csvreader = csv.reader(csvfile, delimiter=',', quotechar='"') # NB - quotechar seems to be OK with e.g. ,WZVCST Broccoli""WN12-95A"",
	for row in csvreader:					### 17.5.2020 - just realied that row here is actually field of row 
		treeTipInfo= '_'.join(row)
		uniqId = row[0]

		# Need  to strip all chars that interfere with Newick format i.e. ( ) : [ ] and ; - any others? Check in Newick definition - see K52.p5 note
		treeTipInfo = re.sub(r'[\]\[\)\(:;]', '_', treeTipInfo)

		if uniqId in csvDict:
			print 'ERROR: id in 1st column of csv file is not unique - it MUST be - offending id:', uniqId
			sys.exit()
		else:
			csvDict[ uniqId ] = treeTipInfo		# NB - don't need to remove the csv file header!


# Prepare to get sample names from file names and test for uniqueness.
# Make copy of sys.argv:
argvList = list(sys.argv)
argvList.pop(0)		# Removes script name
argvList.pop(0)		# Removes sampleTableFile
argvList.pop(0)		# Removes treeTipInfoMapfile --> now only sample files in list
#print argvList
#lenArgv = len(argvList)
#print lenArgv


# Prepare mapfile for Newick Utilites:
fh = open(treeTipInfoMapfile, 'w')

fileNameDict = {}	# Just for determining uniqueness of the filename
for filePathName in argvList:
	print filePathName
	print os.path.split(filePathName)
	fileName = os.path.basename(filePathName)
	print fileName
	# Get filenme prefix (no ending):
	fields = fileName.split('.')
	uniqueSampleId = fields[0]
	print uniqueSampleId

	# First test that filename is unique. 
	# (It might not be unique if filenames do not originate from a single directory.)
	if uniqueSampleId in fileNameDict:
		print 'ERROR: this sample name ID is not unique in the set of filenames - it MUST be:', uniqueSampleId
		sys.exit()
	else:
		if uniqueSampleId not in csvDict:
			print "ERROR: sample name ID in fasta file not present in table - will use filename in the tree tip: ", uniqueSampleId
			# NB - if the filename is not in the csvDict then it MUST still be unique in the csv file sample name IDs, so no real need to check IDs are still unique here.
			fh.write(uniqueSampleId + " " + uniqueSampleId + "\n")
		else:
			fh.write(uniqueSampleId + " " + csvDict[uniqueSampleId] + "\n")
fh.close()		


