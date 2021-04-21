#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import csv
import re


sampleTableFile = sys.argv[1]
treeTipInfoMapfile = sys.argv[2]


csvDict = {}
with open(sampleTableFile, 'r') as csvfile:
	csvreader = csv.reader(csvfile, delimiter=',', quotechar='"') # NB - quotechar seems to be OK with e.g. ,WZVCST Broccoli""WN12-95A"",
	for row in csvreader:
		#print row
		# Removing the sampleId from the list so it doesn't appear in the tree.
		# If still required, user needs to copy this column again into the csv table.
		uniqId = row.pop(0)
		#print uniqId, row
		treeTipInfo = '_'.join(row)
		treeTipInfo = treeTipInfo.rstrip('_')

		# Need  to strip all chars that interfere with Newick format i.e. ( ) : [ ] and ; any others? Check in Newick definition - see K52.p5 note
		# NB - spaces must also be converted to satisfy code in sister script, addTreeTipInfoFromTable_PlusL.py line 26
		treeTipInfo = re.sub(r'[\]\[\)\(:; ]', '_', treeTipInfo)

		if uniqId in csvDict:
			print('ERROR: id in 1st column of csv file is not unique - it MUST be - offending id:', uniqId)
			sys.exit()
		else:
			csvDict[ uniqId ] = treeTipInfo		# NB - don't need to remove the csv file header!


### 29.5.2020 - don't think I need the code just below now

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

# fileNameDict = {}	# Just for determining uniqueness of the filename
# for filePathName in argvList:
# 	#print filePathName
# 	#print os.path.split(filePathName)
# 	fileName = os.path.basename(filePathName)
# 	#print fileName
# 	# Get filenme prefix (no ending):
# 	fields = fileName.split('.')
# 	uniqueSampleId = fields[0]
# 	#print uniqueSampleId

# 	# First test that filename is unique. 
# 	# (It might not be unique if filenames do not originate from a single directory.)
# 	if uniqueSampleId in fileNameDict:
# 		print 'ERROR: this sample name ID is not unique in the set of filenames - it MUST be:', uniqueSampleId
# 		sys.exit()
# 	else:
# 		if uniqueSampleId not in csvDict:
# 			print "WARNING: sample name ID in fasta file not present in table - will use filename in the tree tip: ", uniqueSampleId
# 			# NB - if the filename is not in the csvDict then it MUST still be unique in the csv file sample name IDs, so no real need to check IDs are still unique here.
# 			fh.write(uniqueSampleId + " " + uniqueSampleId + "\n")
# 		else:
# 			fh.write(uniqueSampleId + " " + csvDict[uniqueSampleId] + "\n")

### 27.5.2020 - I've got confused between filename which may or may not comprise of the sample name and the sampleId that's present in the tree.
### It's much easier to just prepare the map file from a user-derived csv file containing the sampleIds on the tree tip instead of also cross-checking with the sample name!
### Just can't check files are unique - but should be doing that earlier.

for key in csvDict:
	fh.write(key + " " + csvDict[key] + "\n")
fh.close()		


